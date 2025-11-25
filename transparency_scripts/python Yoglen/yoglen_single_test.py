#!/usr/bin/env python3
"""
Two-Integral Change Point Model - SINGLE TEST VERSION
Modified from change_point_proper_staged.py to run on individual tests
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
from scipy.optimize import minimize
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import BSpline, make_interp_spline
import warnings
from scipy.stats import linregress
import time
import sys
import os
import argparse
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy import stats

# Import iterative reweighting module
try:
    from chemfarm_analysis.iterative_reweighting import fit_linear_model_irls
except ImportError:
    # For local development/testing
    try:
        from iterative_reweighting import fit_linear_model_irls
    except ImportError:
        print(
            "Warning: iterative_reweighting module not found. Iterative reweighting will not be available."
        )
import scipy.optimize

warnings.filterwarnings("ignore")


def get_test_track_name(test_code):
    """
    Convert test code to test_track format used in bootstrap data.
    
    Args:
        test_code: Test code (e.g., "100900" or "marker.BMI")
    
    Returns:
        Track name (e.g., "lab.100900" or "marker.BMI")
    """
    test_code_str = str(test_code)
    if test_code_str.startswith("marker."):
        return test_code_str
    else:
        return f"lab.{test_code_str}"


def wilson_confidence_interval(y, n, alpha=0.05):
    """
    Calculate Wilson confidence interval for binomial proportion

    Args:
        y: number of successes
        n: number of trials
        alpha: significance level (default 0.05 for 95% CI)

    Returns:
        tuple: (lower_bound, upper_bound)
    """
    z = stats.norm.ppf(1 - alpha / 2)  # 1.96 for 95% CI
    phat = y / n

    # Wilson formula
    den = 1 + (z**2) / n
    ctr = (phat + (z**2) / (2 * n)) / den
    hw = (z / den) * np.sqrt(phat * (1 - phat) / n + (z**2) / (4 * n**2))

    lower = np.maximum(0, ctr - hw)
    upper = np.minimum(1, ctr + hw)

    return lower, upper


def robust_weights_huber(residuals, c=1.345):
    """
    Calculate Huber robust weights

    Args:
        residuals: Pearson residuals
        c: tuning constant (default 1.345)

    Returns:
        robust weights
    """
    return np.minimum(1.0, c / np.abs(residuals))


def robust_weights_tukey(residuals, c=4.685):
    """
    Calculate Tukey biweight robust weights

    Args:
        residuals: Pearson residuals
        c: tuning constant (default 4.685)

    Returns:
        robust weights
    """
    abs_residuals = np.abs(residuals)
    weights = np.zeros_like(residuals)
    mask = abs_residuals <= c
    weights[mask] = (1 - (abs_residuals[mask] / c) ** 2) ** 2
    return weights


def logit_transform_with_robust_weights(test_values, test_ns, max_iterations=3):
    """
    Apply logit transformation with robust binomial GLM weighting

    Args:
        test_values: percentage values (0-100)
        test_ns: sample sizes
        max_iterations: maximum number of robust iterations

    Returns:
        tuple: (transformed_values, robust_weights, iteration_info)
    """
    # Convert percentages to counts and probabilities
    probabilities = test_values / 100.0
    counts = np.round(probabilities * test_ns).astype(int)

    # Clip to avoid log(0) and log(1)
    epsilon = 1e-10
    probabilities_clipped = np.clip(probabilities, epsilon, 1 - epsilon)

    # Initial logit transformation
    logit_values = np.log(probabilities_clipped / (1 - probabilities_clipped))

    # Initial weights (uniform)
    robust_weights = np.ones_like(logit_values)

    iteration_info = []

    for iteration in range(max_iterations):
        # Calculate Pearson residuals for binomial model
        # For logit scale: approximate binomial variance
        p_hat = 1 / (1 + np.exp(-logit_values))
        p_hat = np.clip(p_hat, epsilon, 1 - epsilon)

        # Pearson residuals on logit scale
        # Approximate: SE(logit(p)) ≈ 1/√(n*p*(1-p))
        se_logit = 1 / np.sqrt(test_ns * p_hat * (1 - p_hat))
        se_logit = np.clip(se_logit, epsilon, None)

        # Observed vs expected residuals
        expected_counts = p_hat * test_ns
        observed_counts = counts
        pearson_residuals = (observed_counts - expected_counts) / np.sqrt(
            expected_counts * (1 - p_hat)
        )

        # Calculate robust weights
        new_weights = robust_weights_huber(pearson_residuals, c=1.345)

        # Check for convergence
        weight_change = np.mean(np.abs(new_weights - robust_weights))
        robust_weights = new_weights

        iteration_info.append(
            {
                "iteration": iteration + 1,
                "weight_change": weight_change,
                "mean_weight": np.mean(robust_weights),
                "min_weight": np.min(robust_weights),
                "max_weight": np.max(robust_weights),
            }
        )

        # Early convergence check
        if weight_change < 0.01:
            break

    return logit_values, robust_weights, iteration_info


def calculate_binomial_error_bars(test_values, test_ns, use_logit_transform=False):
    """
    Calculate error bars for plotting - binomial confidence intervals for logit, standard error for regular data

    Args:
        test_values: data values (logit-transformed if use_logit_transform=True)
        test_ns: sample sizes
        use_logit_transform: whether data is logit-transformed

    Returns:
        tuple: (lower_error, upper_error) for plotting
    """
    if test_ns is None or len(test_ns) == 0:
        return None, None

    if use_logit_transform:
        # For logit-transformed data, use binomial confidence intervals
        # Convert logit values back to probabilities for CI calculation
        epsilon = 1e-10
        probabilities = 1 / (1 + np.exp(-test_values))
        probabilities = np.clip(probabilities, epsilon, 1 - epsilon)

        # Convert to counts for Wilson CI
        counts = np.round(probabilities * test_ns).astype(int)

        # Calculate Wilson confidence intervals
        lower_prob, upper_prob = wilson_confidence_interval(counts, test_ns, alpha=0.05)

        # Transform back to logit scale
        lower_prob = np.clip(lower_prob, epsilon, 1 - epsilon)
        upper_prob = np.clip(upper_prob, epsilon, 1 - epsilon)

        lower_logit = np.log(lower_prob / (1 - lower_prob))
        upper_logit = np.log(upper_prob / (1 - upper_prob))

        # Calculate error bars (distance from mean)
        lower_error = test_values - lower_logit
        upper_error = upper_logit - test_values

        return lower_error, upper_error
    else:
        # For regular data, use standard error calculation
        if test_values is None or len(test_values) == 0:
            return None, None

        # Use a simple standard error approximation
        # This is a fallback - in practice, test_sds should be provided
        se = 1.0 / np.sqrt(test_ns + 1)  # Avoid division by zero
        return se, se


##########################################################################
## CONFIGURATION ----------------------------------------------------------
##########################################################################

# Default configuration
DEFAULT_CONFIG = {
    "INFILE": "healthy_menopause_aggregated_data.csv",
    "SURV_FILE": "example_survival_curve.csv",
    # "N_CP_COMBINATIONS": 128,  # Balanced for CPU utilization and memory usage on 32 cores
    "N_CP_COMBINATIONS": 1,  # Testing with 16 combinations
    "OUTPUT_DIR": "change_point_outputs",
    # Grid search parameters
    "first_change_point_range_max": -5,
    "first_change_point_range_min": -20,
    "second_change_point_range_min": 4,
    "second_change_point_range_max": 10,
    # Smoothing parameters
    "ENABLE_GAUSSIAN_SMOOTHING": True,
    "GAUSSIAN_BANDWIDTH": 1.5,
    # Integration accuracy parameters
    # "INTEGRAL_DT": 0.5,
    "INTEGRAL_DT": 0.2,
    # "INTEGRAL_DT": 5,
    # "INTEGRAL_DT": 10,
    # Staged optimization parameters
    # Age filtering parameters
    "MIN_AGE": 25,  # Minimum age to include in analysis (None = no minimum)
    "MAX_AGE": None,  # Maximum age to include in analysis (None = no maximum)
    # Parallel processing parameters
    "N_WORKERS": 32,  # Number of parallel workers
    "THREADS_PER_WORKER": 6,  # Threads per worker
    # Coefficient of variation settings
    "USE_CV": True,  # Whether to use coefficient of variation for error bars and weighted R-squared
    # Logit transformation settings
    "USE_LOGIT_TRANSFORM": False,  # Whether to apply logit transformation to data (for probabilities/percentages)
    # Change point model settings
    "SINGLE_CP_ONLY": False,  # Whether to use only single change point models
    "CP_IMPROVEMENT_THRESHOLD": 0.01,  # Minimum R² improvement to consider two change points
    "SINGLE_CP_BEFORE_ONLY": False,  # Whether to use only single change point models before menopause (default: False, try both)
    # Slope penalization parameters
    "ENABLE_SLOPE_PENALIZATION": True,  # Enable slope penalization for extreme opposite slopes
    "SLOPE_PENALTY_WEIGHT": 0.0,  # Weight for slope penalty in objective function (increased from 1.0)
    "SLOPE_PENALTY_THRESHOLD_RATIO": 0.05,  # Threshold as ratio of jump magnitude (reduced from 0.2 to 0.05 = 5%)
    "SLOPE_PENALTY_MIN_THRESHOLD": 0.01,  # Minimum absolute slope threshold regardless of jump size (reduced from 0.02)
    "SLOPE_PENALTY_DISTANCE_THRESHOLD": 8.0,  # Only penalize if change point is within this distance of menopause (increased from 6.0)
    # Statistical significance testing
    "ENABLE_SIGNIFICANCE_TESTING": True,  # Enable statistical significance testing
    "SIGNIFICANCE_ALPHA": 0.05,  # Alpha level for significance tests
    "MIN_R_SQUARED_IMPROVEMENT": 0.001,  # Minimum R² improvement for considering significant
    "MIN_ABSOLUTE_R_SQUARED": 0.01,  # Minimum absolute R² for the change point model
    # Bootstrap confidence intervals for change points
    "ENABLE_CP_CONFIDENCE_INTERVALS": True,  # Enable confidence interval calculation for change points
    "CP_CONFIDENCE_LEVEL": 0.95,  # Confidence level for change point positions
}

##########################################################################
## HELPER FUNCTIONS -------------------------------------------------------
##########################################################################


def calculate_standard_error(test_sds, test_ns):
    """

    This gives a much more appropriate measure of measurement precision than
    the population SD of z-scores.

    Parameters:
    - test_sds: Standard deviation values from bootstrap data
    - test_ns: Sample size values for each age group

    Returns:
    - standard_errors: Standard error values as SD / √n
    """
    if test_sds is None or test_ns is None or len(test_sds) == 0:
        return test_sds

    # Handle edge cases where n <= 0
    valid_mask = test_ns > 0
    if not np.any(valid_mask):
        # If no valid sample sizes, return small constant errors
        return np.full_like(test_sds, 0.01)

    # Calculate standard error: SE = SD / √n
    standard_errors = np.full_like(test_sds, np.nan)
    standard_errors[valid_mask] = test_sds[valid_mask] / np.sqrt(test_ns[valid_mask])

    # Fill invalid values with small constant
    standard_errors[~valid_mask] = 0.01

    # Ensure SEs are not zero (minimum 0.001 for numerical stability)
    standard_errors = np.maximum(standard_errors, 0.001)

    return standard_errors


def calculate_slope_penalty(slope_to_meno, jump_magnitude, cp_distance_to_meno):
    """
    Calculate penalty for extreme slopes that are opposite to jump direction.

    Parameters:
    - slope_to_meno: slope from change point to menopause
    - jump_magnitude: magnitude of jump at menopause (positive = upward jump)
    - cp_distance_to_meno: distance from change point to menopause

    Returns:
    - penalty: penalty value (0 if no penalty needed)
    """
    # Reduced debug output to prevent I/O contention in parallel processing
    # print(f"         🔍 CHECKING SLOPE PENALTY: slope={slope_to_meno:.3f}, jump={jump_magnitude:.3f}, distance={cp_distance_to_meno:.1f}")

    if not DEFAULT_CONFIG["ENABLE_SLOPE_PENALIZATION"]:
        # print(f"         ❌ SLOPE PENALIZATION DISABLED")
        return 0.0

    # Only apply penalty if change point is close to menopause
    if abs(cp_distance_to_meno) > DEFAULT_CONFIG["SLOPE_PENALTY_DISTANCE_THRESHOLD"]:
        # print(f"         ❌ CP TOO FAR: distance={abs(cp_distance_to_meno):.1f} > threshold={DEFAULT_CONFIG['SLOPE_PENALTY_DISTANCE_THRESHOLD']}")
        return 0.0

    # Calculate threshold relative to jump magnitude
    abs_jump = abs(jump_magnitude)
    if abs_jump > 0:
        # Threshold is a ratio of the jump magnitude
        relative_threshold = DEFAULT_CONFIG["SLOPE_PENALTY_THRESHOLD_RATIO"] * abs_jump
        # But not less than minimum absolute threshold
        slope_threshold = max(
            relative_threshold, DEFAULT_CONFIG["SLOPE_PENALTY_MIN_THRESHOLD"]
        )
    else:
        # If no jump, use minimum threshold
        slope_threshold = DEFAULT_CONFIG["SLOPE_PENALTY_MIN_THRESHOLD"]

    # print(f"         📊 THRESHOLD CALC: abs_jump={abs_jump:.3f}, ratio={DEFAULT_CONFIG['SLOPE_PENALTY_THRESHOLD_RATIO']}, threshold={slope_threshold:.3f}")

    # Check if slope and jump are in opposite directions
    opposite_directions = (slope_to_meno * jump_magnitude) < 0  # Different signs

    if opposite_directions:
        # Calculate the magnitude of the opposite slope relative to jump
        slope_magnitude = abs(slope_to_meno)
        jump_magnitude_abs = abs(jump_magnitude)

        # Calculate how much of the opposite direction occurs before the jump
        # This is the key insight: we want to penalize when the opposite movement
        # is comparable in magnitude to the jump itself
        if jump_magnitude_abs > 0:
            # Calculate the total opposite movement from CP to menopause
            cp_distance_abs = abs(cp_distance_to_meno)
            total_opposite_movement = slope_magnitude * cp_distance_abs

            # Calculate the ratio of opposite movement to jump
            opposite_to_jump_ratio = total_opposite_movement / jump_magnitude_abs

            # PENALIZE when opposite movement is significant compared to jump
            # This prevents extreme "down then up" or "up then down" patterns

            # PENALIZE when opposite movement is significant compared to jump
            # This prevents extreme "down then up" or "up then down" patterns

            # CUTOFF: No penalty if opposite movement is less than 30% of jump
            if opposite_to_jump_ratio < 0.3:
                # Small opposite movements are biologically reasonable - no penalty
                return 0.0

            # Base penalty for significant opposite directions
            base_penalty = 0.0

            # MUCH STRONGER PENALTIES to prevent biologically implausible patterns
            if opposite_to_jump_ratio > 0.3:  # Opposite movement > 30% of jump
                base_penalty += 500.0 * opposite_to_jump_ratio  # Increased from 50.0

            if opposite_to_jump_ratio > 0.5:  # Opposite movement > 50% of jump
                base_penalty += 2000.0 * opposite_to_jump_ratio  # Increased from 200.0

            if opposite_to_jump_ratio > 0.8:  # Opposite movement > 80% of jump
                base_penalty += (
                    10000.0 * opposite_to_jump_ratio
                )  # Increased from 1000.0

            # Additional penalty for very steep opposite slopes regardless of jump size
            if (
                slope_magnitude > 0.05
            ):  # Steep opposite slope (>5% per year) - lowered threshold
                base_penalty += 1000.0 * slope_magnitude  # Increased from 100.0

            if slope_magnitude > 0.1:  # Very steep opposite slope (>10% per year)
                base_penalty += 5000.0 * slope_magnitude  # Increased from 500.0

            if slope_magnitude > 0.15:  # Extremely steep opposite slope (>15% per year)
                base_penalty += 20000.0 * slope_magnitude  # New penalty level

            # Distance factor: stronger penalty for CPs closer to menopause
            distance_factor = max(
                0.1,
                DEFAULT_CONFIG["SLOPE_PENALTY_DISTANCE_THRESHOLD"]
                - abs(cp_distance_to_meno),
            )

            penalty = base_penalty * distance_factor

            # Debug output - commented out to reduce spam
            # print(f"         🚨 BIOLOGICAL PENALTY: slope={slope_to_meno:.3f}, jump={jump_magnitude:.3f}, opposite_ratio={opposite_to_jump_ratio:.3f}, penalty={penalty:.6f}")

            return penalty

        else:
            # No jump, but still penalize steep opposite slopes
            if slope_magnitude > 0.05:  # Steep slope without jump
                return 50.0 * slope_magnitude

    return 0.0

    # print(f"         ✅ NO PENALTY APPLIED")
    return 0.0


def calculate_f_statistic(
    rss_simple, rss_complex, n_params_simple, n_params_complex, n_data
):
    """
    Calculate F-statistic for model comparison.

    Parameters:
    - rss_simple: residual sum of squares for simpler model
    - rss_complex: residual sum of squares for more complex model
    - n_params_simple: number of parameters in simpler model
    - n_params_complex: number of parameters in more complex model
    - n_data: number of data points

    Returns:
    - f_stat: F-statistic
    - p_value: p-value for the F-test
    """
    from scipy.stats import f

    df_complex = n_data - n_params_complex
    df_simple = n_data - n_params_simple
    df_diff = n_params_complex - n_params_simple

    if df_complex <= 0 or rss_complex <= 0 or rss_simple <= rss_complex:
        return np.nan, 1.0  # Invalid or no improvement

    f_stat = ((rss_simple - rss_complex) / df_diff) / (rss_complex / df_complex)
    p_value = 1 - f.cdf(f_stat, df_diff, df_complex)

    return f_stat, p_value


def assess_model_quality(
    test_values, predictions, linear_predictions, n_params, linear_n_params
):
    """
    Assess the quality and significance of a change point model.

    Returns:
    - quality_metrics: dictionary with various quality metrics
    """
    n_data = len(test_values)

    # Calculate residuals and R-squared
    residuals = test_values - predictions
    linear_residuals = test_values - linear_predictions

    rss = np.sum(residuals**2)
    linear_rss = np.sum(linear_residuals**2)
    tss = np.sum((test_values - np.mean(test_values)) ** 2)

    r_squared = 1 - rss / tss if tss > 0 else 0
    linear_r_squared = 1 - linear_rss / tss if tss > 0 else 0
    improvement = r_squared - linear_r_squared

    # Calculate F-statistic for significance test
    f_stat, p_value = calculate_f_statistic(
        linear_rss, rss, linear_n_params, n_params, n_data
    )

    # Calculate AIC and BIC
    aic = n_data * np.log(rss / n_data) + 2 * n_params
    linear_aic = n_data * np.log(linear_rss / n_data) + 2 * linear_n_params

    bic = n_data * np.log(rss / n_data) + np.log(n_data) * n_params
    linear_bic = n_data * np.log(linear_rss / n_data) + np.log(n_data) * linear_n_params

    # Overall quality assessment
    is_significant = (
        p_value < DEFAULT_CONFIG["SIGNIFICANCE_ALPHA"]
        if not np.isnan(p_value)
        else False
    )
    meets_min_improvement = improvement > DEFAULT_CONFIG["MIN_R_SQUARED_IMPROVEMENT"]
    meets_min_r_squared = r_squared > DEFAULT_CONFIG["MIN_ABSOLUTE_R_SQUARED"]

    quality_metrics = {
        "r_squared": r_squared,
        "linear_r_squared": linear_r_squared,
        "improvement": improvement,
        "f_statistic": f_stat,
        "p_value": p_value,
        "is_significant": is_significant,
        "meets_min_improvement": meets_min_improvement,
        "meets_min_r_squared": meets_min_r_squared,
        "aic": aic,
        "linear_aic": linear_aic,
        "bic": bic,
        "linear_bic": linear_bic,
        "aic_improvement": linear_aic - aic,
        "bic_improvement": linear_bic - bic,
        "rss": rss,
        "linear_rss": linear_rss,
        "residuals_std": np.std(residuals),
        "linear_residuals_std": np.std(linear_residuals),
        "is_reasonable": is_significant
        and meets_min_improvement
        and meets_min_r_squared,
    }

    return quality_metrics


class ContinuousGaussianSmoothedSurvival:
    """Truly continuous Gaussian smoothed survival function"""

    def __init__(self, survival_data, bandwidth=1.5):
        self.original_ages = survival_data["age"].values
        self.original_S = survival_data["S"].values
        self.bandwidth = bandwidth

    def __call__(self, age):
        """Evaluate the smoothed survival function at any age(s)"""
        age = np.atleast_1d(age)
        smoothed_S = np.zeros_like(age, dtype=float)

        for i, target_age in enumerate(age):
            weights = np.exp(
                -0.5 * ((self.original_ages - target_age) / self.bandwidth) ** 2
            )
            weights = weights / np.sum(weights)
            smoothed_S[i] = np.sum(weights * self.original_S)

        if len(age) > 1:
            for i in range(1, len(smoothed_S)):
                if age[i] > age[i - 1] and smoothed_S[i] > smoothed_S[i - 1]:
                    smoothed_S[i] = smoothed_S[i - 1]

        smoothed_S = np.clip(smoothed_S, 0, 1)
        return smoothed_S if len(smoothed_S) > 1 else smoothed_S[0]


def apply_gaussian_smoothing_continuous(survival_data, bandwidth=1.5, verbose=True):
    """Apply truly continuous Gaussian kernel smoothing"""
    if not DEFAULT_CONFIG["ENABLE_GAUSSIAN_SMOOTHING"]:
        if verbose:
            print("Gaussian smoothing disabled, returning original data")
        return survival_data

    if verbose:
        print(f"Applying continuous Gaussian smoothing with bandwidth = {bandwidth}")

    smooth_func = ContinuousGaussianSmoothedSurvival(survival_data, bandwidth)

    min_age = survival_data["age"].min()
    max_age = survival_data["age"].max()
    fine_ages = np.arange(
        min_age, max_age + DEFAULT_CONFIG["INTEGRAL_DT"], DEFAULT_CONFIG["INTEGRAL_DT"]
    )
    smoothed_S = smooth_func(fine_ages)

    smoothed_data = pd.DataFrame({"age": fine_ages, "S": smoothed_S})
    smoothed_data.smooth_func = smooth_func

    if verbose:
        print(
            f"Continuous smoothing completed. Data points: {len(survival_data)} -> {len(smoothed_data)} (sampled)"
        )
        print(f"✨ Function is truly continuous - can evaluate at ANY age")

    return smoothed_data


# Global integration grid cache - process-local, no lock needed
_integration_cache = {}


def get_integration_grid(survival_data):
    """Get or create cached integration grid for faster optimization"""
    cache_key = f"{len(survival_data)}_{survival_data['age'].min()}_{survival_data['age'].max()}"

    # Check cache without lock - each process has its own cache after fork
    if cache_key not in _integration_cache:
        # Removed print to eliminate I/O blocking in parallel workers
        start_time = time.time()

        age_grid = np.arange(
            survival_data["age"].min(),
            survival_data["age"].max() + DEFAULT_CONFIG["INTEGRAL_DT"],
            DEFAULT_CONFIG["INTEGRAL_DT"],
        )
        S_interp = np.interp(age_grid, survival_data["age"], survival_data["S"])

        # Calculate event density with proper alignment
        event_density = -np.diff(S_interp)
        event_density = np.concatenate(
            [event_density, [0]]
        )  # Add zero at the end to match age_grid length
        event_density[event_density < 0] = 0

        total_density = np.sum(event_density) * DEFAULT_CONFIG["INTEGRAL_DT"]
        if total_density > 0:
            event_density = event_density / total_density
        else:
            event_density = np.ones_like(event_density) / len(event_density)

        # Safety check: ensure arrays have the same length
        if len(age_grid) != len(event_density):
            # Removed prints to eliminate I/O blocking in parallel workers
            # Fix by truncating the longer array
            min_length = min(len(age_grid), len(event_density))
            age_grid = age_grid[:min_length]
            event_density = event_density[:min_length]

        _integration_cache[cache_key] = {
            "age_grid": age_grid,
            "event_density": event_density,
            "dt": DEFAULT_CONFIG["INTEGRAL_DT"],
        }

    return _integration_cache[cache_key]


def linear_menopause_biological_kernel(t, a, b, c, d, survival_data=None):
    """Linear menopause model biological kernel"""
    if t < 0:
        return c + d * t
    else:
        return a + b * t


def fit_linear_menopause_model(
    test_ages,
    test_values,
    survival_data,
    test_sds=None,
    use_iterative_reweighting=False,
):
    """Fit the linear menopause model for comparison"""
    print("🔬 FITTING LINEAR MENOPAUSE MODEL (for comparison)")

    menopause_params = []
    for age in test_ages:
        survival_prob = np.interp(age, survival_data["age"], survival_data["S"])
        meno_prob = 1 - survival_prob

        early_ages = survival_data["age"][survival_data["age"] <= age]
        early_probs = 1 - survival_data["S"][survival_data["age"] <= age]
        early_probs_diff = np.diff(np.concatenate([[0], early_probs]))

        if np.sum(early_probs_diff) > 0:
            mean_meno_age_earlier = np.sum(early_ages * early_probs_diff) / np.sum(
                early_probs_diff
            )
        else:
            mean_meno_age_earlier = age

        later_ages = survival_data["age"][survival_data["age"] >= age]
        later_probs = 1 - survival_data["S"][survival_data["age"] >= age]
        later_probs_diff = np.diff(np.concatenate([later_probs, [1]]))

        if np.sum(later_probs_diff) > 0:
            mean_meno_age_later = np.sum(later_ages * later_probs_diff) / np.sum(
                later_probs_diff
            )
        else:
            mean_meno_age_later = age

        menopause_params.append(
            {
                "age": age,
                "p": meno_prob,
                "mean_meno_age_earlier": mean_meno_age_earlier,
                "mean_meno_age_later": mean_meno_age_later,
            }
        )

    # Create design matrix
    data_matrix = []
    for i, params in enumerate(menopause_params):
        p = params["p"]
        age = params["age"]
        pm_early = p * (age - params["mean_meno_age_earlier"])
        n1p_late = (1 - p) * (age - params["mean_meno_age_later"])
        data_matrix.append([1, p, pm_early, n1p_late])

    X = np.array(data_matrix)
    y = np.array(test_values)

    # Solve linear regression
    if use_iterative_reweighting:
        # Use iterative reweighted least squares if available
        try:
            print("   🔄 Using iterative reweighted least squares for linear model")
            # Initial weights based on standard deviations if available
            initial_weights = None
            if test_sds is not None and np.any(test_sds > 0):
                initial_weights = 1.0 / (test_sds**2)
                initial_weights = initial_weights / np.mean(
                    initial_weights
                )  # Normalize

            # Fit using IRLS
            coefficients, fitted_values, residuals_vec, final_weights = (
                fit_linear_model_irls(
                    X,
                    y,
                    initial_weights=initial_weights,
                    max_iter=10,
                    weight_func="huber",
                )
            )

            # Calculate weighted metrics
            weighted_rss = np.sum(final_weights * residuals_vec**2)
            weighted_tss = np.sum(
                final_weights * (y - np.average(y, weights=final_weights)) ** 2
            )
            r_squared = 1 - weighted_rss / weighted_tss if weighted_tss > 0 else 0
            rss = weighted_rss

            print(f"   ✅ IRLS completed with {len(final_weights)} weights")

        except (NameError, ImportError) as e:
            print(f"   ⚠️ Error using iterative reweighting: {e}")
            print("   ⚠️ Falling back to standard weighted least squares")
            use_iterative_reweighting = False

    # Standard weighted or unweighted least squares
    if not use_iterative_reweighting:
        if test_sds is not None and np.any(test_sds > 0):
            # Weighted least squares using standard errors
            weights = 1.0 / (test_sds**2)
            W = np.diag(weights)
            # Solve: (X^T W X) beta = X^T W y
            XTWX = X.T @ W @ X
            XTWy = X.T @ W @ y

            try:
                coefficients = np.linalg.solve(XTWX, XTWy)
            except np.linalg.LinAlgError:
                # Singular matrix - fall back to lstsq
                print("   ⚠️ Singular matrix, using lstsq for linear model")
                coefficients, _, _, _ = np.linalg.lstsq(X, y, rcond=None)

            # Calculate weighted predictions and residuals
            fitted_values = X @ coefficients
            residuals_vec = y - fitted_values
            weighted_rss = np.sum(weights * residuals_vec**2)
            weighted_tss = np.sum(weights * (y - np.average(y, weights=weights)) ** 2)
            r_squared = 1 - weighted_rss / weighted_tss if weighted_tss > 0 else 0
            rss = weighted_rss
        else:
            # Unweighted least squares
            coefficients, residuals, rank, s = np.linalg.lstsq(X, y, rcond=None)
            fitted_values = X @ coefficients
            residuals_vec = y - fitted_values
            rss = np.sum(residuals_vec**2)
            tss = np.sum((y - np.mean(y)) ** 2)
            r_squared = 1 - rss / tss if tss > 0 else 0

    # Extract biological parameters
    intercept = coefficients[0]  # c
    coef_p = coefficients[1]  # a - c
    coef_pm_early = coefficients[2]  # b
    coef_n1p_late = coefficients[3]  # d

    c_param = intercept
    a_param = coef_p + c_param
    b_param = coef_pm_early
    d_param = coef_n1p_late

    linear_model = {
        "a": a_param,
        "b": b_param,
        "c": c_param,
        "d": d_param,
        "r_squared": r_squared,
        "rss": rss,
        "fitted_values": fitted_values,
        "residuals": residuals_vec,
        "n_points": len(y),
        "menopause_params": menopause_params,
    }

    print(f"   ✅ Linear menopause model fitted:")
    print(f"   a (post-meno intercept): {a_param:.6f}")
    print(f"   b (post-meno slope): {b_param:.6f}")
    print(f"   c (pre-meno intercept): {c_param:.6f}")
    print(f"   d (pre-meno slope): {d_param:.6f}")
    print(f"   R² = {r_squared:.6f}")
    print(
        f"   Jump at menopause: {a_param - c_param:.6f} (from {c_param:.3f} to {a_param:.3f})"
    )

    return linear_model


def biological_kernel_single_cp(
    t, slope1, val_cp1, val_before_jump, val_after_jump, slope2, cp1
):
    # No parameter validation needed

    if cp1 < 0:
        if t < cp1:
            return val_cp1 + slope1 * (t - cp1)
        elif t < 0:
            if cp1 == 0:
                return val_before_jump
            slope_to_meno = (val_before_jump - val_cp1) / (0 - cp1)
            return val_cp1 + slope_to_meno * (t - cp1)
        elif t == 0:
            return val_after_jump
        else:
            slope_from_meno = slope2
            return val_after_jump + slope_from_meno * (t - 0)
    else:
        if t < 0:
            return val_before_jump + slope1 * (t - 0)
        elif t == 0:
            return val_before_jump
        elif t <= cp1:
            if cp1 == 0:
                return val_after_jump
            slope_from_meno = (val_cp1 - val_after_jump) / (cp1 - 0)
            return val_after_jump + slope_from_meno * t
        else:
            return val_cp1 + slope2 * (t - cp1)


def biological_kernel(
    t, slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, slope2, cp1, cp2
):
    """Biological kernel function: piecewise linear with jump at t=0"""
    # No parameter validation needed

    if t < cp1:
        return val_cp1 + slope1 * (t - cp1)
    elif t < 0:
        if cp1 == 0:
            return val_before_jump
        slope_to_meno = (val_before_jump - val_cp1) / (0 - cp1)
        return val_cp1 + slope_to_meno * (t - cp1)
    elif t == 0:
        return val_after_jump
    elif t <= cp2:
        if cp2 == 0:
            return val_after_jump
        slope_from_meno = (val_cp2 - val_after_jump) / (cp2 - 0)
        return val_after_jump + slope_from_meno * t
    else:
        return val_cp2 + slope2 * (t - cp2)


def biological_kernel_three_cp(
    t, slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, val_cp3, slope2, cp1, cp2, cp3
):
    """Biological kernel function for 3 change points - sorts CPs before use"""
    # Sort change points and track their original indices
    cp_positions = [cp1, cp2, cp3]
    cp_values = [val_cp1, val_cp2, val_cp3]
    
    # Sort CPs and values together
    sorted_indices = np.argsort(cp_positions)
    sorted_cps = [cp_positions[i] for i in sorted_indices]
    sorted_vals = [cp_values[i] for i in sorted_indices]
    
    cp_min, cp_mid, cp_max = sorted_cps
    val_min, val_mid, val_max = sorted_vals
    
    # Handle segments based on sorted CP positions
    if t < cp_min:
        # Before first CP: use slope1
        return val_min + slope1 * (t - cp_min)
    elif t < cp_mid:
        # Between first and middle CP: linear interpolation
        if cp_mid == cp_min:
            return val_min  # CPs collapsed
        slope_seg = (val_mid - val_min) / (cp_mid - cp_min)
        return val_min + slope_seg * (t - cp_min)
    elif t < 0:
        # Between middle CP and menopause: linear to val_before_jump
        if cp_mid == 0:
            return val_before_jump
        slope_to_meno = (val_before_jump - val_mid) / (0 - cp_mid)
        return val_mid + slope_to_meno * (t - cp_mid)
    elif t == 0:
        # At menopause: jump occurs
        return val_after_jump
    elif t <= cp_max:
        # Between menopause and last CP: linear from val_after_jump
        if cp_max == 0:
            return val_after_jump
        slope_from_meno = (val_max - val_after_jump) / (cp_max - 0)
        return val_after_jump + slope_from_meno * t
    else:
        # After last CP: use slope2
        return val_max + slope2 * (t - cp_max)


def pred_age(
    age,
    survival_data,
    slope1,
    val_cp1,
    val_before_jump,
    val_after_jump,
    val_cp2,
    slope2,
    cp1,
    cp2,
):
    """Calculate predicted value using two-integral approach with cached fine integration"""
    grid_cache = get_integration_grid(survival_data)
    age_grid = grid_cache["age_grid"]
    event_density = grid_cache["event_density"]
    dt = grid_cache["dt"]

    # Safety check: ensure arrays have the same length
    if len(age_grid) != len(event_density):
        # Removed print to eliminate I/O blocking
        return 0.0

    integrand = np.zeros_like(age_grid)
    for i, meno_age in enumerate(age_grid):
        t = age - meno_age
        kernel_val = biological_kernel(
            t,
            slope1,
            val_cp1,
            val_before_jump,
            val_after_jump,
            val_cp2,
            slope2,
            cp1,
            cp2,
        )
        integrand[i] = kernel_val * event_density[i]

    prediction = trapezoid(integrand, dx=dt)
    return prediction


def pred_age_three_cp(
    age,
    survival_data,
    slope1,
    val_cp1,
    val_before_jump,
    val_after_jump,
    val_cp2,
    val_cp3,
    slope2,
    cp1,
    cp2,
    cp3,
):
    """Calculate predicted value using two-integral approach with 3 change points"""
    grid_cache = get_integration_grid(survival_data)
    age_grid = grid_cache["age_grid"]
    event_density = grid_cache["event_density"]
    dt = grid_cache["dt"]

    # Safety check: ensure arrays have the same length
    if len(age_grid) != len(event_density):
        return 0.0

    integrand = np.zeros_like(age_grid)
    for i, meno_age in enumerate(age_grid):
        t = age - meno_age
        kernel_val = biological_kernel_three_cp(
            t,
            slope1,
            val_cp1,
            val_before_jump,
            val_after_jump,
            val_cp2,
            val_cp3,
            slope2,
            cp1,
            cp2,
            cp3,
        )
        integrand[i] = kernel_val * event_density[i]

    prediction = trapezoid(integrand, dx=dt)
    return prediction


def biological_kernel_spline(
    t, 
    pre_spline_coeffs,  # Coefficients for pre-menopause spline
    post_spline_coeffs,  # Coefficients for post-menopause spline
    jump_magnitude,  # Jump at menopause
    pre_knots,  # Knot positions for pre-menopause spline
    post_knots,  # Knot positions for post-menopause spline
):
    """Biological kernel function using splines - separate splines before and after menopause"""
    # Create spline evaluators
    # Pre-menopause: t from -20 to 0
    # Post-menopause: t from 0 to 20
    
    if t < 0:
        # Before menopause: use pre-menopause spline
        if len(pre_spline_coeffs) < 2:
            # Fallback to linear if not enough coefficients
            return pre_spline_coeffs[0] + pre_spline_coeffs[1] * t if len(pre_spline_coeffs) >= 2 else pre_spline_coeffs[0]
        
        # Use B-spline evaluation
        # Create B-spline basis
        # For simplicity, use polynomial basis (can be upgraded to B-splines)
        # Using polynomial: y = sum(coeff_i * t^i)
        value = 0.0
        for i, coeff in enumerate(pre_spline_coeffs):
            value += coeff * (t ** i)
        return value
    
    elif t == 0:
        # At menopause: return value from pre-spline + jump
        pre_value = sum(coeff * (0.0 ** i) for i, coeff in enumerate(pre_spline_coeffs))
        return pre_value + jump_magnitude
    
    else:
        # After menopause: use post-menopause spline
        if len(post_spline_coeffs) < 2:
            # Fallback to linear
            return post_spline_coeffs[0] + post_spline_coeffs[1] * t if len(post_spline_coeffs) >= 2 else post_spline_coeffs[0]
        
        # Use polynomial basis for post-menopause
        # Adjust t to be relative to 0 (menopause)
        value = 0.0
        for i, coeff in enumerate(post_spline_coeffs):
            value += coeff * (t ** i)
        return value


def pred_age_spline(
    age,
    survival_data,
    pre_spline_coeffs,
    post_spline_coeffs,
    jump_magnitude,
    pre_knots=None,
    post_knots=None,
):
    """Calculate predicted value using two-integral approach with spline kernel"""
    grid_cache = get_integration_grid(survival_data)
    age_grid = grid_cache["age_grid"]
    event_density = grid_cache["event_density"]
    dt = grid_cache["dt"]

    if len(age_grid) != len(event_density):
        return 0.0

    integrand = np.zeros_like(age_grid)
    for i, meno_age in enumerate(age_grid):
        t = age - meno_age
        kernel_val = biological_kernel_spline(
            t,
            pre_spline_coeffs,
            post_spline_coeffs,
            jump_magnitude,
            pre_knots,
            post_knots,
        )
        integrand[i] = kernel_val * event_density[i]

    prediction = trapezoid(integrand, dx=dt)
    return prediction


def fit_spline_optimization(
    test_ages,
    test_values,
    survival_data,
    linear_model,
    test_sds=None,
    test_ns=None,
    use_iterative_reweighting=False,
    n_coeffs_pre=4,  # Number of coefficients for pre-menopause spline
    n_coeffs_post=4,  # Number of coefficients for post-menopause spline
):
    """SPLINE OPTIMIZATION - separate splines before and after menopause"""
    print(f"   🎯 Starting SPLINE OPTIMIZATION (pre: {n_coeffs_pre} coeffs, post: {n_coeffs_post} coeffs)")
    weights_se = test_sds

    # Estimate initial values from linear model
    data_mean = np.mean(test_values)
    data_std = np.std(test_values)
    
    # Starting strategies
    starting_strategies = [
        {
            "name": "linear_based",
            "pre_coeffs": [linear_model["c"], linear_model["d"], 0.0, 0.0][:n_coeffs_pre],  # [intercept, slope, quad, cubic]
            "post_coeffs": [linear_model["a"], linear_model["b"], 0.0, 0.0][:n_coeffs_post],
            "jump": linear_model["a"] - linear_model["c"],
        },
        {
            "name": "data_based",
            "pre_coeffs": [data_mean, 0.01, 0.0, 0.0][:n_coeffs_pre],
            "post_coeffs": [data_mean, 0.01, 0.0, 0.0][:n_coeffs_post],
            "jump": data_std,
        },
    ]

    best_params = None
    best_loss = np.inf
    best_pure_rss = None

    for strategy in starting_strategies:
        print(f"   🔄 Trying spline strategy: {strategy['name']}")

        # Build initial guess: [pre_coeffs..., post_coeffs..., jump]
        initial_guess = list(strategy["pre_coeffs"]) + list(strategy["post_coeffs"]) + [strategy["jump"]]

        def objective(params):
            # Split parameters
            pre_coeffs = params[:n_coeffs_pre]
            post_coeffs = params[n_coeffs_pre:n_coeffs_pre + n_coeffs_post]
            jump = params[-1]

            # Check for invalid parameters
            if not np.isfinite(params).all():
                return np.inf

            try:
                predictions = []
                for age in test_ages:
                    pred = pred_age_spline(
                        age,
                        survival_data,
                        pre_coeffs,
                        post_coeffs,
                        jump,
                    )
                    if not np.isfinite(pred):
                        return np.inf
                    predictions.append(pred)

                predictions = np.array(predictions)

                if len(predictions) != len(test_values):
                    return np.inf

                residuals = test_values - predictions

                if not np.isfinite(residuals).all():
                    return np.inf

                # Calculate RSS
                if use_iterative_reweighting and "fit_linear_model_irls" in globals():
                    try:
                        X = np.ones((len(test_ages), 1))
                        y = residuals
                        initial_weights = None
                        if (
                            DEFAULT_CONFIG["USE_CV"]
                            and weights_se is not None
                            and np.any(weights_se > 0)
                        ):
                            initial_weights = 1.0 / (weights_se**2)
                        weights, _ = fit_linear_model_irls(X, y, initial_weights)
                        rss = np.sum(weights * residuals**2)
                    except Exception:
                        rss = np.sum(residuals**2)
                else:
                    if DEFAULT_CONFIG["USE_CV"] and weights_se is not None and np.any(weights_se > 0):
                        weights = 1.0 / (weights_se**2)
                        rss = np.sum(weights * residuals**2)
                    else:
                        rss = np.sum(residuals**2)

                # Store pure RSS for BIC calculation
                objective.pure_rss = rss
                return rss

            except Exception as e:
                return np.inf

        # Bounds - allow reasonable ranges for coefficients and jump
        test_range = test_values.max() - test_values.min()
        bounds = []
        
        # Pre-menopause coefficients bounds
        for i in range(n_coeffs_pre):
            if i == 0:  # Intercept
                bounds.append((-test_range, 2 * test_range))
            elif i == 1:  # Linear term
                bounds.append((-2 * test_range, 2 * test_range))
            else:  # Higher order terms
                bounds.append((-test_range / 10, test_range / 10))
        
        # Post-menopause coefficients bounds
        for i in range(n_coeffs_post):
            if i == 0:  # Intercept
                bounds.append((-test_range, 2 * test_range))
            elif i == 1:  # Linear term
                bounds.append((-2 * test_range, 2 * test_range))
            else:  # Higher order terms
                bounds.append((-test_range / 10, test_range / 10))
        
        # Jump magnitude
        bounds.append((-2 * test_range, 2 * test_range))

        try:
            print(f"   🔄 Starting spline optimization with strategy: {strategy['name']}")
            result = minimize(
                objective,
                initial_guess,
                method="L-BFGS-B",
                bounds=bounds,
            )

            if result.success and result.fun < best_loss:
                print(f"   ✅ Spline optimization successful with strategy: {strategy['name']}")
                best_params = result.x
                best_loss = result.fun
                # Recalculate pure RSS
                pre_coeffs = result.x[:n_coeffs_pre]
                post_coeffs = result.x[n_coeffs_pre:n_coeffs_pre + n_coeffs_post]
                jump = result.x[-1]
                
                print(f"   📊 Spline coefficients - Pre: {pre_coeffs[:3]}, Post: {post_coeffs[:3]}, Jump: {jump:.3f}")
                
                # Recalculate predictions for pure RSS
                predictions = []
                for age in test_ages:
                    pred = pred_age_spline(
                        age,
                        survival_data,
                        pre_coeffs,
                        post_coeffs,
                        jump,
                    )
                    predictions.append(pred)
                predictions = np.array(predictions)
                residuals = test_values - predictions
                if DEFAULT_CONFIG["USE_CV"] and weights_se is not None and np.any(weights_se > 0):
                    weights = 1.0 / (weights_se**2)
                    best_pure_rss = np.sum(weights * residuals**2)
                else:
                    best_pure_rss = np.sum(residuals**2)
                print(f"   ✅ RECALC SPLINE: RSS={best_pure_rss:.6f}")
            else:
                print(f"   ❌ Spline optimization failed with strategy: {strategy['name']}")
        except Exception as e:
            print(f"   ❌ Exception in spline optimization: {str(e)}")

    if best_params is None:
        return None, np.inf, None

    return best_params, best_loss, best_pure_rss


def fit_single_cp_optimization(
    cp_position,
    test_ages,
    test_values,
    survival_data,
    linear_model,
    test_sds=None,
    test_ns=None,
    use_iterative_reweighting=False,
):
    """Single CP optimization - cp_position will be optimized as a variable"""
    # SIMPLE FIX: Use test_sds directly (same as linear model) - don't transform to SE!
    # This ensures consistent weighting between linear and change point models
    weights_se = test_sds

    # Reduced debug output to prevent I/O contention
    # Note: cp_position input parameter is ignored - it will be optimized

    # Stage 1: Estimate slopes from data outside menopause transition
    # print(f"      Stage 1: Estimating slopes outside menopause transition")

    # Use slopes from the linear model (disable data-based slope estimation)
    # This matches the approach in fit_staged_optimization to ensure consistency
    pre_slope = linear_model["d"]
    post_slope = linear_model["b"]

    # Stage 2: Estimate jump and change point values
    # print(f"      Stage 2: Estimating jump and change point values")

    data_mean = np.mean(test_values)
    data_std = np.std(test_values)
    data_min = np.min(test_values)
    data_max = np.max(test_values)
    # Two different starting strategies
    # Note: cp_position input is ignored - will be optimized
    starting_strategies = [
        {
            "name": "linear_based",
            "slope1": pre_slope,
            "slope2": post_slope,
            "val_before_jump": linear_model["c"],
            "val_after_jump": linear_model["a"],
            "val_cp": linear_model["c"],
        },
        {
            "name": "data_based",
            "slope1": pre_slope,
            "slope2": post_slope,
            "val_before_jump": data_mean - data_std,
            "val_after_jump": data_mean + data_std,
            "val_cp": data_min,
        },
    ]

    best_params = None
    best_loss = np.inf
    best_pure_rss = None  # Initialize to avoid undefined variable

    for strategy in starting_strategies:
        # Removed prints to eliminate I/O blocking

        initial_guess = [
            strategy["slope1"],
            strategy["val_cp"],
            strategy["val_before_jump"],
            strategy["val_after_jump"],
            strategy["slope2"],
            0.0,  # cp_position initial guess
        ]

        def objective(params):
            slope1, val_cp, val_before_jump, val_after_jump, slope2, opt_cp_position = params

            # Calculate predictions first
            predictions = []
            for age in test_ages:
                pred = pred_age_single_cp(
                    age,
                    survival_data,
                    slope1,
                    val_cp,
                    val_before_jump,
                    val_after_jump,
                    slope2,
                    opt_cp_position,
                )
                predictions.append(pred)

            predictions = np.array(predictions)

            # Debug: Check array shapes
            if len(predictions) != len(test_values):
                # Removed print to eliminate I/O blocking
                return np.inf

            # Calculate residuals
            residuals = test_values - predictions

            # Calculate weighted or unweighted RSS
            if use_iterative_reweighting and "fit_linear_model_irls" in globals():
                # Use IRLS for robust fitting
                try:
                    # Create design matrix for IRLS using residuals instead of raw values
                    X = np.ones((len(test_ages), 1))
                    y = residuals  # Use residuals instead of test_values

                    # Initialize weights if we have standard deviations
                    initial_weights = None
                    if (
                        DEFAULT_CONFIG["USE_CV"]
                        and weights_se is not None
                        and np.any(weights_se > 0)
                    ):
                        initial_weights = 1.0 / (weights_se**2)
                        initial_weights = initial_weights / np.mean(
                            initial_weights
                        )  # Normalize

                    # Fit using IRLS to get weights based on residuals
                    _, _, _, weights = fit_linear_model_irls(
                        X,
                        y,
                        initial_weights=initial_weights,
                        max_iter=10,
                        weight_func="huber",
                    )

                    # Calculate weighted RSS with IRLS weights
                    rss = np.sum(weights * residuals**2)
                except Exception as e:
                    # Fall back to standard weighting
                    if (
                        DEFAULT_CONFIG["USE_CV"]
                        and weights_se is not None
                        and np.any(weights_se > 0)
                    ):
                        weights = 1.0 / (weights_se**2)
                        rss = np.sum(weights * residuals**2)
                    else:
                        rss = np.sum(residuals**2)
            elif (
                DEFAULT_CONFIG["USE_CV"]
                and weights_se is not None
                and np.any(weights_se > 0)
            ):
                if len(weights_se) != len(residuals):
                    # Removed print to eliminate I/O blocking
                    return np.inf
                weights = 1.0 / (weights_se**2)
                rss = np.sum(weights * residuals**2)
            else:
                rss = np.sum(residuals**2)

            # Add slope penalty for single change point model - FIXED VERSION
            total_penalty = 0.0
            if DEFAULT_CONFIG["ENABLE_SLOPE_PENALIZATION"]:
                jump_magnitude = val_after_jump - val_before_jump

                # Calculate slope from change point to menopause - FIXED
                if opt_cp_position != 0:
                    slope_to_meno = (val_before_jump - val_cp) / (0 - opt_cp_position)
                    penalty = calculate_slope_penalty(
                        slope_to_meno, jump_magnitude, opt_cp_position
                    )
                    total_penalty += penalty

            # Store both penalized and unpenalized loss for later BIC calculation
            objective.pure_rss = rss
            return rss + total_penalty

        # Stage 3: Optimize with bounds
        test_range = test_values.max() - test_values.min()
        bounds = [
            (-2 * test_range, 2 * test_range),  # slope1
            (-2 * test_range, 3 * test_range),  # val_cp
            (-test_range, 2 * test_range),  # val_before_jump
            (-test_range, 3 * test_range),  # val_after_jump
            (-2 * test_range, 2 * test_range),  # slope2
            (-20.0, 20.0),  # cp_position
        ]

        try:
            print(
                f"   🔄 Starting one CP optimization with the strategy: {strategy['name']}"
            )

            # Add timeout and more conservative settings
            result = minimize(
                objective,
                initial_guess,
                method="L-BFGS-B",
                bounds=bounds,
                # options={"maxiter": 250, "ftol": 1e-6, "gtol": 1e-6},
            )

            if result.success and result.fun < best_loss:
                print(
                    f"   ✅ One CP optimization successful with the strategy: {strategy['name']}"
                )
                best_params = result.x
                best_loss = result.fun
                # IMPORTANT FIX: Recalculate pure RSS using the BEST parameters from result.x
                # objective.pure_rss is from the LAST evaluation, not necessarily the best!
                slope1, val_cp, val_before_jump, val_after_jump, slope2, opt_cp_position = result.x
                print(f"   🔍 RECALC SINGLE CP: Using best params from result.x")
                print(f"      slope1={slope1:.6f}, slope2={slope2:.6f}, val_cp={val_cp:.3f}, cp_position={opt_cp_position:.3f}")
                _, best_pure_rss = get_final_predictions_for_single_change_point(
                    test_ages,
                    survival_data,
                    slope1,
                    val_cp,
                    val_before_jump,
                    val_after_jump,
                    slope2,
                    opt_cp_position,
                    test_values=test_values,
                    test_sds=test_sds,  # Use same weights as optimization
                    test_ns=None,  # Don't transform - use test_sds directly
                    return_rss=True,
                )
                print(f"   ✅ RECALC SINGLE CP: RSS={best_pure_rss:.6f}")
                # Removed prints to eliminate I/O blocking
                pass
            else:
                print(
                    f"   ❌ One CP optimization failed with the strategy: {strategy['name']}"
                )
                print(f"   ❌ Reason: {result.message}")
                # Removed prints to eliminate I/O blocking
                pass
        except Exception as e:

            print(f"   ❌ Exception in two CP optimization: {e}")
            # Removed prints to eliminate I/O blocking
            pass

    if best_params is None:
        # Removed prints to eliminate I/O blocking
        return None, np.inf, None

    # Calculate final predictions and metrics
    slope1, val_cp, val_before_jump, val_after_jump, slope2, opt_cp_position = best_params

    final_predictions = []
    for age in test_ages:
        pred = pred_age_single_cp(
            age,
            survival_data,
            slope1,
            val_cp,
            val_before_jump,
            val_after_jump,
            slope2,
            opt_cp_position,
        )
        final_predictions.append(pred)

    final_predictions = np.array(final_predictions)
    residuals = test_values - final_predictions

    # Calculate R-squared (weighted if CV is enabled)
    # NOTE: This is just for logging, the actual pure_rss was already calculated above
    if DEFAULT_CONFIG["USE_CV"] and weights_se is not None and np.any(weights_se > 0):
        weights = 1.0 / (weights_se**2)
        rss = np.sum(weights * residuals**2)
        tss = np.sum(
            weights * (test_values - np.average(test_values, weights=weights)) ** 2
        )
    else:
        rss = np.sum(residuals**2)
        tss = np.sum((test_values - np.mean(test_values)) ** 2)

    r_squared = 1 - rss / tss if tss > 0 else 0

    # Removed prints to eliminate I/O blocking

    return best_params, best_loss, best_pure_rss


def pred_age_single_cp(
    age,
    survival_data,
    slope1,
    val_cp,
    val_before_jump,
    val_after_jump,
    slope2,
    cp_position,
):
    """Calculate predicted value using two-integral approach with cached fine integration - EXACT COPY OF TWO CP"""
    grid_cache = get_integration_grid(survival_data)
    age_grid = grid_cache["age_grid"]
    event_density = grid_cache["event_density"]
    dt = grid_cache["dt"]

    integrand = np.zeros_like(age_grid)
    for i, meno_age in enumerate(age_grid):
        t = age - meno_age
        kernel_val = biological_kernel_single_cp(
            t,
            slope1,
            val_cp,
            val_before_jump,
            val_after_jump,
            slope2,
            cp_position,
        )
        integrand[i] = kernel_val * event_density[i]

    prediction = trapezoid(integrand, dx=dt)
    return prediction


def pred_age_linear(age, survival_data, a, b, c, d):
    """Calculate predicted value using two-integral approach for linear model"""
    grid_cache = get_integration_grid(survival_data)
    age_grid = grid_cache["age_grid"]
    event_density = grid_cache["event_density"]
    dt = grid_cache["dt"]

    # Safety check: ensure arrays have the same length
    if len(age_grid) != len(event_density):
        # Removed print to eliminate I/O blocking
        return 0.0

    integrand = np.zeros_like(age_grid)
    for i, meno_age in enumerate(age_grid):
        t = age - meno_age
        kernel_val = biological_kernel_linear(t, a, b, c, d)
        integrand[i] = kernel_val * event_density[i]

    prediction = trapezoid(integrand, dx=dt)
    return prediction


def biological_kernel_linear(t, a, b, c, d):
    """Linear menopause model biological kernel - same as biological_kernel but with explicit parameters"""
    if t < 0:
        return c + d * t
    else:
        return a + b * t


def biological_kernel_polynomial3(t, a0, a1, a2, a3):
    """
    3rd order polynomial biological kernel

    Parameters:
    - a0, a1, a2, a3: polynomial coefficients

    Returns: a0 + a1*t + a2*t^2 + a3*t^3
    """
    return a0 + a1*t + a2*t**2 + a3*t**3


def biological_kernel_polynomial4(t, a0, a1, a2, a3, a4):
    """
    4th order polynomial biological kernel

    Parameters:
    - a0, a1, a2, a3, a4: polynomial coefficients

    Returns: a0 + a1*t + a2*t^2 + a3*t^3 + a4*t^4
    """
    return a0 + a1*t + a2*t**2 + a3*t**3 + a4*t**4


def biological_kernel_piecewise_linear_2dots(t, v_m8_33, v_8_33, slope_before, slope_after):
    """
    Piecewise linear with slope changes ONLY at -8.333 and 8.333 (4 parameters)
    
    The slope is CONSTANT between -8.333 and 8.333 (no change at t=0).

    Parameters:
    - v_m8_33: value at t=-8.333
    - v_8_33: value at t=8.333
    - slope_before: slope for t < -8.333
    - slope_after: slope for t > 8.333

    The function is continuous at -8.333 and 8.333.
    Between -8.333 and 8.333, the slope is constant (calculated from the two values).
    """
    if t < -8.333:
        # Before first transition point: use slope_before
        return v_m8_33 + slope_before * (t - (-8.333))
    elif t <= 8.333:
        # Between -8.333 and 8.333: constant slope (linear interpolation)
        # This ensures NO slope change at t=0
        slope_middle = (v_8_33 - v_m8_33) / (8.333 - (-8.333))
        return v_m8_33 + slope_middle * (t - (-8.333))
    else:
        # After second transition point: use slope_after
        return v_8_33 + slope_after * (t - 8.333)


def biological_kernel_piecewise_linear_continuous(t, v_m12_5, v_0, v_12_5, slope_before, slope_after):
    """
    Piecewise linear with continuous transitions at -12.5, 0, 12.5

    Parameters:
    - v_m12_5: value at t=-12.5
    - v_0: value at t=0
    - v_12_5: value at t=12.5
    - slope_before: slope for t < -12.5
    - slope_after: slope for t > 12.5

    The function is continuous at all transition points.
    """
    if t < -12.5:
        return v_m12_5 + slope_before * (t - (-12.5))
    elif t < 0:
        # Linear interpolation between -12.5 and 0
        alpha = (t - (-12.5)) / (0 - (-12.5))
        return (1 - alpha) * v_m12_5 + alpha * v_0
    elif t < 12.5:
        # Linear interpolation between 0 and 12.5
        alpha = (t - 0) / (12.5 - 0)
        return (1 - alpha) * v_0 + alpha * v_12_5
    else:
        return v_12_5 + slope_after * (t - 12.5)


def biological_kernel_linear_pre_only(t, slope_pre, val_before_jump, val_after_jump):
    """
    Linear model with jump at menopause, slope only before menopause (3 parameters)

    Parameters:
    - slope_pre: slope before menopause
    - val_before_jump: value just before menopause (t=0-)
    - val_after_jump: value just after menopause (t=0+), flat after

    Pre-menopause (t<0): linear with slope
    Post-menopause (t>=0): flat at val_after_jump
    """
    if t < 0:
        return val_before_jump + slope_pre * t
    else:
        return val_after_jump


def biological_kernel_linear_post_only(t, val_before_jump, val_after_jump, slope_post):
    """
    Linear model with jump at menopause, slope only after menopause (3 parameters)

    Parameters:
    - val_before_jump: value just before menopause (t=0-), flat before
    - val_after_jump: value just after menopause (t=0+)
    - slope_post: slope after menopause

    Pre-menopause (t<0): flat at val_before_jump
    Post-menopause (t>=0): linear with slope
    """
    if t < 0:
        return val_before_jump
    else:
        return val_after_jump + slope_post * t


def biological_kernel_exp_linear(t, exp_a, exp_b, val_before_jump, val_after_jump, slope_post):
    """
    Exponential-Linear: exponential before menopause, linear after, with jump at t=0 (5 parameters)

    Parameters:
    - exp_a: exponential coefficient (before menopause)
    - exp_b: exponential rate (before menopause)
    - val_before_jump: value just before menopause (t=0-)
    - val_after_jump: value just after menopause (t=0+)
    - slope_post: linear slope after menopause

    Pre-menopause (t<0): exp_a * exp(exp_b * t)
    Post-menopause (t>=0): val_after_jump + slope_post * t
    Jump at t=0: val_after_jump - val_before_jump
    """
    if t < 0:
        return exp_a * np.exp(exp_b * t)
    else:
        return val_after_jump + slope_post * t


def biological_kernel_linear_exp(t, slope_pre, val_before_jump, val_after_jump, exp_a, exp_b):
    """
    Linear-Exponential: linear before menopause, exponential after, with jump at t=0 (5 parameters)

    Parameters:
    - slope_pre: linear slope before menopause
    - val_before_jump: value just before menopause (t=0-)
    - val_after_jump: value just after menopause (t=0+)
    - exp_a: exponential coefficient (after menopause)
    - exp_b: exponential rate (after menopause)

    Pre-menopause (t<0): val_before_jump + slope_pre * t
    Post-menopause (t>=0): exp_a * exp(exp_b * t)
    Jump at t=0: val_after_jump - val_before_jump
    """
    if t < 0:
        return val_before_jump + slope_pre * t
    else:
        return exp_a * np.exp(exp_b * t)


def biological_kernel_sigmoid(t, c, d, w, h, b):
    """
    Sigmoid transition menopause model biological kernel
    
    Parameters:
    - c: baseline value (pre-menopause intercept at t=0)
    - d: slope before menopause
    - w: width of sigmoid transition (controls how gradual the transition is)
    - h: height of sigmoid (jump magnitude from pre to post menopause)
    - b: slope after menopause
    
    The sigmoid smoothly transitions from pre-menopause line (c + d*t) 
    to post-menopause line ((c + h) + b*t) around t=0.
    """
    # Sigmoid function: transitions from 0 to 1 around t=0
    sigmoid = 1.0 / (1.0 + np.exp(-t / w))
    
    # Pre-menopause value: c + d*t
    pre_value = c + d * t
    
    # Post-menopause value: (c + h) + b*t
    post_value = (c + h) + b * t
    
    # Smooth transition using sigmoid
    return (1.0 - sigmoid) * pre_value + sigmoid * post_value


def pred_age_polynomial3(age, survival_data, a0, a1, a2, a3):
    """Calculate predicted value using 3rd order polynomial with two-integral approach"""
    grid_cache = get_integration_grid(survival_data)
    age_grid = grid_cache["age_grid"]
    event_density = grid_cache["event_density"]
    dt = grid_cache["dt"]

    if len(age_grid) != len(event_density):
        return 0.0

    integrand = np.zeros_like(age_grid)
    for i, meno_age in enumerate(age_grid):
        t = age - meno_age
        kernel_val = biological_kernel_polynomial3(t, a0, a1, a2, a3)
        integrand[i] = kernel_val * event_density[i]

    prediction = trapezoid(integrand, dx=dt)
    return prediction


def pred_age_polynomial4(age, survival_data, a0, a1, a2, a3, a4):
    """Calculate predicted value using 4th order polynomial with two-integral approach"""
    grid_cache = get_integration_grid(survival_data)
    age_grid = grid_cache["age_grid"]
    event_density = grid_cache["event_density"]
    dt = grid_cache["dt"]

    if len(age_grid) != len(event_density):
        return 0.0

    integrand = np.zeros_like(age_grid)
    for i, meno_age in enumerate(age_grid):
        t = age - meno_age
        kernel_val = biological_kernel_polynomial4(t, a0, a1, a2, a3, a4)
        integrand[i] = kernel_val * event_density[i]

    prediction = trapezoid(integrand, dx=dt)
    return prediction


def pred_age_piecewise_linear_2dots(age, survival_data, v_m8_33, v_8_33, slope_before, slope_after):
    """Calculate predicted value using piecewise linear 2-dots model with two-integral approach"""
    grid_cache = get_integration_grid(survival_data)
    age_grid = grid_cache["age_grid"]
    event_density = grid_cache["event_density"]
    dt = grid_cache["dt"]

    if len(age_grid) != len(event_density):
        return 0.0

    integrand = np.zeros_like(age_grid)
    for i, meno_age in enumerate(age_grid):
        t = age - meno_age
        kernel_val = biological_kernel_piecewise_linear_2dots(t, v_m8_33, v_8_33, slope_before, slope_after)
        integrand[i] = kernel_val * event_density[i]

    prediction = trapezoid(integrand, dx=dt)
    return prediction


def pred_age_piecewise_linear_continuous(age, survival_data, v_m12_5, v_0, v_12_5, slope_before, slope_after):
    """Calculate predicted value using piecewise linear continuous model with two-integral approach"""
    grid_cache = get_integration_grid(survival_data)
    age_grid = grid_cache["age_grid"]
    event_density = grid_cache["event_density"]
    dt = grid_cache["dt"]

    if len(age_grid) != len(event_density):
        return 0.0

    integrand = np.zeros_like(age_grid)
    for i, meno_age in enumerate(age_grid):
        t = age - meno_age
        kernel_val = biological_kernel_piecewise_linear_continuous(t, v_m12_5, v_0, v_12_5, slope_before, slope_after)
        integrand[i] = kernel_val * event_density[i]

    prediction = trapezoid(integrand, dx=dt)
    return prediction


def pred_age_linear_pre_only(age, survival_data, slope_pre, val_before_jump, val_after_jump):
    """Calculate predicted value using linear pre-only model with two-integral approach"""
    grid_cache = get_integration_grid(survival_data)
    age_grid = grid_cache["age_grid"]
    event_density = grid_cache["event_density"]
    dt = grid_cache["dt"]

    if len(age_grid) != len(event_density):
        return 0.0

    integrand = np.zeros_like(age_grid)
    for i, meno_age in enumerate(age_grid):
        t = age - meno_age
        kernel_val = biological_kernel_linear_pre_only(t, slope_pre, val_before_jump, val_after_jump)
        integrand[i] = kernel_val * event_density[i]

    prediction = trapezoid(integrand, dx=dt)
    return prediction


def pred_age_linear_post_only(age, survival_data, val_before_jump, val_after_jump, slope_post):
    """Calculate predicted value using linear post-only model with two-integral approach"""
    grid_cache = get_integration_grid(survival_data)
    age_grid = grid_cache["age_grid"]
    event_density = grid_cache["event_density"]
    dt = grid_cache["dt"]

    if len(age_grid) != len(event_density):
        return 0.0

    integrand = np.zeros_like(age_grid)
    for i, meno_age in enumerate(age_grid):
        t = age - meno_age
        kernel_val = biological_kernel_linear_post_only(t, val_before_jump, val_after_jump, slope_post)
        integrand[i] = kernel_val * event_density[i]

    prediction = trapezoid(integrand, dx=dt)
    return prediction


def pred_age_exp_linear(age, survival_data, exp_a, exp_b, val_before_jump, val_after_jump, slope_post):
    """Calculate predicted value using exp-linear model with two-integral approach"""
    grid_cache = get_integration_grid(survival_data)
    age_grid = grid_cache["age_grid"]
    event_density = grid_cache["event_density"]
    dt = grid_cache["dt"]

    if len(age_grid) != len(event_density):
        return 0.0

    integrand = np.zeros_like(age_grid)
    for i, meno_age in enumerate(age_grid):
        t = age - meno_age
        kernel_val = biological_kernel_exp_linear(t, exp_a, exp_b, val_before_jump, val_after_jump, slope_post)
        integrand[i] = kernel_val * event_density[i]

    prediction = trapezoid(integrand, dx=dt)
    return prediction


def pred_age_linear_exp(age, survival_data, slope_pre, val_before_jump, val_after_jump, exp_a, exp_b):
    """Calculate predicted value using linear-exp model with two-integral approach"""
    grid_cache = get_integration_grid(survival_data)
    age_grid = grid_cache["age_grid"]
    event_density = grid_cache["event_density"]
    dt = grid_cache["dt"]

    if len(age_grid) != len(event_density):
        return 0.0

    integrand = np.zeros_like(age_grid)
    for i, meno_age in enumerate(age_grid):
        t = age - meno_age
        kernel_val = biological_kernel_linear_exp(t, slope_pre, val_before_jump, val_after_jump, exp_a, exp_b)
        integrand[i] = kernel_val * event_density[i]

    prediction = trapezoid(integrand, dx=dt)
    return prediction


def pred_age_sigmoid(age, survival_data, c, d, w, h, b):
    """Calculate predicted value using sigmoid transition model with two-integral approach"""
    grid_cache = get_integration_grid(survival_data)
    age_grid = grid_cache["age_grid"]
    event_density = grid_cache["event_density"]
    dt = grid_cache["dt"]

    # Safety check: ensure arrays have the same length
    if len(age_grid) != len(event_density):
        return 0.0

    integrand = np.zeros_like(age_grid)
    for i, meno_age in enumerate(age_grid):
        t = age - meno_age
        kernel_val = biological_kernel_sigmoid(t, c, d, w, h, b)
        integrand[i] = kernel_val * event_density[i]

    prediction = trapezoid(integrand, dx=dt)
    return prediction


def pred_age_biological(age, survival_data):
    """Calculate biological age using survival function"""
    # This is a simplified version - you might want to use your existing biological age calculation
    return age  # For now, just return chronological age


def fit_sigmoid_menopause_model(
    test_ages,
    test_values,
    survival_data,
    linear_model,
    test_sds=None,
    test_ns=None,
    use_iterative_reweighting=False,
):
    """
    Fit the sigmoid transition menopause model
    
    Parameters:
    - test_ages: ages of data points
    - test_values: actual values at those ages
    - survival_data: menopause survival curve
    - linear_model: fitted linear model (used for initialization)
    - test_sds: standard deviations for weighting
    - test_ns: sample sizes for weighting
    
    Returns:
    - sigmoid_model: dict with parameters c, d, w, h, b and fitted values
    """
    print("🔬 FITTING SIGMOID MENOPAUSE MODEL")
    
    # SIMPLE FIX: Use test_sds directly (same as linear model and change point models)
    # This ensures consistent weighting across ALL models
    test_sds_for_weights = test_sds
    
    # Initialize parameters from linear model
    # c: pre-menopause intercept, d: pre-menopause slope
    # h: jump magnitude, b: post-menopause slope
    # w: sigmoid width (start with reasonable default)
    c_init = linear_model["c"]
    d_init = linear_model["d"]
    h_init = linear_model["a"] - linear_model["c"]  # jump magnitude
    b_init = linear_model["b"]
    w_init = 2.0  # Start with 2-year transition width
    
    initial_guess = [c_init, d_init, w_init, h_init, b_init]
    
    def objective(params):
        c, d, w, h, b = params
        
        # Ensure positive width
        if w <= 0:
            return np.inf
        
        # Calculate predictions
        predictions = []
        for age in test_ages:
            pred = pred_age_sigmoid(age, survival_data, c, d, w, h, b)
            if not np.isfinite(pred):
                return np.inf
            predictions.append(pred)
        
        predictions = np.array(predictions)
        residuals = test_values - predictions
        
        # Calculate weighted or unweighted RSS
        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
        else:
            rss = np.sum(residuals**2)
        
        return rss
    
    # Set reasonable bounds
    data_range = test_values.max() - test_values.min()
    bounds = [
        (-2 * data_range, 2 * data_range),  # c: intercept
        (-1.0, 1.0),  # d: pre-menopause slope
        (0.1, 10.0),  # w: sigmoid width (0.1 to 10 years)
        (-2 * data_range, 2 * data_range),  # h: jump magnitude
        (-1.0, 1.0),  # b: post-menopause slope
    ]
    
    # Optimize
    try:
        result = minimize(
            objective,
            initial_guess,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 200, "ftol": 1e-6, "gtol": 1e-6},
        )
        
        if not result.success:
            print(f"   ⚠️ Sigmoid optimization did not fully converge: {result.message}")
        
        c, d, w, h, b = result.x
        
        # Calculate final predictions and metrics
        predictions = []
        for age in test_ages:
            pred = pred_age_sigmoid(age, survival_data, c, d, w, h, b)
            predictions.append(pred)
        predictions = np.array(predictions)
        
        residuals = test_values - predictions
        
        # Calculate R-squared (weighted if CV is enabled)
        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
            tss = np.sum(
                weights * (test_values - np.average(test_values, weights=weights)) ** 2
            )
        else:
            rss = np.sum(residuals**2)
            tss = np.sum((test_values - np.mean(test_values)) ** 2)
        
        r_squared = 1 - rss / tss if tss > 0 else 0
        
        sigmoid_model = {
            "c": c,
            "d": d,
            "w": w,
            "h": h,
            "b": b,
            "r_squared": r_squared,
            "rss": rss,
            "fitted_values": predictions,
            "residuals": residuals,
            "n_points": len(test_values),
        }
        
        print(f"   ✅ Sigmoid menopause model fitted:")
        print(f"   c (pre-meno intercept): {c:.6f}")
        print(f"   d (pre-meno slope): {d:.6f}")
        print(f"   w (sigmoid width): {w:.6f}")
        print(f"   h (jump magnitude): {h:.6f}")
        print(f"   b (post-meno slope): {b:.6f}")
        print(f"   R² = {r_squared:.6f}")
        
        return sigmoid_model
        
    except Exception as e:
        print(f"   ❌ Error fitting sigmoid model: {e}")
        return None


def fit_polynomial3_model(
    test_ages,
    test_values,
    survival_data,
    linear_model,
    test_sds=None,
    test_ns=None,
    use_iterative_reweighting=False,
):
    """Fit 3rd order polynomial model"""
    print("🔬 FITTING 3RD ORDER POLYNOMIAL MODEL")

    test_sds_for_weights = test_sds

    # Initialize with simple values
    a0_init = np.mean(test_values)
    a1_init = 0.0
    a2_init = 0.0
    a3_init = 0.0

    initial_guess = [a0_init, a1_init, a2_init, a3_init]

    def objective(params):
        a0, a1, a2, a3 = params

        predictions = []
        for age in test_ages:
            pred = pred_age_polynomial3(age, survival_data, a0, a1, a2, a3)
            if not np.isfinite(pred):
                return np.inf
            predictions.append(pred)

        predictions = np.array(predictions)
        residuals = test_values - predictions

        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
        else:
            rss = np.sum(residuals**2)

        return rss

    # Set reasonable bounds
    data_range = test_values.max() - test_values.min()
    bounds = [
        (-2 * data_range, 2 * data_range),  # a0
        (-1.0, 1.0),  # a1
        (-0.1, 0.1),  # a2
        (-0.01, 0.01),  # a3
    ]

    try:
        result = minimize(
            objective,
            initial_guess,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 200, "ftol": 1e-6, "gtol": 1e-6},
        )

        if not result.success:
            print(f"   ⚠️ Polynomial3 optimization did not fully converge: {result.message}")

        a0, a1, a2, a3 = result.x

        predictions = []
        for age in test_ages:
            pred = pred_age_polynomial3(age, survival_data, a0, a1, a2, a3)
            predictions.append(pred)
        predictions = np.array(predictions)

        residuals = test_values - predictions

        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
            tss = np.sum(
                weights * (test_values - np.average(test_values, weights=weights)) ** 2
            )
        else:
            rss = np.sum(residuals**2)
            tss = np.sum((test_values - np.mean(test_values)) ** 2)

        r_squared = 1 - rss / tss if tss > 0 else 0

        poly3_model = {
            "a0": a0,
            "a1": a1,
            "a2": a2,
            "a3": a3,
            "r_squared": r_squared,
            "rss": rss,
            "fitted_values": predictions,
            "residuals": residuals,
            "n_points": len(test_values),
        }

        print(f"   ✅ Polynomial3 model fitted: R² = {r_squared:.6f}")

        return poly3_model

    except Exception as e:
        print(f"   ❌ Error fitting polynomial3 model: {e}")
        return None


def fit_piecewise_linear_2dots_model(
    test_ages,
    test_values,
    survival_data,
    linear_model,
    test_sds=None,
    test_ns=None,
    use_iterative_reweighting=False,
):
    """Fit piecewise linear 2-dots model with slope changes ONLY at -8.333 and 8.333"""
    print("🔬 FITTING PIECEWISE LINEAR 2-DOTS MODEL")

    test_sds_for_weights = test_sds

    # Initialize from linear model
    mean_val = np.mean(test_values)
    v_m8_33_init = mean_val
    v_8_33_init = mean_val
    slope_before_init = 0.0
    slope_after_init = 0.0

    initial_guess = [v_m8_33_init, v_8_33_init, slope_before_init, slope_after_init]

    def objective(params):
        v_m8_33, v_8_33, slope_before, slope_after = params

        predictions = []
        for age in test_ages:
            pred = pred_age_piecewise_linear_2dots(age, survival_data, v_m8_33, v_8_33, slope_before, slope_after)
            if not np.isfinite(pred):
                return np.inf
            predictions.append(pred)

        predictions = np.array(predictions)
        residuals = test_values - predictions

        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
        else:
            rss = np.sum(residuals**2)

        return rss

    # Set reasonable bounds
    data_range = test_values.max() - test_values.min()
    bounds = [
        (-2 * data_range, 2 * data_range),  # v_m8_33
        (-2 * data_range, 2 * data_range),  # v_8_33
        (-1.0, 1.0),  # slope_before
        (-1.0, 1.0),  # slope_after
    ]

    try:
        result = minimize(
            objective,
            initial_guess,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 200, "ftol": 1e-6, "gtol": 1e-6},
        )

        if not result.success:
            print(f"   ⚠️ Piecewise linear 2-dots optimization did not fully converge: {result.message}")

        v_m8_33, v_8_33, slope_before, slope_after = result.x

        predictions = []
        for age in test_ages:
            pred = pred_age_piecewise_linear_2dots(age, survival_data, v_m8_33, v_8_33, slope_before, slope_after)
            predictions.append(pred)
        predictions = np.array(predictions)

        residuals = test_values - predictions

        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
            tss = np.sum(
                weights * (test_values - np.average(test_values, weights=weights)) ** 2
            )
        else:
            rss = np.sum(residuals**2)
            tss = np.sum((test_values - np.mean(test_values)) ** 2)

        r_squared = 1 - rss / tss if tss > 0 else 0

        pwl2_model = {
            "v_m8_33": v_m8_33,
            "v_8_33": v_8_33,
            "slope_before": slope_before,
            "slope_after": slope_after,
            "r_squared": r_squared,
            "rss": rss,
            "fitted_values": predictions,
            "residuals": residuals,
            "n_points": len(test_values),
        }

        print(f"   ✅ Piecewise linear 2-dots model fitted: R² = {r_squared:.6f}")

        return pwl2_model

    except Exception as e:
        print(f"   ❌ Error fitting piecewise linear 2-dots model: {e}")
        return None


def fit_linear_pre_only_model(
    test_ages,
    test_values,
    survival_data,
    linear_model,
    test_sds=None,
    test_ns=None,
    use_iterative_reweighting=False,
):
    """Fit linear pre-menopause only model (flat after menopause)"""
    print("🔬 FITTING LINEAR PRE-MENOPAUSE ONLY MODEL")

    test_sds_for_weights = test_sds

    # Initialize from linear model
    mean_val = np.mean(test_values)
    slope_pre_init = 0.0
    val_before_jump_init = mean_val
    val_after_jump_init = mean_val

    initial_guess = [slope_pre_init, val_before_jump_init, val_after_jump_init]

    def objective(params):
        slope_pre, val_before_jump, val_after_jump = params

        predictions = []
        for age in test_ages:
            pred = pred_age_linear_pre_only(age, survival_data, slope_pre, val_before_jump, val_after_jump)
            if not np.isfinite(pred):
                return np.inf
            predictions.append(pred)

        predictions = np.array(predictions)
        residuals = test_values - predictions

        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
        else:
            rss = np.sum(residuals**2)

        return rss

    # Set reasonable bounds
    data_range = test_values.max() - test_values.min()
    bounds = [
        (-1.0, 1.0),  # slope_pre
        (-2 * data_range, 2 * data_range),  # val_before_jump
        (-2 * data_range, 2 * data_range),  # val_after_jump
    ]

    try:
        result = minimize(
            objective,
            initial_guess,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 200, "ftol": 1e-6, "gtol": 1e-6},
        )

        if not result.success:
            print(f"   ⚠️ Linear pre-only optimization did not fully converge: {result.message}")

        slope_pre, val_before_jump, val_after_jump = result.x

        predictions = []
        for age in test_ages:
            pred = pred_age_linear_pre_only(age, survival_data, slope_pre, val_before_jump, val_after_jump)
            predictions.append(pred)
        predictions = np.array(predictions)

        residuals = test_values - predictions

        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
            tss = np.sum(
                weights * (test_values - np.average(test_values, weights=weights)) ** 2
            )
        else:
            rss = np.sum(residuals**2)
            tss = np.sum((test_values - np.mean(test_values)) ** 2)

        r_squared = 1 - rss / tss if tss > 0 else 0

        lin_pre_model = {
            "slope_pre": slope_pre,
            "val_before_jump": val_before_jump,
            "val_after_jump": val_after_jump,
            "r_squared": r_squared,
            "rss": rss,
            "fitted_values": predictions,
            "residuals": residuals,
            "n_points": len(test_values),
        }

        print(f"   ✅ Linear pre-only model fitted: R² = {r_squared:.6f}")

        return lin_pre_model

    except Exception as e:
        print(f"   ❌ Error fitting linear pre-only model: {e}")
        return None


def fit_linear_post_only_model(
    test_ages,
    test_values,
    survival_data,
    linear_model,
    test_sds=None,
    test_ns=None,
    use_iterative_reweighting=False,
):
    """Fit linear post-menopause only model (flat before menopause)"""
    print("🔬 FITTING LINEAR POST-MENOPAUSE ONLY MODEL")

    test_sds_for_weights = test_sds

    # Initialize from linear model
    mean_val = np.mean(test_values)
    val_before_jump_init = mean_val
    val_after_jump_init = mean_val
    slope_post_init = 0.0

    initial_guess = [val_before_jump_init, val_after_jump_init, slope_post_init]

    def objective(params):
        val_before_jump, val_after_jump, slope_post = params

        predictions = []
        for age in test_ages:
            pred = pred_age_linear_post_only(age, survival_data, val_before_jump, val_after_jump, slope_post)
            if not np.isfinite(pred):
                return np.inf
            predictions.append(pred)

        predictions = np.array(predictions)
        residuals = test_values - predictions

        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
        else:
            rss = np.sum(residuals**2)

        return rss

    # Set reasonable bounds
    data_range = test_values.max() - test_values.min()
    bounds = [
        (-2 * data_range, 2 * data_range),  # val_before_jump
        (-2 * data_range, 2 * data_range),  # val_after_jump
        (-1.0, 1.0),  # slope_post
    ]

    try:
        result = minimize(
            objective,
            initial_guess,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 200, "ftol": 1e-6, "gtol": 1e-6},
        )

        if not result.success:
            print(f"   ⚠️ Linear post-only optimization did not fully converge: {result.message}")

        val_before_jump, val_after_jump, slope_post = result.x

        predictions = []
        for age in test_ages:
            pred = pred_age_linear_post_only(age, survival_data, val_before_jump, val_after_jump, slope_post)
            predictions.append(pred)
        predictions = np.array(predictions)

        residuals = test_values - predictions

        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
            tss = np.sum(
                weights * (test_values - np.average(test_values, weights=weights)) ** 2
            )
        else:
            rss = np.sum(residuals**2)
            tss = np.sum((test_values - np.mean(test_values)) ** 2)

        r_squared = 1 - rss / tss if tss > 0 else 0

        lin_post_model = {
            "val_before_jump": val_before_jump,
            "val_after_jump": val_after_jump,
            "slope_post": slope_post,
            "r_squared": r_squared,
            "rss": rss,
            "fitted_values": predictions,
            "residuals": residuals,
            "n_points": len(test_values),
        }

        print(f"   ✅ Linear post-only model fitted: R² = {r_squared:.6f}")

        return lin_post_model

    except Exception as e:
        print(f"   ❌ Error fitting linear post-only model: {e}")
        return None


def fit_polynomial4_model(
    test_ages,
    test_values,
    survival_data,
    linear_model,
    test_sds=None,
    test_ns=None,
    use_iterative_reweighting=False,
):
    """Fit 4th order polynomial model"""
    print("🔬 FITTING 4TH ORDER POLYNOMIAL MODEL")
    
    test_sds_for_weights = test_sds
    
    # Initialize with simple values
    a0_init = np.mean(test_values)
    a1_init = 0.0
    a2_init = 0.0
    a3_init = 0.0
    a4_init = 0.0
    
    initial_guess = [a0_init, a1_init, a2_init, a3_init, a4_init]
    
    def objective(params):
        a0, a1, a2, a3, a4 = params
        
        predictions = []
        for age in test_ages:
            pred = pred_age_polynomial4(age, survival_data, a0, a1, a2, a3, a4)
            if not np.isfinite(pred):
                return np.inf
            predictions.append(pred)
        
        predictions = np.array(predictions)
        residuals = test_values - predictions
        
        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
        else:
            rss = np.sum(residuals**2)
        
        return rss
    
    # Set reasonable bounds
    data_range = test_values.max() - test_values.min()
    bounds = [
        (-2 * data_range, 2 * data_range),  # a0
        (-1.0, 1.0),  # a1
        (-0.1, 0.1),  # a2
        (-0.01, 0.01),  # a3
        (-0.001, 0.001),  # a4
    ]
    
    try:
        result = minimize(
            objective,
            initial_guess,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 200, "ftol": 1e-6, "gtol": 1e-6},
        )
        
        if not result.success:
            print(f"   ⚠️ Polynomial4 optimization did not fully converge: {result.message}")
        
        a0, a1, a2, a3, a4 = result.x
        
        predictions = []
        for age in test_ages:
            pred = pred_age_polynomial4(age, survival_data, a0, a1, a2, a3, a4)
            predictions.append(pred)
        predictions = np.array(predictions)
        
        residuals = test_values - predictions
        
        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
            tss = np.sum(
                weights * (test_values - np.average(test_values, weights=weights)) ** 2
            )
        else:
            rss = np.sum(residuals**2)
            tss = np.sum((test_values - np.mean(test_values)) ** 2)
        
        r_squared = 1 - rss / tss if tss > 0 else 0
        
        poly4_model = {
            "a0": a0,
            "a1": a1,
            "a2": a2,
            "a3": a3,
            "a4": a4,
            "r_squared": r_squared,
            "rss": rss,
            "fitted_values": predictions,
            "residuals": residuals,
            "n_points": len(test_values),
        }
        
        print(f"   ✅ Polynomial4 model fitted: R² = {r_squared:.6f}")
        
        return poly4_model
        
    except Exception as e:
        print(f"   ❌ Error fitting polynomial4 model: {e}")
        return None


def fit_piecewise_linear_continuous_model(
    test_ages,
    test_values,
    survival_data,
    linear_model,
    test_sds=None,
    test_ns=None,
    use_iterative_reweighting=False,
):
    """Fit piecewise linear continuous model with transitions at -12.5, 0, 12.5"""
    print("🔬 FITTING PIECEWISE LINEAR CONTINUOUS MODEL")
    
    test_sds_for_weights = test_sds
    
    # Initialize from linear model
    mean_val = np.mean(test_values)
    v_m12_5_init = mean_val
    v_0_init = mean_val
    v_12_5_init = mean_val
    slope_before_init = 0.0
    slope_after_init = 0.0
    
    initial_guess = [v_m12_5_init, v_0_init, v_12_5_init, slope_before_init, slope_after_init]
    
    def objective(params):
        v_m12_5, v_0, v_12_5, slope_before, slope_after = params
        
        predictions = []
        for age in test_ages:
            pred = pred_age_piecewise_linear_continuous(age, survival_data, v_m12_5, v_0, v_12_5, slope_before, slope_after)
            if not np.isfinite(pred):
                return np.inf
            predictions.append(pred)
        
        predictions = np.array(predictions)
        residuals = test_values - predictions
        
        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
        else:
            rss = np.sum(residuals**2)
        
        return rss
    
    # Set reasonable bounds
    data_range = test_values.max() - test_values.min()
    bounds = [
        (-2 * data_range, 2 * data_range),  # v_m12_5
        (-2 * data_range, 2 * data_range),  # v_0
        (-2 * data_range, 2 * data_range),  # v_12_5
        (-1.0, 1.0),  # slope_before
        (-1.0, 1.0),  # slope_after
    ]
    
    try:
        result = minimize(
            objective,
            initial_guess,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 200, "ftol": 1e-6, "gtol": 1e-6},
        )
        
        if not result.success:
            print(f"   ⚠️ Piecewise linear continuous optimization did not fully converge: {result.message}")
        
        v_m12_5, v_0, v_12_5, slope_before, slope_after = result.x
        
        predictions = []
        for age in test_ages:
            pred = pred_age_piecewise_linear_continuous(age, survival_data, v_m12_5, v_0, v_12_5, slope_before, slope_after)
            predictions.append(pred)
        predictions = np.array(predictions)
        
        residuals = test_values - predictions
        
        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
            tss = np.sum(
                weights * (test_values - np.average(test_values, weights=weights)) ** 2
            )
        else:
            rss = np.sum(residuals**2)
            tss = np.sum((test_values - np.mean(test_values)) ** 2)
        
        r_squared = 1 - rss / tss if tss > 0 else 0
        
        pwl_model = {
            "v_m12_5": v_m12_5,
            "v_0": v_0,
            "v_12_5": v_12_5,
            "slope_before": slope_before,
            "slope_after": slope_after,
            "r_squared": r_squared,
            "rss": rss,
            "fitted_values": predictions,
            "residuals": residuals,
            "n_points": len(test_values),
        }
        
        print(f"   ✅ Piecewise linear continuous model fitted: R² = {r_squared:.6f}")
        
        return pwl_model
        
    except Exception as e:
        print(f"   ❌ Error fitting piecewise linear continuous model: {e}")
        return None


def fit_exp_linear_model(
    test_ages,
    test_values,
    survival_data,
    linear_model,
    test_sds=None,
    test_ns=None,
    use_iterative_reweighting=False,
):
    """Fit exponential-linear model (exp before meno, linear after, with jump)"""
    print("🔬 FITTING EXP-LINEAR MODEL")
    
    test_sds_for_weights = test_sds
    
    # Initialize with reasonable values
    mean_val = np.mean(test_values)
    exp_a_init = mean_val
    exp_b_init = 0.01  # Small exponential rate
    val_before_jump_init = mean_val
    val_after_jump_init = mean_val
    slope_post_init = 0.0
    
    initial_guess = [exp_a_init, exp_b_init, val_before_jump_init, val_after_jump_init, slope_post_init]
    
    def objective(params):
        exp_a, exp_b, val_before_jump, val_after_jump, slope_post = params
        
        predictions = []
        for age in test_ages:
            pred = pred_age_exp_linear(age, survival_data, exp_a, exp_b, val_before_jump, val_after_jump, slope_post)
            if not np.isfinite(pred):
                return np.inf
            predictions.append(pred)
        
        predictions = np.array(predictions)
        residuals = test_values - predictions
        
        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
        else:
            rss = np.sum(residuals**2)
        
        return rss
    
    # Set reasonable bounds
    data_range = test_values.max() - test_values.min()
    bounds = [
        (0.01, 2 * data_range),  # exp_a (must be positive for exponential)
        (-0.5, 0.5),  # exp_b (exponential rate)
        (-2 * data_range, 2 * data_range),  # val_before_jump
        (-2 * data_range, 2 * data_range),  # val_after_jump
        (-1.0, 1.0),  # slope_post
    ]
    
    try:
        result = minimize(
            objective,
            initial_guess,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 200, "ftol": 1e-6, "gtol": 1e-6},
        )
        
        if not result.success:
            print(f"   ⚠️ Exp-linear optimization did not fully converge: {result.message}")
        
        exp_a, exp_b, val_before_jump, val_after_jump, slope_post = result.x
        
        predictions = []
        for age in test_ages:
            pred = pred_age_exp_linear(age, survival_data, exp_a, exp_b, val_before_jump, val_after_jump, slope_post)
            predictions.append(pred)
        predictions = np.array(predictions)
        
        residuals = test_values - predictions
        
        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
            tss = np.sum(
                weights * (test_values - np.average(test_values, weights=weights)) ** 2
            )
        else:
            rss = np.sum(residuals**2)
            tss = np.sum((test_values - np.mean(test_values)) ** 2)
        
        r_squared = 1 - rss / tss if tss > 0 else 0
        
        exp_lin_model = {
            "exp_a": exp_a,
            "exp_b": exp_b,
            "val_before_jump": val_before_jump,
            "val_after_jump": val_after_jump,
            "slope_post": slope_post,
            "r_squared": r_squared,
            "rss": rss,
            "fitted_values": predictions,
            "residuals": residuals,
            "n_points": len(test_values),
        }
        
        print(f"   ✅ Exp-linear model fitted: R² = {r_squared:.6f}")
        
        return exp_lin_model
        
    except Exception as e:
        print(f"   ❌ Error fitting exp-linear model: {e}")
        return None


def fit_linear_exp_model(
    test_ages,
    test_values,
    survival_data,
    linear_model,
    test_sds=None,
    test_ns=None,
    use_iterative_reweighting=False,
):
    """Fit linear-exponential model (linear before meno, exp after, with jump)"""
    print("🔬 FITTING LINEAR-EXP MODEL")
    
    test_sds_for_weights = test_sds
    
    # Initialize with reasonable values
    mean_val = np.mean(test_values)
    slope_pre_init = 0.0
    val_before_jump_init = mean_val
    val_after_jump_init = mean_val
    exp_a_init = mean_val
    exp_b_init = 0.01  # Small exponential rate
    
    initial_guess = [slope_pre_init, val_before_jump_init, val_after_jump_init, exp_a_init, exp_b_init]
    
    def objective(params):
        slope_pre, val_before_jump, val_after_jump, exp_a, exp_b = params
        
        predictions = []
        for age in test_ages:
            pred = pred_age_linear_exp(age, survival_data, slope_pre, val_before_jump, val_after_jump, exp_a, exp_b)
            if not np.isfinite(pred):
                return np.inf
            predictions.append(pred)
        
        predictions = np.array(predictions)
        residuals = test_values - predictions
        
        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
        else:
            rss = np.sum(residuals**2)
        
        return rss
    
    # Set reasonable bounds
    data_range = test_values.max() - test_values.min()
    bounds = [
        (-1.0, 1.0),  # slope_pre
        (-2 * data_range, 2 * data_range),  # val_before_jump
        (-2 * data_range, 2 * data_range),  # val_after_jump
        (0.01, 2 * data_range),  # exp_a (must be positive for exponential)
        (-0.5, 0.5),  # exp_b (exponential rate)
    ]
    
    try:
        result = minimize(
            objective,
            initial_guess,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 200, "ftol": 1e-6, "gtol": 1e-6},
        )
        
        if not result.success:
            print(f"   ⚠️ Linear-exp optimization did not fully converge: {result.message}")
        
        slope_pre, val_before_jump, val_after_jump, exp_a, exp_b = result.x
        
        predictions = []
        for age in test_ages:
            pred = pred_age_linear_exp(age, survival_data, slope_pre, val_before_jump, val_after_jump, exp_a, exp_b)
            predictions.append(pred)
        predictions = np.array(predictions)
        
        residuals = test_values - predictions
        
        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds_for_weights is not None
            and np.any(test_sds_for_weights > 0)
        ):
            weights = 1.0 / (test_sds_for_weights**2)
            rss = np.sum(weights * residuals**2)
            tss = np.sum(
                weights * (test_values - np.average(test_values, weights=weights)) ** 2
            )
        else:
            rss = np.sum(residuals**2)
            tss = np.sum((test_values - np.mean(test_values)) ** 2)
        
        r_squared = 1 - rss / tss if tss > 0 else 0
        
        lin_exp_model = {
            "slope_pre": slope_pre,
            "val_before_jump": val_before_jump,
            "val_after_jump": val_after_jump,
            "exp_a": exp_a,
            "exp_b": exp_b,
            "r_squared": r_squared,
            "rss": rss,
            "fitted_values": predictions,
            "residuals": residuals,
            "n_points": len(test_values),
        }
        
        print(f"   ✅ Linear-exp model fitted: R² = {r_squared:.6f}")
        
        return lin_exp_model
        
    except Exception as e:
        print(f"   ❌ Error fitting linear-exp model: {e}")
        return None





def fit_staged_optimization(
    cp1,
    cp2,
    test_ages,
    test_values,
    survival_data,
    linear_model,
    test_sds=None,
    test_ns=None,
    use_iterative_reweighting=False,
):
    """STAGED OPTIMIZATION APPROACH - cp1 and cp2 will be optimized as variables"""
    print(f"   🎯 Starting TWO CP OPTIMIZATION (cp1 and cp2 will be optimized)")
    # SIMPLE FIX: Use test_sds directly (same as linear model) - don't transform to SE!
    # This ensures consistent weighting between linear and change point models
    weights_se = test_sds

    # Reduced debug output to prevent I/O contention
    print(f"   🎯 STAGED OPTIMIZATION (cp1 and cp2 will be optimized)")

    # Stage 1: Estimate slopes from data outside menopause transition
    # print(f"      Stage 1: Estimating slopes outside menopause transition")

    # Use slopes from the linear model (disable data-based slope estimation)
    pre_slope = linear_model["d"]
    post_slope = linear_model["b"]

    # Stage 2: Estimate jump and change point values
    # print(f"      Stage 2: Estimating jump and change point values")

    data_mean = np.mean(test_values)
    data_std = np.std(test_values)
    data_min = np.min(test_values)
    data_max = np.max(test_values)

    # Two different starting strategies
    # Note: cp1, cp2 input parameters are ignored - they will be optimized
    starting_strategies = [
        {
            "name": "linear_based",
            "slope1": pre_slope,
            "slope2": post_slope,
            "val_before_jump": linear_model["c"],
            "val_after_jump": linear_model["a"],
            "val_cp1": linear_model["c"],
            "val_cp2": linear_model["a"],
        },
        {
            "name": "data_based",
            "slope1": pre_slope,
            "slope2": post_slope,
            "val_before_jump": data_mean - data_std,
            "val_after_jump": data_mean + data_std,
            "val_cp1": data_min,
            "val_cp2": data_max,
        },
    ]

    best_params = None
    best_loss = np.inf
    best_pure_rss = None  # Initialize to avoid undefined variable

    for strategy in starting_strategies:
        print(
            f"   🔄 Trying two CP strategy: {strategy['name']}"
        )

        initial_guess = [
            strategy["slope1"],
            strategy["val_cp1"],
            strategy["val_before_jump"],
            strategy["val_after_jump"],
            strategy["val_cp2"],
            strategy["slope2"],
            -10.0,  # cp1 initial guess
            10.0,   # cp2 initial guess
        ]

        def objective(params):
            slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, slope2, opt_cp1, opt_cp2 = params

            # Check for invalid parameters
            if not np.isfinite(params).all():
                return np.inf

            try:
                predictions = []
                for age in test_ages:
                    pred = pred_age(
                        age,
                        survival_data,
                        slope1,
                        val_cp1,
                        val_before_jump,
                        val_after_jump,
                        val_cp2,
                        slope2,
                        opt_cp1,
                        opt_cp2,
                    )
                    if not np.isfinite(pred):
                        print(
                            f"       ⚠️ Non-finite prediction for age {age}: pred={pred}, params={params}"
                        )
                        return np.inf
                    predictions.append(pred)

                predictions = np.array(predictions)

                # Debug: Check array shapes
                if len(predictions) != len(test_values):
                    print(
                        f"       ⚠️ Shape mismatch: predictions {len(predictions)}, test_values {len(test_values)}"
                    )
                    return np.inf

                residuals = test_values - predictions

                # Check for NaN/inf in residuals
                if not np.isfinite(residuals).all():
                    print(f"       ⚠️ Non-finite residuals detected in two CP objective")
                    return np.inf

            except Exception as e:
                print(
                    f"       ❌ Exception in two CP objective calculation: {type(e).__name__}: {str(e)}"
                )
                return np.inf

            # Check for non-finite predictions
            if len(predictions) > 0:
                pred_min, pred_max = np.min(predictions), np.max(predictions)
                if not np.isfinite(pred_min) or not np.isfinite(pred_max):
                    return np.inf

            # Use weighted RSS if standard deviations are available and CV is enabled
            # Calculate weighted or unweighted RSS
            if use_iterative_reweighting and "fit_linear_model_irls" in globals():
                # Use IRLS for robust fitting
                try:
                    # Create design matrix for IRLS using residuals instead of raw values
                    X = np.ones((len(test_ages), 1))
                    y = residuals  # Use residuals instead of test_values

                    # Initialize weights if we have standard deviations
                    initial_weights = None
                    if (
                        DEFAULT_CONFIG["USE_CV"]
                        and weights_se is not None
                        and np.any(weights_se > 0)
                    ):
                        initial_weights = 1.0 / (weights_se**2)
                        initial_weights = initial_weights / np.mean(
                            initial_weights
                        )  # Normalize

                    # Calculate standard RSS before IRLS for comparison
                    if (
                        DEFAULT_CONFIG["USE_CV"]
                        and weights_se is not None
                        and np.any(weights_se > 0)
                    ):
                        std_weights = 1.0 / (weights_se**2)
                        std_rss = np.sum(std_weights * residuals**2)
                    else:
                        std_rss = np.sum(residuals**2)

                    # Fit using IRLS to get weights based on residuals
                    _, _, _, weights = fit_linear_model_irls(
                        X,
                        y,
                        initial_weights=initial_weights,
                        max_iter=10,
                        weight_func="huber",
                    )

                    # Calculate weighted RSS with IRLS weights
                    rss = np.sum(weights * residuals**2)

                    # Check if any residuals are being significantly down-weighted
                    large_residuals = np.abs(residuals) > 2 * np.std(residuals)
                    if np.any(large_residuals):
                        n_outliers = np.sum(large_residuals)
                        mean_weight_outliers = np.mean(weights[large_residuals])

                except Exception as e:
                    # Fall back to standard weighting
                    if (
                        DEFAULT_CONFIG["USE_CV"]
                        and weights_se is not None
                        and np.any(weights_se > 0)
                    ):
                        weights = 1.0 / (weights_se**2)
                        rss = np.sum(weights * residuals**2)
                    else:
                        rss = np.sum(residuals**2)
            elif (
                DEFAULT_CONFIG["USE_CV"]
                and weights_se is not None
                and np.any(weights_se > 0)
            ):
                weights = 1.0 / (weights_se**2)
                rss = np.sum(weights * residuals**2)
            else:
                rss = np.sum(residuals**2)

            # Add slope penalties if enabled - FIXED VERSION
            total_penalty = 0.0
            if DEFAULT_CONFIG["ENABLE_SLOPE_PENALIZATION"]:
                jump_magnitude = val_after_jump - val_before_jump

                # Calculate slope from cp1 to menopause - FIXED
                if opt_cp1 != 0:
                    slope_cp1_to_meno = (val_before_jump - val_cp1) / (0 - opt_cp1)
                    penalty1 = calculate_slope_penalty(
                        slope_cp1_to_meno, jump_magnitude, opt_cp1
                    )
                    total_penalty += penalty1

                # Calculate slope from menopause to cp2 - FIXED
                if opt_cp2 != 0:
                    slope_meno_to_cp2 = (val_cp2 - val_after_jump) / (opt_cp2 - 0)
                    penalty2 = calculate_slope_penalty(
                        slope_meno_to_cp2, jump_magnitude, opt_cp2
                    )
                    total_penalty += penalty2

            # Store both penalized and unpenalized loss for later BIC calculation
            objective.pure_rss = rss
            return rss + total_penalty

        # Stage 3: Optimize with bounds
        test_range = test_values.max() - test_values.min()
        bounds = [
            (-2 * test_range, 2 * test_range),  # slope1
            (-2 * test_range, 3 * test_range),  # val_cp1
            (-test_range, 2 * test_range),  # val_before_jump
            (-test_range, 3 * test_range),  # val_after_jump
            (
                -2 * test_range,
                3 * test_range,
            ),  # val_cp2 - FIXED: was 2*, now 3* to match older version
            (-2 * test_range, 2 * test_range),  # slope2
            (-20.0, -3.0),  # cp1 position
            (3.0, 20.0),    # cp2 position
        ]

        try:
            print(f"Trying to optimize two CP with the strategy: {strategy['name']}")
            # Add timeout and more conservative settings
            result = minimize(
                objective,
                initial_guess,
                method="L-BFGS-B",
                bounds=bounds,
                # options={"maxiter": 250, "ftol": 1e-6, "gtol": 1e-6},
            )

            if result.success and result.fun < best_loss:
                slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, slope2, opt_cp1, opt_cp2 = (
                    result.x
                )
                print(
                    f"   ⚙️ Two CP optimization result: slope1={slope1:.6f}, slope2={slope2:.6f}"
                )
                print(
                    f"   ⚙️ Two CP values: val_cp1={val_cp1:.2f}, val_cp2={val_cp2:.2f}, positions cp1={opt_cp1:.2f}, cp2={opt_cp2:.2f}"
                )

                # Calculate slopes between change points for comparison
                if opt_cp1 != 0:
                    slope_from_cp1_to_meno = (val_before_jump - val_cp1) / (0 - opt_cp1)
                    print(
                        f"   ⚙️ Calculated slope_from_cp1_to_meno = {slope_from_cp1_to_meno:.6f}"
                    )
                else:
                    slope_from_cp1_to_meno = 0.0

                if opt_cp2 != 0:
                    slope_from_meno_to_cp2 = (val_cp2 - val_after_jump) / (opt_cp2 - 0)
                    print(
                        f"   ⚙️ Calculated slope_from_meno_to_cp2 = {slope_from_meno_to_cp2:.6f}"
                    )
                else:
                    slope_from_meno_to_cp2 = 0.0

                print(
                    f"   ⚙️ Comparing slopes: slope2={slope2:.6f}, slope_from_meno_to_cp2={slope_from_meno_to_cp2:.6f}"
                )

                if abs(slope2) > 1.0:
                    print(
                        f"   ⚠️ after optimization final slope2 is extreme slope2={slope2:.6f}"
                    )
                    print(
                        f"   ⚠️ CP values: val_cp1={val_cp1:.2f}, val_cp2={val_cp2:.2f}, positions cp1={opt_cp1:.2f}, cp2={opt_cp2:.2f}"
                    )

                best_params = result.x
                best_loss = result.fun
                # IMPORTANT FIX: Recalculate pure RSS using the BEST parameters from result.x
                # objective.pure_rss is from the LAST evaluation, not necessarily the best!
                print(f"   🔍 RECALC TWO CP: Using best params from result.x")
                print(f"      slope1={slope1:.6f}, slope2={slope2:.6f}")
                _, best_pure_rss = get_final_predictions_for_change_point(
                    test_ages,
                    survival_data,
                    slope1,
                    val_cp1,
                    val_before_jump,
                    val_after_jump,
                    val_cp2,
                    slope2,
                    opt_cp1,
                    opt_cp2,
                    test_values=test_values,
                    test_sds=test_sds,  # Use same weights as optimization
                    test_ns=None,  # Don't transform - use test_sds directly
                    return_rss=True,
                )
                print(f"   ✅ RECALC TWO CP: RSS={best_pure_rss:.6f}")
            else:
                print(
                    f"   ❌ Two CP optimization failed with the strategy: {strategy['name']}"
                )
                print(f"   ❌ Reason: {result.message}")
        except Exception as e:
            print(f"   ❌ Exception in two CP optimization: {str(e)}")

    if best_params is None:
        # Removed prints to eliminate I/O blocking
        return None, np.inf, None

    # Calculate final predictions and metrics
    slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, slope2, opt_cp1, opt_cp2 = best_params
    # Debug extreme slope2 values
    if abs(slope2) > 0.1:
        print(f"   ⚠️ after optimization final slope2 is extreme slope2={slope2:.6f}")
        print(
            f"   ⚠️ CP values: val_cp1={val_cp1:.2f}, val_cp2={val_cp2:.2f}, positions cp1={opt_cp1:.2f}, cp2={opt_cp2:.2f}"
        )

    final_predictions = []
    for age in test_ages:
        pred = pred_age(
            age,
            survival_data,
            slope1,
            val_cp1,
            val_before_jump,
            val_after_jump,
            val_cp2,
            slope2,
            opt_cp1,
            opt_cp2,
        )
        final_predictions.append(pred)

    final_predictions = np.array(final_predictions)
    residuals = test_values - final_predictions

    # Calculate R-squared (weighted if CV is enabled)
    # NOTE: This is just for logging, the actual pure_rss was already calculated above
    if DEFAULT_CONFIG["USE_CV"] and weights_se is not None and np.any(weights_se > 0):
        weights = 1.0 / (weights_se**2)
        rss = np.sum(weights * residuals**2)
        tss = np.sum(
            weights * (test_values - np.average(test_values, weights=weights)) ** 2
        )
    else:
        rss = np.sum(residuals**2)
        tss = np.sum((test_values - np.mean(test_values)) ** 2)
        r_squared = 1 - rss / tss if tss > 0 else 0

    # Removed prints to eliminate I/O blocking

    return best_params, best_loss, best_pure_rss


def fit_three_cp_optimization(
    test_ages,
    test_values,
    survival_data,
    linear_model,
    test_sds=None,
    test_ns=None,
    use_iterative_reweighting=False,
):
    """THREE CP OPTIMIZATION - all 3 CPs can be in [-20, 20], will be sorted in kernel"""
    print(f"   🎯 Starting THREE CP OPTIMIZATION (all CPs in [-20, 20], will be sorted)")
    weights_se = test_sds

    # Use slopes from the linear model
    pre_slope = linear_model["d"]
    post_slope = linear_model["b"]

    # Estimate values
    data_mean = np.mean(test_values)
    data_std = np.std(test_values)
    data_min = np.min(test_values)
    data_max = np.max(test_values)

    # Starting strategies
    starting_strategies = [
        {
            "name": "linear_based",
            "slope1": pre_slope,
            "slope2": post_slope,
            "val_before_jump": linear_model["c"],
            "val_after_jump": linear_model["a"],
            "val_cp1": linear_model["c"],
            "val_cp2": (linear_model["c"] + linear_model["a"]) / 2,
            "val_cp3": linear_model["a"],
        },
        {
            "name": "data_based",
            "slope1": pre_slope,
            "slope2": post_slope,
            "val_before_jump": data_mean - data_std,
            "val_after_jump": data_mean + data_std,
            "val_cp1": data_min,
            "val_cp2": data_mean,
            "val_cp3": data_max,
        },
    ]

    best_params = None
    best_loss = np.inf
    best_pure_rss = None

    for strategy in starting_strategies:
        print(f"   🔄 Trying three CP strategy: {strategy['name']}")

        initial_guess = [
            strategy["slope1"],
            strategy["val_cp1"],
            strategy["val_before_jump"],
            strategy["val_after_jump"],
            strategy["val_cp2"],
            strategy["val_cp3"],
            strategy["slope2"],
            -10.0,  # cp1 initial guess
            0.0,    # cp2 initial guess (near menopause)
            10.0,   # cp3 initial guess
        ]

        def objective(params):
            slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, val_cp3, slope2, opt_cp1, opt_cp2, opt_cp3 = params

            # Check for invalid parameters
            if not np.isfinite(params).all():
                return np.inf

            try:
                predictions = []
                for age in test_ages:
                    pred = pred_age_three_cp(
                        age,
                        survival_data,
                        slope1,
                        val_cp1,
                        val_before_jump,
                        val_after_jump,
                        val_cp2,
                        val_cp3,
                        slope2,
                        opt_cp1,
                        opt_cp2,
                        opt_cp3,
                    )
                    if not np.isfinite(pred):
                        return np.inf
                    predictions.append(pred)

                predictions = np.array(predictions)

                if len(predictions) != len(test_values):
                    return np.inf

                residuals = test_values - predictions

                if not np.isfinite(residuals).all():
                    return np.inf

                # Calculate RSS
                if use_iterative_reweighting and "fit_linear_model_irls" in globals():
                    try:
                        X = np.ones((len(test_ages), 1))
                        y = residuals
                        initial_weights = None
                        if (
                            DEFAULT_CONFIG["USE_CV"]
                            and weights_se is not None
                            and np.any(weights_se > 0)
                        ):
                            initial_weights = 1.0 / (weights_se**2)
                        weights, _ = fit_linear_model_irls(X, y, initial_weights)
                        rss = np.sum(weights * residuals**2)
                    except Exception:
                        rss = np.sum(residuals**2)
                else:
                    if DEFAULT_CONFIG["USE_CV"] and weights_se is not None and np.any(weights_se > 0):
                        weights = 1.0 / (weights_se**2)
                        rss = np.sum(weights * residuals**2)
                    else:
                        rss = np.sum(residuals**2)

                # Store pure RSS for BIC calculation
                objective.pure_rss = rss
                return rss

            except Exception as e:
                return np.inf

        # Bounds - all CPs can be in [-20, 20]
        test_range = test_values.max() - test_values.min()
        bounds = [
            (-2 * test_range, 2 * test_range),  # slope1
            (-2 * test_range, 3 * test_range),  # val_cp1
            (-test_range, 2 * test_range),  # val_before_jump
            (-test_range, 3 * test_range),  # val_after_jump
            (-2 * test_range, 3 * test_range),  # val_cp2
            (-2 * test_range, 3 * test_range),  # val_cp3
            (-2 * test_range, 2 * test_range),  # slope2
            (-20.0, 20.0),  # cp1 position - full range
            (-20.0, 20.0),  # cp2 position - full range
            (-20.0, 20.0),  # cp3 position - full range
        ]

        try:
            print(f"   🔄 Starting three CP optimization with strategy: {strategy['name']}")
            result = minimize(
                objective,
                initial_guess,
                method="L-BFGS-B",
                bounds=bounds,
            )

            if result.success and result.fun < best_loss:
                print(f"   ✅ Three CP optimization successful with strategy: {strategy['name']}")
                best_params = result.x
                best_loss = result.fun
                # Recalculate pure RSS
                slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, val_cp3, slope2, opt_cp1, opt_cp2, opt_cp3 = result.x
                
                # Sort CPs to show final positions
                sorted_cps = sorted([opt_cp1, opt_cp2, opt_cp3])
                print(f"   📊 Three CP positions (sorted): {sorted_cps[0]:.2f}, {sorted_cps[1]:.2f}, {sorted_cps[2]:.2f}")
                
                # Recalculate predictions for pure RSS
                predictions = []
                for age in test_ages:
                    pred = pred_age_three_cp(
                        age,
                        survival_data,
                        slope1,
                        val_cp1,
                        val_before_jump,
                        val_after_jump,
                        val_cp2,
                        val_cp3,
                        slope2,
                        opt_cp1,
                        opt_cp2,
                        opt_cp3,
                    )
                    predictions.append(pred)
                predictions = np.array(predictions)
                residuals = test_values - predictions
                if DEFAULT_CONFIG["USE_CV"] and weights_se is not None and np.any(weights_se > 0):
                    weights = 1.0 / (weights_se**2)
                    best_pure_rss = np.sum(weights * residuals**2)
                else:
                    best_pure_rss = np.sum(residuals**2)
                print(f"   ✅ RECALC THREE CP: RSS={best_pure_rss:.6f}")
            else:
                print(f"   ❌ Three CP optimization failed with strategy: {strategy['name']}")
        except Exception as e:
            print(f"   ❌ Exception in three CP optimization: {str(e)}")

    if best_params is None:
        return None, np.inf, None

    return best_params, best_loss, best_pure_rss


def generate_additional_single_cp_positions(n_positions):
    """Generate exactly n_positions single change point positions for task expansion"""

    # Check if we should only use positions before menopause
    if DEFAULT_CONFIG.get("SINGLE_CP_BEFORE_ONLY", False):
        # Use all positions for the before range
        n_before = n_positions
        n_after = 0
    else:
        # Split positions between before and after menopause ranges
        # Use 2/3 for before (more important) and 1/3 for after
        n_before = int(n_positions * 2 / 3)
        n_after = n_positions - n_before

    if DEFAULT_CONFIG.get("RANDOMIZE_COMBINATIONS", True):
        # Random sampling with different seed for each run
        import time

        seed = int(time.time() * 1000) % 1000000  # Use milliseconds for seed
        np.random.seed(seed)

        # Before menopause (negative positions) - randomly sample n_before positions
        cp_before = np.random.uniform(
            DEFAULT_CONFIG["first_change_point_range_min"],  # -15
            DEFAULT_CONFIG["first_change_point_range_max"],  # -2
            n_before,
        )

        if n_after > 0:
            # After menopause (positive positions) - randomly sample n_after positions
            cp_after = np.random.uniform(
                DEFAULT_CONFIG["second_change_point_range_min"],  # 2
                DEFAULT_CONFIG["second_change_point_range_max"],  # 15
                n_after,
            )
            # Combine both ranges
            all_positions = np.concatenate([cp_before, cp_after])
        else:
            # Only before positions
            all_positions = cp_before
    else:
        # Grid-based approach
        # Before menopause (negative positions)
        cp_before = np.linspace(
            DEFAULT_CONFIG["first_change_point_range_min"],  # -15
            DEFAULT_CONFIG["first_change_point_range_max"],  # -2
            n_before,
        )

        if n_after > 0:
            # After menopause (positive positions)
            cp_after = np.linspace(
                DEFAULT_CONFIG["second_change_point_range_min"],  # 2
                DEFAULT_CONFIG["second_change_point_range_max"],  # 15
                n_after,
            )
            # Combine both ranges
            all_positions = np.concatenate([cp_before, cp_after])
        else:
            # Only before positions
            all_positions = cp_before

    return all_positions


def generate_single_cp_positions(n_positions, specific_cp_a=None):
    """Generate single change point positions (before and after menopause)"""

    # If a specific CP is provided, use only that for the single CP
    if specific_cp_a is not None:
        print(f"   🎯 Using specific single change point: CP={specific_cp_a}")
        return np.array([specific_cp_a])

    # Check if we should only use positions before menopause
    if DEFAULT_CONFIG.get("SINGLE_CP_BEFORE_ONLY", False):
        print(
            f"   🔒 Single CP: Using ONLY positions BEFORE menopause (before-only mode)"
        )
        # Use all positions for the before range
        n_before = n_positions
        n_after = 0
    else:
        # Use 2 * sqrt(n_positions) for each range to double the coverage
        # This will generate approximately 2x more positions than before
        n_before = int(2 * np.sqrt(n_positions))
        n_after = int(2 * np.sqrt(n_positions))

    if DEFAULT_CONFIG.get("RANDOMIZE_COMBINATIONS", True):
        # Random sampling with different seed for each run
        import time

        seed = int(time.time() * 1000) % 1000000  # Use milliseconds for seed
        np.random.seed(seed)

        # Before menopause (negative positions) - randomly sample n_before positions
        cp_before = np.random.uniform(
            DEFAULT_CONFIG["first_change_point_range_min"],  # -15
            DEFAULT_CONFIG["first_change_point_range_max"],  # -2
            n_before,
        )

        if n_after > 0:
            # After menopause (positive positions) - randomly sample n_after positions
            cp_after = np.random.uniform(
                DEFAULT_CONFIG["second_change_point_range_min"],  # 2
                DEFAULT_CONFIG["second_change_point_range_max"],  # 15
                n_after,
            )
            # Combine both ranges
            all_positions = np.concatenate([cp_before, cp_after])
        else:
            # Only before positions
            all_positions = cp_before

        # If we have more positions than requested, sample randomly
        if len(all_positions) > n_positions:
            all_positions = np.random.choice(all_positions, n_positions, replace=False)
    else:
        # Grid-based approach (original)
        # Before menopause (negative positions)
        cp_before = np.linspace(
            DEFAULT_CONFIG["first_change_point_range_min"],  # -15
            DEFAULT_CONFIG["first_change_point_range_max"],  # -2
            n_before,
        )

        if n_after > 0:
            # After menopause (positive positions)
            cp_after = np.linspace(
                DEFAULT_CONFIG["second_change_point_range_min"],  # 2
                DEFAULT_CONFIG["second_change_point_range_max"],  # 15
                n_after,
            )
            # Combine both ranges
            all_positions = np.concatenate([cp_before, cp_after])
        else:
            # Only before positions
            all_positions = cp_before

        # If we have more positions than requested, sample evenly
        if len(all_positions) > n_positions:
            indices = np.linspace(0, len(all_positions) - 1, n_positions, dtype=int)
            all_positions = all_positions[indices]

    return all_positions


def generate_cp_combinations(
    n_combinations, fixed_seed=None, specific_cp_a=None, specific_cp_b=None
):
    """Generate systematic grid of change point combinations"""
    # If specific CPs are provided, use only that combination
    if specific_cp_a is not None and specific_cp_b is not None:
        print(
            f"   🎯 Using specific change point combination: CP1={specific_cp_a}, CP2={specific_cp_b}"
        )
        return [(specific_cp_a, specific_cp_b)]

    if DEFAULT_CONFIG.get("RANDOMIZE_COMBINATIONS", True):
        # Random sampling with fixed seed if provided, otherwise use time-based seed
        if fixed_seed is not None:
            np.random.seed(fixed_seed)
            print(
                f"   🔒 Using fixed random seed: {fixed_seed} for consistent bootstrap combinations"
            )
        else:
            import time

            seed = int(time.time() * 1000) % 1000000  # Use milliseconds for seed
            np.random.seed(seed)

        # Randomly sample n_combinations from the full range
        cp1_positions = np.random.uniform(
            DEFAULT_CONFIG["first_change_point_range_min"],
            DEFAULT_CONFIG["first_change_point_range_max"],
            n_combinations,
        )
        cp2_positions = np.random.uniform(
            DEFAULT_CONFIG["second_change_point_range_min"],
            DEFAULT_CONFIG["second_change_point_range_max"],
            n_combinations,
        )

        combinations = list(zip(cp1_positions, cp2_positions))
    else:
        # Grid-based approach (original)
        cp1_range = np.linspace(
            DEFAULT_CONFIG["first_change_point_range_min"],
            DEFAULT_CONFIG["first_change_point_range_max"],
            int(np.sqrt(n_combinations)),
        )
        cp2_range = np.linspace(
            DEFAULT_CONFIG["second_change_point_range_min"],
            DEFAULT_CONFIG["second_change_point_range_max"],
            int(np.sqrt(n_combinations)),
        )

        combinations = []
        for cp1 in cp1_range:
            for cp2 in cp2_range:
                combinations.append((cp1, cp2))

        combinations = combinations[:n_combinations]

    return combinations


def optimize_single_cp_position(
    cp_position,
    test_ages,
    test_values,
    survival_data,
    linear_model,
    test_sds=None,
    test_ns=None,
    use_iterative_reweighting=False,
):
    """
    Optimize a single change point position.
    This function is used with multiprocessing for parallel processing.
    """
    # Remove print statement to reduce I/O contention
    # print(f"   🎯 SINGLE CP OPTIMIZATION for cp={cp_position:.1f}")

    params, loss, pure_rss = fit_single_cp_optimization(
        cp_position,
        test_ages,
        test_values,
        survival_data,
        linear_model,
        test_sds,
        test_ns,
        use_iterative_reweighting=use_iterative_reweighting,
    )

    # Extract optimized cp_position from params (last element)
    opt_cp_position = params[-1] if params is not None else None

    return {
        "cp_position": opt_cp_position,  # Return optimized position
        "params": params,
        "loss": loss,
        "pure_rss": pure_rss,
        "success": params is not None,
    }


def optimize_two_change_point(
    cp1,
    cp2,
    test_ages,
    test_values,
    survival_data,
    linear_model,
    test_sds=None,
    test_ns=None,
    use_iterative_reweighting=False,
):
    """
    Optimize a two change pointt change point combination.
    This function is used with multiprocessing for parallel processing.
    """
    # Remove print statement to reduce I/O contention
    # print(f"   🎯 OPTIMIZATION for cp1={cp1:.1f}, cp2={cp2:.1f}")

    params, loss, pure_rss = fit_staged_optimization(
        cp1,
        cp2,
        test_ages,
        test_values,
        survival_data,
        linear_model,
        test_sds,
        test_ns,
        use_iterative_reweighting=use_iterative_reweighting,
    )

    # Extract optimized cp1, cp2 from params (last two elements)
    opt_cp1, opt_cp2 = (params[-2], params[-1]) if params is not None else (None, None)

    return {
        "cp1": opt_cp1,  # Return optimized positions
        "cp2": opt_cp2,
        "params": params,
        "loss": loss,
        "pure_rss": pure_rss,
        "success": params is not None,
    }


def get_final_predictions_for_change_point(
    test_ages,
    survival_data,
    slope1,
    val_cp1,
    val_before_jump,
    val_after_jump,
    val_cp2,
    slope2,
    cp1,
    cp2,
    test_values=None,
    test_sds=None,
    test_ns=None,
    return_rss=False,
):
    """Calculate final predictions for the fitted change point model

    Args:
        return_rss: If True, return (predictions, rss) with RSS calculated using proper weights
        test_values: Actual values (required if return_rss=True)
        test_sds: Original SD values (NOT transformed)
        test_ns: Sample sizes for weight calculation
    """
    predictions = []
    for age in test_ages:
        pred = pred_age(
            age,
            survival_data,
            slope1,
            val_cp1,
            val_before_jump,
            val_after_jump,
            val_cp2,
            slope2,
            cp1,
            cp2,
        )
        predictions.append(pred)
    predictions = np.array(predictions)

    if return_rss:
        if test_values is None:
            raise ValueError("test_values required when return_rss=True")

        residuals = test_values - predictions

        # Calculate RSS with the SAME weights used during optimization
        # SIMPLE FIX: Use test_sds directly (same as linear model)
        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds is not None
            and np.any(test_sds > 0)
        ):
            weights = 1.0 / (test_sds**2)
            rss = np.sum(weights * residuals**2)
        else:
            rss = np.sum(residuals**2)

        return predictions, rss

    return predictions


def get_final_predictions_for_single_change_point(
    test_ages,
    survival_data,
    slope1,
    val_cp1,
    val_before_jump,
    val_after_jump,
    slope2,
    cp_position,
    test_values=None,
    test_sds=None,
    test_ns=None,
    return_rss=False,
):
    """Calculate final predictions for the fitted single change point model

    Args:
        return_rss: If True, return (predictions, rss) with RSS calculated using proper weights
        test_values: Actual values (required if return_rss=True)
        test_sds: Original SD values (NOT transformed)
        test_ns: Sample sizes for weight calculation
    """
    predictions = []
    for age in test_ages:
        pred = pred_age_single_cp(
            age,
            survival_data,
            slope1,
            val_cp1,
            val_before_jump,
            val_after_jump,
            slope2,
            cp_position,
        )
        predictions.append(pred)
    predictions = np.array(predictions)

    if return_rss:
        if test_values is None:
            raise ValueError("test_values required when return_rss=True")

        residuals = test_values - predictions

        # Calculate RSS with the SAME weights used during optimization
        # SIMPLE FIX: Use test_sds directly (same as linear model)
        if (
            DEFAULT_CONFIG["USE_CV"]
            and test_sds is not None
            and np.any(test_sds > 0)
        ):
            weights = 1.0 / (test_sds**2)
            rss = np.sum(weights * residuals**2)
        else:
            rss = np.sum(residuals**2)

        return predictions, rss

    return predictions


def get_final_predictions_for_linear(test_ages, survival_data, a, b, c, d):
    """Calculate final predictions for the fitted linear model using integral approach"""
    predictions = []
    for age in test_ages:
        pred = pred_age_linear(age, survival_data, a, b, c, d)
        predictions.append(pred)
    return np.array(predictions)


def calculate_pure_rss(
    test_ages,
    test_values,
    survival_data,
    slope1=None,
    val_cp1=None,
    val_before_jump=None,
    val_after_jump=None,
    val_cp2=None,
    slope2=None,
    cp1=None,
    cp2=None,
    linear_a=None,
    linear_b=None,
    linear_c=None,
    linear_d=None,
    test_sds=None,
    test_ns=None,
):
    """
    Calculate pure weighted RSS (without penalties) from fitted parameters.

    This recalculates the predictions and computes weighted RSS without any slope penalties
    that may have been added during optimization. Uses the same weighting scheme as optimization.

    For linear models, pass linear_a, linear_b, linear_c, linear_d
    For change point models, pass the slope and cp parameters.
    
    IMPORTANT: Pass test_sds (SD) + test_ns for proper SE transformation!
    """
    # Check if this is a linear model
    if (
        linear_a is not None
        and linear_b is not None
        and linear_c is not None
        and linear_d is not None
    ):
        predictions = get_final_predictions_for_linear(
            test_ages,
            survival_data,
            linear_a,
            linear_b,
            linear_c,
            linear_d,
        )
    # Check if this is a single change point model (cp1 == cp2)
    elif cp1 == cp2:
        # Use single change point prediction function
        predictions = get_final_predictions_for_single_change_point(
            test_ages,
            survival_data,
            slope1,
            val_cp1,
            val_before_jump,
            val_after_jump,
            slope2,
            cp1,  # cp_position
        )
    else:
        # Use two change point prediction function
        predictions = get_final_predictions_for_change_point(
            test_ages,
            survival_data,
            slope1,
            val_cp1,
            val_before_jump,
            val_after_jump,
            val_cp2,
            slope2,
            cp1,
            cp2,
        )

    residuals = test_values - predictions

    # Use the same weighting scheme as optimization for BIC calculation
    # SIMPLE FIX: Use test_sds directly (same as linear model)
    if DEFAULT_CONFIG["USE_CV"] and test_sds is not None and np.any(test_sds > 0):
        weights = 1.0 / (test_sds**2)
        rss = np.sum(weights * residuals**2)
    else:
        rss = np.sum(residuals**2)

    return rss


def create_plots(
    test_ages,
    test_values,
    predictions,
    linear_model,
    slope1,
    val_cp1,
    val_before_jump,
    val_after_jump,
    val_cp2,
    slope2,
    best_cp1,
    best_cp2,
    r_squared,
    test_name,
    output_dir,
    test_sds=None,
    test_ns=None,
    single_cp_params=None,
    single_cp_predictions=None,
    single_cp_r_squared=None,
    single_cp_position=None,
    best_two_cp_params=None,
    survival_data=None,
    best_two_cp_cp1=None,
    best_two_cp_cp2=None,
    data_is_logit_transformed=None,
    sigmoid_model=None,
    poly3_model=None,
    poly4_model=None,
    pwl_model=None,
    pwl2_model=None,
    exp_lin_model=None,
    lin_exp_model=None,
    lin_pre_model=None,
    lin_post_model=None,
):
    """Create comprehensive plots"""

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))

    # Determine if data is logit-transformed
    # Use the parameter if provided, otherwise fall back to the global config
    is_logit_transformed = (
        data_is_logit_transformed
        if data_is_logit_transformed is not None
        else DEFAULT_CONFIG["USE_LOGIT_TRANSFORM"]
    )

    # Plot 1: Actual vs Predicted comparisons
    if (
        DEFAULT_CONFIG["USE_CV"]
        and test_sds is not None
        and test_ns is not None
        and np.any(test_sds > 0)
    ):
        # Calculate error bars - binomial CI for logit, standard error for regular data
        if is_logit_transformed:
            lower_error, upper_error = calculate_binomial_error_bars(
                test_values, test_ns, use_logit_transform=True
            )
            ax1.errorbar(
                test_ages,
                test_values,
                yerr=[lower_error, upper_error],
                fmt="o",
                alpha=0.7,
                label=f"Actual {test_name} (Logit)",
                color="blue",
                capsize=3,
                capthick=1,
            )
        else:
            # Calculate standard errors for regular data
            standard_errors = calculate_standard_error(test_sds, test_ns)
            ax1.errorbar(
                test_ages,
                test_values,
                yerr=standard_errors,
                fmt="o",
                alpha=0.7,
                label=f"Actual {test_name}",
                color="blue",
                capsize=3,
                capthick=1,
            )
    else:
        ax1.scatter(
            test_ages,
            test_values,
            alpha=0.7,
            label=f"Actual {test_name}",
            color="blue",
            s=30,
        )
    # Determine which model was selected
    is_single_cp_selected = best_cp1 == best_cp2

    # Calculate two CP model predictions - ALWAYS show both models
    if (
        best_two_cp_params is not None
        and best_two_cp_cp1 is not None
        and best_two_cp_cp2 is not None
    ):
        (
            two_cp_slope1,
            two_cp_val_cp1,
            two_cp_val_before_jump,
            two_cp_val_after_jump,
            two_cp_val_cp2,
            two_cp_slope2,
            two_cp_opt_cp1,
            two_cp_opt_cp2,
        ) = best_two_cp_params
        # Debug extreme slope2 values
        if abs(slope2) > 0.1:
            print(f"   ⚠️ create_plots called with extreme slope2={slope2:.6f}")
            print(
                f"   ⚠️ CP positions: cp1={two_cp_val_cp1:.2f}, cp2={two_cp_val_cp2:.2f}. for cp1={best_two_cp_cp1:.2f}, cp2={best_two_cp_cp2:.2f}"
            )

        two_cp_predictions = []
        for age in test_ages:
            pred = pred_age(
                age,
                survival_data,
                two_cp_slope1,
                two_cp_val_cp1,
                two_cp_val_before_jump,
                two_cp_val_after_jump,
                two_cp_val_cp2,
                two_cp_slope2,
                two_cp_opt_cp1,  # Use optimized positions from params
                two_cp_opt_cp2,
            )
            two_cp_predictions.append(pred)
        two_cp_predictions = np.array(two_cp_predictions)

        # Plot two CP model (check --only-basic-models flag)
        if not DEFAULT_CONFIG.get("ONLY_BASIC_MODELS", False):
            ax1.plot(
                test_ages,
                two_cp_predictions,
                (
                    "r-" if not is_single_cp_selected else "r:"
                ),  # Solid if selected, dotted if not
                linewidth=2,
                label=f"Two CP Model\n(R² = {r_squared:.4f})",
            )

    # Plot linear model (always dotted)
    ax1.plot(
        test_ages,
        linear_model["fitted_values"],
        "g--",
        linewidth=2,
        label=f'Linear Menopause\n(R² = {linear_model["r_squared"]:.4f})',
    )

    # Plot sigmoid model if available
    if sigmoid_model is not None:
        ax1.plot(
            test_ages,
            sigmoid_model["fitted_values"],
            "m-.",
            linewidth=2,
            label=f'Sigmoid Menopause\n(R² = {sigmoid_model["r_squared"]:.4f})',
        )

    # Plot polynomial4 model if available
    if poly4_model is not None:
        ax1.plot(
            test_ages,
            poly4_model["fitted_values"],
            "c-.",
            linewidth=2,
            label=f'Polynomial4\n(R² = {poly4_model["r_squared"]:.4f})',
        )

    # Plot piecewise linear continuous model if available
    if pwl_model is not None:
        ax1.plot(
            test_ages,
            pwl_model["fitted_values"],
            "y-.",
            linewidth=2,
            label=f'Piecewise Linear\n(R² = {pwl_model["r_squared"]:.4f})',
        )

    # Plot exp-linear model if available
    if exp_lin_model is not None:
        ax1.plot(
            test_ages,
            exp_lin_model["fitted_values"],
            color="darkgreen",
            linestyle="-.",
            linewidth=2,
            label=f'Exp-Linear\n(R² = {exp_lin_model["r_squared"]:.4f})',
        )

    # Plot linear-exp model if available
    if lin_exp_model is not None:
        ax1.plot(
            test_ages,
            lin_exp_model["fitted_values"],
            color="darkblue",
            linestyle="-.",
            linewidth=2,
            label=f'Linear-Exp\n(R² = {lin_exp_model["r_squared"]:.4f})',
        )

    # Plot polynomial3 model if available
    if poly3_model is not None:
        ax1.plot(
            test_ages,
            poly3_model["fitted_values"],
            color="orange",
            linestyle="-.",
            linewidth=2,
            label=f'Polynomial3\n(R² = {poly3_model["r_squared"]:.4f})',
        )

    # Plot piecewise linear 2-dots model if available
    if pwl2_model is not None:
        ax1.plot(
            test_ages,
            pwl2_model["fitted_values"],
            color="brown",
            linestyle="-.",
            linewidth=2,
            label=f'PWL-2dots\n(R² = {pwl2_model["r_squared"]:.4f})',
        )

    # Plot linear pre-only model if available
    if lin_pre_model is not None:
        ax1.plot(
            test_ages,
            lin_pre_model["fitted_values"],
            color="pink",
            linestyle="-.",
            linewidth=2,
            label=f'Linear Pre\n(R² = {lin_pre_model["r_squared"]:.4f})',
        )

    # Plot linear post-only model if available
    if lin_post_model is not None:
        ax1.plot(
            test_ages,
            lin_post_model["fitted_values"],
            color="purple",
            linestyle="-.",
            linewidth=2,
            label=f'Linear Post\n(R² = {lin_post_model["r_squared"]:.4f})',
        )

    # Plot single CP model if available (check --only-basic-models flag)
    if single_cp_predictions is not None and not DEFAULT_CONFIG.get("ONLY_BASIC_MODELS", False):
        ax1.plot(
            test_ages,
            single_cp_predictions,
            "b-" if is_single_cp_selected else "b:",  # Solid if selected, dotted if not
            linewidth=2,
            label=f"Single CP Model\n(R² = {single_cp_r_squared:.4f})",
        )
    ax1.set_xlabel("Age (years)")
    ax1.set_ylabel(f"{test_name} Level")
    ax1.set_title(f"Staged Change Point Model vs Linear Menopause Model - {test_name}")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Biological kernel functions comparison
    t_range = np.linspace(-20, 20, 4000)

    cp_kernel_values = [
        biological_kernel(
            t,
            slope1,
            val_cp1,
            val_before_jump,
            val_after_jump,
            val_cp2,
            slope2,
            best_cp1,
            best_cp2,
        )
        for t in t_range
    ]
    linear_kernel_values = [
        linear_menopause_biological_kernel(
            t,
            linear_model["a"],
            linear_model["b"],
            linear_model["c"],
            linear_model["d"],
            None,
        )
        for t in t_range
    ]

    # Determine which model was selected
    is_single_cp_selected = best_cp1 == best_cp2

    # Get best two CP model parameters - ALWAYS show both models
    if (
        best_two_cp_params is not None
        and best_two_cp_cp1 is not None
        and best_two_cp_cp2 is not None
    ):
        (
            two_cp_slope1,
            two_cp_val_cp1,
            two_cp_val_before_jump,
            two_cp_val_after_jump,
            two_cp_val_cp2,
            two_cp_slope2,
            two_cp_opt_cp1,
            two_cp_opt_cp2,
        ) = best_two_cp_params

        # Calculate two CP model kernel values
        two_cp_kernel_values = [
            biological_kernel(
                t,
                two_cp_slope1,
                two_cp_val_cp1,
                two_cp_val_before_jump,
                two_cp_val_after_jump,
                two_cp_val_cp2,
                two_cp_slope2,
                two_cp_opt_cp1,  # Use optimized positions from params
                two_cp_opt_cp2,
            )
            for t in t_range
        ]

        # Plot two CP model (always show, solid if selected, dotted if not)
        ax2.plot(
            t_range,
            two_cp_kernel_values,
            (
                "r-" if not is_single_cp_selected else "r:"
            ),  # Solid if selected, dotted if not
            linewidth=3,
            label=f"Two CP Model\n(R²={r_squared:.4f})",
            alpha=0.8,
        )

    # Plot linear model (always dotted)
    ax2.plot(
        t_range,
        linear_kernel_values,
        "g--",
        linewidth=3,
        label=f'Linear Menopause Model\n(R²={linear_model["r_squared"]:.4f})',
        alpha=0.8,
    )

    # Plot sigmoid model if available
    if sigmoid_model is not None:
        sigmoid_kernel_values = [
            biological_kernel_sigmoid(
                t,
                sigmoid_model["c"],
                sigmoid_model["d"],
                sigmoid_model["w"],
                sigmoid_model["h"],
                sigmoid_model["b"],
            )
            for t in t_range
        ]
        ax2.plot(
            t_range,
            sigmoid_kernel_values,
            "m-.",
            linewidth=3,
            label=f'Sigmoid Model\n(R²={sigmoid_model["r_squared"]:.4f}, w={sigmoid_model["w"]:.1f}y)',
            alpha=0.8,
        )

    # Plot polynomial4 model if available
    if poly4_model is not None:
        poly4_kernel_values = [
            biological_kernel_polynomial4(t, poly4_model.get("a0", 0), poly4_model.get("a1", 0), 
                                         poly4_model.get("a2", 0), poly4_model.get("a3", 0), 
                                         poly4_model.get("a4", 0))
            for t in t_range
        ]
        ax2.plot(
            t_range,
            poly4_kernel_values,
            "c-.",
            linewidth=3,
            label=f'Polynomial4\n(R²={poly4_model["r_squared"]:.4f})',
            alpha=0.8,
        )

    # Plot piecewise linear continuous model if available
    if pwl_model is not None:
        pwl_kernel_values = [
            biological_kernel_piecewise_linear_continuous(t, pwl_model.get("v_m12_5", 0), 
                                                         pwl_model.get("v_0", 0), 
                                                         pwl_model.get("v_12_5", 0),
                                                         pwl_model.get("slope_before", 0), 
                                                         pwl_model.get("slope_after", 0))
            for t in t_range
        ]
        ax2.plot(
            t_range,
            pwl_kernel_values,
            "y-.",
            linewidth=3,
            label=f'Piecewise Linear\n(R²={pwl_model["r_squared"]:.4f})',
            alpha=0.8,
        )

    # Plot exp-linear model if available
    if exp_lin_model is not None:
        exp_lin_kernel_values = [
            biological_kernel_exp_linear(t, exp_lin_model.get("exp_a", 0), 
                                        exp_lin_model.get("exp_b", 0),
                                        exp_lin_model.get("val_before_jump", 0),
                                        exp_lin_model.get("val_after_jump", 0),
                                        exp_lin_model.get("slope_post", 0))
            for t in t_range
        ]
        ax2.plot(
            t_range,
            exp_lin_kernel_values,
            color="darkgreen",
            linestyle="-.",
            linewidth=3,
            label=f'Exp-Linear\n(R²={exp_lin_model["r_squared"]:.4f})',
            alpha=0.8,
        )

    # Plot linear-exp model if available
    if lin_exp_model is not None:
        lin_exp_kernel_values = [
            biological_kernel_linear_exp(t, lin_exp_model.get("slope_pre", 0), 
                                        lin_exp_model.get("val_before_jump", 0),
                                        lin_exp_model.get("val_after_jump", 0),
                                        lin_exp_model.get("exp_a", 0),
                                        lin_exp_model.get("exp_b", 0))
            for t in t_range
        ]
        ax2.plot(
            t_range,
            lin_exp_kernel_values,
            color="darkblue",
            linestyle="-.",
            linewidth=3,
            label=f'Linear-Exp\n(R²={lin_exp_model["r_squared"]:.4f})',
            alpha=0.8,
        )

    # Plot polynomial3 model if available
    if poly3_model is not None:
        poly3_kernel_values = [
            biological_kernel_polynomial3(t, poly3_model.get("a0", 0), poly3_model.get("a1", 0), 
                                         poly3_model.get("a2", 0), poly3_model.get("a3", 0))
            for t in t_range
        ]
        ax2.plot(
            t_range,
            poly3_kernel_values,
            color="orange",
            linestyle="-.",
            linewidth=3,
            label=f'Polynomial3\n(R²={poly3_model["r_squared"]:.4f})',
            alpha=0.8,
        )

    # Plot piecewise linear 2-dots model if available
    if pwl2_model is not None:
        pwl2_kernel_values = [
            biological_kernel_piecewise_linear_2dots(t, pwl2_model.get("v_m8_33", 0), 
                                                     pwl2_model.get("v_8_33", 0),
                                                     pwl2_model.get("slope_before", 0),
                                                     pwl2_model.get("slope_after", 0))
            for t in t_range
        ]
        ax2.plot(
            t_range,
            pwl2_kernel_values,
            color="brown",
            linestyle="-.",
            linewidth=3,
            label=f'PWL-2dots\n(R²={pwl2_model["r_squared"]:.4f})',
            alpha=0.8,
        )

    # Plot linear pre-only model if available
    if lin_pre_model is not None:
        lin_pre_kernel_values = [
            biological_kernel_linear_pre_only(t, lin_pre_model.get("slope_pre", 0), 
                                             lin_pre_model.get("val_before_jump", 0),
                                             lin_pre_model.get("val_after_jump", 0))
            for t in t_range
        ]
        ax2.plot(
            t_range,
            lin_pre_kernel_values,
            color="pink",
            linestyle="-.",
            linewidth=3,
            label=f'Linear Pre\n(R²={lin_pre_model["r_squared"]:.4f})',
            alpha=0.8,
        )

    # Plot linear post-only model if available
    if lin_post_model is not None:
        lin_post_kernel_values = [
            biological_kernel_linear_post_only(t, lin_post_model.get("val_before_jump", 0), 
                                              lin_post_model.get("val_after_jump", 0),
                                              lin_post_model.get("slope_post", 0))
            for t in t_range
        ]
        ax2.plot(
            t_range,
            lin_post_kernel_values,
            color="purple",
            linestyle="-.",
            linewidth=3,
            label=f'Linear Post\n(R²={lin_post_model["r_squared"]:.4f})',
            alpha=0.8,
        )

        # Plot single CP model if available
    if single_cp_params is not None:
        (
            single_slope1,
            single_val_cp1,
            single_val_before_jump,
            single_val_after_jump,
            single_slope2,
            single_cp_opt_position,  # Include optimized cp_position (6th parameter)
        ) = single_cp_params
        single_cp_kernel_values = [
            biological_kernel_single_cp(
                t,
                single_slope1,
                single_val_cp1,
                single_val_before_jump,
                single_val_after_jump,
                single_slope2,
                single_cp_position,  # cp1
            )
            for t in t_range
        ]
        ax2.plot(
            t_range,
            single_cp_kernel_values,
            "b-" if is_single_cp_selected else "b:",  # Solid if selected, dotted if not
            linewidth=3,
            label=f"Single CP Model\n(R²={single_cp_r_squared:.4f})",
            alpha=0.8,
        )

    # Mark key points for two CP model
    if (
        best_two_cp_params is not None
        and best_two_cp_cp1 is not None
        and best_two_cp_cp2 is not None
    ):
        (
            two_cp_slope1,
            two_cp_val_cp1,
            two_cp_val_before_jump,
            two_cp_val_after_jump,
            two_cp_val_cp2,
            two_cp_slope2,
            two_cp_opt_cp1,
            two_cp_opt_cp2,
        ) = best_two_cp_params

        ax2.scatter(
            [two_cp_opt_cp1],  # Use optimized position from params
            [two_cp_val_cp1],
            color="red",
            s=120,
            zorder=5,
            label=f"Two CP Model CP1 ({two_cp_opt_cp1:.1f}, {two_cp_val_cp1:.3f})",
        )
        ax2.scatter(
            [0],
            [two_cp_val_before_jump],
            color="orange",
            s=120,
            zorder=5,
            label=f"Two CP Model Pre-jump (0, {two_cp_val_before_jump:.3f})",
        )
        ax2.scatter(
            [0],
            [two_cp_val_after_jump],
            color="red",
            s=120,
            zorder=5,
            label=f"Two CP Model Post-jump (0, {two_cp_val_after_jump:.3f})",
        )
        # Only plot CP2 if it's a two CP model (different from CP1)
        if best_two_cp_cp1 != best_two_cp_cp2:
            ax2.scatter(
                [best_two_cp_cp2],
                [two_cp_val_cp2],
                color="green",
                s=120,
                zorder=5,
                label=f"Two CP Model CP2 ({best_two_cp_cp2:.1f}, {two_cp_val_cp2:.3f})",
            )

    # Mark key points for best single CP model if available
    if single_cp_params is not None and single_cp_position is not None:
        (
            single_slope1,
            single_val_cp1,
            single_val_before_jump,
            single_val_after_jump,
            single_slope2,
            single_cp_opt_position,  # Include optimized cp_position (6th parameter)
        ) = single_cp_params

        # Add single CP change points with different markers
        ax2.scatter(
            [single_cp_position],
            [single_val_cp1],
            color="blue",
            s=100,
            marker="s",  # Square marker
            zorder=6,
            label=f"Single CP Model CP ({single_cp_position:.1f}, {single_val_cp1:.3f})",
        )
        ax2.scatter(
            [0],
            [single_val_before_jump],
            color="brown",
            s=100,
            marker="s",  # Square marker
            zorder=6,
            label=f"Single CP Model Pre-jump (0, {single_val_before_jump:.3f})",
        )
        ax2.scatter(
            [0],
            [single_val_after_jump],
            color="pink",
            s=100,
            marker="s",  # Square marker
            zorder=6,
            label=f"Single CP Model Post-jump (0, {single_val_after_jump:.3f})",
        )

        # Add vertical lines for single CP
        ax2.axvline(
            x=single_cp_position, color="purple", linestyle=":", alpha=0.5, linewidth=1
        )

    ax2.axvline(x=0, color="black", linestyle="--", alpha=0.7, linewidth=2)
    # Plot vertical lines for the best two CP model
    if best_two_cp_cp1 is not None:
        ax2.axvline(
            x=best_two_cp_cp1, color="red", linestyle=":", alpha=0.7, linewidth=2
        )
    if best_two_cp_cp2 is not None and best_two_cp_cp1 != best_two_cp_cp2:
        ax2.axvline(
            x=best_two_cp_cp2, color="green", linestyle=":", alpha=0.7, linewidth=2
        )

    ax2.set_xlabel("Time relative to menopause (years)")
    ax2.set_ylabel("Kernel Value")
    ax2.set_title(f"Biological Kernel Functions - {test_name}")
    ax2.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(-20, 20)

    # Plot 3: RSS (Residual Sum of Squares) Comparison
    model_names = []
    rss_values = []
    model_descriptions = []
    
    # Linear model (4 params: 2 slopes + 2 jump values)
    model_names.append("Linear\n(4p)")
    rss_values.append(linear_model.get("rss", np.sum(linear_model["residuals"]**2)))
    model_descriptions.append("Linear with jump at menopause")
    
    # Sigmoid model (5 params)
    if sigmoid_model is not None:
        model_names.append("Sigmoid\n(5p)")
        rss_values.append(sigmoid_model.get("rss", np.sum(sigmoid_model["residuals"]**2)))
        model_descriptions.append("Smooth S-curve transition")
    
    # Polynomial4 model (5 params)
    if poly4_model is not None:
        model_names.append("Poly4\n(5p)")
        rss_values.append(poly4_model.get("rss", 0))
        model_descriptions.append("4th order polynomial")
    
    # Piecewise Linear Continuous model (5 params)
    if pwl_model is not None:
        model_names.append("PWL\n(5p)")
        rss_values.append(pwl_model.get("rss", 0))
        model_descriptions.append("Piecewise linear, continuous")
    
    # Exp-Linear model (5 params)
    if exp_lin_model is not None:
        model_names.append("Exp-Lin\\n(5p)")
        rss_values.append(exp_lin_model.get("rss", 0))
        model_descriptions.append("Exp before, linear after")
    
    # Linear-Exp model (5 params)
    if lin_exp_model is not None:
        model_names.append("Lin-Exp\\n(5p)")
        rss_values.append(lin_exp_model.get("rss", 0))
        model_descriptions.append("Linear before, exp after")
    
    # Polynomial3 model (4 params)
    if poly3_model is not None:
        model_names.append("Poly3\n(4p)")
        rss_values.append(poly3_model.get("rss", 0))
        model_descriptions.append("3rd order polynomial")
    
    # Piecewise Linear 2-dots model (4 params)
    if pwl2_model is not None:
        model_names.append("PWL-2\n(4p)")
        rss_values.append(pwl2_model.get("rss", 0))
        model_descriptions.append("Piecewise linear, 2 dots")
    
    # Linear Pre-only model (3 params)
    if lin_pre_model is not None:
        model_names.append("Lin-Pre\n(3p)")
        rss_values.append(lin_pre_model.get("rss", 0))
        model_descriptions.append("Linear before meno only")
    
    # Linear Post-only model (3 params)
    if lin_post_model is not None:
        model_names.append("Lin-Post\n(3p)")
        rss_values.append(lin_post_model.get("rss", 0))
        model_descriptions.append("Linear after meno only")
    
    
    # Single CP model (6 params) - check --only-basic-models flag
    if single_cp_predictions is not None and single_cp_r_squared is not None and not DEFAULT_CONFIG.get("ONLY_BASIC_MODELS", False):
        model_names.append("Single CP\n(6p)")
        single_cp_residuals = test_values - single_cp_predictions
        rss_values.append(np.sum(single_cp_residuals**2))
        model_descriptions.append("One change point + jump")
    
    # Two CP model (8 params) - check --only-basic-models flag
    if best_two_cp_params is not None and not DEFAULT_CONFIG.get("ONLY_BASIC_MODELS", False):
        model_names.append("Two CP\n(8p)")
        two_cp_residuals = test_values - predictions
        rss_values.append(np.sum(two_cp_residuals**2))
        model_descriptions.append("Two change points + jump")
    
    # Find best model (lowest RSS)
    best_rss_idx = np.argmin(rss_values)
    
    # Create colors: best model in color, others in grey
    colors = ['lightgrey'] * len(model_names)
    colors[best_rss_idx] = 'green'  # Highlight best model in green
    
    # Create bar plot
    x_pos = np.arange(len(model_names))
    ax3.bar(x_pos, rss_values, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(model_names, rotation=45, ha='right', fontsize=9)
    ax3.set_ylabel("RSS (Residual Sum of Squares)")
    ax3.set_title(f"Model Comparison: RSS (Lower is Better)\nBest: {model_names[best_rss_idx].split(chr(10))[0]} - {test_name}")
    ax3.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for i, (name, val) in enumerate(zip(model_names, rss_values)):
        ax3.text(i, val, f'{val:.2f}', ha='center', va='bottom', fontsize=8, fontweight='bold' if i == best_rss_idx else 'normal')

    # Plot 4: BIC (Bayesian Information Criterion) Comparison
    bic_names = []
    bic_values = []
    bic_n_params = []
    
    n_data = len(test_values)
    
    # Linear model (4 params)
    if "rss" in linear_model:
        linear_rss = linear_model["rss"]
        linear_n_params = 4
        linear_bic = n_data * np.log(linear_rss / n_data) + np.log(n_data) * linear_n_params
        bic_names.append("Linear\n(4p)")
        bic_values.append(linear_bic)
        bic_n_params.append(linear_n_params)
    
    # Sigmoid model (5 params)
    if sigmoid_model is not None and "rss" in sigmoid_model:
        sigmoid_rss = sigmoid_model["rss"]
        sigmoid_n_params = 5
        sigmoid_bic = n_data * np.log(sigmoid_rss / n_data) + np.log(n_data) * sigmoid_n_params
        bic_names.append("Sigmoid\n(5p)")
        bic_values.append(sigmoid_bic)
        bic_n_params.append(sigmoid_n_params)
    
    # Polynomial4 model (5 params)
    if poly4_model is not None and "rss" in poly4_model:
        poly4_rss = poly4_model["rss"]
        poly4_n_params = 5
        poly4_bic = n_data * np.log(poly4_rss / n_data) + np.log(n_data) * poly4_n_params
        bic_names.append("Poly4\n(5p)")
        bic_values.append(poly4_bic)
        bic_n_params.append(poly4_n_params)
    
    # Piecewise Linear Continuous model (5 params)
    if pwl_model is not None and "rss" in pwl_model:
        pwl_rss = pwl_model["rss"]
        pwl_n_params = 5
        pwl_bic = n_data * np.log(pwl_rss / n_data) + np.log(n_data) * pwl_n_params
        bic_names.append("PWL\n(5p)")
        bic_values.append(pwl_bic)
        bic_n_params.append(pwl_n_params)
    
    # Exp-Linear model (5 params)
    if exp_lin_model is not None and "rss" in exp_lin_model:
        exp_lin_rss = exp_lin_model["rss"]
        exp_lin_n_params = 5
        exp_lin_bic = n_data * np.log(exp_lin_rss / n_data) + np.log(n_data) * exp_lin_n_params
        bic_names.append("Exp-Lin\\n(5p)")
        bic_values.append(exp_lin_bic)
        bic_n_params.append(exp_lin_n_params)
    
    # Linear-Exp model (5 params)
    if lin_exp_model is not None and "rss" in lin_exp_model:
        lin_exp_rss = lin_exp_model["rss"]
        lin_exp_n_params = 5
        lin_exp_bic = n_data * np.log(lin_exp_rss / n_data) + np.log(n_data) * lin_exp_n_params
        bic_names.append("Lin-Exp\\n(5p)")
        bic_values.append(lin_exp_bic)
        bic_n_params.append(lin_exp_n_params)
    
    # Polynomial3 model (4 params)
    if poly3_model is not None and "rss" in poly3_model:
        poly3_rss = poly3_model["rss"]
        poly3_n_params = 4
        poly3_bic = n_data * np.log(poly3_rss / n_data) + np.log(n_data) * poly3_n_params
        bic_names.append("Poly3\n(4p)")
        bic_values.append(poly3_bic)
        bic_n_params.append(poly3_n_params)
    
    # Piecewise Linear 2-dots model (4 params)
    if pwl2_model is not None and "rss" in pwl2_model:
        pwl2_rss = pwl2_model["rss"]
        pwl2_n_params = 4
        pwl2_bic = n_data * np.log(pwl2_rss / n_data) + np.log(n_data) * pwl2_n_params
        bic_names.append("PWL-2\n(4p)")
        bic_values.append(pwl2_bic)
        bic_n_params.append(pwl2_n_params)
    
    # Linear Pre-only model (3 params)
    if lin_pre_model is not None and "rss" in lin_pre_model:
        lin_pre_rss = lin_pre_model["rss"]
        lin_pre_n_params = 3
        lin_pre_bic = n_data * np.log(lin_pre_rss / n_data) + np.log(n_data) * lin_pre_n_params
        bic_names.append("Lin-Pre\n(3p)")
        bic_values.append(lin_pre_bic)
        bic_n_params.append(lin_pre_n_params)
    
    # Linear Post-only model (3 params)
    if lin_post_model is not None and "rss" in lin_post_model:
        lin_post_rss = lin_post_model["rss"]
        lin_post_n_params = 3
        lin_post_bic = n_data * np.log(lin_post_rss / n_data) + np.log(n_data) * lin_post_n_params
        bic_names.append("Lin-Post\n(3p)")
        bic_values.append(lin_post_bic)
        bic_n_params.append(lin_post_n_params)
    
    # Single CP model (6 params) - check --only-basic-models flag
    if single_cp_predictions is not None and not DEFAULT_CONFIG.get("ONLY_BASIC_MODELS", False):
        single_cp_residuals = test_values - single_cp_predictions
        single_cp_rss = np.sum(single_cp_residuals**2)
        single_cp_n_params = 6
        single_cp_bic = n_data * np.log(single_cp_rss / n_data) + np.log(n_data) * single_cp_n_params
        bic_names.append("Single CP\n(6p)")
        bic_values.append(single_cp_bic)
        bic_n_params.append(single_cp_n_params)
    
    # Two CP model (8 params) - check --only-basic-models flag
    if best_two_cp_params is not None and not DEFAULT_CONFIG.get("ONLY_BASIC_MODELS", False):
        two_cp_residuals = test_values - predictions
        two_cp_rss = np.sum(two_cp_residuals**2)
        two_cp_n_params = 8
        two_cp_bic = n_data * np.log(two_cp_rss / n_data) + np.log(n_data) * two_cp_n_params
        bic_names.append("Two CP\n(8p)")
        bic_values.append(two_cp_bic)
        bic_n_params.append(two_cp_n_params)
    
    # Create bar plot
    if bic_names:
        # Find best model (lowest BIC)
        best_bic_idx = np.argmin(bic_values)
        
        # Create colors: best model in color, others in grey
        bic_colors = ['lightgrey'] * len(bic_names)
        bic_colors[best_bic_idx] = 'blue'  # Highlight best model in blue
        
        x_pos = np.arange(len(bic_names))
        ax4.bar(x_pos, bic_values, color=bic_colors, alpha=0.7, edgecolor='black', linewidth=1.5)
        ax4.set_xticks(x_pos)
        ax4.set_xticklabels(bic_names, rotation=45, ha='right', fontsize=9)
        ax4.set_ylabel("BIC (Bayesian Information Criterion)")
        ax4.set_title(f"Model Comparison: BIC (Lower is Better)\nBest: {bic_names[best_bic_idx].split(chr(10))[0]} - {test_name}")
        ax4.grid(True, alpha=0.3, axis='y')
        
        # Add value labels on bars
        for i, (name, val) in enumerate(zip(bic_names, bic_values)):
            ax4.text(i, val, f'{val:.1f}', ha='center', va='bottom', fontsize=8, fontweight='bold' if i == best_bic_idx else 'normal')
    
    # Fallback if no BIC values
    if not bic_names:
        ax4.text(0.5, 0.5, "No BIC values available", ha='center', va='center', transform=ax4.transAxes)
        ax4.set_title(f"Model Comparison: BIC - {test_name}")

    plt.tight_layout()

    # Save plot
    plot_path = os.path.join(output_dir, f"{test_name}_change_point_analysis.pdf")
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()

    return plot_path


def run_single_test_analysis(
    test_code,
    test_name,
    system_name,
    output_base_dir,
    bootstrap_folder=None,
    specific_cp_a=None,
    specific_cp_b=None,
    use_iterative_reweighting=False,
    override_bootstrap=False,
):
    """Run change point analysis for a single test"""

    print(f"\n🚀 CHANGE POINT ANALYSIS FOR {test_name} (Code: {test_code})")
    print("=" * 80)

    # Create output directories
    system_dir = os.path.join(output_base_dir, system_name)
    test_dir = os.path.join(system_dir, test_name.replace(" ", "_").replace("/", "_"))
    os.makedirs(test_dir, exist_ok=True)

    try:
        # 1. Load data
        print("1. Loading data...")

        if bootstrap_folder:
            # Load bootstrap data
            print(f"   Loading bootstrap data from: {bootstrap_folder}")
            bootstrap_files = [
                f for f in os.listdir(bootstrap_folder) if f.endswith(".csv")
            ]
            bootstrap_files.sort()  # Ensure consistent ordering

            print(
                f"   🔍 DEBUG: Limited to {len(bootstrap_files)} bootstrap file(s) for debugging"
            )

            if not bootstrap_files:
                print(
                    f"❌ ERROR: No CSV files found in bootstrap folder {bootstrap_folder}"
                )
                return None

            print(f"   Found {len(bootstrap_files)} bootstrap files")

            # Load first file to get structure
            first_data = pd.read_csv(os.path.join(bootstrap_folder, bootstrap_files[0]))

            # Determine the correct track name for the test code
            track_name = get_test_track_name(test_code)

            test_data = first_data[first_data["test_track"] == track_name]

            if len(test_data) == 0:
                print(
                    f"❌ ERROR: No data found for test code {test_code} in bootstrap files"
                )
                return None

            # Get age columns from first file
            age_columns = [
                col
                for col in first_data.columns[9:]
                if not col.endswith("_sd") and not col.endswith("_n")
            ]
            sd_columns = [col for col in first_data.columns[9:] if col.endswith("_sd")]
            n_columns = [col for col in first_data.columns[9:] if col.endswith("_n")]
            test_ages = np.array([float(col) for col in age_columns])

            # Load all bootstrap values
            all_bootstrap_values = []
            all_bootstrap_sds = []
            all_bootstrap_ns = []

            # First, determine the common valid mask across all bootstrap files
            common_valid_mask = None

            for bootstrap_file in bootstrap_files:
                data = pd.read_csv(os.path.join(bootstrap_folder, bootstrap_file))

                # Determine the correct track name for the test code
                track_name = get_test_track_name(test_code)

                test_data = data[data["test_track"] == track_name]

                if len(test_data) > 0:
                    test_values = np.array(
                        [test_data.iloc[0][col] for col in age_columns]
                    )
                    valid_mask = ~np.isnan(test_values)

                    if common_valid_mask is None:
                        common_valid_mask = valid_mask
                    else:
                        # Ensure all bootstrap files have the same valid mask
                        if not np.array_equal(common_valid_mask, valid_mask):
                            print(
                                f"   ⚠️  Warning: Different valid masks in bootstrap files, using intersection"
                            )
                            common_valid_mask = common_valid_mask & valid_mask

            if common_valid_mask is None:
                print(f"❌ ERROR: No valid data found in any bootstrap file")
                return None

            # Apply common mask to ages
            test_ages = test_ages[common_valid_mask]

            # Apply age limits if specified (for bootstrap data)
            if (
                "MIN_AGE" in DEFAULT_CONFIG and DEFAULT_CONFIG["MIN_AGE"] is not None
            ) or (
                "MAX_AGE" in DEFAULT_CONFIG and DEFAULT_CONFIG["MAX_AGE"] is not None
            ):
                age_mask = np.ones_like(test_ages, dtype=bool)

                if (
                    "MIN_AGE" in DEFAULT_CONFIG
                    and DEFAULT_CONFIG["MIN_AGE"] is not None
                ):
                    min_age = DEFAULT_CONFIG["MIN_AGE"]
                    age_mask &= test_ages >= min_age
                    print(f"   Filtering bootstrap ages >= {min_age}")

                if (
                    "MAX_AGE" in DEFAULT_CONFIG
                    and DEFAULT_CONFIG["MAX_AGE"] is not None
                ):
                    max_age = DEFAULT_CONFIG["MAX_AGE"]
                    age_mask &= test_ages <= max_age
                    print(f"   Filtering bootstrap ages <= {max_age}")

                # Apply age filter to ages
                test_ages = test_ages[age_mask]
                print(f"   After age filtering: {len(test_ages)} age points")

            # Now load all bootstrap values with the common mask
            for bootstrap_file in bootstrap_files:
                data = pd.read_csv(os.path.join(bootstrap_folder, bootstrap_file))

                # Determine the correct track name for the test code
                track_name = get_test_track_name(test_code)

                test_data = data[data["test_track"] == track_name]

                if len(test_data) > 0:
                    test_values = np.array(
                        [test_data.iloc[0][col] for col in age_columns]
                    )
                    # Ensure SD columns are properly aligned with age columns
                    if sd_columns:
                        # Create SD array with same length as age_columns
                        test_sds = np.zeros_like(test_values)
                        for i, age_col in enumerate(age_columns):
                            sd_col = age_col + "_sd"
                            if sd_col in test_data.columns:
                                test_sds[i] = test_data.iloc[0][sd_col]
                    else:
                        test_sds = np.zeros_like(test_values)

                    # Ensure n columns are properly aligned with age columns
                    if n_columns:
                        # Create n array with same length as age_columns
                        test_ns = np.zeros_like(test_values, dtype=int)
                        for i, age_col in enumerate(age_columns):
                            n_col = age_col + "_n"
                            if n_col in test_data.columns:
                                test_ns[i] = test_data.iloc[0][n_col]
                    else:
                        test_ns = np.ones_like(
                            test_values, dtype=int
                        )  # Default to 1 if no n data

                    # Apply the common mask
                    test_values = test_values[common_valid_mask]
                    test_sds = test_sds[common_valid_mask]
                    test_ns = test_ns[common_valid_mask]

                    # Apply age limits if specified
                    if (
                        "MIN_AGE" in DEFAULT_CONFIG
                        and DEFAULT_CONFIG["MIN_AGE"] is not None
                    ) or (
                        "MAX_AGE" in DEFAULT_CONFIG
                        and DEFAULT_CONFIG["MAX_AGE"] is not None
                    ):
                        # Create age mask for the original test_values array
                        age_mask = np.ones_like(test_values, dtype=bool)

                        # Get the original ages before filtering (same length as test_values)
                        original_ages = np.array([float(col) for col in age_columns])
                        original_ages = original_ages[
                            common_valid_mask
                        ]  # Apply common mask first

                        if (
                            "MIN_AGE" in DEFAULT_CONFIG
                            and DEFAULT_CONFIG["MIN_AGE"] is not None
                        ):
                            min_age = DEFAULT_CONFIG["MIN_AGE"]
                            age_mask &= original_ages >= min_age

                        if (
                            "MAX_AGE" in DEFAULT_CONFIG
                            and DEFAULT_CONFIG["MAX_AGE"] is not None
                        ):
                            max_age = DEFAULT_CONFIG["MAX_AGE"]
                            age_mask &= original_ages <= max_age

                        # Apply age filter to values, standard deviations, and sample sizes
                        test_values = test_values[age_mask]
                        test_sds = test_sds[age_mask]
                        test_ns = test_ns[age_mask]

                    all_bootstrap_values.append(test_values)
                    all_bootstrap_sds.append(test_sds)
                    all_bootstrap_ns.append(test_ns)

            if not all_bootstrap_values:
                print(
                    f"❌ ERROR: No valid bootstrap data found for test code {test_code}"
                )
                return None

            # Calculate mean and standard deviation across all bootstrap samples
            all_bootstrap_values_array = np.array(all_bootstrap_values)
            all_bootstrap_ns_array = np.array(all_bootstrap_ns)
            test_values = np.mean(all_bootstrap_values_array, axis=0)
            test_sds = np.std(all_bootstrap_values_array, axis=0)
            test_ns = np.mean(all_bootstrap_ns_array, axis=0).astype(
                int
            )  # Mean sample size, convert to int

            # Apply logit transformation if enabled
            if DEFAULT_CONFIG["USE_LOGIT_TRANSFORM"]:
                print("   🔄 Applying robust logit transformation to bootstrap data...")
                # Use robust binomial GLM approach
                logit_values, robust_weights, iteration_info = (
                    logit_transform_with_robust_weights(
                        test_values, test_ns, max_iterations=3
                    )
                )

                test_values = logit_values

                # Update test_sds to use robust weights instead of standard deviations
                robust_test_sds = (
                    robust_weights  # Store robust weights in test_sds for later use
                )
                test_sds = robust_test_sds

                print(
                    f"   📊 Bootstrap logit range: {test_values.min():.3f} to {test_values.max():.3f}"
                )
                print(f"   🔧 Robust weighting iterations: {len(iteration_info)}")
                for info in iteration_info:
                    print(
                        f"      Iteration {info['iteration']}: weight change={info['weight_change']:.4f}, mean weight={info['mean_weight']:.3f}"
                    )

            # Store which bootstrap file was used for plotting
            first_bootstrap_file = bootstrap_files[0]

            print(
                f"Found {test_name} data: {len(test_values)} data points from {len(all_bootstrap_values)} bootstrap runs"
            )
            print(f"Age range: {test_ages.min():.3f} to {test_ages.max():.3f}")
            print(
                f"{test_name} range: {test_values.min():.3f} to {test_values.max():.3f}"
            )
            print(
                f"Sample size range: {test_ns.min()} to {test_ns.max()} (mean: {test_ns.mean():.0f})"
            )

            # Store bootstrap data for later use
            bootstrap_data = {
                "ages": test_ages,
                "values": all_bootstrap_values,
                "sds": all_bootstrap_sds,
                "ns": all_bootstrap_ns,
                "files": bootstrap_files,
                "first_file": first_bootstrap_file,
            }
        else:
            # Original single file loading
            data = pd.read_csv(DEFAULT_CONFIG["INFILE"])
            test_data = data[data["test_code"] == test_code]

            if len(test_data) == 0:
                print(f"❌ ERROR: No data found for test code {test_code}")
                return None

            test_data = test_data.iloc[0]

            # Filter out _sd columns and get only age columns
            age_columns = [col for col in data.columns[9:] if not col.endswith("_sd")]
            sd_columns = [col for col in data.columns[9:] if col.endswith("_sd")]

            test_ages = np.array([float(col) for col in age_columns])
            test_values = np.array([test_data[col] for col in age_columns])

            # Get standard deviations for error bars
            if sd_columns:
                # Ensure SD columns are properly aligned with age columns
                test_sds = np.zeros_like(test_values)
                for i, age_col in enumerate(age_columns):
                    sd_col = age_col + "_sd"
                    if sd_col in data.columns:
                        test_sds[i] = test_data[sd_col]
            else:
                test_sds = np.zeros_like(test_values)

            # Get sample sizes for error bars
            n_columns = [col for col in data.columns[9:] if col.endswith("_n")]
            if n_columns:
                # Create n array with same length as age_columns
                test_ns = np.zeros_like(test_values, dtype=int)
                for i, age_col in enumerate(age_columns):
                    n_col = age_col + "_n"
                    if n_col in data.columns:
                        test_ns[i] = test_data[n_col]
            else:
                test_ns = np.ones_like(
                    test_values, dtype=int
                )  # Default to 1 if no n data

            valid_mask = ~np.isnan(test_values)
            test_ages = test_ages[valid_mask]
            test_values = test_values[valid_mask]
            test_sds = test_sds[valid_mask]  # Apply same mask to standard deviations
            test_ns = test_ns[valid_mask]  # Apply same mask to sample sizes

            # Apply logit transformation if enabled
            if DEFAULT_CONFIG["USE_LOGIT_TRANSFORM"]:
                print("   🔄 Applying robust logit transformation to data...")
                # Use robust binomial GLM approach
                logit_values, robust_weights, iteration_info = (
                    logit_transform_with_robust_weights(
                        test_values, test_ns, max_iterations=3
                    )
                )

                # Store original values for plotting
                original_test_values = test_values.copy()
                test_values = logit_values

                # Update test_sds to use robust weights instead of standard deviations
                # For binomial data, we'll calculate Wilson confidence intervals later
                robust_test_sds = (
                    robust_weights  # Store robust weights in test_sds for later use
                )

                print(
                    f"   📊 Original range: {original_test_values.min():.1f}% to {original_test_values.max():.1f}%"
                )
                print(
                    f"   📊 Logit range: {test_values.min():.3f} to {test_values.max():.3f}"
                )
                print(f"   🔧 Robust weighting iterations: {len(iteration_info)}")
                for info in iteration_info:
                    print(
                        f"      Iteration {info['iteration']}: weight change={info['weight_change']:.4f}, mean weight={info['mean_weight']:.3f}"
                    )

                # Store robust weights for later use in fitting
                test_sds = robust_test_sds

            # Apply age limits if specified
            if (
                "MIN_AGE" in DEFAULT_CONFIG and DEFAULT_CONFIG["MIN_AGE"] is not None
            ) or (
                "MAX_AGE" in DEFAULT_CONFIG and DEFAULT_CONFIG["MAX_AGE"] is not None
            ):
                age_mask = np.ones_like(test_ages, dtype=bool)

                if (
                    "MIN_AGE" in DEFAULT_CONFIG
                    and DEFAULT_CONFIG["MIN_AGE"] is not None
                ):
                    min_age = DEFAULT_CONFIG["MIN_AGE"]
                    age_mask &= test_ages >= min_age
                    print(f"   Filtering ages >= {min_age}")

                if (
                    "MAX_AGE" in DEFAULT_CONFIG
                    and DEFAULT_CONFIG["MAX_AGE"] is not None
                ):
                    max_age = DEFAULT_CONFIG["MAX_AGE"]
                    age_mask &= test_ages <= max_age
                    print(f"   Filtering ages <= {max_age}")

                # Apply age filter
                test_ages = test_ages[age_mask]
                test_values = test_values[age_mask]
                test_sds = test_sds[age_mask]

                print(f"   After age filtering: {len(test_values)} data points")

            print(f"Found {test_name} data: {len(test_values)} data points")
            print(f"Age range: {test_ages.min():.3f} to {test_ages.max():.3f}")
            print(
                f"{test_name} range: {test_values.min():.3f} to {test_values.max():.3f}"
            )

            bootstrap_data = None

        # 2. Load and process survival curve
        print("2. Loading survival curve...")
        survival_data_orig = pd.read_csv(DEFAULT_CONFIG["SURV_FILE"])
        print(f"Original survival data: {len(survival_data_orig)} points")

        survival_data = apply_gaussian_smoothing_continuous(
            survival_data_orig,
            bandwidth=DEFAULT_CONFIG["GAUSSIAN_BANDWIDTH"],
        )

        print(f"Smoothed survival data: {len(survival_data)} points")

        # 3. Fit linear menopause model for comparison
        print("3. Fitting linear menopause model for comparison...")
        linear_model = fit_linear_menopause_model(
            test_ages,
            test_values,
            survival_data,
            test_sds,
            use_iterative_reweighting=use_iterative_reweighting,
        )

        # Extract linear model R-squared for model selection (must be done early!)
        linear_r_squared = linear_model["r_squared"]

        # 3.5 Fit sigmoid menopause model for comparison
        print("3.5 Fitting sigmoid menopause model for comparison...")
        sigmoid_model = fit_sigmoid_menopause_model(
            test_ages,
            test_values,
            survival_data,
            linear_model,
            test_sds,
            test_ns,
            use_iterative_reweighting=use_iterative_reweighting,
        )

        # 3.6 Fit polynomial4 model for comparison
        print("3.6 Fitting 4th order polynomial model for comparison...")
        poly4_model = fit_polynomial4_model(
            test_ages,
            test_values,
            survival_data,
            linear_model,
            test_sds,
            test_ns,
            use_iterative_reweighting=use_iterative_reweighting,
        )

        # 3.7 Fit piecewise linear continuous model for comparison
        print("3.7 Fitting piecewise linear continuous model for comparison...")
        pwl_model = fit_piecewise_linear_continuous_model(
            test_ages,
            test_values,
            survival_data,
            linear_model,
            test_sds,
            test_ns,
            use_iterative_reweighting=use_iterative_reweighting,
        )

        # 3.8 Fit exp-linear model for comparison
        print("3.8 Fitting exp-linear model for comparison...")
        exp_lin_model = fit_exp_linear_model(
            test_ages,
            test_values,
            survival_data,
            linear_model,
            test_sds,
            test_ns,
            use_iterative_reweighting=use_iterative_reweighting,
        )

        # 3.9 Fit linear-exp model for comparison
        print("3.9 Fitting linear-exp model for comparison...")
        lin_exp_model = fit_linear_exp_model(
            test_ages,
            test_values,
            survival_data,
            linear_model,
            test_sds,
            test_ns,
            use_iterative_reweighting=use_iterative_reweighting,
        )

        # 3.10 Fit polynomial3 model for comparison
        print("3.10 Fitting 3rd order polynomial model for comparison...")
        poly3_model = fit_polynomial3_model(
            test_ages,
            test_values,
            survival_data,
            linear_model,
            test_sds,
            test_ns,
            use_iterative_reweighting=use_iterative_reweighting,
        )

        # 3.11 Fit piecewise linear 2-dots model for comparison
        print("3.11 Fitting piecewise linear 2-dots model for comparison...")
        pwl2_model = fit_piecewise_linear_2dots_model(
            test_ages,
            test_values,
            survival_data,
            linear_model,
            test_sds,
            test_ns,
            use_iterative_reweighting=use_iterative_reweighting,
        )

        # 3.12 Fit linear pre-menopause only model for comparison
        print("3.12 Fitting linear pre-menopause only model for comparison...")
        lin_pre_model = fit_linear_pre_only_model(
            test_ages,
            test_values,
            survival_data,
            linear_model,
            test_sds,
            test_ns,
            use_iterative_reweighting=use_iterative_reweighting,
        )

        # 3.13 Fit linear post-menopause only model for comparison
        print("3.13 Fitting linear post-menopause only model for comparison...")
        lin_post_model = fit_linear_post_only_model(
            test_ages,
            test_values,
            survival_data,
            linear_model,
            test_sds,
            test_ns,
            use_iterative_reweighting=use_iterative_reweighting,
        )

        # 4. Two-stage optimization approach
        print("4. Starting TWO-STAGE optimization...")

        # Stage 1: Single change point optimization
        print("   Stage 1: Single change point optimization")
        # Use same number of positions as two CP model (10 before + 10 after = 20 total)
        n_single_positions = DEFAULT_CONFIG[
            "N_CP_COMBINATIONS"
        ]  # Try same number as two CP combinations
        single_cp_positions = generate_single_cp_positions(
            n_single_positions, specific_cp_a=specific_cp_a
        )
        if DEFAULT_CONFIG.get("SINGLE_CP_BEFORE_ONLY", False):
            print(
                f"   Generated {len(single_cp_positions)} single change point positions (before menopause only)"
            )
        else:
            print(
                f"   Generated {len(single_cp_positions)} single change point positions (before and after menopause)"
            )

        # Use multiprocessing instead of Dask for LSF compatibility
        print(f"🚀 Using multiprocessing with {DEFAULT_CONFIG['N_WORKERS']} workers...")
        # Get the number of available CPUs from LSF environment
        n_cpus = int(os.environ.get("LSB_DJOB_NUMPROC", DEFAULT_CONFIG["N_WORKERS"]))
        print(f"   Available CPUs from LSF: {n_cpus}")

        # Balance CPU utilization with memory usage
        total_tasks = len(single_cp_positions) + DEFAULT_CONFIG["N_CP_COMBINATIONS"]
        if len(single_cp_positions) < n_cpus:
            print(
                f"   ⚠️  Warning: Only {len(single_cp_positions)} single CP tasks for {n_cpus} cores"
            )
            print(
                f"   💡 Consider increasing --n-combinations for better CPU utilization"
            )

        print(
            f"   🎯 Target CPU utilization: {n_cpus} cores with {total_tasks} total tasks"
        )
        print(
            f"   🧠 Memory-efficient processing: using 'fork' context and balanced task size"
        )

        try:
            # Skip optimization when doing bootstrap-only analysis or only basic models
            if not bootstrap_folder and not DEFAULT_CONFIG.get("ONLY_BASIC_MODELS", False):
                # Stage 1: Single change point optimization
                print("🔥 Creating single CP optimization tasks...")
                print(
                    f"⚡ Executing {len(single_cp_positions)} single CP tasks in parallel..."
                )
                start_time = time.time()

                # Force utilization of ALL available cores by ensuring enough tasks
                min_tasks_per_core = (
                    4  # Minimum 4 tasks per core to ensure full utilization
                )
                if len(single_cp_positions) < n_cpus * min_tasks_per_core:
                    additional_positions = generate_additional_single_cp_positions(
                        n_cpus * min_tasks_per_core - len(single_cp_positions)
                    )
                    single_cp_positions = np.concatenate(
                        [single_cp_positions, additional_positions]
                    )
                    print(
                        f"   📈 Increased single CP tasks to {len(single_cp_positions)} for better CPU utilization"
                    )

                    # Use ProcessPoolExecutor for parallel processing with aggressive settings
                # Use all available cores and fork context for better LSF compatibility
                with ProcessPoolExecutor(
                    max_workers=n_cpus, mp_context=multiprocessing.get_context("fork")
                ) as executor:
                    # Submit all tasks with immediate execution
                    # Submitting tasks to workers without print to reduce I/O contention
                    future_to_cp = {
                        executor.submit(
                            optimize_single_cp_position,
                            cp_pos,
                            test_ages,
                            test_values,
                            survival_data,
                            linear_model,
                            test_sds,
                            test_ns,
                            use_iterative_reweighting=use_iterative_reweighting,
                        ): cp_pos
                        for cp_pos in single_cp_positions
                    }

                    # Collect results as they complete
                    single_cp_results = []
                    for future in as_completed(future_to_cp):
                        try:
                            result = future.result()
                            single_cp_results.append(result)
                        except Exception as exc:
                            cp_pos = future_to_cp[future]
                            print(
                                f"❌ Task for CP position {cp_pos} generated an exception: {exc}"
                            )
                            raise

                single_cp_time = time.time() - start_time

                print(
                    f"✅ Single CP optimization completed in {single_cp_time:.2f} seconds!"
                )

                # Initialize variables for result processing
                best_single_cp_params = None
                best_single_cp_loss = np.inf
                best_single_cp_position = None
                single_cp_successful = 0

                # Skip result processing when doing bootstrap-only analysis
                if not bootstrap_folder:
                    # Find best single change point result
                    for result in single_cp_results:
                        if result["success"] and result["loss"] < best_single_cp_loss:
                            best_single_cp_params = result["params"]
                            best_single_cp_loss = result["loss"]
                            best_single_cp_position = result["cp_position"]
                            # Found new best single CP result (removed print to reduce I/O contention)
                        if result["success"]:
                            single_cp_successful += 1

                # Skip printing summary when doing bootstrap-only analysis
                if not bootstrap_folder:
                    print(f"📊 Single CP optimization summary:")
                    print(f"   Total positions: {len(single_cp_positions)}")
                    print(f"   Successful optimizations: {single_cp_successful}")
                    print(
                        f"   Success rate: {100*single_cp_successful/len(single_cp_positions):.1f}%"
                    )

            # Skip two CP optimization when doing bootstrap-only analysis or only basic models
            if not bootstrap_folder and not DEFAULT_CONFIG.get("ONLY_BASIC_MODELS", False):
                # Stage 2: Two change point optimization (always run both models)
                print("   Stage 2: Two change point optimization")
                # Use fixed seed for bootstrap consistency
                fixed_seed = (
                    42 if bootstrap_folder else None
                )  # Use seed 42 for bootstrap consistency
                combinations = generate_cp_combinations(
                    DEFAULT_CONFIG["N_CP_COMBINATIONS"],
                    fixed_seed,
                    specific_cp_a=specific_cp_a,
                    specific_cp_b=specific_cp_b,
                )
                print(f"   Generated {len(combinations)} two CP combinations")

                print(f"⚡ Executing {len(combinations)} two CP tasks in parallel...")
                start_time = time.time()

                # Force utilization of ALL available cores by ensuring enough tasks
                min_tasks_per_core = (
                    4  # Minimum 4 tasks per core to ensure full utilization
                )
                if len(combinations) < n_cpus * min_tasks_per_core:
                    # When using specific CPs, don't add additional combinations
                    if specific_cp_a is not None and specific_cp_b is not None:
                        additional_combinations = []
                    else:
                        additional_combinations = generate_cp_combinations(
                            n_cpus * min_tasks_per_core - len(combinations), fixed_seed
                        )
                    combinations.extend(additional_combinations)
                    print(
                        f"   📈 Increased two CP tasks to {len(combinations)} for better CPU utilization"
                    )

                # Use ProcessPoolExecutor for parallel processing with aggressive settings
                # Use all available cores and fork context for better LSF compatibility
                with ProcessPoolExecutor(
                    max_workers=n_cpus, mp_context=multiprocessing.get_context("fork")
                ) as executor:
                    # Submit all tasks with immediate execution
                    # Submitting two CP tasks to workers without print to reduce I/O contention
                    future_to_cp = {
                        executor.submit(
                            optimize_two_change_point,
                            cp1,
                            cp2,
                            test_ages,
                            test_values,
                            survival_data,
                            linear_model,
                            test_sds,
                            test_ns,
                            use_iterative_reweighting=use_iterative_reweighting,
                        ): (cp1, cp2)
                        for cp1, cp2 in combinations
                    }

                    # Collect results as they complete
                    two_cp_results = []
                    for future in as_completed(future_to_cp):
                        try:
                            result = future.result()
                            two_cp_results.append(result)
                        except Exception as exc:
                            cp1, cp2 = future_to_cp[future]
                            print(
                                f"❌ Task for CP positions ({cp1}, {cp2}) generated an exception: {exc}"
                            )
                            raise

                two_cp_time = time.time() - start_time

                print(f"✅ Two CP optimization completed in {two_cp_time:.2f} seconds!")

                # Initialize variables for result processing
                best_two_cp_params = None
                best_two_cp_loss = np.inf
                best_two_cp_cp1, best_two_cp_cp2 = None, None
                two_cp_successful = 0

                # Skip two CP result processing when doing bootstrap-only analysis
                if not bootstrap_folder:
                    # Find best two change point result
                    for result in two_cp_results:
                        if result["success"] and result["loss"] < best_two_cp_loss:
                            best_two_cp_params = result["params"]
                            best_two_cp_loss = result["loss"]
                            best_two_cp_cp1, best_two_cp_cp2 = (
                                result["cp1"],
                                result["cp2"],
                            )
                            # Found new best two CP result (removed print to reduce I/O contention)
                        if result["success"]:
                            two_cp_successful += 1

                    print(f"📊 Two CP optimization summary:")
                    print(f"   Total combinations: {len(combinations)}")
                    print(f"   Successful optimizations: {two_cp_successful}")
                    print(
                        f"   Success rate: {100*two_cp_successful/len(combinations):.1f}%"
                    )

            # Fit 3 CP model (only for non-bootstrap analysis, skip if only basic models)
            best_three_cp_params = None
            best_three_cp_loss = np.inf
            best_three_cp_pure_rss = None
            if not bootstrap_folder and not DEFAULT_CONFIG.get("ONLY_BASIC_MODELS", False):
                print("🔬 FITTING THREE CP MODEL...")
                best_three_cp_params, best_three_cp_loss, best_three_cp_pure_rss = fit_three_cp_optimization(
                    test_ages,
                    test_values,
                    survival_data,
                    linear_model,
                    test_sds=test_sds,
                    test_ns=test_ns,
                    use_iterative_reweighting=use_iterative_reweighting,
                )
                if best_three_cp_params is not None:
                    slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, val_cp3, slope2, opt_cp1, opt_cp2, opt_cp3 = best_three_cp_params
                    sorted_cps = sorted([opt_cp1, opt_cp2, opt_cp3])
                    print(f"✅ Three CP model fitted: CPs at {sorted_cps[0]:.2f}, {sorted_cps[1]:.2f}, {sorted_cps[2]:.2f}")
                else:
                    print("❌ Three CP model fitting failed")

            # Fit spline model (only for non-bootstrap analysis, skip if only basic models)
            best_spline_params = None
            best_spline_loss = np.inf
            best_spline_pure_rss = None
            if not bootstrap_folder and not DEFAULT_CONFIG.get("ONLY_BASIC_MODELS", False):
                print("🔬 FITTING SPLINE MODEL...")
                best_spline_params, best_spline_loss, best_spline_pure_rss = fit_spline_optimization(
                    test_ages,
                    test_values,
                    survival_data,
                    linear_model,
                    test_sds=test_sds,
                    test_ns=test_ns,
                    use_iterative_reweighting=use_iterative_reweighting,
                    n_coeffs_pre=4,  # 4 coefficients for pre-menopause (cubic polynomial)
                    n_coeffs_post=4,  # 4 coefficients for post-menopause (cubic polynomial)
                )
                if best_spline_params is not None:
                    n_coeffs_pre = 4
                    n_coeffs_post = 4
                    pre_coeffs = best_spline_params[:n_coeffs_pre]
                    post_coeffs = best_spline_params[n_coeffs_pre:n_coeffs_pre + n_coeffs_post]
                    jump = best_spline_params[-1]
                    print(f"✅ Spline model fitted: Jump={jump:.3f}, Pre coeffs={pre_coeffs[:2]}, Post coeffs={post_coeffs[:2]}")
                else:
                    print("❌ Spline model fitting failed")

            # Initialize variables for model selection
            best_params = None
            best_loss = np.inf
            best_cp1 = 0
            best_cp2 = 0
            successful_results = 0
            r_squared = 0
            predictions = test_values  # dummy predictions

            # Simple model selection: Choose between single CP, two CP, three CP, and spline using BIC
            if not bootstrap_folder and best_single_cp_params is not None:
                print(f"📊 MODEL COMPARISON:")
                print(f"   Single CP Loss: {best_single_cp_loss:.6f}")
                print(f"   Two CP Loss: {best_two_cp_loss:.6f}")
                if best_three_cp_pure_rss is not None:
                    print(f"   Three CP Loss: {best_three_cp_pure_rss:.6f}")
                if best_spline_pure_rss is not None:
                    print(f"   Spline Loss: {best_spline_pure_rss:.6f}")

                # Calculate BIC for all models with stronger penalty for complexity
                n_params_single = 6  # slope1, val_cp, val_before_jump, val_after_jump, slope2, cp_position
                n_params_two = 8  # slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, slope2, cp1_pos, cp2_pos
                n_params_three = 10  # slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, val_cp3, slope2, cp1_pos, cp2_pos, cp3_pos
                n_params_spline = 9  # 4 pre_coeffs + 4 post_coeffs + 1 jump
                n_data = len(test_values)

                # Use the pure RSS for BIC calculation
                pure_rss_single = best_single_cp_loss
                pure_rss_two = best_two_cp_loss
                pure_rss_three = best_three_cp_pure_rss if best_three_cp_pure_rss is not None else np.inf
                pure_rss_spline = best_spline_pure_rss if best_spline_pure_rss is not None else np.inf

                # Calculate BIC for all models using pure RSS (no penalties)
                bic_single = (
                    n_data * np.log(pure_rss_single / n_data)
                    + np.log(n_data) * n_params_single
                )
                bic_two = (
                    n_data * np.log(pure_rss_two / n_data)
                    + np.log(n_data) * n_params_two
                )
                bic_three = (
                    n_data * np.log(pure_rss_three / n_data)
                    + np.log(n_data) * n_params_three
                ) if best_three_cp_pure_rss is not None else np.inf
                bic_spline = (
                    n_data * np.log(pure_rss_spline / n_data)
                    + np.log(n_data) * n_params_spline
                ) if best_spline_pure_rss is not None else np.inf

                # Calculate likelihood ratio test statistic using pure RSS
                # LRT = -2 * log(L_simple / L_complex) = -2 * log(RSS_simple / RSS_complex)
                # Under null hypothesis (simple model is true), LRT ~ chi-squared(df_complex - df_simple)
                lrt_statistic = n_data * np.log(pure_rss_single / pure_rss_two)
                df_diff = n_params_two - n_params_single  # 6 - 5 = 1

                # Critical value for chi-squared with 1 df at 5% significance
                critical_value = 3.84  # chi-squared(1, 0.05)

                print(f"   Likelihood Ratio Test:")
                print(f"      LRT statistic: {lrt_statistic:.2f}")
                print(f"      Critical value (5%): {critical_value:.2f}")
                print(
                    f"      Significant improvement: {lrt_statistic > critical_value}"
                )

                print(f"   Single CP BIC: {bic_single:.2f}")
                print(f"   Two CP BIC: {bic_two:.2f}")
                if best_three_cp_pure_rss is not None:
                    print(f"   Three CP BIC: {bic_three:.2f}")
                if best_spline_pure_rss is not None:
                    print(f"   Spline BIC: {bic_spline:.2f}")
                
                # Calculate BIC for new basic models
                if poly4_model is not None:
                    poly4_n_params = 5
                    poly4_rss = poly4_model["rss"]
                    bic_poly4 = (
                        n_data * np.log(poly4_rss / n_data)
                        + np.log(n_data) * poly4_n_params
                    )
                    print(f"   Polynomial4 BIC: {bic_poly4:.2f}")
                else:
                    bic_poly4 = np.inf
                
                if pwl_model is not None:
                    pwl_n_params = 5
                    pwl_rss = pwl_model["rss"]
                    bic_pwl = (
                        n_data * np.log(pwl_rss / n_data)
                        + np.log(n_data) * pwl_n_params
                    )
                    print(f"   Piecewise Linear Continuous BIC: {bic_pwl:.2f}")
                else:
                    bic_pwl = np.inf
                
                if exp_lin_model is not None:
                    exp_lin_n_params = 5
                    exp_lin_rss = exp_lin_model["rss"]
                    bic_exp_lin = (
                        n_data * np.log(exp_lin_rss / n_data)
                        + np.log(n_data) * exp_lin_n_params
                    )
                    print(f"   Exp-Linear BIC: {bic_exp_lin:.2f}")
                else:
                    bic_exp_lin = np.inf
                
                if lin_exp_model is not None:
                    lin_exp_n_params = 5
                    lin_exp_rss = lin_exp_model["rss"]
                    bic_lin_exp = (
                        n_data * np.log(lin_exp_rss / n_data)
                        + np.log(n_data) * lin_exp_n_params
                    )
                    print(f"   Linear-Exp BIC: {bic_lin_exp:.2f}")
                else:
                    bic_lin_exp = np.inf
                
                if poly3_model is not None:
                    poly3_n_params = 4
                    poly3_rss = poly3_model["rss"]
                    bic_poly3 = (
                        n_data * np.log(poly3_rss / n_data)
                        + np.log(n_data) * poly3_n_params
                    )
                    print(f"   Polynomial3 BIC: {bic_poly3:.2f}")
                else:
                    bic_poly3 = np.inf
                
                if pwl2_model is not None:
                    pwl2_n_params = 4
                    pwl2_rss = pwl2_model["rss"]
                    bic_pwl2 = (
                        n_data * np.log(pwl2_rss / n_data)
                        + np.log(n_data) * pwl2_n_params
                    )
                    print(f"   Piecewise Linear 2-dots BIC: {bic_pwl2:.2f}")
                else:
                    bic_pwl2 = np.inf
                
                if lin_pre_model is not None:
                    lin_pre_n_params = 3
                    lin_pre_rss = lin_pre_model["rss"]
                    bic_lin_pre = (
                        n_data * np.log(lin_pre_rss / n_data)
                        + np.log(n_data) * lin_pre_n_params
                    )
                    print(f"   Linear Pre-only BIC: {bic_lin_pre:.2f}")
                else:
                    bic_lin_pre = np.inf
                
                if lin_post_model is not None:
                    lin_post_n_params = 3
                    lin_post_rss = lin_post_model["rss"]
                    bic_lin_post = (
                        n_data * np.log(lin_post_rss / n_data)
                        + np.log(n_data) * lin_post_n_params
                    )
                    print(f"   Linear Post-only BIC: {bic_lin_post:.2f}")
                else:
                    bic_lin_post = np.inf
                
                # Find best model by BIC (including new basic models)
                bic_values = {
                    "single_cp": bic_single,
                    "two_cp": bic_two,
                    "poly3": bic_poly3,
                    "poly4": bic_poly4,
                    "pwl": bic_pwl,
                    "pwl2": bic_pwl2,
                    "exp_lin": bic_exp_lin,
                    "lin_exp": bic_lin_exp,
                    "lin_pre": bic_lin_pre,
                    "lin_post": bic_lin_post,
                }
                if best_three_cp_pure_rss is not None:
                    bic_values["three_cp"] = bic_three
                if best_spline_pure_rss is not None:
                    bic_values["spline"] = bic_spline
                
                best_model_type = min(bic_values, key=bic_values.get)
                print(f"   🏆 Best model by BIC: {best_model_type} (BIC={bic_values[best_model_type]:.2f})")

                # Calculate R-squared for both models to use in model selection
                # Calculate single CP R-squared
                if best_single_cp_params is not None:
                    (
                        single_cp_slope1,
                        single_cp_val_cp,
                        single_cp_val_before_jump,
                        single_cp_val_after_jump,
                        single_cp_slope2,
                        single_cp_opt_position,
                    ) = best_single_cp_params
                    single_cp_predictions = []
                    for age in test_ages:
                        pred = pred_age_single_cp(
                            age,
                            survival_data,
                            single_cp_slope1,
                            single_cp_val_cp,
                            single_cp_val_before_jump,
                            single_cp_val_after_jump,
                            single_cp_slope2,
                            single_cp_opt_position,  # Use optimized position from params
                        )
                        single_cp_predictions.append(pred)
                    single_cp_predictions = np.array(single_cp_predictions)

                    single_cp_residuals = test_values - single_cp_predictions
                    if (
                        DEFAULT_CONFIG["USE_CV"]
                        and test_sds is not None
                        and np.any(test_sds > 0)
                    ):
                        weights = 1.0 / (test_sds**2)
                        single_cp_rss = np.sum(weights * single_cp_residuals**2)
                        single_cp_tss = np.sum(
                            weights
                            * (test_values - np.average(test_values, weights=weights))
                            ** 2
                        )
                    else:
                        single_cp_rss = np.sum(single_cp_residuals**2)
                        single_cp_tss = np.sum(
                            (test_values - np.mean(test_values)) ** 2
                        )
                    single_cp_r_squared = (
                        1 - single_cp_rss / single_cp_tss if single_cp_tss > 0 else 0
                    )
                else:
                    single_cp_r_squared = 0

                # Calculate two CP R-squared
                if best_two_cp_params is not None:
                    (
                        two_cp_slope1,
                        two_cp_val_cp1,
                        two_cp_val_before_jump,
                        two_cp_val_after_jump,
                        two_cp_val_cp2,
                        two_cp_slope2,
                        two_cp_opt_cp1,
                        two_cp_opt_cp2,
                    ) = best_two_cp_params
                    two_cp_predictions = []
                    for age in test_ages:
                        pred = pred_age(
                            age,
                            survival_data,
                            two_cp_slope1,
                            two_cp_val_cp1,
                            two_cp_val_before_jump,
                            two_cp_val_after_jump,
                            two_cp_val_cp2,
                            two_cp_slope2,
                            two_cp_opt_cp1,  # Use optimized positions from params
                            two_cp_opt_cp2,
                        )
                        two_cp_predictions.append(pred)
                    two_cp_predictions = np.array(two_cp_predictions)

                    two_cp_residuals = test_values - two_cp_predictions
                    if (
                        DEFAULT_CONFIG["USE_CV"]
                        and test_sds is not None
                        and np.any(test_sds > 0)
                    ):
                        weights = 1.0 / (test_sds**2)
                        two_cp_rss = np.sum(weights * two_cp_residuals**2)
                        two_cp_tss = np.sum(
                            weights
                            * (test_values - np.average(test_values, weights=weights))
                            ** 2
                        )
                    else:
                        two_cp_rss = np.sum(two_cp_residuals**2)
                        two_cp_tss = np.sum((test_values - np.mean(test_values)) ** 2)
                    two_cp_r_squared = (
                        1 - two_cp_rss / two_cp_tss if two_cp_tss > 0 else 0
                    )
                else:
                    two_cp_r_squared = 0

                # Skip model selection when doing bootstrap-only analysis
                if not bootstrap_folder:
                    # Choose model based on likelihood ratio test AND biological plausibility

                    # Enhanced model selection: Choose between single CP and two CP using BIC and biological plausibility
                    # Apply biological veto for single CP models when change point is before menopause

                    # Check if single CP model has change point before menopause
                    single_cp_before_menopause = best_single_cp_position < 0

                    # Calculate improvement for single CP model
                    single_cp_improvement = single_cp_r_squared - linear_r_squared

                    # For single CP models with change point before menopause, be more conservative
                    # They often perform worse than linear models due to slope limitations
                    if single_cp_before_menopause and single_cp_improvement < 0:
                        # Single CP model performs worse than linear - prefer two CP or linear
                        if lrt_statistic > critical_value:
                            print(
                                f"   🏆 SELECTING TWO CP MODEL (LRT significant: {lrt_statistic:.2f} > {critical_value:.2f})"
                            )
                            print(
                                f"   💡 Single CP model rejected: performs worse than linear (improvement={single_cp_improvement:.6f})"
                            )
                            best_params = best_two_cp_params
                            best_loss = best_two_cp_loss
                            best_cp1, best_cp2 = best_two_cp_cp1, best_two_cp_cp2
                            successful_results = two_cp_successful
                            r_squared = two_cp_r_squared  # Use two CP R-squared
                            predictions = two_cp_predictions  # Use two CP predictions
                        else:
                            print(
                                f"   ⚠️  SELECTING SINGLE CP MODEL despite poor performance (no significant two CP)"
                            )
                            print(
                                f"   💡 Single CP improvement: {single_cp_improvement:.6f} (worse than linear)"
                            )
                            best_params = best_single_cp_params
                            best_loss = best_single_cp_loss
                            best_cp1, best_cp2 = (
                                best_single_cp_position,
                                best_single_cp_position,
                            )
                            successful_results = single_cp_successful
                            r_squared = single_cp_r_squared  # Use single CP R-squared
                            predictions = (
                                single_cp_predictions  # Use single CP predictions
                            )
                    else:
                        # Standard selection based on LRT
                        if lrt_statistic > critical_value:
                            print(
                                f"   🏆 SELECTING TWO CP MODEL (LRT significant: {lrt_statistic:.2f} > {critical_value:.2f})"
                            )
                            best_params = best_two_cp_params
                            best_loss = best_two_cp_loss
                            best_cp1, best_cp2 = best_two_cp_cp1, best_two_cp_cp2
                            successful_results = two_cp_successful
                            r_squared = two_cp_r_squared  # Use two CP R-squared
                            predictions = two_cp_predictions  # Use two CP predictions
                        else:
                            print(
                                f"   🏆 SELECTING SINGLE CP MODEL (LRT not significant)"
                            )
                            if single_cp_before_menopause:
                                print(
                                    f"   💡 Single CP change point before menopause - limited slope modeling"
                                )
                            best_params = best_single_cp_params
                            best_loss = best_single_cp_loss
                            best_cp1, best_cp2 = (
                                best_single_cp_position,
                                best_single_cp_position,
                            )
                        successful_results = single_cp_successful
                        r_squared = single_cp_r_squared  # Use single CP R-squared
                        predictions = single_cp_predictions  # Use single CP predictions
            elif not bootstrap_folder and best_two_cp_params is not None:
                print("   ❌ No successful single CP optimization, using two CP result")
                best_params = best_two_cp_params
                best_loss = best_two_cp_loss
                best_cp1, best_cp2 = best_two_cp_cp1, best_two_cp_cp2

        finally:
            # Cleanup completed
            pass

        # Set default values for bootstrap-only mode
        if bootstrap_folder and best_params is None:
            # For bootstrap mode, we don't need real optimization results
            # Just set dummy values to allow bootstrap analysis to proceed
            # Use 6 parameters for single CP: [slope1, val_cp, val_before_jump, val_after_jump, slope2, cp_position]
            best_params = [0, 0, 0, 0, 0, 0]  # dummy parameters (6 for single CP)
            best_loss = 0
            best_cp1 = 0
            best_cp2 = 0
            successful_results = 0
            r_squared = 0
            predictions = test_values  # dummy predictions

        if best_params is None:
            print("🚨 ERROR: No valid staged model found!")
            return None

        print("5. Computing final predictions...")

        # Determine if we're using single or two change point model
        if best_cp1 == best_cp2:  # Single change point model
            slope1, val_cp, val_before_jump, val_after_jump, slope2, opt_cp_position = best_params
            predictions = []
            for age in test_ages:
                pred = pred_age_single_cp(
                    age,
                    survival_data,
                    slope1,
                    val_cp,
                    val_before_jump,
                    val_after_jump,
                    slope2,
                    opt_cp_position,  # Use optimized cp_position
                )
                predictions.append(pred)
            predictions = np.array(predictions)
        else:  # Two change point model
            slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, slope2, opt_cp1, opt_cp2 = (
                best_params
            )
            predictions = get_final_predictions_for_change_point(
                test_ages,
                survival_data,
                slope1,
                val_cp1,
                val_before_jump,
                val_after_jump,
                val_cp2,
                slope2,
                opt_cp1,  # Use optimized positions from params
                opt_cp2,
            )

        # Assess model quality with statistical significance testing
        n_params_cp = (
            8 if best_cp1 != best_cp2 else 6
        )  # Two CP (8) vs Single CP (6) - includes change point positions
        n_params_linear = 4  # Linear menopause model has 4 parameters (a, b, c, d)

        quality_metrics = assess_model_quality(
            test_values,
            predictions,
            linear_model["fitted_values"],
            n_params_cp,
            n_params_linear,
        )

        print(f"\n🎉 FINAL RESULTS FOR {test_name}:")
        print(f"Staged R-squared: {r_squared:.9f}")
        print(f"Linear R-squared: {linear_model['r_squared']:.9f}")
        print(f"Improvement: {r_squared - linear_model['r_squared']:.9f}")
        print(f"\n📊 MODEL QUALITY ASSESSMENT:")
        print(f"F-statistic: {quality_metrics['f_statistic']:.4f}")
        print(f"P-value: {quality_metrics['p_value']:.6f}")
        print(f"Is significant: {quality_metrics['is_significant']}")
        print(f"AIC improvement: {quality_metrics['aic_improvement']:.2f}")
        print(f"BIC improvement: {quality_metrics['bic_improvement']:.2f}")
        print(f"Model is reasonable: {quality_metrics['is_reasonable']}")

        # Initialize result variables for both modes
        single_cp_predictions = None
        single_cp_r_squared = None
        single_cp_aic = None
        single_cp_bic = None
        two_cp_predictions = None
        two_cp_r_squared = None
        two_cp_aic = None
        two_cp_bic = None
        final_predictions = test_values  # dummy
        final_r_squared = 0
        single_cp_params = None
        linear_aic = 0
        linear_bic = 0

        # Initialize quality metrics for bootstrap mode
        quality_metrics = {
            "f_statistic": 0,
            "p_value": 1,
            "is_significant": False,
            "meets_min_improvement": False,
            "meets_min_r_squared": False,
            "is_reasonable": False,
            "aic": 0,
            "linear_aic": 0,
            "aic_improvement": 0,
            "bic": 0,
            "linear_bic": 0,
            "bic_improvement": 0,
            "residuals_std": 0,
            "linear_residuals_std": 0,
        }

        # Skip result processing when doing bootstrap-only analysis
        if not bootstrap_folder:
            # Calculate predictions, R², AIC, and BIC for ALL THREE models consistently
            n_data = len(test_values)

            # 1. Linear model (already calculated)
            linear_predictions = linear_model["fitted_values"]
            linear_residuals = test_values - linear_predictions
            linear_n_params = 4  # a, b, c, d
            linear_r_squared = linear_model[
                "r_squared"
            ]  # Use R² from model (correctly weighted)

            # Calculate linear RSS consistently with change point models
            # Use the same approach as change point models: calculate BIC using integral-based RSS
            # This ensures exactly the same calculation method as CP models
            linear_rss_for_bic = calculate_pure_rss(
                test_ages,
                test_values,
                survival_data,
                linear_a=linear_model["a"],
                linear_b=linear_model["b"],
                linear_c=linear_model["c"],
                linear_d=linear_model["d"],
                test_sds=test_sds,
                test_ns=None,  # Don't transform - use test_sds directly
            )

            # Calculate BIC using the same approach as change point models
            # Use n_data for both AIC and BIC (not n_eff) to ensure exact consistency
            linear_aic = (
                n_data * np.log(linear_rss_for_bic / n_data) + 2 * linear_n_params
            )
            linear_bic = (
                n_data * np.log(linear_rss_for_bic / n_data)
                + np.log(n_data) * linear_n_params
            )

            # 2. Single CP model - use values calculated earlier
            if best_single_cp_params is not None:
                single_cp_predictions = (
                    single_cp_predictions  # Already calculated above
                )
                single_cp_r_squared = single_cp_r_squared  # Already calculated above

                # Calculate AIC and BIC for single CP model
                single_cp_n_params = 6  # slope1, val_cp1, val_before_jump, val_after_jump, slope2, cp_position (search adds complexity)
                single_cp_residuals = test_values - single_cp_predictions
                if test_sds is not None and np.any(test_sds > 0):
                    weights = 1.0 / (test_sds**2)
                    single_cp_rss = np.sum(weights * single_cp_residuals**2)
                    # Calculate effective sample size (only used for information, not for BIC)
                    n_eff = (np.sum(weights)) ** 2 / np.sum(weights**2)

                    # Always use n_data for BIC calculation to ensure consistency
                    single_cp_aic = (
                        n_data * np.log(single_cp_rss / n_data) + 2 * single_cp_n_params
                    )
                    single_cp_bic = (
                        n_data * np.log(single_cp_rss / n_data)
                        + np.log(n_data) * single_cp_n_params
                    )
                else:
                    single_cp_rss = np.sum(single_cp_residuals**2)
                    single_cp_aic = (
                        n_data * np.log(single_cp_rss / n_data) + 2 * single_cp_n_params
                    )
                    single_cp_bic = (
                        n_data * np.log(single_cp_rss / n_data)
                        + np.log(n_data) * single_cp_n_params
                    )
            else:
                single_cp_predictions = None
                single_cp_r_squared = None
                single_cp_aic = None
                single_cp_bic = None

            # 3. Two CP model - use values calculated earlier
            if best_two_cp_params is not None:
                two_cp_predictions = two_cp_predictions  # Already calculated above
                two_cp_r_squared = two_cp_r_squared  # Already calculated above

                # Calculate AIC and BIC for two CP model
                two_cp_n_params = 8  # slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, slope2, cp1_pos, cp2_pos (search adds complexity)
                two_cp_residuals = test_values - two_cp_predictions
                if test_sds is not None and np.any(test_sds > 0):
                    weights = 1.0 / (test_sds**2)
                    two_cp_rss = np.sum(weights * two_cp_residuals**2)
                    # Calculate effective sample size (only used for information, not for BIC)
                    n_eff = (np.sum(weights)) ** 2 / np.sum(weights**2)

                    # Always use n_data for BIC calculation to ensure consistency
                    two_cp_aic = (
                        n_data * np.log(two_cp_rss / n_data) + 2 * two_cp_n_params
                    )
                    two_cp_bic = (
                        n_data * np.log(two_cp_rss / n_data)
                        + np.log(n_data) * two_cp_n_params
                    )
                else:
                    two_cp_rss = np.sum(two_cp_residuals**2)
                    two_cp_aic = (
                        n_data * np.log(two_cp_rss / n_data) + 2 * two_cp_n_params
                    )
                    two_cp_bic = (
                        n_data * np.log(two_cp_rss / n_data)
                        + np.log(n_data) * two_cp_n_params
                    )
            else:
                two_cp_predictions = None
                two_cp_r_squared = None
                two_cp_aic = None
                two_cp_bic = None

            # Extract parameters based on model type for compatibility
            if best_cp1 == best_cp2:  # Single change point model selected
                slope1, val_cp, val_before_jump, val_after_jump, slope2, _ = best_params
                # Use the selected model's predictions and R²
                final_predictions = (
                    single_cp_predictions
                    if single_cp_predictions is not None
                    else predictions
                )
                final_r_squared = (
                    single_cp_r_squared
                    if single_cp_r_squared is not None
                    else r_squared
                )
                # Set single_cp_params for plotting
                single_cp_params = best_params
            else:  # Two change point model selected
                slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, slope2, _, _ = (
                    best_params
                )
                # Use the selected model's predictions and R²
                final_predictions = (
                    two_cp_predictions
                    if two_cp_predictions is not None
                    else predictions
                )
                final_r_squared = (
                    two_cp_r_squared if two_cp_r_squared is not None else r_squared
                )
                # Set single_cp_params for plotting (use best single CP if available)
                single_cp_params = (
                    best_single_cp_params if best_single_cp_params is not None else None
                )

            # 6. Save detailed results
            print("6. Saving results...")

            # Save fitted values CSV with ALL THREE models
            csv_data = {
                "age": test_ages,
                "actual": test_values,
                "predicted_linear": linear_predictions,
                "residual_linear": test_values - linear_predictions,
            }

            # Add single CP model predictions and residuals
            if single_cp_predictions is not None:
                csv_data["predicted_single_cp"] = single_cp_predictions
                csv_data["residual_single_cp"] = test_values - single_cp_predictions
            else:
                # Fill with NaN if not available
                csv_data["predicted_single_cp"] = np.full_like(test_values, np.nan)
                csv_data["residual_single_cp"] = np.full_like(test_values, np.nan)

            # Add two CP model predictions and residuals
            if two_cp_predictions is not None:
                csv_data["predicted_two_cp"] = two_cp_predictions
                csv_data["residual_two_cp"] = test_values - two_cp_predictions
            else:
                # Fill with NaN if not available
                csv_data["predicted_two_cp"] = np.full_like(test_values, np.nan)
                csv_data["residual_two_cp"] = np.full_like(test_values, np.nan)

            # Add the selected model for backward compatibility
            csv_data["predicted_staged"] = final_predictions
            csv_data["residual_staged"] = test_values - final_predictions

            results_df = pd.DataFrame(csv_data)

        # Initialize variables for bootstrap mode
        if bootstrap_folder:
            fitted_values_path = None
            cp_params_path = None
            plot_path = None
            results_df = None

        # Skip saving fitted values when doing bootstrap-only analysis
        if not bootstrap_folder:
            fitted_values_path = os.path.join(
                test_dir, f"{test_name}_fitted_values.csv"
            )
            results_df.to_csv(fitted_values_path, index=False)
        else:
            fitted_values_path = None

        # Save change point parameters CSV
        # Calculate the actual slopes from change points to menopause for penalty analysis
        if best_cp1 == best_cp2:  # Single CP model
            slope_cp1_to_meno = (
                (val_before_jump - val_cp) / (0 - best_cp1) if best_cp1 != 0 else 0
            )
            slope_meno_to_cp2 = 0  # No second change point
            val_cp1_value = val_cp
            val_cp2_value = val_cp  # Same as cp1 for single CP
        else:  # Two CP model
            slope_cp1_to_meno = (
                (val_before_jump - val_cp1) / (0 - best_cp1) if best_cp1 != 0 else 0
            )
            slope_meno_to_cp2 = (
                (val_cp2 - val_after_jump) / (best_cp2 - 0) if best_cp2 != 0 else 0
            )
            val_cp1_value = val_cp1
            val_cp2_value = val_cp2

        cp_params_df = pd.DataFrame(
            {
                "parameter": [
                    "cp1_position",
                    "cp2_position",
                    "slope1",
                    "slope1_cp_to_meno",
                    "val_cp1",
                    "val_before_jump",
                    "val_after_jump",
                    "val_cp2",
                    "slope2",
                    "slope2_meno_to_cp",
                    "r_squared_linear",
                    "r_squared_single_cp",
                    "r_squared_two_cp",
                    "r_squared_staged",
                    "improvement",
                    "jump_magnitude",
                    "aic_linear",
                    "aic_single_cp",
                    "aic_two_cp",
                    "bic_linear",
                    "bic_single_cp",
                    "bic_two_cp",
                    "aic_improvement_single_cp",
                    "aic_improvement_two_cp",
                    "bic_improvement_single_cp",
                    "bic_improvement_two_cp",
                ],
                "value": [
                    best_cp1,
                    best_cp2,
                    slope1,
                    slope_cp1_to_meno,
                    val_cp1_value,
                    val_before_jump,
                    val_after_jump,
                    val_cp2_value,
                    slope2,
                    slope_meno_to_cp2,
                    linear_r_squared,
                    single_cp_r_squared if single_cp_r_squared is not None else np.nan,
                    two_cp_r_squared if two_cp_r_squared is not None else np.nan,
                    final_r_squared,
                    final_r_squared - linear_r_squared,
                    val_after_jump - val_before_jump,
                    linear_aic,
                    single_cp_aic if single_cp_aic is not None else np.nan,
                    two_cp_aic if two_cp_aic is not None else np.nan,
                    linear_bic,
                    single_cp_bic if single_cp_bic is not None else np.nan,
                    two_cp_bic if two_cp_bic is not None else np.nan,
                    (
                        (linear_aic - single_cp_aic)
                        if single_cp_aic is not None
                        else np.nan
                    ),
                    (linear_aic - two_cp_aic) if two_cp_aic is not None else np.nan,
                    (
                        (linear_bic - single_cp_bic)
                        if single_cp_bic is not None
                        else np.nan
                    ),
                    (linear_bic - two_cp_bic) if two_cp_bic is not None else np.nan,
                ],
            }
        )
        # Skip saving parameters when doing bootstrap-only analysis
        if not bootstrap_folder:
            cp_params_path = os.path.join(
                test_dir, f"{test_name}_change_point_parameters.csv"
            )
            cp_params_df.to_csv(cp_params_path, index=False)
        else:
            cp_params_path = None

        # 7. Create and save plots
        print("7. Generating plots...")

        # Get the best two CP model parameters (always available)
        # Note: best_two_cp_params was calculated earlier and should be preserved

        # Add bootstrap info to plot title if using bootstrap data
        plot_title = test_name
        if bootstrap_data:
            plot_title = f"{test_name} (Bootstrap: {bootstrap_data['first_file']})"

        # Skip plotting and file creation when doing bootstrap-only analysis
        if not bootstrap_folder:
            plot_path = create_plots(
                test_ages,
                test_values,
                predictions,
                linear_model,
                slope1,
                val_cp1_value,
                val_before_jump,
                val_after_jump,
                val_cp2_value,
                slope2,
                best_cp1,
                best_cp2,
                r_squared,
                plot_title,  # Use modified title
                test_dir,
                test_sds,
                (
                    test_ns if "test_ns" in locals() else None
                ),  # Add sample sizes for SE calculation
                single_cp_params,
                single_cp_predictions,
                single_cp_r_squared,
                (
                    best_single_cp_position
                    if "best_single_cp_position" in locals()
                    else None
                ),
                best_two_cp_params,  # Add best two CP model parameters
                survival_data,  # Add survival_data parameter
                best_two_cp_cp1,  # Add best two CP model change point positions
                best_two_cp_cp2,
                data_is_logit_transformed=DEFAULT_CONFIG["USE_LOGIT_TRANSFORM"],
                sigmoid_model=sigmoid_model if "sigmoid_model" in locals() else None,  # Add sigmoid model
                poly3_model=poly3_model if "poly3_model" in locals() else None,  # Add polynomial3 model
                poly4_model=poly4_model if "poly4_model" in locals() else None,  # Add polynomial4 model
                pwl_model=pwl_model if "pwl_model" in locals() else None,  # Add piecewise linear model
                pwl2_model=pwl2_model if "pwl2_model" in locals() else None,  # Add piecewise linear 2-dots model
                exp_lin_model=exp_lin_model if "exp_lin_model" in locals() else None,  # Add exp-linear model
                lin_exp_model=lin_exp_model if "lin_exp_model" in locals() else None,  # Add linear-exp model
                lin_pre_model=lin_pre_model if "lin_pre_model" in locals() else None,  # Add linear pre-only model
                lin_post_model=lin_post_model if "lin_post_model" in locals() else None,  # Add linear post-only model
            )
        else:
            plot_path = None

        print(f"✅ Analysis completed for {test_name}")
        print(f"📁 Results saved in: {test_dir}")

        # Return summary for all_data.csv
        base_result = {
            "test_code": test_code,
            "test_name": test_name,
            "system": system_name,
            "n_data_points": len(test_values),
            "age_range_min": test_ages.min(),
            "age_range_max": test_ages.max(),
            "test_range_min": test_values.min(),
            "test_range_max": test_values.max(),
            "output_dir": test_dir,
            "fitted_values_file": fitted_values_path,
            "parameters_file": cp_params_path,
            "plot_file": plot_path,
            # Add quality metrics to results
            "f_statistic": quality_metrics["f_statistic"],
            "p_value": quality_metrics["p_value"],
            "is_significant": quality_metrics["is_significant"],
            "meets_min_improvement": quality_metrics["meets_min_improvement"],
            "meets_min_r_squared": quality_metrics["meets_min_r_squared"],
            "is_reasonable": quality_metrics["is_reasonable"],
            "aic": quality_metrics["aic"],
            "linear_aic": quality_metrics["linear_aic"],
            "aic_improvement": quality_metrics["aic_improvement"],
            "bic": quality_metrics["bic"],
            "linear_bic": quality_metrics["linear_bic"],
            "bic_improvement": quality_metrics["bic_improvement"],
            "residuals_std": quality_metrics["residuals_std"],
            "linear_residuals_std": quality_metrics["linear_residuals_std"],
        }

        if best_cp1 == best_cp2:  # Single change point model
            slope1, val_cp1, val_before_jump, val_after_jump, slope2, opt_cp_position = best_params
            base_result.update(
                {
                    "model_type": "single_cp",
                    "cp1_position": opt_cp_position,  # Use optimized position
                    "cp2_position": opt_cp_position,  # Same as cp1 for single CP
                    "slope1": slope1,
                    "val_cp1": val_cp1,
                    "val_before_jump": val_before_jump,
                    "val_after_jump": val_after_jump,
                    "slope2": slope2,
                    "jump_magnitude": val_after_jump - val_before_jump,
                    "r_squared_staged": r_squared,
                    "r_squared_linear": linear_model["r_squared"],
                    "improvement": r_squared - linear_model["r_squared"],
                }
            )
        else:  # Two change point model
            slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, slope2, opt_cp1, opt_cp2 = (
                best_params
            )

            # Calculate the actual slopes from change points to menopause for penalty analysis
            slope_cp1_to_meno = (
                (val_before_jump - val_cp1) / (0 - opt_cp1) if opt_cp1 != 0 else 0
            )
            slope_meno_to_cp2 = (
                (val_cp2 - val_after_jump) / (opt_cp2 - 0) if opt_cp2 != 0 else 0
            )

            base_result.update(
                {
                    "model_type": "two_cp",
                    "cp1_position": opt_cp1,  # Use optimized position from params
                    "cp2_position": opt_cp2,  # Use optimized position from params
                    "slope1": slope1,  # Original optimization parameter
                    "slope1_cp_to_meno": slope_cp1_to_meno,  # Actual slope from CP1 to menopause
                    "val_cp1": val_cp1,
                    "val_before_jump": val_before_jump,
                    "val_after_jump": val_after_jump,
                    "val_cp2": val_cp2,
                    "slope2": slope2,  # Original optimization parameter
                    "slope2_meno_to_cp": slope_meno_to_cp2,  # Actual slope from menopause to CP2
                    "jump_magnitude": val_after_jump - val_before_jump,
                    "r_squared_staged": r_squared,
                    "r_squared_linear": linear_model["r_squared"],
                    "improvement": r_squared - linear_model["r_squared"],
                }
            )

        # If bootstrap data is available, run analysis on all bootstrap files
        if bootstrap_data:
            print(
                f"\n🔄 Running bootstrap analysis on {len(bootstrap_data['files'])} bootstrap files..."
            )
            bootstrap_results = []

            # Check for existing bootstrap results
            bootstrap_dir = os.path.join(test_dir, "bootstrap_results")
            os.makedirs(bootstrap_dir, exist_ok=True)

            # Load existing results if summary file exists
            existing_results = []
            summary_bootstrap_file = os.path.join(
                output_base_dir, system_name, f"{test_name}_bootstrap_results.csv"
            )
            if os.path.exists(summary_bootstrap_file):
                try:
                    existing_df = pd.read_csv(summary_bootstrap_file)
                    existing_results = existing_df.to_dict("records")
                    print(
                        f"   📂 Found existing summary file with {len(existing_results)} results"
                    )
                except Exception as e:
                    print(f"   ⚠️  Could not load existing summary file: {e}")

            # Create a set of already processed bootstrap files
            processed_files = set()
            for result in existing_results:
                if "bootstrap_file" in result:
                    processed_files.add(result["bootstrap_file"])

            print(
                f"   🔍 Found {len(processed_files)} already processed bootstrap files"
            )

            # Prepare bootstrap tasks for parallel processing
            bootstrap_tasks = []
            for i, (bootstrap_file, bootstrap_values, bootstrap_sds) in enumerate(
                zip(
                    bootstrap_data["files"],
                    bootstrap_data["values"],
                    bootstrap_data["sds"],
                )
            ):
                # Check if this bootstrap file has already been processed
                individual_bootstrap_file = os.path.join(
                    bootstrap_dir, f"{bootstrap_file.replace('.csv', '_result.csv')}"
                )
                override_bootstrap = os.environ.get("OVERRIDE_BOOTSTRAP", "0") == "1"

                if (
                    os.path.exists(individual_bootstrap_file)
                    and bootstrap_file in processed_files
                    and not override_bootstrap
                ):
                    print(f"   ⏭️  Skipping {bootstrap_file} - already processed")
                    # Load the existing result
                    try:
                        existing_df = pd.read_csv(individual_bootstrap_file)
                        if len(existing_df) > 0:
                            existing_result = existing_df.iloc[0].to_dict()
                            bootstrap_results.append(existing_result)
                            print(f"   ✅ Loaded existing result for {bootstrap_file}")
                        continue
                    except Exception as e:
                        print(
                            f"   ⚠️  Could not load existing result for {bootstrap_file}: {e}"
                        )
                        # Continue with processing if loading fails
                elif override_bootstrap and os.path.exists(individual_bootstrap_file):
                    print(
                        f"   🔄 Overriding existing bootstrap result for {bootstrap_file}"
                    )

                # OLD CODE (BROKEN) - Same linear model for ALL bootstrap samples:
                # bootstrap_tasks.append((
                #     test_code, test_name, system_name,
                #     bootstrap_data['ages'], bootstrap_values, bootstrap_sds,
                #     survival_data, linear_model, test_dir, bootstrap_file  # <-- SAME linear_model!
                # ))

                # NEW CODE (FIXED) - Unique linear model for EACH bootstrap sample:
                bootstrap_linear_model = fit_linear_menopause_model(
                    bootstrap_data["ages"],
                    bootstrap_values,
                    survival_data,
                    bootstrap_sds,
                    use_iterative_reweighting=use_iterative_reweighting,
                )

                # Add to tasks for parallel processing
                bootstrap_tasks.append(
                    (
                        test_code,
                        test_name,
                        system_name,
                        bootstrap_data["ages"],
                        bootstrap_values,
                        bootstrap_sds,
                        test_ns,
                        survival_data,
                        bootstrap_linear_model,
                        test_dir,
                        bootstrap_file,
                        specific_cp_a,
                        specific_cp_b,
                        use_iterative_reweighting,  # Pass the new parameters
                    )
                )

            # Process bootstrap tasks in parallel
            if bootstrap_tasks:
                print(
                    f"   🚀 Processing {len(bootstrap_tasks)} bootstrap files in parallel..."
                )

                # Use more aggressive worker count for bootstrap processing
                bootstrap_workers = min(
                    max(1, n_cpus // 2), len(bootstrap_tasks), 16
                )  # Use half the cores (minimum 1), max 16 workers
                print(
                    f"   🔧 Using {bootstrap_workers} workers for bootstrap processing"
                )

                # Ensure we have enough bootstrap tasks to keep all workers busy
                if len(bootstrap_tasks) < bootstrap_workers:
                    bootstrap_workers = len(bootstrap_tasks)
                    print(
                        f"   🔧 Reduced to {bootstrap_workers} workers to match task count"
                    )

                with ProcessPoolExecutor(
                    max_workers=bootstrap_workers,
                    mp_context=multiprocessing.get_context("fork"),
                ) as executor:
                    # Submit all bootstrap tasks
                    bootstrap_futures = {
                        executor.submit(
                            run_single_bootstrap_analysis,
                            task[0],  # test_code
                            task[1],  # test_name
                            task[2],  # system_name
                            task[3],  # bootstrap_data["ages"]
                            task[4],  # bootstrap_values
                            task[5],  # bootstrap_sds
                            task[6],  # test_ns
                            task[7],  # survival_data
                            task[8],  # bootstrap_linear_model
                            task[9],  # test_dir
                            task[10],  # bootstrap_file
                            task[11] if len(task) > 11 else None,  # specific_cp_a
                            task[12] if len(task) > 12 else None,  # specific_cp_b
                            (
                                task[13] if len(task) > 13 else False
                            ),  # use_iterative_reweighting
                        ): task[
                            10
                        ]  # task[10] is bootstrap_file
                        for task in bootstrap_tasks
                    }

                    # Collect results as they complete
                    for future in as_completed(bootstrap_futures):
                        bootstrap_file = bootstrap_futures[future]
                        try:
                            bootstrap_result = future.result()
                            if bootstrap_result:
                                bootstrap_results.append(bootstrap_result)

                                # Save individual bootstrap result immediately
                                individual_bootstrap_file = os.path.join(
                                    test_dir,
                                    "bootstrap_results",
                                    f"{bootstrap_file.replace('.csv', '_result.csv')}",
                                )
                                individual_df = pd.DataFrame([bootstrap_result])
                                individual_df.to_csv(
                                    individual_bootstrap_file, index=False
                                )
                                print(
                                    f"   💾 Individual bootstrap result saved: {individual_bootstrap_file}"
                                )

                                # Also update the summary file after each bootstrap
                                if bootstrap_results:
                                    summary_bootstrap_file = os.path.join(
                                        output_base_dir,
                                        system_name,
                                        f"{test_name}_bootstrap_results.csv",
                                    )
                                    summary_df = pd.DataFrame(bootstrap_results)
                                    summary_df.to_csv(
                                        summary_bootstrap_file, index=False
                                    )
                                    print(
                                        f"   📊 Updated summary file: {summary_bootstrap_file} ({len(bootstrap_results)} results)"
                                    )
                        except Exception as exc:
                            print(
                                f"   ❌ Bootstrap task for {bootstrap_file} generated an exception: {exc}"
                            )

                print(f"   ✅ Completed parallel bootstrap processing")

            # Add bootstrap results to the main result
            if bootstrap_results:
                base_result["bootstrap_results"] = bootstrap_results
                base_result["n_bootstrap_runs"] = len(bootstrap_results)
                print(f"   ✅ Completed {len(bootstrap_results)} bootstrap analyses")

                # Final summary of all saved files
                print(f"📊 All bootstrap results saved:")
                print(
                    f"   - Summary: {output_base_dir}/{system_name}/{test_name}_bootstrap_results.csv"
                )
                print(f"   - Individual files: {test_dir}/bootstrap_results/")

            else:
                print(f"   ❌ No successful bootstrap analyses")

        return base_result

    except Exception as e:
        print(f"❌ ERROR analyzing {test_name}: {str(e)}")
        import traceback

        traceback.print_exc()
        return None


def run_single_bootstrap_analysis(
    test_code,
    test_name,
    system_name,
    test_ages,
    test_values,
    test_sds,
    test_ns,
    survival_data,
    linear_model,
    test_dir,
    bootstrap_file,
    specific_cp_a=None,
    specific_cp_b=None,
    use_iterative_reweighting=False,
):
    """Run change point analysis for a single bootstrap file"""
    
    # Enable profiling for this worker if requested
    if DEFAULT_CONFIG.get("PROFILE_WORKERS", False):
        import cProfile
        import pstats
        
        profiler = cProfile.Profile()
        profiler.enable()
        
        # Call the actual worker function
        result = _run_single_bootstrap_analysis_impl(
            test_code, test_name, system_name, test_ages, test_values,
            test_sds, test_ns, survival_data, linear_model, test_dir,
            bootstrap_file, specific_cp_a, specific_cp_b, use_iterative_reweighting
        )
        
        profiler.disable()
        
        # Save worker profile
        bootstrap_name = bootstrap_file.replace('.csv', '')
        profile_file = os.path.join(test_dir, "bootstrap_results", 
                                    f"profile_worker_{bootstrap_name}.prof")
        profiler.dump_stats(profile_file)
        print(f"   💾 Worker profile saved: {profile_file}")
        
        return result
    else:
        # No profiling - call directly
        return _run_single_bootstrap_analysis_impl(
            test_code, test_name, system_name, test_ages, test_values,
            test_sds, test_ns, survival_data, linear_model, test_dir,
            bootstrap_file, specific_cp_a, specific_cp_b, use_iterative_reweighting
        )


def _run_single_bootstrap_analysis_impl(
    test_code,
    test_name,
    system_name,
    test_ages,
    test_values,
    test_sds,
    test_ns,
    survival_data,
    linear_model,
    test_dir,
    bootstrap_file,
    specific_cp_a=None,
    specific_cp_b=None,
    use_iterative_reweighting=False,
):
    """Implementation of single bootstrap analysis (separated for profiling)"""

    # Fix: Convert all input arrays to numpy arrays to avoid pandas Series indexing errors
    test_ages = np.asarray(test_ages)
    test_values = np.asarray(test_values)
    if test_sds is not None:
        test_sds = np.asarray(test_sds)
    if test_ns is not None:
        test_ns = np.asarray(test_ns)
    
    # SIMPLE FIX: Use test_sds directly for weights (same as linear model)
    # Don't transform to SE - this ensures consistent weighting
    print(f"   🔍 BOOTSTRAP: test_sds range: {test_sds.min():.6f} to {test_sds.max():.6f}")
    print(f"   🔍 BOOTSTRAP: test_ns range: {test_ns.min()} to {test_ns.max()}")
    print(f"   🔍 BOOTSTRAP: Using test_sds DIRECTLY for weights (no SE transformation)")
    print(f"   🔍 About to enter try block for {bootstrap_file}")
    try:
        # Create output directories
        bootstrap_dir = os.path.join(test_dir, "bootstrap_results")
        os.makedirs(bootstrap_dir, exist_ok=True)

        # Run the same analysis as the main function but with bootstrap data
        # Stage 1: Single change point optimization
        n_single_positions = DEFAULT_CONFIG["N_CP_COMBINATIONS"]
        single_cp_positions = generate_single_cp_positions(
            n_single_positions, specific_cp_a=specific_cp_a
        )

        # Stage 1: Single change point optimization using multiprocessing
        # For bootstrap, use fewer CPUs to avoid resource contention since multiple bootstrap processes run in parallel
        bootstrap_cpus = max(
            2, DEFAULT_CONFIG["N_WORKERS"] // 4
        )  # Use 1/4 of available CPUs per bootstrap task
        print(f"   🔧 Bootstrap using {bootstrap_cpus} CPUs per task")

        # Ensure we have enough tasks to keep all cores busy
        if len(single_cp_positions) < bootstrap_cpus * 4:
            additional_positions = generate_additional_single_cp_positions(
                bootstrap_cpus * 4 - len(single_cp_positions)
            )
            single_cp_positions = np.concatenate(
                [single_cp_positions, additional_positions]
            )
            print(
                f"   📈 Increased single CP tasks to {len(single_cp_positions)} for better CPU utilization"
            )

        with ProcessPoolExecutor(
            max_workers=bootstrap_cpus, mp_context=multiprocessing.get_context("fork")
        ) as executor:
            # Submit single CP optimization tasks
            single_cp_futures = []
            for cp_pos in single_cp_positions:
                future = executor.submit(
                    optimize_single_cp_position,
                    cp_pos,
                    test_ages,
                    test_values,
                    survival_data,
                    linear_model,
                    test_sds,
                    test_ns,
                    use_iterative_reweighting=use_iterative_reweighting,
                )
                single_cp_futures.append(future)

            # Collect results
            single_cp_results = []
            for future in as_completed(single_cp_futures):
                single_cp_results.append(future.result())

            # Find best single change point result
            best_single_cp_params = None
            best_single_cp_loss = np.inf
            best_single_cp_position = None

            for result in single_cp_results:
                if result["success"] and result["loss"] < best_single_cp_loss:
                    best_single_cp_params = result["params"]
                    best_single_cp_loss = result["loss"]
                    best_single_cp_position = result["cp_position"]

            # Stage 2: Two change point optimization
            # Use fixed seed for bootstrap consistency
            fixed_seed = 42  # Use seed 42 for bootstrap consistency
            combinations = generate_cp_combinations(
                DEFAULT_CONFIG["N_CP_COMBINATIONS"],
                fixed_seed,
                specific_cp_a=specific_cp_a,
                specific_cp_b=specific_cp_b,
            )

            # Ensure we have enough two CP tasks to keep all cores busy
            if len(combinations) < bootstrap_cpus * 4:
                additional_combinations = generate_cp_combinations(
                    bootstrap_cpus * 4 - len(combinations), fixed_seed
                )
                combinations.extend(additional_combinations)
                print(
                    f"   📈 Increased two CP tasks to {len(combinations)} for better CPU utilization"
                )

            # Ensure we have enough single CP tasks to keep all cores busy (for the case where we need more tasks)
            if len(single_cp_positions) < bootstrap_cpus * 4:
                additional_positions = generate_additional_single_cp_positions(
                    bootstrap_cpus * 4 - len(single_cp_positions)
                )
                single_cp_positions = np.concatenate(
                    [single_cp_positions, additional_positions]
                )
                print(
                    f"   📈 Increased single CP tasks to {len(single_cp_positions)} for better CPU utilization"
                )

            # Submit two CP optimization tasks with aggressive parallel processing
            print(
                f"   🚀 Submitting {len(combinations)} two CP combinations for {bootstrap_file}"
            )
            two_cp_futures = []
            for cp1, cp2 in combinations:
                future = executor.submit(
                    optimize_two_change_point,
                    cp1,
                    cp2,
                    test_ages,
                    test_values,
                    survival_data,
                    linear_model,
                    test_sds,
                    test_ns,
                    use_iterative_reweighting=use_iterative_reweighting,
                )
                two_cp_futures.append(future)
            print(f"   ✅ Submitted {len(two_cp_futures)} two CP futures")

            # Collect results
            two_cp_results = []
            for future in as_completed(two_cp_futures):
                two_cp_results.append(future.result())
            print(
                f"   📊 Collected {len(two_cp_results)} two CP optimization results for {bootstrap_file}"
            )

            # print the amount of successful two CP results
            print(
                f"   ✅ {len([succ_result for succ_result in two_cp_results if succ_result['success']])} successful two CP results"
            )
            # Find best two change point result
            best_two_cp_params = None
            best_two_cp_loss = np.inf
            best_two_cp_cp1, best_two_cp_cp2 = None, None

            for result in two_cp_results:
                if result["success"] and result["loss"] < best_two_cp_loss:
                    best_two_cp_params = result["params"]
                    best_two_cp_loss = result["loss"]
                    best_two_cp_cp1, best_two_cp_cp2 = result["cp1"], result["cp2"]

            # Save all optimization results for CP1 heatmap analysis
            import pandas as pd

            # Save all single CP results with BIC calculation
            single_cp_all_results = []
            for result in single_cp_results:
                if result["success"]:
                    # Calculate pure RSS (without penalties) for BIC calculation
                    n_params = 6  # slope1, val_cp1, val_before_jump, val_after_jump, slope2, cp_position (search adds complexity)

                    # Extract parameters from result (now includes cp_position as last parameter)
                    (
                        slope1,
                        val_cp1,
                        val_before_jump,
                        val_after_jump,
                        slope2,
                        cp_position,
                    ) = result["params"]

                    _, pure_rss = get_final_predictions_for_single_change_point(
                        test_ages,
                        survival_data,
                        slope1,
                        val_cp1,
                        val_before_jump,
                        val_after_jump,
                        slope2,
                        cp_position,
                        test_values=test_values,
                        test_sds=test_sds,  # Use same weights as optimization
                        test_ns=None,  # Don't transform - use test_sds directly
                        return_rss=True,
                    )
                    
                    # VERIFICATION: Compare with optimization's pure_rss to detect issues
                    print(f"   ✅ SINGLE CP: Recalc RSS={pure_rss:.6f} for cp={cp_position:.2f} using result['params']")
                    if "pure_rss" in result and result["pure_rss"] is not None:
                        opt_pure_rss = result["pure_rss"]
                        rss_diff_pct = 100*abs(pure_rss - opt_pure_rss)/max(abs(opt_pure_rss), 1e-10)
                        if rss_diff_pct > 1.0:
                            print(f"   ⚠️  SINGLE CP RSS MISMATCH:")
                            print(f"      Optimization RSS:  {opt_pure_rss:.6f}")
                            print(f"      Recalculated RSS:  {pure_rss:.6f}")
                            print(f"      Difference:        {abs(pure_rss - opt_pure_rss):.6f} ({rss_diff_pct:.1f}%)")
                            print(f"      This means objective.pure_rss was from WRONG parameters!")
                        else:
                            print(f"   ✓ RSS match: {rss_diff_pct:.2f}% diff")

                    # Always use n_data for BIC calculation to ensure consistency
                    n_data = len(test_values)
                    bic = n_data * np.log(pure_rss / n_data) + np.log(n_data) * n_params

                    # Calculate log-likelihood for normal distribution
                    log_likelihood = (
                        -n_data
                        / 2
                        * (np.log(2 * np.pi) + np.log(pure_rss / n_data) + 1)
                    )

                    # Calculate slope information for single CP model
                    cp1_pos = float(result["cp_position"])

                    # For single CP: slope_before_cp1 = slope1, slope_after_cp1 = slope2
                    # The slope change at cp1 is slope2 - slope1
                    slope_before_cp1 = slope1
                    slope_after_cp1 = slope2
                    slope_change_cp1 = slope2 - slope1

                    # Log the values before storing in results
                    print(f"   📊 Storing single CP result: cp1={cp1_pos:.2f}")
                    print(f"   📊 Slopes: slope1={slope1:.6f}, slope2={slope2:.6f}")
                    print(
                        f"   📊 Values: val_cp1={val_cp1:.2f}, val_before_jump={val_before_jump:.2f}, val_after_jump={val_after_jump:.2f}"
                    )
                    print(f"   📊 RSS: pure_rss={pure_rss:.2f}, BIC={bic:.2f}")

                    # Check if slope2 is extreme
                    if abs(slope2) > 0.1:
                        print(f"   ⚠️ EXTREME SLOPE2 in single CP: slope2={slope2:.6f}")

                    single_cp_all_results.append(
                        {
                            "cp1_position": cp1_pos,
                            "cp2_position": cp1_pos,  # Same as cp1 for single CP
                            "loss": float(
                                result["loss"]
                            ),  # Keep penalized loss for optimization records
                            "pure_rss": float(pure_rss),  # Add pure RSS
                            "bic": float(bic),
                            "log_likelihood": float(log_likelihood),
                            "success": bool(result["success"]),
                            "model_type": "single_cp",
                            "slope_before_cp1": float(slope_before_cp1),
                            "slope_after_cp1": float(slope_after_cp1),
                            "slope_change_cp1": float(slope_change_cp1),
                            "val_cp1": float(val_cp1),
                            "val_before_jump": float(val_before_jump),
                            "val_after_jump": float(val_after_jump),
                            "slope2": float(
                                slope2
                            ),  # Slope after menopause (after change point)
                        }
                    )

            if single_cp_all_results:
                single_cp_df = pd.DataFrame(single_cp_all_results)
                single_cp_csv_path = os.path.join(
                    test_dir,
                    "bootstrap_results",
                    f"{bootstrap_file.replace('.csv', '')}_single_cp_all_combinations.csv",
                )
                single_cp_df.to_csv(single_cp_csv_path, index=False)

            # Save all two CP results with BIC calculation
            print(
                f"   📊 Processing {len(two_cp_results)} two CP results for {bootstrap_file}"
            )
            two_cp_all_results = []
            successful_two_cp_count = 0
            for result in two_cp_results:
                if result["success"]:
                    successful_two_cp_count += 1
                    # Calculate pure RSS (without penalties) for BIC calculation
                    n_params = 8  # slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, slope2, cp1_pos, cp2_pos (search adds complexity)

                    # Extract parameters from result (now includes cp1, cp2 positions as last two parameters)
                    (
                        slope1,
                        val_cp1,
                        val_before_jump,
                        val_after_jump,
                        val_cp2,
                        slope2,
                        cp1,
                        cp2,
                    ) = result["params"]

                    # Debug print for extreme slope2 values
                    if abs(slope2) > 0.1:
                        print(f"   🚨 EXTREME slope2 detected: {slope2:.6f}")
                        print(
                            f"   🔍 CP positions: cp1={result['cp1']:.2f}, cp2={result['cp2']:.2f}"
                        )
                        print(
                            f"   📊 Values: val_cp1={val_cp1:.2f}, val_before_jump={val_before_jump:.2f}, val_after_jump={val_after_jump:.2f}, val_cp2={val_cp2:.2f}"
                        )

                    # Calculate RSS using the helper function with transformed SE values
                    # This ensures we use the same weight calculation as during optimization
                    # Note: cp1 and cp2 are now extracted from params above

                    _, pure_rss = get_final_predictions_for_change_point(
                        test_ages,
                        survival_data,
                        slope1,
                        val_cp1,
                        val_before_jump,
                        val_after_jump,
                        val_cp2,
                        slope2,
                        cp1,
                        cp2,
                        test_values=test_values,
                        test_sds=test_sds,  # Use same weights as optimization
                        test_ns=None,  # Don't transform - use test_sds directly
                        return_rss=True,
                    )
                    
                    # VERIFICATION: Compare with optimization's pure_rss to detect issues
                    print(f"   ✅ TWO CP: Recalc RSS={pure_rss:.6f} for cp1={cp1:.2f}, cp2={cp2:.2f} using result['params']")
                    if "pure_rss" in result and result["pure_rss"] is not None:
                        opt_pure_rss = result["pure_rss"]
                        rss_diff_pct = 100*abs(pure_rss - opt_pure_rss)/max(abs(opt_pure_rss), 1e-10)
                        if rss_diff_pct > 1.0:
                            print(f"   ⚠️  TWO CP RSS MISMATCH:")
                            print(f"      Optimization RSS:  {opt_pure_rss:.6f}")
                            print(f"      Recalculated RSS:  {pure_rss:.6f}")
                            print(f"      Difference:        {abs(pure_rss - opt_pure_rss):.6f} ({rss_diff_pct:.1f}%)")
                            print(f"      This means objective.pure_rss was from WRONG parameters!")
                        else:
                            print(f"   ✓ RSS match: {rss_diff_pct:.2f}% diff")

                    # Debug print for extreme RSS values
                    if pure_rss > 1000:
                        print(f"   🚨 EXTREME RSS detected: {pure_rss:.2f}")
                        print(
                            f"   📊 Actual range: {test_values.min():.2f} to {test_values.max():.2f}"
                        )

                    # Always use n_data for BIC calculation to ensure consistency
                    n_data = len(test_values)
                    bic = n_data * np.log(pure_rss / n_data) + np.log(n_data) * n_params

                    # Calculate log-likelihood for normal distribution
                    log_likelihood = (
                        -n_data
                        / 2
                        * (np.log(2 * np.pi) + np.log(pure_rss / n_data) + 1)
                    )

                    # Calculate slope information for two CP model
                    cp1_pos = float(result["cp1"])
                    cp2_pos = float(result["cp2"])

                    # For two CP: slope_before_cp1 = slope1, slope_after_cp1 = slope from cp1 to menopause
                    # Calculate the slope from cp1 to menopause (t=0)
                    if cp1_pos != 0:
                        slope_from_cp1_to_meno = (val_before_jump - val_cp1) / (
                            0 - cp1_pos
                        )
                    else:
                        slope_from_cp1_to_meno = 0.0

                    # The slope change at cp1 is slope_from_cp1_to_meno - slope1
                    slope_change_cp1 = slope_from_cp1_to_meno - slope1

                    # For the second change point: slope_before_cp2 = slope from menopause to cp2
                    # Calculate the slope from menopause (t=0) to cp2
                    if cp2_pos != 0:
                        slope_from_meno_to_cp2 = (val_cp2 - val_after_jump) / (
                            cp2_pos - 0
                        )
                    else:
                        slope_from_meno_to_cp2 = 0.0

                    # The slope change at cp2 is slope2 - slope_from_meno_to_cp2
                    slope_change_cp2 = slope2 - slope_from_meno_to_cp2

                    # Log the values before storing in results
                    print(
                        f"   📊 Storing two CP result: cp1={cp1_pos:.2f}, cp2={cp2_pos:.2f}"
                    )
                    print(f"   📊 Slopes: slope1={slope1:.6f}, slope2={slope2:.6f}")
                    print(
                        f"   📊 Calculated slopes: slope_from_cp1_to_meno={slope_from_cp1_to_meno:.6f}, slope_from_meno_to_cp2={slope_from_meno_to_cp2:.6f}"
                    )
                    print(
                        f"   📊 Values: val_cp1={val_cp1:.2f}, val_cp2={val_cp2:.2f}, val_before_jump={val_before_jump:.2f}, val_after_jump={val_after_jump:.2f}"
                    )
                    print(f"   📊 RSS: pure_rss={pure_rss:.2f}, BIC={bic:.2f}")

                    # Check if slope2 is very different from slope_from_meno_to_cp2
                    if abs(slope2) > 1.0 and abs(slope2 - slope_from_meno_to_cp2) > 0.1:
                        print(
                            f"   ⚠️ EXTREME SLOPE DIFFERENCE: slope2={slope2:.6f}, slope_from_meno_to_cp2={slope_from_meno_to_cp2:.6f}"
                        )
                        print(
                            f"   ⚠️ Difference: {abs(slope2 - slope_from_meno_to_cp2):.6f}"
                        )

                    two_cp_all_results.append(
                        {
                            "cp1_position": cp1_pos,
                            "cp2_position": cp2_pos,
                            "loss": float(
                                result["loss"]
                            ),  # Keep penalized loss for optimization records
                            "pure_rss": float(pure_rss),  # Add pure RSS
                            "bic": float(bic),
                            "log_likelihood": float(log_likelihood),
                            "success": bool(result["success"]),
                            "model_type": "two_cp",
                            "slope_before_cp1": float(slope1),
                            "slope_after_cp1": float(slope_from_cp1_to_meno),
                            "slope_change_cp1": float(slope_change_cp1),
                            "slope_before_cp2": float(slope_from_meno_to_cp2),
                            "slope_after_cp2": float(slope2),
                            "slope_change_cp2": float(slope_change_cp2),
                            "val_cp1": float(val_cp1),
                            "val_before_jump": float(val_before_jump),
                            "val_after_jump": float(val_after_jump),
                            "val_cp2": float(val_cp2),
                        }
                    )

            print(
                f"   ✅ Processed {successful_two_cp_count} successful two CP results out of {len(two_cp_results)} total"
            )
            if two_cp_all_results:
                print(
                    f"   💾 Creating two CP DataFrame with {len(two_cp_all_results)} rows"
                )
                two_cp_df = pd.DataFrame(two_cp_all_results)
                two_cp_csv_path = os.path.join(
                    test_dir,
                    "bootstrap_results",
                    f"{bootstrap_file.replace('.csv', '')}_two_cp_all_combinations.csv",
                )
                print(f"   📝 Saving two CP results to {two_cp_csv_path}")
                two_cp_df.to_csv(two_cp_csv_path, index=False)
                print(f"   ✅ Two CP file saved successfully")
            else:
                print(
                    f"   ❌ No successful two CP results to save for {bootstrap_file}"
                )

            # Fit sigmoid model
            print(f"   🔄 Fitting sigmoid model for bootstrap {bootstrap_file}")
            sigmoid_model = fit_sigmoid_menopause_model(
                test_ages,
                test_values,
                survival_data,
                linear_model,
                test_sds,
                test_ns,
                use_iterative_reweighting=use_iterative_reweighting,
            )

            # Fit polynomial4 model
            print(f"   🔄 Fitting polynomial4 model for bootstrap {bootstrap_file}")
            poly4_model = fit_polynomial4_model(
                test_ages,
                test_values,
                survival_data,
                linear_model,
                test_sds,
                test_ns,
                use_iterative_reweighting=use_iterative_reweighting,
            )

            # Fit piecewise linear continuous model
            print(f"   🔄 Fitting piecewise linear continuous model for bootstrap {bootstrap_file}")
            pwl_model = fit_piecewise_linear_continuous_model(
                test_ages,
                test_values,
                survival_data,
                linear_model,
                test_sds,
                test_ns,
                use_iterative_reweighting=use_iterative_reweighting,
            )

            # Fit exp-linear model
            print(f"   🔄 Fitting exp-linear model for bootstrap {bootstrap_file}")
            exp_lin_model = fit_exp_linear_model(
                test_ages,
                test_values,
                survival_data,
                linear_model,
                test_sds,
                test_ns,
                use_iterative_reweighting=use_iterative_reweighting,
            )

            # Fit linear-exp model
            print(f"   🔄 Fitting linear-exp model for bootstrap {bootstrap_file}")
            lin_exp_model = fit_linear_exp_model(
                test_ages,
                test_values,
                survival_data,
                linear_model,
                test_sds,
                test_ns,
                use_iterative_reweighting=use_iterative_reweighting,
            )

            # Fit polynomial3 model
            print(f"   🔄 Fitting polynomial3 model for bootstrap {bootstrap_file}")
            poly3_model = fit_polynomial3_model(
                test_ages,
                test_values,
                survival_data,
                linear_model,
                test_sds,
                test_ns,
                use_iterative_reweighting=use_iterative_reweighting,
            )

            # Fit piecewise linear 2-dots model
            print(f"   🔄 Fitting piecewise linear 2-dots model for bootstrap {bootstrap_file}")
            pwl2_model = fit_piecewise_linear_2dots_model(
                test_ages,
                test_values,
                survival_data,
                linear_model,
                test_sds,
                test_ns,
                use_iterative_reweighting=use_iterative_reweighting,
            )

            # Calculate results for ALL models (linear, sigmoid, poly3, poly4, pwl2, pwl3, exp_lin, lin_exp, single CP, two CP)
            bootstrap_result = {
                "bootstrap_file": bootstrap_file,
            }

            # Linear model results (always available)
            # Use the bootstrap-specific linear model fitted values
            linear_predictions = linear_model["fitted_values"]
            linear_residuals = test_values - linear_predictions
            n_data = len(test_values)
            linear_n_params = 4  # a, b, c, d

            # Calculate linear RSS using integral approach for consistency with CP models
            linear_rss_for_bic = calculate_pure_rss(
                test_ages,
                test_values,
                survival_data,
                linear_a=linear_model["a"],
                linear_b=linear_model["b"],
                linear_c=linear_model["c"],
                linear_d=linear_model["d"],
                test_sds=test_sds,
                test_ns=None,  # Don't transform - use test_sds directly
            )

            # Calculate BIC using the same approach as change point models
            # Use n_data for both AIC and BIC (not n_eff) to ensure exact consistency
            # SIMPLE FIX: Use test_sds directly (same as linear model)
            if test_sds is not None and np.any(test_sds > 0):
                weights = 1.0 / (test_sds**2)
                # For R-squared, calculate TSS with proper weighting
                linear_tss = np.sum(
                    weights
                    * (test_values - np.average(test_values, weights=weights)) ** 2
                )
                # Use the RSS from residuals for both AIC and BIC (same as CP models use their loss)
                linear_aic = (
                    n_data * np.log(linear_rss_for_bic / n_data) + 2 * linear_n_params
                )
                linear_bic = (
                    n_data * np.log(linear_rss_for_bic / n_data)
                    + np.log(n_data) * linear_n_params
                )
                linear_log_likelihood = (
                    -n_data
                    / 2
                    * (np.log(2 * np.pi) + np.log(linear_rss_for_bic / n_data) + 1)
                )
                linear_r_squared = (
                    1 - linear_rss_for_bic / linear_tss if linear_tss > 0 else 0
                )
            else:
                # For unweighted case
                linear_tss = np.sum((test_values - np.mean(test_values)) ** 2)
                # Use the RSS from residuals for both AIC and BIC (same as CP models use their loss)
                linear_aic = (
                    n_data * np.log(linear_rss_for_bic / n_data) + 2 * linear_n_params
                )
                linear_bic = (
                    n_data * np.log(linear_rss_for_bic / n_data)
                    + np.log(n_data) * linear_n_params
                )
                linear_log_likelihood = (
                    -n_data
                    / 2
                    * (np.log(2 * np.pi) + np.log(linear_rss_for_bic / n_data) + 1)
                )
                linear_r_squared = (
                    1 - linear_rss_for_bic / linear_tss if linear_tss > 0 else 0
                )

            # Add linear model parameters for jump magnitude calculation
            bootstrap_result.update(
                {
                    "linear_r_squared": linear_r_squared,
                    "linear_aic": linear_aic,
                    "linear_bic": linear_bic,
                    "linear_log_likelihood": linear_log_likelihood,
                    "linear_pure_rss": linear_rss_for_bic,  # Add pure_rss to the result
                    "linear_a": linear_model["a"],  # Post-menopause intercept
                    "linear_b": linear_model["b"],  # Post-menopause slope
                    "linear_c": linear_model["c"],  # Pre-menopause intercept
                    "linear_d": linear_model["d"],  # Pre-menopause slope
                    "linear_jump_magnitude": linear_model["a"]
                    - linear_model["c"],  # Jump = a - c
                }
            )

            # Sigmoid model results
            if sigmoid_model is not None:
                # Calculate AIC and BIC for sigmoid model
                sigmoid_n_params = 5  # c, d, w, h, b
                sigmoid_rss = sigmoid_model["rss"]
                
                n_data = len(test_values)
                sigmoid_aic = n_data * np.log(sigmoid_rss / n_data) + 2 * sigmoid_n_params
                sigmoid_bic = (
                    n_data * np.log(sigmoid_rss / n_data)
                    + np.log(n_data) * sigmoid_n_params
                )
                
                sigmoid_log_likelihood = (
                    -n_data / 2 * (np.log(2 * np.pi) + np.log(sigmoid_rss / n_data) + 1)
                )
                
                bootstrap_result.update(
                    {
                        "sigmoid_c": sigmoid_model["c"],
                        "sigmoid_d": sigmoid_model["d"],
                        "sigmoid_w": sigmoid_model["w"],
                        "sigmoid_h": sigmoid_model["h"],
                        "sigmoid_b": sigmoid_model["b"],
                        "sigmoid_r_squared": sigmoid_model["r_squared"],
                        "sigmoid_improvement": sigmoid_model["r_squared"] - linear_r_squared,
                        "sigmoid_aic": sigmoid_aic,
                        "sigmoid_bic": sigmoid_bic,
                        "sigmoid_log_likelihood": sigmoid_log_likelihood,
                        "sigmoid_pure_rss": sigmoid_rss,
                        "sigmoid_aic_improvement": linear_aic - sigmoid_aic,
                        "sigmoid_bic_improvement": linear_bic - sigmoid_bic,
                    }
                )
                print(f"   ✅ Sigmoid model results added to bootstrap: R²={sigmoid_model['r_squared']:.4f}, w={sigmoid_model['w']:.2f}")
            else:
                # Set default values when sigmoid model is not available
                bootstrap_result.update(
                    {
                        "sigmoid_c": np.nan,
                        "sigmoid_d": np.nan,
                        "sigmoid_w": np.nan,
                        "sigmoid_h": np.nan,
                        "sigmoid_b": np.nan,
                        "sigmoid_r_squared": np.nan,
                        "sigmoid_improvement": np.nan,
                        "sigmoid_aic": np.nan,
                        "sigmoid_bic": np.nan,
                        "sigmoid_log_likelihood": np.nan,
                        "sigmoid_pure_rss": np.nan,
                        "sigmoid_aic_improvement": np.nan,
                        "sigmoid_bic_improvement": np.nan,
                    }
                )
                print(f"   ⚠️  Sigmoid model fitting failed for bootstrap {bootstrap_file}")

            # Polynomial4 model results
            if poly4_model is not None:
                poly4_n_params = 5  # a0, a1, a2, a3, a4
                poly4_rss = poly4_model["rss"]
                
                n_data = len(test_values)
                poly4_aic = n_data * np.log(poly4_rss / n_data) + 2 * poly4_n_params
                poly4_bic = (
                    n_data * np.log(poly4_rss / n_data)
                    + np.log(n_data) * poly4_n_params
                )
                
                bootstrap_result.update(
                    {
                        "poly4_r_squared": poly4_model["r_squared"],
                        "poly4_aic": poly4_aic,
                        "poly4_bic": poly4_bic,
                        "poly4_pure_rss": poly4_rss,
                        "poly4_a0": poly4_model["a0"],
                        "poly4_a1": poly4_model["a1"],
                        "poly4_a2": poly4_model["a2"],
                        "poly4_a3": poly4_model["a3"],
                        "poly4_a4": poly4_model["a4"],
                    }
                )
                print(f"   ✅ Polynomial4 model results added: R²={poly4_model['r_squared']:.4f}")
            else:
                bootstrap_result.update(
                    {
                        "poly4_r_squared": np.nan,
                        "poly4_aic": np.nan,
                        "poly4_bic": np.nan,
                        "poly4_pure_rss": np.nan,
                    }
                )

            # Piecewise linear continuous model results
            if pwl_model is not None:
                pwl_n_params = 5  # v_m12_5, v_0, v_12_5, slope_before, slope_after
                pwl_rss = pwl_model["rss"]
                
                n_data = len(test_values)
                pwl_aic = n_data * np.log(pwl_rss / n_data) + 2 * pwl_n_params
                pwl_bic = (
                    n_data * np.log(pwl_rss / n_data)
                    + np.log(n_data) * pwl_n_params
                )
                
                bootstrap_result.update(
                    {
                        "pwl_r_squared": pwl_model["r_squared"],
                        "pwl_aic": pwl_aic,
                        "pwl_bic": pwl_bic,
                        "pwl_pure_rss": pwl_rss,
                        "pwl_v_m12_5": pwl_model["v_m12_5"],
                        "pwl_v_0": pwl_model["v_0"],
                        "pwl_v_12_5": pwl_model["v_12_5"],
                        "pwl_slope_before": pwl_model["slope_before"],
                        "pwl_slope_after": pwl_model["slope_after"],
                    }
                )
                print(f"   ✅ Piecewise linear continuous model results added: R²={pwl_model['r_squared']:.4f}")
            else:
                bootstrap_result.update(
                    {
                        "pwl_r_squared": np.nan,
                        "pwl_aic": np.nan,
                        "pwl_bic": np.nan,
                        "pwl_pure_rss": np.nan,
                    }
                )

            # Exp-linear model results
            if exp_lin_model is not None:
                exp_lin_n_params = 5
                exp_lin_rss = exp_lin_model["rss"]
                
                n_data = len(test_values)
                exp_lin_aic = n_data * np.log(exp_lin_rss / n_data) + 2 * exp_lin_n_params
                exp_lin_bic = (
                    n_data * np.log(exp_lin_rss / n_data)
                    + np.log(n_data) * exp_lin_n_params
                )
                
                bootstrap_result.update(
                    {
                        "exp_lin_r_squared": exp_lin_model["r_squared"],
                        "exp_lin_aic": exp_lin_aic,
                        "exp_lin_bic": exp_lin_bic,
                        "exp_lin_pure_rss": exp_lin_rss,
                        "exp_lin_exp_a": exp_lin_model["exp_a"],
                        "exp_lin_exp_b": exp_lin_model["exp_b"],
                        "exp_lin_val_before_jump": exp_lin_model["val_before_jump"],
                        "exp_lin_val_after_jump": exp_lin_model["val_after_jump"],
                        "exp_lin_slope_post": exp_lin_model["slope_post"],
                    }
                )
                print(f"   ✅ Exp-linear model results added: R²={exp_lin_model['r_squared']:.4f}")
            else:
                bootstrap_result.update(
                    {
                        "exp_lin_r_squared": np.nan,
                        "exp_lin_aic": np.nan,
                        "exp_lin_bic": np.nan,
                        "exp_lin_pure_rss": np.nan,
                    }
                )

            # Linear-exp model results
            if lin_exp_model is not None:
                lin_exp_n_params = 5
                lin_exp_rss = lin_exp_model["rss"]
                
                n_data = len(test_values)
                lin_exp_aic = n_data * np.log(lin_exp_rss / n_data) + 2 * lin_exp_n_params
                lin_exp_bic = (
                    n_data * np.log(lin_exp_rss / n_data)
                    + np.log(n_data) * lin_exp_n_params
                )
                
                bootstrap_result.update(
                    {
                        "lin_exp_r_squared": lin_exp_model["r_squared"],
                        "lin_exp_aic": lin_exp_aic,
                        "lin_exp_bic": lin_exp_bic,
                        "lin_exp_pure_rss": lin_exp_rss,
                        "lin_exp_slope_pre": lin_exp_model["slope_pre"],
                        "lin_exp_val_before_jump": lin_exp_model["val_before_jump"],
                        "lin_exp_val_after_jump": lin_exp_model["val_after_jump"],
                        "lin_exp_exp_a": lin_exp_model["exp_a"],
                        "lin_exp_exp_b": lin_exp_model["exp_b"],
                    }
                )
                print(f"   ✅ Linear-exp model results added: R²={lin_exp_model['r_squared']:.4f}")
            else:
                bootstrap_result.update(
                    {
                        "lin_exp_r_squared": np.nan,
                        "lin_exp_aic": np.nan,
                        "lin_exp_bic": np.nan,
                        "lin_exp_pure_rss": np.nan,
                    }
                )

            # Polynomial3 model results
            if poly3_model is not None:
                poly3_n_params = 4
                poly3_rss = poly3_model["rss"]
                
                n_data = len(test_values)
                poly3_aic = n_data * np.log(poly3_rss / n_data) + 2 * poly3_n_params
                poly3_bic = (
                    n_data * np.log(poly3_rss / n_data)
                    + np.log(n_data) * poly3_n_params
                )
                
                bootstrap_result.update(
                    {
                        "poly3_r_squared": poly3_model["r_squared"],
                        "poly3_aic": poly3_aic,
                        "poly3_bic": poly3_bic,
                        "poly3_pure_rss": poly3_rss,
                        "poly3_a0": poly3_model["a0"],
                        "poly3_a1": poly3_model["a1"],
                        "poly3_a2": poly3_model["a2"],
                        "poly3_a3": poly3_model["a3"],
                    }
                )
                print(f"   ✅ Polynomial3 model results added: R²={poly3_model['r_squared']:.4f}")
            else:
                bootstrap_result.update(
                    {
                        "poly3_r_squared": np.nan,
                        "poly3_aic": np.nan,
                        "poly3_bic": np.nan,
                        "poly3_pure_rss": np.nan,
                    }
                )

            # Piecewise linear 2-dots model results
            if pwl2_model is not None:
                pwl2_n_params = 4
                pwl2_rss = pwl2_model["rss"]
                
                n_data = len(test_values)
                pwl2_aic = n_data * np.log(pwl2_rss / n_data) + 2 * pwl2_n_params
                pwl2_bic = (
                    n_data * np.log(pwl2_rss / n_data)
                    + np.log(n_data) * pwl2_n_params
                )
                
                bootstrap_result.update(
                    {
                        "pwl2_r_squared": pwl2_model["r_squared"],
                        "pwl2_aic": pwl2_aic,
                        "pwl2_bic": pwl2_bic,
                        "pwl2_pure_rss": pwl2_rss,
                        "pwl2_v_m8_33": pwl2_model["v_m8_33"],
                        "pwl2_v_8_33": pwl2_model["v_8_33"],
                        "pwl2_slope_before": pwl2_model["slope_before"],
                        "pwl2_slope_after": pwl2_model["slope_after"],
                    }
                )
                print(f"   ✅ Piecewise linear 2-dots model results added: R²={pwl2_model['r_squared']:.4f}")
            else:
                bootstrap_result.update(
                    {
                        "pwl2_r_squared": np.nan,
                        "pwl2_aic": np.nan,
                        "pwl2_bic": np.nan,
                        "pwl2_pure_rss": np.nan,
                    }
                )

            # Single CP model results
            if best_single_cp_params is not None:
                slope1, val_cp1, val_before_jump, val_after_jump, slope2, opt_cp_position = (
                    best_single_cp_params
                )

                # Calculate predictions and metrics for single CP
                single_cp_predictions = []
                for age in test_ages:
                    pred = pred_age_single_cp(
                        age,
                        survival_data,
                        slope1,
                        val_cp1,
                        val_before_jump,
                        val_after_jump,
                        slope2,
                        opt_cp_position,  # Use optimized position from params
                    )
                    single_cp_predictions.append(pred)
                single_cp_predictions = np.array(single_cp_predictions)

                # Calculate R-squared for single CP using the same weighting scheme as linear model
                single_cp_n_params = 6  # slope1, val_cp1, val_before_jump, val_after_jump, slope2, cp_position (search adds complexity)
                single_cp_residuals = test_values - single_cp_predictions

                # Calculate BIC using pure RSS (without penalties)
                # Use the pure RSS from optimization if available, otherwise recalculate
                if (
                    hasattr(best_single_cp_params, "pure_rss")
                    and best_single_cp_params.pure_rss is not None
                ):
                    pure_rss_single = best_single_cp_params.pure_rss
                else:
                    # Fallback: recalculate pure RSS from fitted parameters
                    # Use pure_rss directly from optimization result instead of recalculating
                    if (
                        hasattr(best_single_cp_params, "pure_rss")
                        and best_single_cp_params.pure_rss is not None
                    ):
                        pure_rss_single = best_single_cp_params.pure_rss
                        print(
                            f"   📊 Using pure_rss from optimization for single CP: {pure_rss_single:.6f}"
                        )
                    else:
                        # Fallback to calculating it if not available
                        print(
                            "   ⚠️ pure_rss not found in optimization result, recalculating"
                        )
                        pure_rss_single = calculate_pure_rss(
                            test_ages,
                            test_values,
                            survival_data,
                            slope1,
                            val_cp1,
                            val_before_jump,
                            val_after_jump,
                            val_cp1,
                            slope2,
                            best_single_cp_position,
                            best_single_cp_position,  # cp2 = cp1 for single CP
                            test_sds=test_sds,
                            test_ns=None,  # Don't transform - use test_sds directly
                        )

                # SIMPLE FIX: Use test_sds directly (same as linear model)
                if test_sds is not None and np.any(test_sds > 0):
                    weights = 1.0 / (test_sds**2)
                    # For R-squared calculation, use weighted RSS
                    single_cp_rss = np.sum(weights * single_cp_residuals**2)
                    single_cp_tss = np.sum(
                        weights
                        * (test_values - np.average(test_values, weights=weights)) ** 2
                    )
                    # Note: We use n_data instead of n_eff for BIC calculation to ensure consistency
                    # This ensures BIC values are comparable across different models and files
                    n_data = len(test_values)
                    single_cp_aic = (
                        n_data * np.log(pure_rss_single / n_data)
                        + 2 * single_cp_n_params
                    )
                    single_cp_bic = (
                        n_data * np.log(pure_rss_single / n_data)
                        + np.log(n_data) * single_cp_n_params
                    )
                else:
                    single_cp_rss = np.sum(single_cp_residuals**2)
                    n_data = len(test_values)
                    single_cp_tss = np.sum((test_values - np.mean(test_values)) ** 2)
                    single_cp_aic = (
                        n_data * np.log(pure_rss_single / n_data)
                        + 2 * single_cp_n_params
                    )
                    single_cp_bic = (
                        n_data * np.log(pure_rss_single / n_data)
                        + np.log(n_data) * single_cp_n_params
                    )

                single_cp_r_squared = (
                    1 - single_cp_rss / single_cp_tss if single_cp_tss > 0 else 0
                )

                # Add single CP results
                bootstrap_result.update(
                    {
                        "single_cp_model_type": "single_cp",
                        "single_cp_position": best_single_cp_position,
                        "single_cp_slope1": slope1,
                        "single_cp_val_cp": val_cp1,
                        "single_cp_val_before_jump": val_before_jump,
                        "single_cp_val_after_jump": val_after_jump,
                        "single_cp_slope2": slope2,
                        "single_cp_jump_magnitude": val_after_jump - val_before_jump,
                        "single_cp_r_squared": single_cp_r_squared,
                        "single_cp_improvement": single_cp_r_squared - linear_r_squared,
                        "single_cp_loss": best_single_cp_loss,
                        "single_cp_pure_rss": pure_rss_single,  # Add pure_rss to the result
                        "single_cp_aic": single_cp_aic,
                        "single_cp_bic": single_cp_bic,
                        "single_cp_aic_improvement": linear_aic - single_cp_aic,
                        "single_cp_bic_improvement": linear_bic - single_cp_bic,
                    }
                )
            else:
                # Set default values when single CP model is not available
                bootstrap_result.update(
                    {
                        "single_cp_model_type": "none",
                        "single_cp_position": np.nan,
                        "single_cp_slope1": np.nan,
                        "single_cp_val_cp": np.nan,
                        "single_cp_val_before_jump": np.nan,
                        "single_cp_val_after_jump": np.nan,
                        "single_cp_slope2": np.nan,
                        "single_cp_jump_magnitude": np.nan,
                        "single_cp_r_squared": np.nan,
                        "single_cp_improvement": np.nan,
                        "single_cp_loss": np.nan,
                        "single_cp_pure_rss": np.nan,
                        "single_cp_aic": np.nan,
                        "single_cp_bic": np.nan,
                        "single_cp_aic_improvement": np.nan,
                        "single_cp_bic_improvement": np.nan,
                    }
                )

            # Two CP model results
            if best_two_cp_params is not None:
                slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, slope2, opt_cp1, opt_cp2 = (
                    best_two_cp_params
                )

                # Calculate predictions and metrics for two CP
                two_cp_predictions = []
                for age in test_ages:
                    pred = pred_age(
                        age,
                        survival_data,
                        slope1,
                        val_cp1,
                        val_before_jump,
                        val_after_jump,
                        val_cp2,
                        slope2,
                        opt_cp1,  # Use optimized positions from params
                        opt_cp2,
                    )
                    two_cp_predictions.append(pred)
                two_cp_predictions = np.array(two_cp_predictions)
                print()
                print(
                    f"   📊 DEBUG: Prediction range: {two_cp_predictions.min():.6f} to {two_cp_predictions.max():.6f}"
                )
                if not np.all(np.isfinite(two_cp_predictions)):
                    print(f"   🚨 DEBUG: Non-finite predictions detected!")
                    print(
                        f"   📊 DEBUG: First 10 predictions: {two_cp_predictions[:10]}"
                    )

                # Calculate R-squared for two CP using the same weighting scheme as linear model
                two_cp_n_params = 8  # slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, slope2, cp1_pos, cp2_pos (search adds complexity)
                two_cp_residuals = test_values - two_cp_predictions

                # Calculate BIC using pure RSS (without penalties)
                # Use pure_rss directly from optimization result instead of recalculating
                if (
                    hasattr(best_two_cp_params, "pure_rss")
                    and best_two_cp_params.pure_rss is not None
                ):
                    pure_rss_two = best_two_cp_params.pure_rss
                    print(
                        f"   📊 Using pure_rss from optimization for two CP: {pure_rss_two:.6f}"
                    )
                else:
                    # Fallback to calculating it if not available
                    print(
                        "   ⚠️ pure_rss not found in optimization result for two CP, recalculating"
                    )
                    pure_rss_two = calculate_pure_rss(
                        test_ages,
                        test_values,
                        survival_data,
                        slope1,
                        val_cp1,
                        val_before_jump,
                        val_after_jump,
                        val_cp2,
                        slope2,
                        best_two_cp_cp1,
                        best_two_cp_cp2,
                        test_sds=test_sds,
                        test_ns=None,  # Don't transform - use test_sds directly
                    )

                # SIMPLE FIX: Use test_sds directly (same as linear model)
                if test_sds is not None and np.any(test_sds > 0):
                    weights = 1.0 / (test_sds**2)
                    # For R-squared calculation, use weighted RSS
                    two_cp_rss = np.sum(weights * two_cp_residuals**2)
                    two_cp_tss = np.sum(
                        weights
                        * (test_values - np.average(test_values, weights=weights)) ** 2
                    )
                    # Note: We use n_data instead of n_eff for BIC calculation to ensure consistency
                    # This ensures BIC values are comparable across different models and files
                    n_data = len(test_values)
                    two_cp_aic = (
                        n_data * np.log(pure_rss_two / n_data) + 2 * two_cp_n_params
                    )
                    two_cp_bic = (
                        n_data * np.log(pure_rss_two / n_data)
                        + np.log(n_data) * two_cp_n_params
                    )
                else:
                    two_cp_rss = np.sum(two_cp_residuals**2)
                    n_data = len(test_values)
                    two_cp_tss = np.sum((test_values - np.mean(test_values)) ** 2)
                    two_cp_aic = (
                        n_data * np.log(pure_rss_two / n_data) + 2 * two_cp_n_params
                    )
                    two_cp_bic = (
                        n_data * np.log(pure_rss_two / n_data)
                        + np.log(n_data) * two_cp_n_params
                    )

                two_cp_r_squared = 1 - two_cp_rss / two_cp_tss if two_cp_tss > 0 else 0

                # Add two CP results
                bootstrap_result.update(
                    {
                        "two_cp_model_type": "two_cp",
                        "two_cp_cp1_position": best_two_cp_cp1,
                        "two_cp_cp2_position": best_two_cp_cp2,
                        "two_cp_slope1": slope1,
                        "two_cp_val_cp1": val_cp1,
                        "two_cp_val_before_jump": val_before_jump,
                        "two_cp_val_after_jump": val_after_jump,
                        "two_cp_val_cp2": val_cp2,
                        "two_cp_slope2": slope2,
                        "two_cp_jump_magnitude": val_after_jump - val_before_jump,
                        "two_cp_r_squared": two_cp_r_squared,
                        "two_cp_improvement": two_cp_r_squared - linear_r_squared,
                        "two_cp_loss": best_two_cp_loss,
                        "two_cp_pure_rss": pure_rss_two,  # Add pure_rss to the result
                        "two_cp_aic": two_cp_aic,
                        "two_cp_bic": two_cp_bic,
                        "two_cp_aic_improvement": linear_aic - two_cp_aic,
                        "two_cp_bic_improvement": linear_bic - two_cp_bic,
                    }
                )
            else:
                # Set default values when two CP model is not available
                bootstrap_result.update(
                    {
                        "two_cp_model_type": "none",
                        "two_cp_cp1_position": np.nan,
                        "two_cp_cp2_position": np.nan,
                        "two_cp_slope1": np.nan,
                        "two_cp_val_cp1": np.nan,
                        "two_cp_val_before_jump": np.nan,
                        "two_cp_val_after_jump": np.nan,
                        "two_cp_val_cp2": np.nan,
                        "two_cp_slope2": np.nan,
                        "two_cp_jump_magnitude": np.nan,
                        "two_cp_r_squared": np.nan,
                        "two_cp_improvement": np.nan,
                        "two_cp_loss": np.nan,
                        "two_cp_pure_rss": np.nan,
                        "two_cp_aic": np.nan,
                        "two_cp_bic": np.nan,
                        "two_cp_aic_improvement": np.nan,
                        "two_cp_bic_improvement": np.nan,
                    }
                )

            # Model selection for backward compatibility
            if best_single_cp_params is not None and best_two_cp_params is not None:
                n_params_single = 6  # slope1, val_cp, val_before_jump, val_after_jump, slope2, cp_position
                n_params_two = 8  # slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, slope2, cp1_pos, cp2_pos
                n_data = len(test_values)

                lrt_statistic = n_data * np.log(best_single_cp_loss / best_two_cp_loss)
                critical_value = 3.84

                if lrt_statistic > critical_value:
                    best_params = best_two_cp_params
                    best_loss = best_two_cp_loss
                    best_cp1, best_cp2 = best_two_cp_cp1, best_two_cp_cp2
                    model_type = "two_cp"
                else:
                    best_params = best_single_cp_params
                    best_loss = best_single_cp_loss
                    best_cp1, best_cp2 = (
                        best_single_cp_position,
                        best_single_cp_position,
                    )
                    model_type = "single_cp"
            elif best_single_cp_params is not None:
                best_params = best_single_cp_params
                best_loss = best_single_cp_loss
                best_cp1, best_cp2 = best_single_cp_position, best_single_cp_position
                model_type = "single_cp"
            elif best_two_cp_params is not None:
                best_params = best_two_cp_params
                best_loss = best_two_cp_loss
                best_cp1, best_cp2 = best_two_cp_cp1, best_two_cp_cp2
                model_type = "two_cp"
            else:
                return None
            if model_type == "two_cp":
                # Add backward compatibility fields with all three R², AIC, and BIC values
                slope1, val_cp1, val_before_jump, val_after_jump, val_cp2, slope2, _, _ = (
                    best_params
                )
            elif model_type == "single_cp":
                slope1, val_cp1, val_before_jump, val_after_jump, slope2, _ = best_params
                val_cp2 = val_cp1

            # Log the final selected model parameters
            print(f"\n   🏆 SELECTED MODEL: {model_type}")
            print(f"   🏆 CP positions: cp1={best_cp1:.2f}, cp2={best_cp2:.2f}")
            print(f"   🏆 Slopes: slope1={slope1:.6f}, slope2={slope2:.6f}")
            print(
                f"   🏆 Values: val_cp1={val_cp1:.2f}, val_cp2={val_cp2:.2f}, val_before_jump={val_before_jump:.2f}, val_after_jump={val_after_jump:.2f}"
            )
            print(f"   🏆 Jump magnitude: {val_after_jump - val_before_jump:.2f}")

            if abs(slope2) > 1.0:
                print(
                    f"   ⚠️ EXTREME SLOPE2 in FINAL SELECTED MODEL: slope2={slope2:.6f}"
                )

            bootstrap_result.update(
                {
                    "model_type": model_type,
                    "cp1_position": best_cp1,
                    "cp2_position": best_cp2,
                    "slope1": slope1,
                    "val_cp1": val_cp1,
                    "val_before_jump": val_before_jump,
                    "val_after_jump": val_after_jump,
                    "val_cp2": val_cp2,
                    "slope2": slope2,
                    "jump_magnitude": val_after_jump - val_before_jump,
                    "r_squared_linear": linear_r_squared,
                    "r_squared_single_cp": bootstrap_result.get(
                        "single_cp_r_squared", np.nan
                    ),
                    "r_squared_two_cp": bootstrap_result.get(
                        "two_cp_r_squared", np.nan
                    ),
                    "r_squared_staged": bootstrap_result.get(
                        f"{model_type}_r_squared", 0
                    ),
                    "improvement": bootstrap_result.get(f"{model_type}_improvement", 0),
                    "loss": best_loss,
                    "aic_linear": linear_aic,
                    "aic_single_cp": bootstrap_result.get("single_cp_aic", np.nan),
                    "aic_two_cp": bootstrap_result.get("two_cp_aic", np.nan),
                    "bic_linear": linear_bic,
                    "bic_single_cp": bootstrap_result.get("single_cp_bic", np.nan),
                    "bic_two_cp": bootstrap_result.get("two_cp_bic", np.nan),
                    "aic_improvement_single_cp": bootstrap_result.get(
                        "single_cp_aic_improvement", np.nan
                    ),
                    "aic_improvement_two_cp": bootstrap_result.get(
                        "two_cp_aic_improvement", np.nan
                    ),
                    "bic_improvement_single_cp": bootstrap_result.get(
                        "single_cp_bic_improvement", np.nan
                    ),
                    "bic_improvement_two_cp": bootstrap_result.get(
                        "two_cp_bic_improvement", np.nan
                    ),
                }
            )

            # Create individual plot for this bootstrap
            try:
                # Create plots folder
                plots_folder = os.path.join(test_dir, "plots_folder")
                os.makedirs(plots_folder, exist_ok=True)

                # Extract best model parameters for plotting
                model_type = bootstrap_result.get("model_type", "single_cp")
                slope1 = bootstrap_result.get("slope1", 0)
                val_cp1 = bootstrap_result.get("val_cp1", 0)
                val_before_jump = bootstrap_result.get("val_before_jump", 0)
                val_after_jump = bootstrap_result.get("val_after_jump", 0)
                val_cp2 = bootstrap_result.get("val_cp2", val_cp1)
                slope2 = bootstrap_result.get("slope2", 0)
                cp1 = bootstrap_result.get("cp1_position", 0)
                cp2 = bootstrap_result.get("cp2_position", cp1)
                r_squared = bootstrap_result.get("r_squared_staged", 0)

                print(f"\n   🎨 PLOTTING MODEL: {model_type}")
                print(f"   🎨 CP positions: cp1={cp1:.2f}, cp2={cp2:.2f}")
                print(f"   🎨 Slopes: slope1={slope1:.6f}, slope2={slope2:.6f}")
                print(
                    f"   🎨 Values: val_cp1={val_cp1:.2f}, val_cp2={val_cp2:.2f}, val_before_jump={val_before_jump:.2f}, val_after_jump={val_after_jump:.2f}"
                )

                # Generate predictions for best model using proper integration
                if model_type == "single_cp":
                    predictions = get_final_predictions_for_single_change_point(
                        test_ages,
                        survival_data,
                        slope1,
                        val_cp1,
                        val_before_jump,
                        val_after_jump,
                        slope2,
                        cp1,
                    )
                else:  # two_cp
                    predictions = get_final_predictions_for_change_point(
                        test_ages,
                        survival_data,
                        slope1,
                        val_cp1,
                        val_before_jump,
                        val_after_jump,
                        val_cp2,
                        slope2,
                        cp1,
                        cp2,
                    )

                # Extract single CP model parameters for plotting
                single_cp_params = None
                single_cp_predictions = None
                single_cp_r_squared = None
                single_cp_position = None

                if bootstrap_result.get("single_cp_slope1") is not None:
                    single_cp_slope1 = bootstrap_result.get("single_cp_slope1", 0)
                    single_cp_val_cp = bootstrap_result.get("single_cp_val_cp", 0)
                    single_cp_val_before_jump = bootstrap_result.get(
                        "single_cp_val_before_jump", 0
                    )
                    single_cp_val_after_jump = bootstrap_result.get(
                        "single_cp_val_after_jump", 0
                    )
                    single_cp_slope2 = bootstrap_result.get("single_cp_slope2", 0)
                    single_cp_position = bootstrap_result.get("single_cp_position", 0)
                    single_cp_r_squared = bootstrap_result.get("single_cp_r_squared", 0)

                    single_cp_params = [
                        single_cp_slope1,
                        single_cp_val_cp,
                        single_cp_val_before_jump,
                        single_cp_val_after_jump,
                        single_cp_slope2,
                        single_cp_position,  # Include optimized cp_position
                    ]

                    # Generate single CP predictions using proper integration
                    single_cp_predictions = (
                        get_final_predictions_for_single_change_point(
                            test_ages,
                            survival_data,
                            single_cp_slope1,
                            single_cp_val_cp,
                            single_cp_val_before_jump,
                            single_cp_val_after_jump,
                            single_cp_slope2,
                            single_cp_position,
                        )
                    )

                # Extract two CP model parameters for plotting
                best_two_cp_params = None
                two_cp_cp1 = cp1  # Default to main model cp1
                two_cp_cp2 = cp2  # Default to main model cp2

                if bootstrap_result.get("two_cp_slope1") is not None:
                    two_cp_slope1 = bootstrap_result.get("two_cp_slope1", 0)
                    two_cp_val_cp1 = bootstrap_result.get("two_cp_val_cp1", 0)
                    two_cp_val_before_jump = bootstrap_result.get(
                        "two_cp_val_before_jump", 0
                    )
                    two_cp_val_after_jump = bootstrap_result.get(
                        "two_cp_val_after_jump", 0
                    )
                    two_cp_val_cp2 = bootstrap_result.get("two_cp_val_cp2", 0)
                    two_cp_slope2 = bootstrap_result.get("two_cp_slope2", 0)
                    two_cp_cp1 = bootstrap_result.get("two_cp_cp1_position", 0)
                    two_cp_cp2 = bootstrap_result.get("two_cp_cp2_position", 0)

                    best_two_cp_params = [
                        two_cp_slope1,
                        two_cp_val_cp1,
                        two_cp_val_before_jump,
                        two_cp_val_after_jump,
                        two_cp_val_cp2,
                        two_cp_slope2,
                        two_cp_cp1,  # Include optimized cp1 position
                        two_cp_cp2,  # Include optimized cp2 position
                    ]

                # Create plot title with bootstrap info
                plot_title = (
                    f"{test_name} (Bootstrap: {bootstrap_file.replace('.csv', '')})"
                )

                # Create the plot
                plot_filename = f"{bootstrap_file.replace('.csv', '')}.pdf"
                plot_path = os.path.join(plots_folder, plot_filename)

                # Generate linear model predictions using proper integration
                linear_predictions = get_final_predictions_for_linear(
                    test_ages,
                    survival_data,
                    linear_model["a"],
                    linear_model["b"],
                    linear_model["c"],
                    linear_model["d"],
                )

                # Update linear model with integrated predictions for plotting
                linear_model_for_plot = linear_model.copy()
                linear_model_for_plot["fitted_values"] = linear_predictions

                # Extract new model parameters from bootstrap_result
                poly4_model_for_plot = None
                if bootstrap_result.get("poly4_r_squared") is not None:
                    poly4_model_for_plot = {
                        "a0": bootstrap_result.get("poly4_a0", 0),
                        "a1": bootstrap_result.get("poly4_a1", 0),
                        "a2": bootstrap_result.get("poly4_a2", 0),
                        "a3": bootstrap_result.get("poly4_a3", 0),
                        "a4": bootstrap_result.get("poly4_a4", 0),
                        "r_squared": bootstrap_result.get("poly4_r_squared", 0),
                        "rss": bootstrap_result.get("poly4_pure_rss", 0),
                        "fitted_values": np.array([
                            pred_age_polynomial4(age, survival_data, 
                                               bootstrap_result.get("poly4_a0", 0),
                                               bootstrap_result.get("poly4_a1", 0),
                                               bootstrap_result.get("poly4_a2", 0),
                                               bootstrap_result.get("poly4_a3", 0),
                                               bootstrap_result.get("poly4_a4", 0))
                            for age in test_ages
                        ]),
                    }

                pwl_model_for_plot = None
                if bootstrap_result.get("pwl_r_squared") is not None:
                    pwl_model_for_plot = {
                        "v_m12_5": bootstrap_result.get("pwl_v_m12_5", 0),
                        "v_0": bootstrap_result.get("pwl_v_0", 0),
                        "v_12_5": bootstrap_result.get("pwl_v_12_5", 0),
                        "slope_before": bootstrap_result.get("pwl_slope_before", 0),
                        "slope_after": bootstrap_result.get("pwl_slope_after", 0),
                        "r_squared": bootstrap_result.get("pwl_r_squared", 0),
                        "rss": bootstrap_result.get("pwl_pure_rss", 0),
                        "fitted_values": np.array([
                            pred_age_piecewise_linear_continuous(age, survival_data,
                                                                bootstrap_result.get("pwl_v_m12_5", 0),
                                                                bootstrap_result.get("pwl_v_0", 0),
                                                                bootstrap_result.get("pwl_v_12_5", 0),
                                                                bootstrap_result.get("pwl_slope_before", 0),
                                                                bootstrap_result.get("pwl_slope_after", 0))
                            for age in test_ages
                        ]),
                    }

                exp_lin_model_for_plot = None
                if bootstrap_result.get("exp_lin_r_squared") is not None:
                    exp_lin_model_for_plot = {
                        "exp_a": bootstrap_result.get("exp_lin_exp_a", 0),
                        "exp_b": bootstrap_result.get("exp_lin_exp_b", 0),
                        "val_before_jump": bootstrap_result.get("exp_lin_val_before_jump", 0),
                        "val_after_jump": bootstrap_result.get("exp_lin_val_after_jump", 0),
                        "slope_post": bootstrap_result.get("exp_lin_slope_post", 0),
                        "r_squared": bootstrap_result.get("exp_lin_r_squared", 0),
                        "rss": bootstrap_result.get("exp_lin_pure_rss", 0),
                        "fitted_values": np.array([
                            pred_age_exp_linear(age, survival_data,
                                               bootstrap_result.get("exp_lin_exp_a", 0),
                                               bootstrap_result.get("exp_lin_exp_b", 0),
                                               bootstrap_result.get("exp_lin_val_before_jump", 0),
                                               bootstrap_result.get("exp_lin_val_after_jump", 0),
                                               bootstrap_result.get("exp_lin_slope_post", 0))
                            for age in test_ages
                        ]),
                    }

                lin_exp_model_for_plot = None
                if bootstrap_result.get("lin_exp_r_squared") is not None:
                    lin_exp_model_for_plot = {
                        "slope_pre": bootstrap_result.get("lin_exp_slope_pre", 0),
                        "val_before_jump": bootstrap_result.get("lin_exp_val_before_jump", 0),
                        "val_after_jump": bootstrap_result.get("lin_exp_val_after_jump", 0),
                        "exp_a": bootstrap_result.get("lin_exp_exp_a", 0),
                        "exp_b": bootstrap_result.get("lin_exp_exp_b", 0),
                        "r_squared": bootstrap_result.get("lin_exp_r_squared", 0),
                        "rss": bootstrap_result.get("lin_exp_pure_rss", 0),
                        "fitted_values": np.array([
                            pred_age_linear_exp(age, survival_data,
                                               bootstrap_result.get("lin_exp_slope_pre", 0),
                                               bootstrap_result.get("lin_exp_val_before_jump", 0),
                                               bootstrap_result.get("lin_exp_val_after_jump", 0),
                                               bootstrap_result.get("lin_exp_exp_a", 0),
                                               bootstrap_result.get("lin_exp_exp_b", 0))
                            for age in test_ages
                        ]),
                    }

                poly3_model_for_plot = None
                if bootstrap_result.get("poly3_r_squared") is not None:
                    poly3_model_for_plot = {
                        "a0": bootstrap_result.get("poly3_a0", 0),
                        "a1": bootstrap_result.get("poly3_a1", 0),
                        "a2": bootstrap_result.get("poly3_a2", 0),
                        "a3": bootstrap_result.get("poly3_a3", 0),
                        "r_squared": bootstrap_result.get("poly3_r_squared", 0),
                        "rss": bootstrap_result.get("poly3_pure_rss", 0),
                        "fitted_values": np.array([
                            pred_age_polynomial3(age, survival_data, 
                                               bootstrap_result.get("poly3_a0", 0),
                                               bootstrap_result.get("poly3_a1", 0),
                                               bootstrap_result.get("poly3_a2", 0),
                                               bootstrap_result.get("poly3_a3", 0))
                            for age in test_ages
                        ]),
                    }

                pwl2_model_for_plot = None
                if bootstrap_result.get("pwl2_r_squared") is not None:
                    pwl2_model_for_plot = {
                        "v_m8_33": bootstrap_result.get("pwl2_v_m8_33", 0),
                        "v_8_33": bootstrap_result.get("pwl2_v_8_33", 0),
                        "slope_before": bootstrap_result.get("pwl2_slope_before", 0),
                        "slope_after": bootstrap_result.get("pwl2_slope_after", 0),
                        "r_squared": bootstrap_result.get("pwl2_r_squared", 0),
                        "rss": bootstrap_result.get("pwl2_pure_rss", 0),
                        "fitted_values": np.array([
                            pred_age_piecewise_linear_2dots(age, survival_data,
                                                           bootstrap_result.get("pwl2_v_m8_33", 0),
                                                           bootstrap_result.get("pwl2_v_8_33", 0),
                                                           bootstrap_result.get("pwl2_slope_before", 0),
                                                           bootstrap_result.get("pwl2_slope_after", 0))
                            for age in test_ages
                        ]),
                    }

                # Call create_plots with all model parameters
                create_plots(
                    test_ages,
                    test_values,
                    predictions,
                    linear_model_for_plot,
                    slope1,
                    val_cp1,
                    val_before_jump,
                    val_after_jump,
                    val_cp2,
                    slope2,
                    cp1,
                    cp2,
                    r_squared,
                    plot_title,
                    plots_folder,  # output_dir
                    test_sds,
                    test_ns,
                    single_cp_params,  # single_cp_params
                    single_cp_predictions,  # single_cp_predictions
                    single_cp_r_squared,  # single_cp_r_squared
                    single_cp_position,  # single_cp_position
                    best_two_cp_params,  # best_two_cp_params
                    survival_data,
                    two_cp_cp1,  # best_two_cp_cp1 (optimized positions from two CP model if available)
                    two_cp_cp2,  # best_two_cp_cp2
                    data_is_logit_transformed=DEFAULT_CONFIG[
                        "USE_LOGIT_TRANSFORM"
                    ],  # Use global flag to determine error bar type
                    sigmoid_model=sigmoid_model,  # Add sigmoid model
                    poly3_model=poly3_model_for_plot,  # Add polynomial3 model
                    poly4_model=poly4_model_for_plot,  # Add polynomial4 model
                    pwl2_model=pwl2_model_for_plot,  # Add piecewise linear 2-dots model
                    pwl_model=pwl_model_for_plot,  # Add piecewise linear 3-dots model (PWL3)
                    exp_lin_model=exp_lin_model_for_plot,  # Add exp-linear model
                    lin_exp_model=lin_exp_model_for_plot,  # Add linear-exp model
                )

                print(f"   📊 Created plot for bootstrap {bootstrap_file}")

            except Exception as plot_error:
                print(
                    f"   ⚠️  Failed to create plot for bootstrap {bootstrap_file}: {plot_error}"
                )

            return bootstrap_result

    except Exception as e:
        print(f"   ❌ Error in bootstrap {bootstrap_file}: {str(e)}")
        import traceback

        print("   📋 Full traceback:")
        traceback.print_exc()

        # Return minimal result with linear model data even on failure
        bootstrap_result = {
            "bootstrap_file": bootstrap_file,
            "analysis_status": "failed",
            "error": str(e),
            # Linear model data (should be available)
            "linear_a": linear_model.get("a", None) if linear_model else None,
            "linear_b": linear_model.get("b", None) if linear_model else None,
            "linear_c": linear_model.get("c", None) if linear_model else None,
            "linear_d": linear_model.get("d", None) if linear_model else None,
            "linear_bic": None,
            "linear_aic": None,
            "linear_pure_rss": None,
            # Failed model results
            "single_cp_bic": None,
            "single_cp_aic": None,
            "single_cp_pure_rss": None,
            "two_cp_bic": None,
            "two_cp_aic": None,
            "two_cp_pure_rss": None,
        }

        # Try to calculate linear model metrics even on failure
        if linear_model and all(k in linear_model for k in ["a", "b", "c", "d"]):
            try:
                linear_rss_for_bic = calculate_pure_rss(
                    test_ages,
                    test_values,
                    survival_data,
                    linear_a=linear_model["a"],
                    linear_b=linear_model["b"],
                    linear_c=linear_model["c"],
                    linear_d=linear_model["d"],
                    test_sds=test_sds,
                    test_ns=None,  # Don't transform - use test_sds directly
                )
                n_data = len(test_values)
                linear_n_params = 4
                bootstrap_result["linear_pure_rss"] = linear_rss_for_bic
                bootstrap_result["linear_aic"] = (
                    n_data * np.log(linear_rss_for_bic / n_data) + 2 * linear_n_params
                )
                bootstrap_result["linear_bic"] = (
                    n_data * np.log(linear_rss_for_bic / n_data)
                    + np.log(n_data) * linear_n_params
                )
                print(f"   ✅ Calculated linear model metrics for failed bootstrap")
            except Exception as calc_e:
                print(f"   ⚠️ Could not calculate linear metrics: {calc_e}")

        return bootstrap_result


def main():
    """Main function with command line argument parsing"""
    parser = argparse.ArgumentParser(
        description="Run change point analysis for a single test"
    )
    parser.add_argument("test_code", help="Test code to analyze")
    parser.add_argument("test_name", help="Test name")
    parser.add_argument("system_name", help="System name")
    parser.add_argument(
        "--output-dir", default="change_point_outputs", help="Output directory"
    )
    parser.add_argument(
        "--data-file",
        default="healthy_menopause_aggregated_data.csv",
        help="Input data file",
    )
    parser.add_argument(
        "--survival-file",
        default="example_survival_curve.csv",
        help="Survival curve file",
    )
    parser.add_argument(
        "--n-workers",
        type=int,
        default=32,
        help="Number of parallel workers for optimization (default: 32)",
    )
    parser.add_argument(
        "--n-combinations",
        type=int,
        default=100,
        help="Number of change point combinations to test (default: 100)",
    )
    parser.add_argument(
        "--use-cv",
        action="store_false",
        help="Use coefficient of variation (CV) for error bars and weighted R-squared (default: True)",
    )
    parser.add_argument(
        "--single-cp-only",
        action="store_true",
        help="Use only single change point models (faster, default: False)",
    )
    parser.add_argument(
        "--cp-improvement-threshold",
        type=float,
        default=0.01,
        help="Minimum R² improvement to consider two change points (default: 0.01)",
    )
    parser.add_argument(
        "--dont-randomize-combinations",
        action="store_true",
        help="Use grid-based combinations instead of random sampling (default: random sampling)",
    )
    parser.add_argument(
        "--bootstrap-folder",
        help="Folder containing bootstrap CSV files (default: None)",
    )
    parser.add_argument(
        "--min-age",
        type=float,
        help="Minimum age to include in analysis (default: None, use all ages)",
    )
    parser.add_argument(
        "--max-age",
        type=float,
        help="Maximum age to include in analysis (default: None, use all ages)",
    )
    parser.add_argument(
        "--gaussian-bandwidth",
        type=float,
        default=1.5,
        help="Gaussian smoothing bandwidth (sigma) for survival curve (default: 1.5)",
    )
    parser.add_argument(
        "--integral-dt",
        type=float,
        default=0.2,
        help="Integration step size for survival curve calculations (default: 0.2)",
    )
    parser.add_argument(
        "--single-cp-before-only",
        action="store_true",
        help="For single change point models, only try positions before menopause (default: False, try both before and after)",
    )
    parser.add_argument(
        "--disable-slope-penalization",
        action="store_true",
        help="Disable slope penalization for extreme opposite slopes (default: enabled)",
    )
    parser.add_argument(
        "--specific-cp-a",
        type=float,
        help="Test only a specific first change point position",
    )
    parser.add_argument(
        "--specific-cp-b",
        type=float,
        help="Test only a specific second change point position",
    )
    parser.add_argument(
        "--use-iterative-reweighting",
        action="store_true",
        help="Use iterative reweighted least squares instead of simple weighted least squares",
    )
    parser.add_argument(
        "--slope-penalty-weight",
        type=float,
        default=10.0,
        help="Weight for slope penalty in objective function (default: 10.0)",
    )
    parser.add_argument(
        "--slope-penalty-threshold-ratio",
        type=float,
        default=0.05,
        help="Threshold as ratio of jump magnitude for considering slope extreme (default: 0.05 = 5% of jump)",
    )
    parser.add_argument(
        "--slope-penalty-min-threshold",
        type=float,
        default=0.01,
        help="Minimum absolute slope threshold regardless of jump size (default: 0.01)",
    )
    parser.add_argument(
        "--slope-penalty-distance-threshold",
        type=float,
        default=8.0,
        help="Only penalize slopes if change point is within this distance of menopause (default: 5.0 years)",
    )
    parser.add_argument(
        "--disable-significance-testing",
        action="store_true",
        help="Disable statistical significance testing (default: enabled)",
    )
    parser.add_argument(
        "--significance-alpha",
        type=float,
        default=0.05,
        help="Alpha level for significance tests (default: 0.05)",
    )
    parser.add_argument(
        "--min-r-squared-improvement",
        type=float,
        default=0.001,
        help="Minimum R² improvement for considering significant (default: 0.001)",
    )
    parser.add_argument(
        "--min-absolute-r-squared",
        type=float,
        default=0.01,
        help="Minimum absolute R² for the change point model (default: 0.01)",
    )
    parser.add_argument(
        "--override-bootstrap",
        action="store_true",
        help="Override existing bootstrap results",
    )
    parser.add_argument(
        "--use-logit-transform",
        action="store_true",
        help="Apply logit transformation to data (for probabilities/percentages 0-100)",
    )
    parser.add_argument(
        "--only-basic-models",
        action="store_true",
        help="Only fit linear and sigmoid models (skip all change point models and spline)",
    )
    parser.add_argument(
        "--profile",
        action="store_true",
        help="Enable profiling and save results to file"
    )
    parser.add_argument(
        "--profile-workers",
        action="store_true",
        help="Enable profiling for each worker process (bootstrap analysis)"
    )

    args = parser.parse_args()

    # Update configuration
    DEFAULT_CONFIG["INFILE"] = args.data_file
    DEFAULT_CONFIG["SURV_FILE"] = args.survival_file
    DEFAULT_CONFIG["OUTPUT_DIR"] = args.output_dir
    DEFAULT_CONFIG["GAUSSIAN_BANDWIDTH"] = args.gaussian_bandwidth
    DEFAULT_CONFIG["INTEGRAL_DT"] = args.integral_dt
    DEFAULT_CONFIG["N_WORKERS"] = args.n_workers
    DEFAULT_CONFIG["N_CP_COMBINATIONS"] = args.n_combinations
    DEFAULT_CONFIG["USE_CV"] = args.use_cv
    DEFAULT_CONFIG["SINGLE_CP_ONLY"] = args.single_cp_only
    DEFAULT_CONFIG["CP_IMPROVEMENT_THRESHOLD"] = args.cp_improvement_threshold
    DEFAULT_CONFIG["RANDOMIZE_COMBINATIONS"] = not args.dont_randomize_combinations
    DEFAULT_CONFIG["SINGLE_CP_BEFORE_ONLY"] = args.single_cp_before_only

    # Check for only-basic-models flag (from args or environment variable)
    only_basic_models = args.only_basic_models or os.environ.get("ONLY_BASIC_MODELS", "").lower() in ("1", "true", "yes")
    DEFAULT_CONFIG["ONLY_BASIC_MODELS"] = only_basic_models
    
    # Update slope penalization configuration
    DEFAULT_CONFIG["ENABLE_SLOPE_PENALIZATION"] = not args.disable_slope_penalization
    DEFAULT_CONFIG["SLOPE_PENALTY_WEIGHT"] = args.slope_penalty_weight
    DEFAULT_CONFIG["SLOPE_PENALTY_THRESHOLD_RATIO"] = args.slope_penalty_threshold_ratio
    DEFAULT_CONFIG["SLOPE_PENALTY_MIN_THRESHOLD"] = args.slope_penalty_min_threshold
    DEFAULT_CONFIG["SLOPE_PENALTY_DISTANCE_THRESHOLD"] = (
        args.slope_penalty_distance_threshold
    )

    # Update significance testing configuration
    DEFAULT_CONFIG["ENABLE_SIGNIFICANCE_TESTING"] = (
        not args.disable_significance_testing
    )
    DEFAULT_CONFIG["SIGNIFICANCE_ALPHA"] = args.significance_alpha
    DEFAULT_CONFIG["MIN_R_SQUARED_IMPROVEMENT"] = args.min_r_squared_improvement
    DEFAULT_CONFIG["MIN_ABSOLUTE_R_SQUARED"] = args.min_absolute_r_squared

    # Add age limits to configuration
    if args.min_age is not None:
        DEFAULT_CONFIG["MIN_AGE"] = args.min_age
    if args.max_age is not None:
        DEFAULT_CONFIG["MAX_AGE"] = args.max_age

    # Update logit transformation setting
    DEFAULT_CONFIG["USE_LOGIT_TRANSFORM"] = args.use_logit_transform
    
    # Update worker profiling setting
    DEFAULT_CONFIG["PROFILE_WORKERS"] = args.profile_workers

    # Store specific change points if provided
    specific_cp_a = args.specific_cp_a
    specific_cp_b = args.specific_cp_b
    if specific_cp_a is not None or specific_cp_b is not None:
        print(
            f"⚠️  Using specific change points: CP1={specific_cp_a}, CP2={specific_cp_b}"
        )

    # Store iterative reweighting flag
    use_iterative_reweighting = args.use_iterative_reweighting
    if use_iterative_reweighting:
        print("⚠️  Using iterative reweighted least squares!")

    # Run analysis
    result = run_single_test_analysis(
        args.test_code,
        args.test_name,
        args.system_name,
        args.output_dir,
        args.bootstrap_folder,
        specific_cp_a=specific_cp_a,
        specific_cp_b=specific_cp_b,
        use_iterative_reweighting=use_iterative_reweighting,
        override_bootstrap=args.override_bootstrap,
    )

    if result:
        # Prepare result for CSV (handle bootstrap results properly)
        csv_result = result.copy()

        # If bootstrap results exist, save them separately and add summary stats
        if "bootstrap_results" in result and result["bootstrap_results"]:
            bootstrap_results = result["bootstrap_results"]

            # Calculate summary statistics for bootstrap results
            bootstrap_summary = {}

            # Collect all change point positions from both models
            all_cp1_positions = []
            all_cp2_positions = []
            all_slope1_values = []
            all_slope2_values = []
            all_jump_magnitudes = []
            all_r_squared_values = []
            all_improvement_values = []

            # Determine which model type was selected for this test
            selected_model_type = result.get(
                "model_type", "two_cp"
            )  # Default to two_cp

            for br in bootstrap_results:
                # Only collect statistics from the selected model type
                if selected_model_type == "single_cp":
                    # Single CP model
                    if "single_cp_position" in br and not np.isnan(
                        br["single_cp_position"]
                    ):
                        single_cp_pos = br["single_cp_position"]
                        if single_cp_pos < 0:
                            all_cp1_positions.append(single_cp_pos)
                        else:
                            all_cp2_positions.append(single_cp_pos)

                    if "single_cp_slope1" in br:
                        all_slope1_values.append(br["single_cp_slope1"])
                    if "single_cp_slope2" in br:
                        all_slope2_values.append(br["single_cp_slope2"])
                    if "single_cp_jump_magnitude" in br:
                        all_jump_magnitudes.append(br["single_cp_jump_magnitude"])
                    if "single_cp_r_squared" in br:
                        all_r_squared_values.append(br["single_cp_r_squared"])
                    if "single_cp_improvement" in br:
                        all_improvement_values.append(br["single_cp_improvement"])

                elif selected_model_type == "two_cp":
                    # Two CP model
                    if "two_cp_cp1_position" in br and not np.isnan(
                        br["two_cp_cp1_position"]
                    ):
                        two_cp_cp1_pos = br["two_cp_cp1_position"]
                        if two_cp_cp1_pos < 0:
                            all_cp1_positions.append(two_cp_cp1_pos)

                    if "two_cp_cp2_position" in br and not np.isnan(
                        br["two_cp_cp2_position"]
                    ):
                        two_cp_cp2_pos = br["two_cp_cp2_position"]
                        if two_cp_cp2_pos > 0:
                            all_cp2_positions.append(two_cp_cp2_pos)

                    if "two_cp_slope1" in br:
                        all_slope1_values.append(br["two_cp_slope1"])
                    if "two_cp_slope2" in br:
                        all_slope2_values.append(br["two_cp_slope2"])
                    if "two_cp_jump_magnitude" in br:
                        all_jump_magnitudes.append(br["two_cp_jump_magnitude"])
                    if "two_cp_r_squared" in br:
                        all_r_squared_values.append(br["two_cp_r_squared"])
                    if "two_cp_improvement" in br:
                        all_improvement_values.append(br["two_cp_improvement"])

                # Backward compatibility: also collect from old format if no specific model type
                if "cp1_position" in br and not np.isnan(br["cp1_position"]):
                    cp1_pos = br["cp1_position"]
                    if cp1_pos < 0:
                        all_cp1_positions.append(cp1_pos)
                    else:
                        all_cp2_positions.append(cp1_pos)

                if "cp2_position" in br and not np.isnan(br["cp2_position"]):
                    cp2_pos = br["cp2_position"]
                    if cp2_pos > 0:
                        all_cp2_positions.append(cp2_pos)

                if "slope1" in br:
                    all_slope1_values.append(br["slope1"])
                if "slope2" in br:
                    all_slope2_values.append(br["slope2"])
                if "jump_magnitude" in br:
                    all_jump_magnitudes.append(br["jump_magnitude"])
                if "r_squared_staged" in br:
                    all_r_squared_values.append(br["r_squared_staged"])
                if "improvement" in br:
                    all_improvement_values.append(br["improvement"])

            # Calculate combined statistics
            if all_cp1_positions:
                bootstrap_summary["cp1_position_median"] = np.median(all_cp1_positions)
                bootstrap_summary["cp1_position_min"] = np.min(all_cp1_positions)
                bootstrap_summary["cp1_position_max"] = np.max(all_cp1_positions)
                bootstrap_summary["cp1_position_std"] = np.std(all_cp1_positions)

            if all_cp2_positions:
                bootstrap_summary["cp2_position_median"] = np.median(all_cp2_positions)
                bootstrap_summary["cp2_position_min"] = np.min(all_cp2_positions)
                bootstrap_summary["cp2_position_max"] = np.max(all_cp2_positions)
                bootstrap_summary["cp2_position_std"] = np.std(all_cp2_positions)

            if all_slope1_values:
                bootstrap_summary["slope1_median"] = np.median(all_slope1_values)
                bootstrap_summary["slope1_min"] = np.min(all_slope1_values)
                bootstrap_summary["slope1_max"] = np.max(all_slope1_values)
                bootstrap_summary["slope1_std"] = np.std(all_slope1_values)

            if all_slope2_values:
                bootstrap_summary["slope2_median"] = np.median(all_slope2_values)
                bootstrap_summary["slope2_min"] = np.min(all_slope2_values)
                bootstrap_summary["slope2_max"] = np.max(all_slope2_values)
                bootstrap_summary["slope2_std"] = np.std(all_slope2_values)

            if all_jump_magnitudes:
                bootstrap_summary["jump_magnitude_median"] = np.median(
                    all_jump_magnitudes
                )
                bootstrap_summary["jump_magnitude_min"] = np.min(all_jump_magnitudes)
                bootstrap_summary["jump_magnitude_max"] = np.max(all_jump_magnitudes)
                bootstrap_summary["jump_magnitude_std"] = np.std(all_jump_magnitudes)

            if all_r_squared_values:
                bootstrap_summary["r_squared_staged_median"] = np.median(
                    all_r_squared_values
                )
                bootstrap_summary["r_squared_staged_min"] = np.min(all_r_squared_values)
                bootstrap_summary["r_squared_staged_max"] = np.max(all_r_squared_values)
                bootstrap_summary["r_squared_staged_std"] = np.std(all_r_squared_values)

            if all_improvement_values:
                bootstrap_summary["improvement_median"] = np.median(
                    all_improvement_values
                )
                bootstrap_summary["improvement_min"] = np.min(all_improvement_values)
                bootstrap_summary["improvement_max"] = np.max(all_improvement_values)
                bootstrap_summary["improvement_std"] = np.std(all_improvement_values)

            # Also calculate separate statistics for each model type for backward compatibility
            for key in [
                "cp1_position",
                "cp2_position",
                "slope1",
                "val_cp1",
                "val_before_jump",
                "val_after_jump",
                "val_cp2",
                "slope2",
                "jump_magnitude",
                "r_squared_staged",
                "improvement",
                "loss",
            ]:
                values = [br[key] for br in bootstrap_results if key in br]
                if values:
                    bootstrap_summary[f"{key}_median"] = np.median(values)
                    bootstrap_summary[f"{key}_min"] = np.min(values)
                    bootstrap_summary[f"{key}_max"] = np.max(values)
                    bootstrap_summary[f"{key}_std"] = np.std(values)

            # Add bootstrap summary to CSV result
            csv_result.update(bootstrap_summary)

        # Remove the bootstrap_results from CSV (too large for CSV)
        if "bootstrap_results" in csv_result:
            del csv_result["bootstrap_results"]

        # Write to all_data.csv with proper file locking (only once at the end)
        all_data_path = os.path.join(args.output_dir, "all_data.csv")

        # Use file locking to prevent multiple jobs from writing to all_data.csv simultaneously
        import fcntl
        import time

        lock_file_path = all_data_path + ".lock"
        max_retries = 10
        retry_delay = 2

        for attempt in range(max_retries):
            try:
                with open(lock_file_path, "w") as lock_file:
                    # Try to acquire an exclusive lock
                    fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)

                    try:
                        if os.path.exists(all_data_path):
                            # Load existing data
                            existing_df = pd.read_csv(all_data_path)

                            # Check if this test already exists
                            test_mask = (
                                (existing_df["test_code"] == csv_result["test_code"])
                                & (existing_df["test_name"] == csv_result["test_name"])
                                & (existing_df["system"] == csv_result["system"])
                            )

                            if test_mask.any():
                                # Replace existing entry
                                existing_df.loc[test_mask] = pd.Series(csv_result)
                                print(
                                    f"📊 Updated existing entry for {csv_result['test_name']}"
                                )
                            else:
                                # Append new entry
                                new_df = pd.DataFrame([csv_result])
                                existing_df = pd.concat(
                                    [existing_df, new_df], ignore_index=True
                                )
                                print(
                                    f"📊 Added new entry for {csv_result['test_name']}"
                                )

                            combined_df = existing_df
                        else:
                            # Create new file
                            combined_df = pd.DataFrame([csv_result])
                            print(
                                f"📊 Created new all_data.csv with {csv_result['test_name']}"
                            )

                        combined_df.to_csv(all_data_path, index=False)
                        print(
                            f"📊 Results saved to {all_data_path} ({len(combined_df)} total tests)"
                        )
                        break  # Success, exit retry loop

                    finally:
                        # Release the lock
                        fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)

            except (IOError, OSError) as e:
                if attempt < max_retries - 1:
                    print(
                        f"⚠️  File lock busy, retrying in {retry_delay} seconds... (attempt {attempt + 1}/{max_retries})"
                    )
                    time.sleep(retry_delay)
                else:
                    print(
                        f"❌ Failed to acquire file lock after {max_retries} attempts: {e}"
                    )
                    # Fallback: write to a temporary file and let the user merge later
                    temp_path = all_data_path + f".{os.getpid()}.tmp"
                    temp_df = pd.DataFrame([csv_result])
                    temp_df.to_csv(temp_path, index=False)
                    print(f"📊 Results saved to temporary file: {temp_path}")
                    return 0

        return 0
    else:
        return 1


if __name__ == "__main__":
    # Check if profiling is enabled
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", action="store_true")
    
    # Quick parse just to check for --profile flag
    profile_args, _ = parser.parse_known_args()
    
    if profile_args.profile:
        import cProfile
        import pstats
        
        print("🔍 Profiling enabled - performance data will be saved")
        
        profiler = cProfile.Profile()
        profiler.enable()
        exit_code = main()
        profiler.disable()
        
        # Get output directory from the actual arguments
        full_parser = argparse.ArgumentParser()
        full_parser.add_argument("test_code")
        full_parser.add_argument("test_name")
        full_parser.add_argument("system_name")
        full_parser.add_argument("--output-dir", default="change_point_outputs")
        full_args, _ = full_parser.parse_known_args()
        
        # Save profile to file
        profile_file = os.path.join(full_args.output_dir, f"profile_{full_args.test_code}.prof")
        profiler.dump_stats(profile_file)
        print(f"📊 Profile data saved to: {profile_file}")
        
        # Print top 20 time-consuming functions
        print("\n" + "="*80)
        print("TOP 20 TIME-CONSUMING FUNCTIONS (by cumulative time)")
        print("="*80)
        stats = pstats.Stats(profiler)
        stats.sort_stats('cumulative')
        stats.print_stats(20)
        
        sys.exit(exit_code)
    else:
        sys.exit(main())
