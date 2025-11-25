#!/usr/bin/env python3
"""
Cross-validation script for linear and sigmoid models
Tests models fitted on one bootstrap against different bootstrap raw data
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import time
from pathlib import Path

# Import functions from change_point_single_test.py
from src.chemfarm_analysis.change_point_single_test import (
    get_final_predictions_for_linear,
    pred_age_sigmoid,
    apply_gaussian_smoothing_continuous,
    DEFAULT_CONFIG,
)


def load_bootstrap_raw_data(bootstrap_folder, test_code, bootstrap_file):
    """Load raw data from a specific bootstrap file"""
    bootstrap_path = os.path.join(bootstrap_folder, bootstrap_file)
    
    if not os.path.exists(bootstrap_path):
        raise FileNotFoundError(f"Bootstrap file not found: {bootstrap_path}")
    
    data = pd.read_csv(bootstrap_path)
    
    # Determine the correct prefix for the test code
    if test_code in ["marker.BMI", "marker.bp.high.no_meds", "marker.bp.low.no_meds"]:
        track_name = test_code
    else:
        track_name = f"lab.{test_code}"
    
    test_data = data[data["test_track"] == track_name]
    
    if len(test_data) == 0:
        raise ValueError(f"No data found for test code {test_code} in {bootstrap_file}")
    
    # Get age columns (exclude _sd and _n columns)
    age_columns = [col for col in data.columns[9:] 
                   if not col.endswith("_sd") and not col.endswith("_n")]
    
    test_ages = np.array([float(col) for col in age_columns])
    test_values = np.array([test_data.iloc[0][col] for col in age_columns])
    
    # Get standard deviations
    test_sds = np.zeros_like(test_values)
    for i, age_col in enumerate(age_columns):
        sd_col = age_col + "_sd"
        if sd_col in data.columns:
            test_sds[i] = test_data.iloc[0][sd_col]
    
    # Get sample sizes
    test_ns = np.ones_like(test_values, dtype=int)
    for i, age_col in enumerate(age_columns):
        n_col = age_col + "_n"
        if n_col in data.columns:
            test_ns[i] = test_data.iloc[0][n_col]
    
    # Remove NaN values
    valid_mask = ~np.isnan(test_values)
    test_ages = test_ages[valid_mask]
    test_values = test_values[valid_mask]
    test_sds = test_sds[valid_mask]
    test_ns = test_ns[valid_mask]
    
    return test_ages, test_values, test_sds, test_ns


def calculate_cross_validation_metrics(predictions, actual_values, weights=None):
    """Calculate RSS and RMSE"""
    residuals = actual_values - predictions
    
    if weights is not None and np.any(weights > 0):
        # Weighted RSS
        rss = np.sum(weights * residuals**2)
        # Weighted RMSE
        rmse = np.sqrt(np.sum(weights * residuals**2) / np.sum(weights))
    else:
        # Unweighted RSS
        rss = np.sum(residuals**2)
        # Unweighted RMSE
        rmse = np.sqrt(np.mean(residuals**2))
    
    return rss, rmse


def cross_validate_single_test(test_code, test_name, system_name, output_base_dir, 
                                bootstrap_folder, survival_file):
    """
    Cross-validate models for a single test
    """
    print(f"\n{'='*80}")
    print(f"🔬 CROSS-VALIDATING {test_name} (Code: {test_code})")
    print(f"{'='*80}")
    
    # Find test directory
    system_dir = os.path.join(output_base_dir, system_name)
    test_dir = os.path.join(system_dir, test_name.replace(" ", "_").replace("/", "_"))
    
    if not os.path.exists(test_dir):
        print(f"❌ Test directory not found: {test_dir}")
        return None
    
    # Find bootstrap results directory
    bootstrap_results_dir = os.path.join(test_dir, "bootstrap_results")
    
    if not os.path.exists(bootstrap_results_dir):
        print(f"❌ Bootstrap results directory not found: {bootstrap_results_dir}")
        return None
    
    # Create cross-validation output directory
    cv_output_dir = os.path.join(test_dir, "cross_validation_results")
    os.makedirs(cv_output_dir, exist_ok=True)
    
    print(f"📁 Output directory: {cv_output_dir}")
    
    # Load survival curve with smoothing
    print("📊 Loading survival curve...")
    survival_data_orig = pd.read_csv(survival_file)
    survival_data = apply_gaussian_smoothing_continuous(
        survival_data_orig,
        bandwidth=DEFAULT_CONFIG["GAUSSIAN_BANDWIDTH"],
        verbose=False
    )
    
    # Find all bootstrap result files
    result_files = sorted([f for f in os.listdir(bootstrap_results_dir) 
                          if f.endswith("_result.csv")])
    
    if len(result_files) == 0:
        print(f"❌ No bootstrap result files found in {bootstrap_results_dir}")
        return None
    
    print(f"📋 Found {len(result_files)} bootstrap result files")
    
    # Cross-validation results
    cv_results = []
    
    # Process each bootstrap (train on i, test on i+1, wrapping around)
    for i, result_file in enumerate(result_files):
        # Get next bootstrap index (wrap around)
        test_idx = (i + 1) % len(result_files)
        test_result_file = result_files[test_idx]
        
        train_bootstrap = result_file.replace("_result.csv", ".csv")
        test_bootstrap = test_result_file.replace("_result.csv", ".csv")
        
        print(f"\n🔄 Train: {train_bootstrap} → Test: {test_bootstrap}")
        
        # Load training model parameters
        result_path = os.path.join(bootstrap_results_dir, result_file)
        result_df = pd.read_csv(result_path)
        
        if len(result_df) == 0:
            print(f"   ⚠️  Empty result file, skipping...")
            continue
        
        result_row = result_df.iloc[0]
        
        # Check if models exist
        has_linear = all(k in result_row and not pd.isna(result_row[k]) 
                        for k in ["linear_a", "linear_b", "linear_c", "linear_d"])
        has_sigmoid = all(k in result_row and not pd.isna(result_row[k])
                         for k in ["sigmoid_c", "sigmoid_d", "sigmoid_w", "sigmoid_h", "sigmoid_b"])
        
        if not has_linear and not has_sigmoid:
            print(f"   ⚠️  No valid models found, skipping...")
            continue
        
        # Load test data
        try:
            test_ages, test_values, test_sds, test_ns = load_bootstrap_raw_data(
                bootstrap_folder, test_code, test_bootstrap
            )
        except Exception as e:
            print(f"   ❌ Error loading test data: {e}")
            continue
        
        print(f"   📊 Test data: {len(test_values)} points")
        
        # Calculate weights
        weights = None
        if DEFAULT_CONFIG["USE_CV"] and test_sds is not None and np.any(test_sds > 0):
            weights = 1.0 / (test_sds**2)
        
        cv_row = {
            "train_bootstrap": train_bootstrap,
            "test_bootstrap": test_bootstrap,
            "n_test_points": len(test_values),
        }
        
        # Cross-validate LINEAR model
        if has_linear:
            linear_a = result_row["linear_a"]
            linear_b = result_row["linear_b"]
            linear_c = result_row["linear_c"]
            linear_d = result_row["linear_d"]
            
            # Get predictions
            linear_predictions = get_final_predictions_for_linear(
                test_ages, survival_data, linear_a, linear_b, linear_c, linear_d
            )
            
            # Calculate metrics
            linear_rss, linear_rmse = calculate_cross_validation_metrics(
                linear_predictions, test_values, weights
            )
            
            cv_row.update({
                "linear_cv_rss": linear_rss,
                "linear_cv_rmse": linear_rmse,
                "linear_a": linear_a,
                "linear_b": linear_b,
                "linear_c": linear_c,
                "linear_d": linear_d,
            })
            
            print(f"   ✅ Linear: RSS={linear_rss:.4f}, RMSE={linear_rmse:.4f}")
        else:
            cv_row.update({
                "linear_cv_rss": np.nan,
                "linear_cv_rmse": np.nan,
                "linear_a": np.nan,
                "linear_b": np.nan,
                "linear_c": np.nan,
                "linear_d": np.nan,
            })
            linear_predictions = None
        
        # Cross-validate SIGMOID model
        if has_sigmoid:
            sigmoid_c = result_row["sigmoid_c"]
            sigmoid_d = result_row["sigmoid_d"]
            sigmoid_w = result_row["sigmoid_w"]
            sigmoid_h = result_row["sigmoid_h"]
            sigmoid_b = result_row["sigmoid_b"]
            
            # Get predictions
            sigmoid_predictions = []
            for age in test_ages:
                pred = pred_age_sigmoid(age, survival_data, sigmoid_c, sigmoid_d, 
                                       sigmoid_w, sigmoid_h, sigmoid_b)
                sigmoid_predictions.append(pred)
            sigmoid_predictions = np.array(sigmoid_predictions)
            
            # Calculate metrics
            sigmoid_rss, sigmoid_rmse = calculate_cross_validation_metrics(
                sigmoid_predictions, test_values, weights
            )
            
            cv_row.update({
                "sigmoid_cv_rss": sigmoid_rss,
                "sigmoid_cv_rmse": sigmoid_rmse,
                "sigmoid_c": sigmoid_c,
                "sigmoid_d": sigmoid_d,
                "sigmoid_w": sigmoid_w,
                "sigmoid_h": sigmoid_h,
                "sigmoid_b": sigmoid_b,
            })
            
            print(f"   ✅ Sigmoid: RSS={sigmoid_rss:.4f}, RMSE={sigmoid_rmse:.4f}")
        else:
            cv_row.update({
                "sigmoid_cv_rss": np.nan,
                "sigmoid_cv_rmse": np.nan,
                "sigmoid_c": np.nan,
                "sigmoid_d": np.nan,
                "sigmoid_w": np.nan,
                "sigmoid_h": np.nan,
                "sigmoid_b": np.nan,
            })
            sigmoid_predictions = None
        
        # Store predictions for plotting
        cv_row["test_ages"] = test_ages
        cv_row["test_values"] = test_values
        cv_row["linear_predictions"] = linear_predictions
        cv_row["sigmoid_predictions"] = sigmoid_predictions
        cv_row["test_sds"] = test_sds
        
        cv_results.append(cv_row)
    
    if len(cv_results) == 0:
        print(f"❌ No cross-validation results generated")
        return None
    
    # Save CSV results (without array columns)
    csv_results = []
    for row in cv_results:
        csv_row = {k: v for k, v in row.items() 
                   if not isinstance(v, np.ndarray)}
        csv_results.append(csv_row)
    
    cv_df = pd.DataFrame(csv_results)
    csv_path = os.path.join(cv_output_dir, "cross_validation_metrics.csv")
    cv_df.to_csv(csv_path, index=False)
    print(f"\n💾 Saved metrics: {csv_path}")
    
    # Create plots
    print(f"\n📊 Creating plots...")
    create_cv_plots(cv_results, test_name, cv_output_dir)
    
    # Calculate summary statistics
    print(f"\n📈 SUMMARY STATISTICS:")
    if "linear_cv_rss" in cv_df.columns:
        linear_mean_rss = cv_df["linear_cv_rss"].mean()
        linear_mean_rmse = cv_df["linear_cv_rmse"].mean()
        print(f"   Linear:  Mean RSS={linear_mean_rss:.4f}, Mean RMSE={linear_mean_rmse:.4f}")
    
    if "sigmoid_cv_rss" in cv_df.columns:
        sigmoid_mean_rss = cv_df["sigmoid_cv_rss"].mean()
        sigmoid_mean_rmse = cv_df["sigmoid_cv_rmse"].mean()
        print(f"   Sigmoid: Mean RSS={sigmoid_mean_rss:.4f}, Mean RMSE={sigmoid_mean_rmse:.4f}")
    
    print(f"\n✅ Cross-validation complete for {test_name}")
    
    return cv_results


def create_cv_plots(cv_results, test_name, output_dir):
    """Create plots showing cross-validation fits"""
    n_results = len(cv_results)
    
    # Determine grid size (5x5 for 25 bootstraps)
    n_cols = 5
    n_rows = (n_results + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 4*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    
    for idx, result in enumerate(cv_results):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]
        
        test_ages = result["test_ages"]
        test_values = result["test_values"]
        linear_preds = result["linear_predictions"]
        sigmoid_preds = result["sigmoid_predictions"]
        test_sds = result["test_sds"]
        
        # Plot data with error bars
        if test_sds is not None:
            se = test_sds / np.sqrt(len(test_values))  # Approximate standard error
            ax.errorbar(test_ages, test_values, yerr=se, fmt='o', 
                       alpha=0.5, color='gray', capsize=3, label='Test Data')
        else:
            ax.scatter(test_ages, test_values, alpha=0.5, color='gray', 
                      s=20, label='Test Data')
        
        # Plot predictions
        if linear_preds is not None:
            ax.plot(test_ages, linear_preds, 'b-', linewidth=2, 
                   label=f'Linear (RMSE={result["linear_cv_rmse"]:.3f})')
        
        if sigmoid_preds is not None:
            ax.plot(test_ages, sigmoid_preds, 'r-', linewidth=2,
                   label=f'Sigmoid (RMSE={result["sigmoid_cv_rmse"]:.3f})')
        
        # Formatting
        train_num = result["train_bootstrap"].split("_")[1].split(".")[0]
        test_num = result["test_bootstrap"].split("_")[1].split(".")[0]
        ax.set_title(f'Train: {train_num} → Test: {test_num}', fontsize=10)
        ax.set_xlabel('Age', fontsize=8)
        ax.set_ylabel('Value', fontsize=8)
        ax.legend(fontsize=7, loc='best')
        ax.grid(True, alpha=0.3)
    
    # Hide empty subplots
    for idx in range(n_results, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].axis('off')
    
    plt.suptitle(f'Cross-Validation: {test_name}\n(Train on one bootstrap, test on next)', 
                fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    plot_path = os.path.join(output_dir, "cross_validation_fits.pdf")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"   💾 Saved plot: {plot_path}")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Cross-validate linear and sigmoid models across bootstrap samples"
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory containing test results"
    )
    parser.add_argument(
        "--bootstrap-folder",
        required=True,
        help="Folder containing bootstrap CSV files"
    )
    parser.add_argument(
        "--survival-file",
        default="/home/labs/alon/yoavh/Chackouts/LongitudonalTrajectories/tmp/meno_survival_data.csv",
        help="Survival curve file"
    )
    parser.add_argument(
        "--test-name",
        help="Specific test name to analyze (optional, analyzes all if not provided)"
    )
    parser.add_argument(
        "--system-name",
        help="System name (required if --test-name is provided)"
    )
    parser.add_argument(
        "--test-code",
        help="Test code (required if --test-name is provided)"
    )
    
    args = parser.parse_args()
    
    if not os.path.exists(args.output_dir):
        print(f"❌ Output directory not found: {args.output_dir}")
        return 1
    
    if not os.path.exists(args.bootstrap_folder):
        print(f"❌ Bootstrap folder not found: {args.bootstrap_folder}")
        return 1
    
    if not os.path.exists(args.survival_file):
        print(f"❌ Survival file not found: {args.survival_file}")
        return 1
    
    print(f"🔍 Output directory: {args.output_dir}")
    print(f"🔍 Bootstrap folder: {args.bootstrap_folder}")
    print(f"🔍 Survival file: {args.survival_file}")
    
    # If specific test provided, analyze only that test
    if args.test_name:
        if not args.system_name or not args.test_code:
            print("❌ --system-name and --test-code required when --test-name is provided")
            return 1
        
        start_time = time.time()
        cross_validate_single_test(
            args.test_code,
            args.test_name,
            args.system_name,
            args.output_dir,
            args.bootstrap_folder,
            args.survival_file
        )
        elapsed = time.time() - start_time
        print(f"\n⏱️  Total time: {elapsed:.2f} seconds")
        
        return 0
    
    # Otherwise, analyze all tests
    print("\n🔍 Searching for tests to analyze...")
    
    # Find all test directories
    tests_to_analyze = []
    for system_name in os.listdir(args.output_dir):
        system_dir = os.path.join(args.output_dir, system_name)
        if not os.path.isdir(system_dir):
            continue
        
        for test_dir_name in os.listdir(system_dir):
            test_dir = os.path.join(system_dir, test_dir_name)
            if not os.path.isdir(test_dir):
                continue
            
            # Check if bootstrap_results exists
            bootstrap_results_dir = os.path.join(test_dir, "bootstrap_results")
            if not os.path.exists(bootstrap_results_dir):
                continue
            
            # Try to extract test info from results
            result_files = [f for f in os.listdir(bootstrap_results_dir) 
                          if f.endswith("_result.csv")]
            if len(result_files) == 0:
                continue
            
            tests_to_analyze.append({
                "system_name": system_name,
                "test_dir_name": test_dir_name,
                "test_dir": test_dir
            })
    
    print(f"📋 Found {len(tests_to_analyze)} tests to analyze")
    
    if len(tests_to_analyze) == 0:
        print("❌ No tests found to analyze")
        return 1
    
    # Analyze each test
    start_time = time.time()
    for i, test_info in enumerate(tests_to_analyze):
        print(f"\n{'='*80}")
        print(f"Progress: {i+1}/{len(tests_to_analyze)}")
        print(f"{'='*80}")
        
        # For now, skip tests where we can't determine test_code
        # In production, you'd need to map test_dir_name back to test_code
        print(f"⚠️  Skipping {test_info['test_dir_name']} - need test_code mapping")
    
    elapsed = time.time() - start_time
    print(f"\n⏱️  Total time: {elapsed:.2f} seconds")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

