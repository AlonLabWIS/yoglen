# background
yoglen was developed using R version 4.5.1. yoglen contains scripts for systems analysis of menopause from large scale studies where menopause timing is unknown. yoglen uses latent variable modelling to infer the time to menopause changes of lab tests and outcomes. it includes two parallel analyses using two related but unique models and two different datasets (NHANES and Clalit).

# installation
1. Clone the package


# yoglen latent variable ("deconvolution") algorithm for inferring unknown menopause timing
**Background**
yoglen was developed using R version 4.5.1. yoglen contains scripts for systems analysis of menopause from large scale studies where menopause timing is unknown. yoglen uses latent variable modelling to infer the time to menopause changes of lab tests and outcomes. it includes two parallel analyses using two related but unique models and two different datasets (NHANES and Clalit).

---

# Citation info
Please cite as:


---

## 📘 Description
Contains:  
-essential analysis steps (i.e. yoglen model fitting)  
-plot scripts for each figure  
-aggregated data needed for each figure and summary stats  
-simulated NHANES data  
-unedited analysis scripts for transparency (transparency folder)  

Does not include individualized data:  
-NHANES data (download here: https://wwwn.cdc.gov/nchs/nhanes/)  
-Clalit electronic medical records

There is an associated Zenodo repository that is recommended but not strictly necessary:


---

## Dependencies
1. Rstudio is needed to run .rmd files (https://posit.co/downloads/)  
2. R packages (dependencies.R): c("survival","mgcv","survPen","dplyr, ggplot2","cowplot","scico","gridExtra","ggrepel")  
3. Some files are too large for GitHub and are on the associated Zenodo repository. Two options: (i) don't download the Zenodo files, continue as normal and use simulated data. (ii) download the Zenodo files that contain the exact analysis files.

---

## ⚙️ Installation
1. Clone repository  
2. Install dependencies  
3. Install and open Rstudio   
4. Ready to go!  

```bash
git clone git@github.com:AlonLabWIS/yoglen.git
cd yoglen
R dependencies.R
```

---

## Usage
### Making plots
Run plot.Rmd 

-This will use existing aggregated data from the original analysis

### Analysis
-This will default to using simulated data, unless the user provides real data.  
The general steps are:  
1. Estimate the final menstrual period (FMP) distribution and impute a distribution of values for each person (age_of_menopause_imputation.Rmd),  
2. Model selection (model_selection.Rmd),  
3. Plot results (plot.Rmd).  

### Simulated data
-yoglen_simulation_study.Rmd will use the fitted models to simulate realistic data.  
-yoglen_validation_step_function.Rmd simulates a simple step function to show that the approaches are insensitive to shifting and scaling of the FMP distribution.

## yoglen algorithm
The yoglen algorithm uses expectation maximization to maximize the expected log-likelihood using a discrete grid for the latent variable. The latent variable is the unknown age of menopause. Regression models as functions of time to menopause can be optionally supplied assuming conditional independence (i.e. ignoring direct correlations between variables but allowing them via time to menopause).


(Please forgive the use of $\pi_{ij}$ to denote discrete probability parameter set.)

The full $p(a_m|y,a,m)p(a_m|a,m)jp(y_j|a_m,a,m)$ used a regression model for each $y_j$, $p(y_j|a_m,a,m)$. This model permits a smooth time to menopause dependence $t\equiv a-a_m$ with the potential for a jump at menopause $y(t)f(t)+\delta(t)$. Assuming Gaussian error the full log-likelihood including the latent distribution is
$$
  l=\sum_i \big[ ln(\pi_{ij})-\frac{1}{2\sigma(t_i)^2}(y_i-f(t_i))^2-\frac{1}{2}ln(2\pi\sigma(t_i)^2) \big]
$$
where $\pi_{ij}$ is the probability of $t_i$ equalling the $j$th grid value given the age of the individual $a_i$ (normalized to $\sum_j \pi_{ij} = 1$).
Then we use expectation maximization on the log-likelihood. The conditional expectation over all possible t given age is
$$
  E[l|a]=\sum_{ij} \pi_{ij} \big[ ln(\pi_{ij})-\frac{1}{2\sigma(t_i)^2}(y_i-f(t_i))^2-\frac{1}{2}ln(2\pi\sigma(t_i)^2) \big].
$$
In the expectation maximization approach we iteratively estimate the expected latent variable distribution, $\pi_{ij}$, and use it to maximize the log-likelihood, $max_fE(l|a)$. Since the parameters of $f$ do not mix with the unknown $t$ the likelihood can be separately maximized using the partial expected likelihood
$$
  max_f E[l|a]=max_f\big[ -\frac{1}{2}\sum_{ij} \pi_{ij} \big[ \frac{1}{\sigma(t_i)^2}(y_i-f(t_i))^2+ln(2\pi\sigma(t_i)^2) \big] \big],
$$
this is exactly the weighted log-likelihood for regression with weights $\pi_{ij}$. We can therefore substitute any regression model into the maximization step (formally the errors must be Gaussian but informally everything works.)

### YoGlen Algorithm Steps
1. Estimate the starting distribution, $\pi_{ij}$, for each individual based on the menopause vs age curve.
2. While expected log-likelihood $E[l|a]$ is decreasing:
    2a. Estimate a weighted regression model using the vectorized ij i.e. $f(t_ij)$ with weight $\pi_{ij}$
    2b. Update $\pi_{ij}=p(t_i|a_i,y_i)$.
For the linear model we solve the likelihood exactly. For generalized linear models we used a spline model with default settings and weighted with the ($t_j$,$\pi_{ij}$) grid. For multivariate adaptive regression splines we used the $\pi_{ij}$ as weights to resample and proceeded normally with default settings.

There are two methods used in the original paper:

Method 1. (used on NHANES)
$y_i$ is menopause status. If menopause age is known then it is set to a normal random variable with standard deviation 1 (can be changed as function argument).

Method 2. (used on Clalit)
Piecewise linear model for $f(t)$ with no menopause status or menopause age. Solved exactly in the paper.
