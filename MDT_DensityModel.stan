data {
  
  // Dimensions of data
	int<lower=1> n_transects;  // total number of transects
  int<lower=1> n_segments;   // total number of segments
  int<lower=1> n_segments1;  // number of segments with at least one detection
  int<lower=1> n_years;      // number of years included in analyses
  int<lower=1> n_recov;      // number of recovery units (RUs)
  int<lower=1> n_telem_obs;  // number of telemetry observations
  int<lower=1> n_rmtorts;    // number of radiomarked tortoises  
  
  // Number of regression coefficients
  int<lower=1> n_cov_lam;    // number of (segment-level) fixed effects in model of abundance
  int<lower=1> n_cov_detect; // number of fixed effects in detection model
  int<lower=1> n_cov_telem;  // number of fixed effects in availability model
  
  // Covariate data
  matrix[n_segments, n_cov_lam] cov_lam;       // covariate values for each segment in model of abundance
  matrix[n_segments, n_cov_detect] cov_detect; // covariate values for each segment in detection model
  matrix[n_segments, n_cov_telem] cov_g0;      // covariate values for each segment in availability model
  matrix[n_telem_obs, n_cov_telem] cov_telem;  // covariate values for each observation of a radiomarked tortoise (to model availability) 
  
  // Indices
	int<lower=1> transectID[n_segments];                   // transect ID (index) associated with each segment
  int<lower=1> recovID[n_segments];                      // recovery unit ID (index) associated with each segment
  int<lower=1> year[n_segments];                         // year (index) associated with each segment
  vector[n_years] yeartrend;                             // vector of year index, minYear(1):maxYear(18)
	int<lower=0,upper=1> lateyrs[n_segments];              // indicator for survey years 2004-2018
  int<lower=1,upper=n_rmtorts> rmtortID[n_telem_obs];    // tortoise ID (index) associated with each telemetry observation
  
  // Response data
  int<lower=0,upper=1> y_telem[n_telem_obs];  // telemetered tortoise visible? (1/0)
  int<lower=0> n[n_segments];                 // number of segments surveyed
  matrix<lower=0>[n_segments1,max(n)] y;      // distances of observed tortoises for segments with >=1 detection 
										                          // [replaced NAs with 0s, as we loop over the number of observations on each segment]
  vector<lower=0>[n_segments] L;              // length of each segment (m)
  real<lower=0> W;                            // max observation distance (when data truncated, as they are here, W = T)
  real<lower=0> T;                            // truncation distance
  
  // Data for predictions
  matrix[n_recov,n_cov_lam] RUmeans;  // mean covariate values in each recovery unit
  int<lower=1> n_preds;               // number of surveyed year-RU combinations (= number of Density predictions)
  int<lower=1> recov_pred[n_preds];   // recovery unit ID (index) for density predictions
  int<lower=1> yr_pred[n_preds];      // year (index) for density predictions

}

parameters {

  // Std normal variates (unscaled random intercepts)
  matrix[n_recov,n_years] recov_l_raw;     // trend random effects     
	vector[n_transects]     tran_l_raw;      // transect-level random effects in abundance model
  vector[n_segments]      seg_l_raw;       // segment-level random effects in abundance model (only needed for years 2004-2018)
  vector[n_segments]      detect_raw;      // segment-level random effects in detection model
  vector[n_rmtorts]       g0_raw;          // random effects for radio-marked tortoises

  // Regression coefficients
  vector[n_cov_lam]    beta_lam;         // coefficients for segment-level covariates in abundance model
  vector[n_cov_detect] beta_detect;      // coefficients in detection model
  vector[n_cov_telem]  beta_g0;          // coefficients in telemetry model

  // Temporal trends (different in each recovery unit)
  vector[n_recov] trend;

  // Intercepts
  vector<lower=0,upper=20>[n_recov] mu_recov; // mean densities in year before surveys started (implicit uniform prior)
  real<lower=0,upper=20> mu_sigma;            // mean sigma in detection model (implicit uniform prior)
  real<lower=0,upper=1> mu_g0;                // mean availability (implicit uniform prior)
  
  // Variance components (as std dev)
  real<lower=0> sd_recov;
	real<lower=0> sd_tran;
  real<lower=0> sd_seg;
  real<lower=0> sd_detect;
  real<lower=0> sd_g0;
  
} 

transformed parameters {

  vector<lower=0,upper=1>[n_telem_obs] telem_p;
  vector<lower=0,upper=1>[n_segments] g0;
  vector<lower=0>[n_segments] sigma; 
  vector<lower=0,upper=1>[n_segments] Pc;
  matrix[n_recov,n_years] recov_l; 
  vector[n_segments] seg_l_mu;
  vector[n_segments] seg_l;
  vector<lower=0>[n_segments] lambda;
  vector<lower=0>[n_segments] lambda_p;
  vector[n_rmtorts] randtort;
   
  // Random effect specification (throughout): using SD * raw standard-normal variate 
	// This greatly improves mixing

  // Availability (estimated from telemetry data; g0 are predicted values for each segment)
  randtort = sd_g0 * g0_raw;
  for(i in 1:n_telem_obs)  
		telem_p[i] = inv_logit(logit(mu_g0) + cov_telem[i] * beta_g0 + randtort[rmtortID[i]]);
  g0 = inv_logit(logit(mu_g0) + cov_g0 * beta_g0);

  // SD in half-normal detection function
  sigma = exp(log(mu_sigma) + sd_detect * detect_raw + cov_detect * beta_detect);

  // Unconditional detection probability (includes g0)
  Pc = (sqrt(2 * pi()) / W * g0) .* sigma .* (Phi_approx(W ./ sigma) - 0.5);

  // Tortoise density: trend and random yearly effects in recovery units
  recov_l = log(mu_recov) * rep_row_vector(1, n_years) + sd_recov * recov_l_raw + trend * yeartrend'; 

  // Tortoise density: covariate and random effects in abundance model
  for (i in 1:n_segments) 
    seg_l_mu[i] = recov_l[recovID[i],year[i]];
	for (i in 1:n_segments)
		seg_l[i] = seg_l_mu[i] + cov_lam[i] * beta_lam + sd_tran * tran_l_raw[transectID[i]] + lateyrs[i] * sd_seg * seg_l_raw[i]; 

  // Poisson intensity for true abundance 
  lambda = (2 * W * L) .* exp(seg_l - log(1000000));  //need factor of 1,000,000 to change units from m2 to km2
  
  // Mean of Poisson for observed counts (ie, assuming a thinned Poisson process)
  lambda_p = lambda .* Pc;
  
}

model {

  // Priors on variance components
  target += cauchy_lpdf(sd_recov  | 0, 1);
	target += cauchy_lpdf(sd_tran   | 0, 1);
  target += cauchy_lpdf(sd_seg    | 0, 1);
  target += cauchy_lpdf(sd_detect | 0, 1);
  target += cauchy_lpdf(sd_g0     | 0, 1);

  // Priors on trend and regression coefficients
  target += normal_lpdf(trend       | 0, sqrt(10));
  target += normal_lpdf(beta_lam    | 0, sqrt(10));
  target += normal_lpdf(beta_detect | 0, sqrt(10));
  target += normal_lpdf(beta_g0     | 0, sqrt(10));

  // Unscaled random intercepts
  target += std_normal_lpdf(to_vector(recov_l_raw));
	target += std_normal_lpdf(tran_l_raw);
  target += std_normal_lpdf(seg_l_raw);
  target += std_normal_lpdf(detect_raw);
  target += std_normal_lpdf(g0_raw);

  // Telemetry submodel
  target += bernoulli_lpmf(y_telem | telem_p);

  // Abundance submodel
  // Replacing Poisson-binomial mixture with thinned Poisson process (mean = lambda*p)
  for (i in 1:n_segments)
	target += poisson_lpmf(n[i] | lambda_p[i]);

  // Distance detection submodel
  for (i in 1:n_segments1)
  {
    target += -n[i] * log(Phi_approx(T/sigma[i]) - 0.5) - n[i] * log(sigma[i]) - n[i] * 0.5 * log(2 * pi());
    for (j in 1:n[i])
      target += -0.5 * pow(y[i,j]/sigma[i], 2);
  }
}

generated quantities {

  matrix<lower=0>[n_recov,n_years] predTrend;
  vector<lower=0>[n_preds] D_survyrs;
  vector<lower=0>[n_segments] E_obs; // fit of observed data
  vector<lower=0>[n_segments] n_new; // replicate values
  vector<lower=0>[n_segments] E_new; // fit of replicate values
  real<lower=0> fit_obs;
  real<lower=0> fit_new;

  // Calculate trends in each recovery unit
  predTrend = exp((log(mu_recov) + RUmeans * beta_lam) * rep_row_vector(1, n_years) + trend * yeartrend');

  // Calculate mean annual density in each recovery unit (only in years unit was surveyed)
  for(i in 1:n_preds)
	D_survyrs[i] = exp(log(mu_recov[recov_pred[i]]) + RUmeans[recov_pred[i]] * beta_lam
		           + trend[recov_pred[i]] * yr_pred[i] + sd_recov * recov_l_raw[recov_pred[i],yr_pred[i]]); 

  // Fit statistics
  for(i in 1:n_segments)
  {
	E_obs[i] = pow((n[i] - lambda_p[i]), 2) / (lambda_p[i] + 0.5);
	n_new[i] = poisson_rng(lambda_p[i]);
	E_new[i] = pow((n_new[i] - lambda_p[i]), 2) / (lambda_p[i] + 0.5);
  }

  fit_obs = sum(E_obs[]);
  fit_new = sum(E_new[]);

}


