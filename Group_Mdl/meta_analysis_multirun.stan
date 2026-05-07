data {
  int<lower=1> N;           // total number of runs
  int<lower=1> K;           // number of subjects
  int<lower=1> p;           // number of parameters per run
  int<lower=1> q_subj;      // number of subject-level covariates
  int<lower=1> q_run;       // number of run-level covariates

  array[N] vector[p] y;          // run-level posterior means
  array[N] matrix[p,p] S;        // run-level posterior covariances

  matrix[K,q_subj] X_subj;       // subject-level covariates
  matrix[N,q_run]  X_run;        // run-level covariates

  array[N] int<lower=1,upper=K> subj;  // subject index for each run
}

parameters {
  vector[p] alpha;                  // group-level intercepts
  matrix[p, q_subj] B_subj;         // slopes for subject-level covariates
  matrix[p, q_run]  B_run;          // slopes for run-level covariates
  vector<lower=0>[p] tau;           // between-subject SDs
  cholesky_factor_corr[p] Lcorr;    // Cholesky factor of correlation
}

transformed parameters {
  matrix[p,p] Ltau = diag_pre_multiply(tau, Lcorr);
  matrix[p,p] Tau = multiply_lower_tri_self_transpose(Ltau);  // Tau = Ltau * Ltau'
}

model {
  // Priors
  alpha ~ normal(0, 1);
  to_vector(B_subj) ~ normal(0, 1);
  to_vector(B_run)  ~ normal(0, 1);
  tau ~ normal(0, 1);
  Lcorr ~ lkj_corr_cholesky(2);

  for (i in 1:N) {
    // Linear predictor for run i
    vector[p] mu_i = alpha
                     + B_subj * to_vector(X_subj[subj[i],])  // subject-level
                     + B_run  * to_vector(X_run[i,]);        // run-level

    // Cholesky factor for total covariance
    matrix[p,p] L_i = cholesky_decompose(Tau + S[i]);

    // Likelihood
    target += multi_normal_cholesky_lpdf(y[i] | mu_i, L_i);
  }
}
