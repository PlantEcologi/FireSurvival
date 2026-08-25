data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  array[N] int<lower=0, upper=1> y;
  array[N] int<lower=0, upper=1> include_row;
}

parameters {
  vector[K] beta;
}

model {
  // Match legacy diffuse priors (R implementation used Normal(0, 1000)).
  beta ~ normal(0, 1000);

  for (n in 1:N) {
    if (include_row[n] == 1) {
      if (y[n] == 1) {
        target += normal_lcdf(dot_product(X[n], beta) | 0, 1);
      } else {
        target += normal_lccdf(dot_product(X[n], beta) | 0, 1);
      }
    }
  }
}

generated quantities {
  vector[N] log_lik;
  vector[N] p_fire;

  for (n in 1:N) {
    real eta = dot_product(X[n], beta);
    p_fire[n] = Phi(eta);

    if (include_row[n] == 1) {
      if (y[n] == 1) {
        log_lik[n] = normal_lcdf(eta | 0, 1);
      } else {
        log_lik[n] = normal_lccdf(eta | 0, 1);
      }
    } else {
      log_lik[n] = 0;
    }
  }
}
