  /* Maginal GP */

  functions{
     matrix cov_exp_L(real sigmasq, real tausq,
                                    real phi, vector vdist, int N){

          int h = 0;
          matrix[N, N] K;

          for (j in 1:(N - 1)){
              K[j, j] = sigmasq + tausq;
              for (k in (j + 1):N){
                  h = h + 1;
                  K[j, k] = sigmasq * exp(- phi * vdist[h]);
                  K[k, j] = K[j, k];
              }
          }
          K[N, N] = sigmasq + tausq;
          return cholesky_decompose(K);
      }

      vector get_vdist(matrix coords){
          int h = 0;
          int N = dims(coords)[1];
          vector[N * (N - 1) / 2] vdist;

          for(j in 1:(N - 1)){
              for(k in (j + 1):N){
                  h = h + 1;
                  vdist[h] = distance(coords[j, ], coords[k, ]);
              }
          }
          return vdist;
      }
  }

  data {
      int<lower=1> N;
      int<lower=1> P;
      vector[N] Y;
      matrix[N, P + 1] X;
      matrix[N, 2] coords;
      vector[P + 1] uB;
      matrix[P + 1, P + 1] VB;
      real ss;
      real st;
      real ap;
      real bp;
  }

  transformed data {
      cholesky_factor_cov[P + 1] L_VB = cholesky_decompose(VB);
      vector[N * (N - 1) / 2] vdist = get_vdist(coords);
  }

 parameters {
      vector[P + 1] beta;
      real<lower = 0> sigma;
      real<lower = 0> tau;
      real<lower = 0> phi;
  }

  transformed parameters{
      real sigmasq = square(sigma);
      real tausq = square(tau);
  }

  model{
      beta ~ multi_normal_cholesky(uB, L_VB);
      phi ~ gamma(ap, bp);
      sigma ~ normal(0, ss);
      tau ~ normal(0, st);

      Y ~ multi_normal_cholesky(X * beta,
                                cov_exp_L(sigmasq, tausq, phi, vdist, N));
  }















