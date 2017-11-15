  /* Random-effects NNGP model with spatial random effect centered at intercept */

  functions{
      real nngp_w_lpdf(vector w_b1, real sigmasq, real phi, matrix NN_dist,
                       matrix NN_distM,int[,] NN_ind, int N, int M,
                       real intercept){

          vector[N] V;
          vector[N] I_Aw;
          vector[N] w;
          int dim;
          int h;
          real out;
          w = w_b1 - intercept;
          I_Aw = w;

          for (i in 2:N) {
              matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M]
              iNNdistM;
              matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M]
              iNNCholL;
              vector[ i < (M + 1)? (i - 1) : M] iNNcorr;
              vector[ i < (M + 1)? (i - 1) : M] v;
              row_vector[i < (M + 1)? (i - 1) : M] v2;

              dim = (i < (M + 1))? (i - 1) : M;

              // get exp(-phi * NN_distM)

              if(dim == 1){iNNdistM[1, 1] = 1;}
              else{
                  h = 0;
                  for (j in 1:(dim - 1)){
                      for (k in (j + 1):dim){
                          h = h + 1;
                          iNNdistM[j, k] = exp(- phi * NN_distM[(i - 1), h]);
                          iNNdistM[k, j] = iNNdistM[j, k];
                      }
                  }
                  for(j in 1:dim){
                      iNNdistM[j, j] = 1;
                  }
              }

              iNNCholL = cholesky_decompose(iNNdistM);
              for (j in 1: dim){
                  iNNcorr[j] = exp(- phi * NN_dist[(i - 1), j]);
              }
            
              //vector[dim] v;
              v = mdivide_left_tri_low(iNNCholL, iNNcorr);

              V[i] = 1 - dot_self(v);

              v2 = mdivide_right_tri_low(v', iNNCholL);

              for (j in 1:dim){
                  I_Aw[i] = I_Aw[i] - v2[j] * w[NN_ind[(i - 1), j]];
              }
          }
          V[1] = 1;
          out = - 0.5 * ( 1 / sigmasq * dot_product(I_Aw, (I_Aw ./ V)) +
                          sum(log(V)) + N * log(sigmasq));
          return out;
      }
  }


  data {
      int<lower=1> N;
      int<lower=1> M;
      int<lower=1> P;
      vector[N] Y;
      matrix[N, P + 1] X;
      int NN_ind[N - 1, M];
      matrix[N - 1, M] NN_dist;
      matrix[N - 1, (M * (M - 1) / 2)] NN_distM;
      vector[P + 1] uB;
      matrix[P + 1, P + 1] VB;
      real ss;
      real st;
      real ap;
      real bp;
  }

  transformed data {
      cholesky_factor_cov[P + 1] L_VB;
      L_VB = cholesky_decompose(VB);
  }

  parameters{
      vector[P + 1] beta;
      real<lower = 0> sigma;
      real<lower = 0> tau;
      real<lower = ap, upper = bp> phi;
      vector[N] w_b1;
  }

  transformed parameters {
      real sigmasq = sigma^2;
      real tausq = tau^2;
  }

  model{
      beta ~ multi_normal_cholesky(uB, L_VB);
      sigma ~ normal(0, ss);
      tau ~ normal(0, st);
      w_b1 ~ nngp_w(sigmasq, phi, NN_dist, NN_distM, NN_ind, N, M, beta[1]);
      Y ~ normal(block(X, 1, 2, N, P) * tail(beta, P) + w_b1, tau);
  }
















