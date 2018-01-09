  /* Response NNGP model */

  functions{
      real nngp_lpdf(vector Y, vector X_beta, real sigmasq, real tausq,
                     real phi, matrix NN_dist, matrix NN_distM, int[,] NN_ind,
                     int N, int M){

          vector[N] V;
          vector[N] YXb = Y - X_beta;
          vector[N] U = YXb;
          real kappa_p_1 = tausq / sigmasq + 1;
          int dim;
          int h;

          for (i in 2:N) {
              matrix[ i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1): M]
              iNNdistM;
              matrix[ i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1): M]
              iNNCholL;
              vector[ i < (M + 1) ? (i - 1) : M] iNNcorr;
              vector[ i < (M + 1) ? (i - 1) : M] v;
              row_vector[i < (M + 1) ? (i - 1) : M] v2;
              dim = (i < (M + 1))? (i - 1) : M;

              if(dim == 1){iNNdistM[1, 1] = kappa_p_1;}
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
                      iNNdistM[j, j] = kappa_p_1;
                  }
              }

              iNNCholL = cholesky_decompose(iNNdistM);
              iNNcorr = to_vector(exp(- phi * NN_dist[(i - 1), 1: dim]));

             v = mdivide_left_tri_low(iNNCholL, iNNcorr);

             V[i] = kappa_p_1 - dot_self(v);

             v2 = mdivide_right_tri_low(v', iNNCholL);

             U[i] = U[i] - v2 * YXb[NN_ind[(i - 1), 1:dim]];
          }
          V[1] = kappa_p_1;
          return - 0.5 * ( 1 / sigmasq * dot_product(U, (U ./ V)) +
                          sum(log(V)) + N * log(sigmasq));
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
      real<lower = 0> phi;
  }

  transformed parameters {
      real sigmasq = square(sigma);
      real tausq = square(tau);
  }

  model{
      beta ~ multi_normal_cholesky(uB, L_VB);
      phi ~ gamma(ap, bp);
      sigma ~ normal(0, ss);
      tau ~ normal(0, st);
      Y ~ nngp(X * beta, sigmasq, tausq, phi, NN_dist, NN_distM, NN_ind, N, M);
  }













