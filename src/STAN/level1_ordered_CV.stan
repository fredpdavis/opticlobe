/*
 * Models log(O_gs) as mixture of on and off gaussian components
 *
 * 1. P(E_g) in bimodal-on v bimodal-off: mix of 2 normals
 *    1. 'on' nuclei that express gene
 *    2. 'off' nuclei that do not express gene

 * given GS knowns: O_gs, observed expression matrix (genes x samples)
 *
 * infer 5 unknowns:
 * - P(E_g|on) mu and sigma
 * - P(E_g|off) mu and sigma
 * - pi_on mixture weight
 *
 */


data {

   int<lower=2> nSamples_t ;
   int<lower=2> nSamples_h ;

   vector<lower=0,upper=14> [nSamples_t] logE_t;  // observed tx abundance
   vector<lower=0,upper=14> [nSamples_h] logE_h;  // observed tx abundance

}


parameters {

   positive_ordered[2]        mu ;
   real<lower=0.001,upper=3>  sd1 ;
   real<lower=0,upper=1>      pi_on ;    // on/off mixing weight

}

model {

   mu ~ normal(7,5) ;

   for (s in 1:nSamples_t) {
      target += log_mix(pi_on,
                        normal_lpdf(logE_t[s] | mu[2], sd1),
                        normal_lpdf(logE_t[s] | mu[1], sd1));
   }

}


generated quantities {

   vector [nSamples_t] log_lik_t ;
   vector [nSamples_h] log_lik_h ;

   for (s in 1:nSamples_t) {
      log_lik_t[s] = log_mix(pi_on,
                        normal_lpdf(logE_t[s] | mu[2], sd1),
                        normal_lpdf(logE_t[s] | mu[1], sd1)); }

   for (s in 1:nSamples_h) {
      log_lik_h[s] = log_mix(pi_on,
                        normal_lpdf(logE_h[s] | mu[2], sd1),
                        normal_lpdf(logE_h[s] | mu[1], sd1)); }


}
