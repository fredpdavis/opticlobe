/*
 * Models log(O_gs) as arising from a single gaussian component
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

   vector<lower=-12,upper=20> [nSamples_t] logE_t;  // observed tx abundance
   vector<lower=-12,upper=20> [nSamples_h] logE_h;  // observed tx abundance

}


parameters {

// gene-specific estimates of P(E_g) component 1 and 2
// GOAL: assign gene,samples  to on/off components
//   real<lower=-10,upper=17>     mu1 ;
   real<lower=0,upper=14>     mu1 ;
   real<lower=0.001,upper=3>  sd1 ;

}


model {

   for (s in 1:nSamples_t) {
      target += normal_lpdf(logE_t[s] | mu1, sd1) ;
   }

}


generated quantities {

   vector [nSamples_t] log_lik_t ;
   vector [nSamples_h] log_lik_h ;

   for (s in 1:nSamples_t) {
      log_lik_t[s] = normal_lpdf(logE_t[s] | mu1, sd1) ; }

   for (s in 1:nSamples_h) {
      log_lik_h[s] = normal_lpdf(logE_h[s] | mu1, sd1) ; }

}
