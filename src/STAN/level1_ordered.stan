/*
 * 170729_0103 -- can we improve mixing / sampling efficiency by
 *                using ordered mu's?
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

   int<lower=2> nSamples ;
   int<lower=1> nDrivers ;
   int<lower=1> nCells ;

   vector<lower=0,upper=14>    [nSamples] logE;       // observed tx abundance
   int                         driver[nSamples] ;     // driver
   int                         cell[nSamples] ;       // cell type

}


parameters {

// gene-specific estimates of P(E_g) component 1 and 2
// GOAL: assign gene,samples  to on/off components

   positive_ordered[2]        mu ;       // component means
   real<lower=0.001,upper=3>  sd1 ;      // equal variance
   real<lower=0,upper=1>      pi_on ;    // on/off mixing weight

}

transformed parameters {

   real mu_on ;
   real sd_on ;
   real mu_off ;
   real sd_off ;

   mu_on  = mu[2] ;
   mu_off = mu[1] ;
   sd_on  = sd1 ;
   sd_off = sd1 ;

}


model {

   mu ~ normal(7,5) ;

   for (s in 1:nSamples) {
      target += log_mix(pi_on,
                        normal_lpdf(logE[s] | mu_on, sd_on),
                        normal_lpdf(logE[s] | mu_off, sd_off));
   }

}


generated quantities {

   vector [nSamples]    pon_gs ;   // posterior likelihood of expression
   vector [nDrivers]    pon_gd ;   // cell-level likelihood of expression
   vector [nCells]      pon_gc ;   // cell-level likelihood of expression

   vector [nDrivers]    pon_gd_onmass ;
   vector [nDrivers]    pon_gd_offmass ;

   vector [nCells]      pon_gc_onmass ;
   vector [nCells]      pon_gc_offmass ;


   for (d in 1:nDrivers) {
      pon_gd_onmass[d] = 0 ;
      pon_gd_offmass[d] = 0 ;
   }

   for (c in 1:nCells) {
      pon_gc_onmass[c] = 0 ;
      pon_gc_offmass[c] = 0 ;
   }

   for (s in 1:nSamples) {

      real ps[2];  // temp for log component P
      real sumps;  // temp for sum of log component P

      ps[1] = normal_lpdf(logE[s] | mu_on, sd_on) ;
      ps[2] = normal_lpdf(logE[s] | mu_off, sd_off) ;

      sumps = log_mix(pi_on, ps[1], ps[2]) ;

      pon_gs[s] = log(pi_on) + ps[1] - sumps ;

      pon_gd_onmass[driver[s]]  = pon_gd_onmass[driver[s]] + ps[1] ;
      pon_gd_offmass[driver[s]] = pon_gd_offmass[driver[s]] + ps[2] ;

      pon_gc_onmass[cell[s]]  = pon_gc_onmass[cell[s]] + ps[1] ;
      pon_gc_offmass[cell[s]] = pon_gc_offmass[cell[s]] + ps[2] ;

   }

   for (c in 1:nCells) {
      pon_gc[c] = log(pi_on) + pon_gc_onmass[c] -
                  log_mix(pi_on, pon_gc_onmass[c], pon_gc_offmass[c]) ;
   }

   for (d in 1:nDrivers) {
      pon_gd[d] = log(pi_on) + pon_gd_onmass[d] -
                  log_mix(pi_on, pon_gd_onmass[d], pon_gd_offmass[d]) ;
   }

}
