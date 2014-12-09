//************************************************************************
//    sd-opacity.c - 
//------------------------------------------------------------------------
//    Written by:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute
//    Bocni II 1401/1a, 141-31 Praha 4, Czech Republic
//------------------------------------------------------------------------
//    Created: 10.12.2012
//************************************************************************


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "fluxmodel.h"
#include "flux-sd.h"
#include "flux-nt.h"
#include "sim5kerr.h"
#include "sim5const.h"
#include "sim5math.h"



int main(int argc, char *argv[])
//***************************************************
{
    if (argc!=6) {
        fprintf(stderr, "Use: %s spin L alpha rmax nrad\n", argv[0]);
        exit(-1);
    }

    double M     = 10.0;              // BH mass
    double a     = atof(argv[1]);     // BH spin
    double lum   = atof(argv[2]);     // BH luminosity
    double alpha = atof(argv[3]);     // alpha parameter
    double rmax  = atof(argv[4]);     // maximal radius to consider
    double rmin  = r_horizon(a)*1.01; // minimal radius
    int    nrad  = atoi(argv[5]);     // number of rings

    // initialize flux model
    flux_model_init("./flux-nt.so", M, a, lum, alpha, NULL, SD_INIT_OPT_LUMINOSITY);

    int ir;
    for (ir=0; ir<nrad; ir++) {
        double R = exp(log(rmin) + (log(rmax)-log(rmin))*(double)(ir)/(double)(nrad-1));
        double T = flux_model_eval_temp(R);                              // [K]
        double S = flux_model_eval_sigma(R);                             // [g cm^-2]
        double H = flux_model_eval_h(R)*M*grav_radius_cgs;               // [cm]
        double l = flux_model_eval_l(R);                                 // [?]
        double V = flux_model_eval_vr(R);                                // [(c)]
        double C = flux_model_eval_custom(R, SD_EVAL_TEMP_CENTRAL);      // [K]
        // central density
        // density profile: // rho = rho0 * (1 - z^2/H^2)^N;   N=3...polytropic index
        // Sigma = \int_0^H rho dz = 16/35*H*rho0 (for N=3) ---> rho0 = 35/16*Sigma/H
        double D = 35./16.*S/H;                                                 // [g cm^-3]

        double tau_ff = 6.4e22 * 5.*M_PI/32. * pow(D,2.) * pow(C,-3.5) * H;
        double tau_es = 0.34 * 16./35. * D * H;
        double tau_eff = sqrt(tau_ff*(tau_ff+tau_es));
        double tau_tot = tau_ff+tau_es;
        double t_thermal  = (H/tau_eff)/speed_of_light_cgs;
        double t_freefall = (sqrt(2.*sqr3(R))/3. - 4./3.)*M*grav_radius_cgs/speed_of_light_cgs;

        fprintf(stdout, "%.4e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n", R, C, T, D, H, tau_tot, tau_eff, tau_ff, tau_es, t_thermal/t_freefall);
    }

    flux_model_done();

    return 0;
}




