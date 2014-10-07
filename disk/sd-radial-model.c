//************************************************************************
//    sd-radial-model.c - print out radial model of disk
//------------------------------------------------------------------------
//    Written by:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute
//    Bocni II 1401/1a, 141-31 Praha 4, Czech Republic
//------------------------------------------------------------------------
//    Created: 18.10.2012
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


typedef struct {
    double R;
    double T;
    double S;
    double H;
    double l;
    double V;
    double Q;
    double totL;
} radial_structure;

radial_structure* RS = NULL;
int nRS = 0;



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
    double rmin  = r_horizon(a)*1.05; // minimal radius
    int    nrad  = atoi(argv[5]);     // number of rings

    // initialize flux model
    flux_model_init("./flux-sd.so", M, a, lum, alpha, NULL, SD_INIT_OPT_LUMINOSITY);

    nRS = nrad*100;
    RS  = (radial_structure*) calloc(nRS, sizeof(radial_structure));
    double R0 = r_horizon(a);

    // evaluate local model parameters at r=r0
    int ir;
    for (ir=0; ir<nRS; ir++) {
        RS[ir].R = exp(log(rmin) + (log(rmax)-log(rmin))*(double)(ir)/(double)(nRS-1));
        RS[ir].T = flux_model_eval_temp(RS[ir].R);
        RS[ir].S = flux_model_eval_sigma(RS[ir].R);
        RS[ir].H = flux_model_eval_h(RS[ir].R);
        RS[ir].l = flux_model_eval_l(RS[ir].R);
        RS[ir].V = flux_model_eval_vr(RS[ir].R);
        RS[ir].totL = (ir>0?RS[ir-1].totL:0.0) + pow(RS[ir].T,4.)*sqr3(RS[ir].R)*(RS[ir].R-R0);     // dL = F*r*dr

        // vertical gravity (see Zhu et al. 2012)
        double g00, g11, g22, g33, g03;
        kerr_metric(a, RS[ir].R, 0.0, &g00, &g11, &g22, &g33, &g03);
    	double Omega = Omega_from_ell(RS[ir].l, g00, g33, g03);
    	double u[4];
      	u[0] = sqrt(-1.0/(g00 + 2.*Omega*g03 + sqr(Omega)*g33));
       	u[1] = 0.0;
       	u[2] = 0.0;
       	u[3] = u[0]*Omega;
        if (fabs(RS[ir].V)>0.0) {
            // radial velocity (V) is measured in fluid co-rotating frame
            double e0[4], e1[4], e2[4], e3[4];
            ortho_tetrad_U_azimuthal_motion(u, g00, g11, g22, g33, g03, e0, e1, e2, e3);
            double u_loc[4];
            u_loc[0] = sqrt(1.0+sqr(RS[ir].V));     // normalization to u.u=-1
            u_loc[1] = RS[ir].V;                    // r-vector of the ortonorm frame points outwards & vr is negative
            u_loc[2] = 0.0;
            u_loc[3] = 0.0;
            on2bl(u_loc, u, g00, g11, g22, g33, g03, e0, e1, e2, e3);
        }
        double u_t = u[0]*g00 + u[3]*g03;
        double u_f = u[0]*g03 + u[3]*g33;
        RS[ir].Q = M*solar_mass*grav_const/pow(RS[ir].R*M*grav_radius,3.) * (sqr(u_f) + a*a*(u_t-1.))/RS[ir].R;

        R0 = RS[ir].R;
    }

    double dL = (log(RS[nRS-1].totL) - log(RS[0].totL))/(nrad-1);
    int ID = 0;
    for (ir=0; ir<nRS; ir++) {
        if ((ir==nRS-1) || ((log(RS[ir].totL)-log(RS[0].totL)) > ID*dL)) { 
            fprintf(
                stdout, "%04d   %.4e    %e    %e    %e    %e    %e    %e\n", 
                ID+1, RS[ir].R, RS[ir].T, RS[ir].S, RS[ir].Q, RS[ir].H, RS[ir].l, RS[ir].V
            );
            ID++;
        }
    }

    free(RS);
    // close flux model
    flux_model_done();

    return 0;
}



