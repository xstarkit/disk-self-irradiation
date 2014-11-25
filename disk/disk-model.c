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
#include "sim5lib.h"


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
    if (argc!=7) {
        fprintf(stderr, "Use: %s disk-model spin L alpha rmax nrad\n", argv[0]);
        exit(-1);
    }

    double M     = 10.0;              // BH mass
    char*  model = argv[1];           // disk model
    double a     = atof(argv[2]);     // BH spin
    double lum   = atof(argv[3]);     // BH luminosity
    double alpha = atof(argv[4]);     // alpha parameter
    double rmax  = atof(argv[5]);     // maximal radius to consider
    double rmin  = 0.0;               // minimal radius
    int    nrad  = atoi(argv[6]);     // number of rings

    nRS = nrad*100;
    RS  = (radial_structure*) calloc(nRS, sizeof(radial_structure));

    disk_nt_setup(M, a, lum, alpha, 0);
    rmin = disk_nt_r_min();

    fprintf(stdout, "# radial disk structure model\n");
    fprintf(stdout, "# model: %s\n", model);
    fprintf(stdout, "# mass:  %.4f\n", M);
    fprintf(stdout, "# spin:  %.4f\n", a);
    fprintf(stdout, "# lum:   %.4f\n", lum);
    fprintf(stdout, "# alpha: %.4f\n", alpha);
    fprintf(stdout, "# r_min: %.2f\n", rmin);
    fprintf(stdout, "# r_max: %.2f\n", rmax);
    fprintf(stdout, "# N_rad: %d\n", nrad);
    fprintf(stdout, "# model mdot: %e\n", disk_nt_mdot());
    fprintf(stdout, "# model lumi: %e\n", disk_nt_lum());

    // evaluate local model parameters at r=r0
    int ir;
    double Ltot = 0.0;
    for (ir=0; ir<nRS; ir++) {
        RS[ir].R = exp(log(rmin) + (log(rmax)-log(rmin))*(double)(ir)/(double)(nRS-1));
        float R2 = exp(log(rmin) + (log(rmax)-log(rmin))*(double)(ir+1)/(double)(nRS-1));

        RS[ir].T = disk_nt_flux(RS[ir].R);
        RS[ir].S = disk_nt_sigma(RS[ir].R);
        RS[ir].H = disk_nt_h(RS[ir].R);
        RS[ir].l = disk_nt_ell(RS[ir].R);
        RS[ir].V = disk_nt_vr(RS[ir].R);
        
        Ltot += pow(RS[ir].T,4.)*RS[ir].R*(R2-RS[ir].R);     // dL = F*r*dr
        RS[ir].totL = Ltot>0 ? (Ltot) : 0.0; 

        // vertical gravity (see Zhu et al. 2012)
        sim5metric m;;
        kerr_metric(a, RS[ir].R, 0.0, &m);
    	double Omega = Omega_from_ell(RS[ir].l, &m);
    	double u[4];
        if (fabs(RS[ir].V)>0.0) {
            // radial velocity (V) is measured in fluid co-rotating frame
            sim5tetrad t;
            tetrad_azimuthal(&m, Omega, &t);
            double u_loc[4];
            u_loc[0] = sqrt(1.0+sqr(RS[ir].V));     // normalization to u.u=-1
            u_loc[1] = RS[ir].V;                    // r-vector of the ortonorm frame points outwards & vr is negative
            u_loc[2] = 0.0;
            u_loc[3] = 0.0;
            on2bl(u_loc, u, &t);
        } else {
          	u[0] = sqrt(-1.0/(m.g00 + 2.*Omega*m.g03 + sqr(Omega)*m.g33));
           	u[1] = 0.0;
           	u[2] = 0.0;
           	u[3] = u[0]*Omega;
        }
        double u_t = u[0]*m.g00 + u[3]*m.g03;
        double u_f = u[0]*m.g03 + u[3]*m.g33;
        RS[ir].Q = M*solar_mass*grav_const/pow(RS[ir].R*M*grav_radius,3.) * (sqr(u_f) + a*a*(u_t-1.))/RS[ir].R;
    }

    fprintf(stderr, "L0: %e / %.3f\n", RS[nRS-1].totL, log(RS[nRS-1].totL));
    fprintf(stderr, "L1: %e / %.3f\n", RS[1].totL, log(RS[1].totL));
    double dL = (RS[nRS-2].totL-RS[1].totL);
    fprintf(stderr, "Lmax=%.3f dL: %.3f\n", (RS[nRS-2].totL-RS[1].totL), dL);
    int ID;
    double prevL, prevR;


    do {
        dL /= 1.1;
        ID = 0;
        prevL = prevR = 0.0;
        for (ir=0; ir<nRS; ir++) {
        if ((ir==0) || (ir==nRS-1) || ((RS[ir].totL-prevL >= dL)&&(RS[ir].R-prevR>0.2)) || ((RS[ir].R-prevR)/RS[ir].R>0.5)) { 
                ID++;
                prevL = RS[ir].totL;
                prevR = RS[ir].R;
            }
        }
        fprintf(stderr, "dL= %e  %d/%d\n", ID, nrad);
    } while (ID!=nrad);

    ID = 0;
    prevL = prevR = 0.0;
    for (ir=0; ir<nRS; ir++) {
        if ((ir==0) || (ir==nRS-1) || ((RS[ir].totL-prevL >= dL)&&(RS[ir].R-prevR>0.2)) || ((RS[ir].R-prevR)/RS[ir].R>0.3)) { 
            fprintf(
                stdout, "%04d/%04d  %.2e  %.4e    %e    %e    %e    %e    %e    %e\n", 
                ID+1, ir, RS[ir].totL, RS[ir].R, RS[ir].T, RS[ir].S, RS[ir].Q, RS[ir].H, RS[ir].l, RS[ir].V
            );
            ID++;
            prevL = RS[ir].totL;
            prevR = RS[ir].R;
        }
    }

    free(RS);

    // close flux model
    disk_nt_finish();

    return 0;
}



