#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include <fftw3.h>

#include "../cosmolike_core/cfftlog/cfftlog.h"
#include "../cosmolike_core/cfftlog/utils.h"
#include "../cosmolike_core/cfastpt/cfastpt.h"
#include "../cosmolike_core/cfastpt/utils.h"

#include "../cosmolike_core/theory/baryons.h"
#include "../cosmolike_core/theory/basics.c"
#include "../cosmolike_core/theory/structs.c"
#include "../cosmolike_core/theory/parameters.c"
#include "../cosmolike_core/emu17/P_cb/emu.c"
#include "../cosmolike_core/theory/recompute.c"
#include "../cosmolike_core/theory/cosmo3D.c"
#include "../cosmolike_core/theory/redshift_spline.c"
#include "../cosmolike_core/theory/halo.c"
#include "../cosmolike_core/theory/pt_cfastpt.c"
#include "../cosmolike_core/theory/radial_weights.c"
#include "init_y3.c"
#include "../cosmolike_core/theory/init_baryon.c"
#include "../cosmolike_core/theory/priors_mpp.c"
#include "../cosmolike_core/theory/cosmo2D_exact_fft.c"
#include "../cosmolike_core/theory/reduced_shear.c"
#include "../cosmolike_core/theory/cosmo2D_fullsky_TATT.c"
#include "../cosmolike_core/cfftlog/utils_complex.h"
#include "../cosmolike_core/cfastpt/utils_complex.h"
double X_lens = 0.87;
//double point_mass[10] = {3.23, 3.04, 2.85, 2.71, 2.54,0.,0.,0.,0.,0.};
double point_mass[10] = {3.,3.,3.,3.,3.,0.,0.,0.,0.,0.};
typedef struct input_nuisance_params_y3 {
    double bias[10];
    double bias2[10];
    double b_mag[10];
    double A_z[10];
    double A2_z[10];
    double b_ta[10];
    double lens_z_bias[10];
    double source_z_bias[10];
    double shear_m[10];
    double pm[10];
} input_nuisance_params_y3;

typedef struct input_cosmo_params_y3 {
    double omega_m;
    double sigma_8;
    double A_s;
    double n_s;
    double w0;
    double wa;
    double omega_b;
    double omega_nuh2;
    double h0;
    double MGSigma;
    double MGmu;
    double theta_s;
} input_cosmo_params_y3;

void set_params_from_python(input_cosmo_params_y3 ic, input_nuisance_params_y3 in);
void compute_w_1tomobin(double *data, int z1);
void compute_gamma_1tomobin(double *data, int zl,int zs);
void compute_xip_1tomobin(double *data, int z1,int z2);
void compute_xim_1tomobin(double *data, int z1,int z2);

double log_like_wrapper(input_cosmo_params_y3 ic, input_nuisance_params_y3 in);
void write_datavector_wrapper(char *filename, input_cosmo_params_y3 ic, input_nuisance_params_y3 in);

double get_sigma_8(input_cosmo_params_y3 ic);
double get_h0(input_cosmo_params_y3 ic);

double get_h0(input_cosmo_params_y3 ic){
  return cosmology.h0;
}

double get_sigma_8(input_cosmo_params_y3 ic){
  if (ic.A_s != cosmology.A_s){
    printf("cosmology changed before calling get_sigma_8\n");
    return -1.;
  }
  return cosmology.sigma_8;
}

double get_thetas(input_cosmo_params_y3 ic);
double get_thetas(input_cosmo_params_y3 ic){
  return cosmology.theta_s;
}

double xi_shear_tomo_sys(int pm, double theta, int nt, int z1, int z2)
{
  double xi;
  xi= xi_pm_TATT(pm,nt,z1,z2); //cosmo2D_real now includes NLA IA terms

  if(like.shearcalib==1) xi *=(1.0+nuisance.shear_calibration_m[z1])*(1.0+nuisance.shear_calibration_m[z2]);
return xi;
}


void compute_xip_1tomobin(double *data, int z1, int z2){
  static double *theta;
  if(theta==0){
    theta= create_double_vector(0, like.Ntheta-1);
    double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
    for(int i=0; i<like.Ntheta ; i++){
      double thetamin=exp(log(like.vtmin)+(i+0.0)*logdt);
      double thetamax=exp(log(like.vtmin)+(i+1.0)*logdt);
      theta[i] = 2./3.*(pow(thetamax,3.)-pow(thetamin,3.))/(pow(thetamax,2.)-pow(thetamin,2.));
    }
  }
  for (int i = 0; i < like.Ntheta; i++){
    data[i] = xi_shear_tomo_sys(1,theta[i],i,z1,z2);
  }
}

void compute_xim_1tomobin(double *data, int z1, int z2){
  static double *theta;
  if(theta==0){
    theta= create_double_vector(0, like.Ntheta-1);
    double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
    for(int i=0; i<like.Ntheta ; i++){
      double thetamin=exp(log(like.vtmin)+(i+0.0)*logdt);
      double thetamax=exp(log(like.vtmin)+(i+1.0)*logdt);
      theta[i] = 2./3.*(pow(thetamax,3.)-pow(thetamin,3.))/(pow(thetamax,2.)-pow(thetamin,2.));
    }
  }
  for (int i = 0; i < like.Ntheta; i++){
    data[i] = xi_shear_tomo_sys(-1,theta[i],i,z1,z2);
  }
}

void set_data_shear(double *theta, double *data, int start)
{
  int i,z1,z2,nz,j;
  for (nz = 0; nz < tomo.shear_Npowerspectra; nz++){
    z1 = Z1(nz); z2 = Z2(nz);
    for (i = 0; i < like.Ntheta; i++){
      if (mask(like.Ntheta*nz+i)){
        data[like.Ntheta*nz+i] = xi_shear_tomo_sys(1,theta[i],i,z1,z2);
      }
      if (mask(like.Ntheta*(tomo.shear_Npowerspectra+nz)+i)){
        data[like.Ntheta*(tomo.shear_Npowerspectra+nz)+i] = xi_shear_tomo_sys(-1,theta[i],i,z1,z2);
      }
    }
  }
}



void set_cosmology_params(double OMM, double NORM, double NS, double W0,double WA, double OMB, double OMNUh2, double H0, double MGSigma, double MGmu, double THETA_S)
{
  cosmology.Omega_m=OMM;
  cosmology.Omega_v= 1.0-cosmology.Omega_m;
  if (NORM < 1.e-7){
    cosmology.A_s = NORM;
    cosmology.sigma_8 = 0.;
  }
  else{
    cosmology.sigma_8=NORM;
    cosmology.A_s = 0.;
  }
  cosmology.theta_s = THETA_S;
  cosmology.n_spec= NS;
  cosmology.w0=W0;
  cosmology.wa=WA;
  cosmology.omb=OMB;
  if (H0> 0){
    cosmology.Omega_nu=OMNUh2/H0/H0;
  }
  else{cosmology.Omega_nu =0.0;}
  cosmology.h0=H0;
  cosmology.MGSigma =  MGSigma;
  cosmology.MGmu =  MGmu;
//  printf("cosmo: %e %e   %e %e   %e %e   %e %e\n",cosmology.Omega_m,cosmology.Omega_v,cosmology.A_s,cosmology.sigma_8,cosmology.theta_s,cosmology.h0,cosmology.omb,cosmology.n_spec);
//  printf(" %e %e   %e    %e %e\n",cosmology.w0,cosmology.wa, cosmology.Omega_nu,cosmology.MGSigma,cosmology.MGmu);
}



void set_nuisance_shear_calib(double *M){
  for (int i =0; i < tomo.shear_Nbin; i++){ nuisance.shear_calibration_m[i] = M[i];}
}

void set_nuisance_shear_photoz(double *SP){
  for (int i =0; i < tomo.shear_Nbin; i++){nuisance.bias_zphot_shear[i]=SP[i];}
}



void set_nuisance_ia_mpp(double *A1,double *A2,double *B_TA){
  //update to c1rhocrit
  //this should be udated in structs.c, but I didn't want to confuse other projects
  nuisance.c1rhocrit_ia = 0.01389;
  if (like.IA ==3 || like.IA ==5){
    for (int i = 0; i < tomo.shear_Nbin; i ++){
      nuisance.A_z[i] = A1[i];
      nuisance.A2_z[i] = A2[i];
      nuisance.b_ta_z[i] = B_TA[i];
      if (B_TA[i] == 0.){nuisance.b_ta_z[i] = B_TA[0];}
    }
  }
  if (like.IA==4 || like.IA==6){
    nuisance.A_ia = A1[0];
    nuisance.eta_ia = A1[1];
    nuisance.oneplusz0_ia = 1.62;
    nuisance.A2_ia = A2[0];
    nuisance.eta_ia_tt = A2[1];
    nuisance.b_ta_z[0] = B_TA[0];
  //  printf("set eta_ia = %.2f, eta_ia_tt=%.2f\n",nuisance.eta_ia,nuisance.eta_ia_tt);
    for (int i = 2; i < tomo.shear_Nbin; i ++){
      if (A1[i] != 0. || A2[i] != 0. || B_TA[i] != 0){
        printf("like_real_y3.c: nuisance.A_z[%d]=%e, nuisance.A2_z[%d]=%e, nuisance.b_ta[%d]=%e specified, while power-law evolution specified\nDid you forget to call init_IA_mpp()?\nEXIT\n",i, nuisance.A_z[i],i, nuisance.A2_z[i],i, nuisance.b_ta_z[i]);
        exit(1);}
    }
  }
}



double log_multi_like(double OMM, double NORM, double NS, double W0,double WA, double OMB, double OMNUh2, double H0, double MGSigma, double MGmu, double THETA_S,
  double *B1, double *B2, double *B_MAG,
  double *A1, double *A2, double *B_TA,
  double *SP, double *CP, double *M, double *PM)
{

  int i,j,k,m=0,l;
  static double *pred;
  static double *theta;
  static double *thetamin;
  static double *thetamax;
  static double *dtheta;
  static double logdt;
  double chisqr,a,log_L_prior=0.0;

  if(theta==0){
    pred= create_double_vector(0, like.Ndata-1);
    theta= create_double_vector(0, like.Ntheta-1);
    thetamin= create_double_vector(0, like.Ntheta-1);
    thetamax= create_double_vector(0, like.Ntheta-1);
    dtheta= create_double_vector(0, like.Ntheta-1);
    logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
    for(i=0; i<like.Ntheta ; i++){
      thetamin[i]=exp(log(like.vtmin)+(i+0.0)*logdt);
      thetamax[i]=exp(log(like.vtmin)+(i+1.0)*logdt);
      theta[i] = 2./3.*(pow(thetamax[i],3.)-pow(thetamin[i],3.))/(pow(thetamax[i],2.)-pow(thetamin[i],2.));
      dtheta[i]=thetamax[i]-thetamin[i];
    }
    like.theta=theta;
  }
  //flat prior on Omega_b*h^2 required by BBN constraints used in CLASS
  log_L_prior=0.0;
  set_cosmology_params(OMM,NORM,NS,W0,WA,OMB,OMNUh2,H0, MGSigma, MGmu, THETA_S);
  if (strcmp(pdeltaparams.runmode,"class")==0||strcmp(pdeltaparams.runmode,"CLASS")==0) {
    int status = 0;
    if (H0> 0 &&(OMB*H0*H0 >= 0.04 || OMB*H0*H0 <= 0.005)){printf("BBN\n"); return -1.e+15;}
    p_class(1.,1.,0,&status);
    if (status){printf("CLASS error\n"); return -1.e+15;}
  }
  set_nuisance_shear_calib(M);
  set_nuisance_shear_photoz(SP);
  set_nuisance_ia_mpp(A1,A2,B_TA);
  //gaussian priors on nuisance parameters
  if(like.wlphotoz!=0) log_L_prior+=log_L_wlphotoz();
  if(like.shearcalib==1) log_L_prior+=log_L_shear_calib();

  for (i=0; i<like.Ndata; i++){
    pred[i] = 0.;
  }

  set_data_shear(theta, pred, 0);

  chisqr=0.0;
  for (i=0; i<like.Ndata; i++){
      //printf("%d %d %e %e\n",i,mask(i),pred[i],data_read(1,i));
    for (j=0; j<like.Ndata; j++){
      if((mask(i)>0.1) && (mask(j)>0.1)){ //mask is either 0.0 or 1.0
        a=(pred[i]-data_read(1,i))*invcov_mask(1,i,j)*(pred[j]-data_read(1,j));
        //printf("%d %d %le %le\n",i,j,a,invcov_mask(1,i,j));
        chisqr=chisqr+a;
      }
    }
  //printf("%d %le %le\n",i,pred[i],data_read(1,i));
  }
  if (chisqr<0.0){
    printf("errror: chisqr < 0\n");
  }
  if (chisqr<-1.0) exit(EXIT_FAILURE);
  if (isnan(chisqr)){return -1.e+16;}
//  printf("%le %le\n",-0.5*chisqr,log_L_prior);
  return -0.5*chisqr+log_L_prior;
}


void compute_data_vector(char *filename,
  double OMM, double NORM, double NS, double W0,double WA, double OMB, double OMNUh2, double H0, double MGSigma, double MGmu, double THETA_S,
  double *B1, double *B2, double *B_MAG,
  double *A1, double *A2, double *B_TA,
  double *SP, double *CP, double *M, double *PM)
{
  if (like.Ndata == 0 || like.shear_shear + like.shear_pos +like.pos_pos == 0){
    printf("like_real_mpp_init.c: compute_data_vector called with like.Ndata = %d, like.shear_shear =%d, like.shear_pos = %d, like.pos_pos =%d\n",like.Ndata,like.shear_shear,like.shear_pos,like.pos_pos);
    printf("Nothing to be done.\nCosmoLike needs to be initialized first.\nEXIT\n");
    exit(1);
  }
  int i,j,k,m=0,l;
  static double *pred;
  static double *theta;
  static double *thetamin;
  static double *thetamax;
  static double *dtheta;
  static double logdt;
  double chisqr,a,log_L_prior=0.0;

  if(theta==0){
    pred= create_double_vector(0, like.Ndata-1);
    theta= create_double_vector(0, like.Ntheta-1);
    thetamin= create_double_vector(0, like.Ntheta-1);
    thetamax= create_double_vector(0, like.Ntheta-1);
    dtheta= create_double_vector(0, like.Ntheta-1);
    logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
    for(i=0; i<like.Ntheta ; i++){
      thetamin[i]=exp(log(like.vtmin)+(i+0.0)*logdt);
      thetamax[i]=exp(log(like.vtmin)+(i+1.0)*logdt);
      theta[i] = 2./3.*(pow(thetamax[i],3.)-pow(thetamin[i],3.))/(pow(thetamax[i],2.)-pow(thetamin[i],2.));
      dtheta[i]=thetamax[i]-thetamin[i];
      // printf("theta_%d: %le\n",i,theta[i]/constants.arcmin);
    }
    like.theta=theta; // like.theta ponts to memory of theta
  }

  /*for (l=0;l<like.Ntheta;l++){
    printf("%d %le\n",l,theta[l]/constants.arcmin);
  }*/
  set_cosmology_params(OMM,NORM,NS,W0,WA,OMB,OMNUh2,H0, MGSigma, MGmu, THETA_S);
  set_nuisance_shear_calib(M);
  set_nuisance_shear_photoz(SP);
  set_nuisance_ia_mpp(A1,A2,B_TA);
  set_data_shear(theta, pred, 0);

 // printf("%d\n",start);
  FILE *F;
//  char filename[300];
//  sprintf(filename,"datav/%s_%s_%s",survey.name,like.probes,details);
  F=fopen(filename,"w");
  for (i=0;i<like.Ndata; i++){
    fprintf(F,"%d %le\n",i,pred[i]);
    //printf("%d %le\n",i,pred[i]);
  }
  fclose(F);
}

void set_params(double OMM, double NORM, double NS, double W0,double WA, double OMB, double OMNUh2, double H0, double MGSigma, double MGmu, double THETA_S,
  double *B1, double *B2, double *B_MAG,
  double *A1, double *A2, double *B_TA,
  double *SP, double *CP, double *M, double *PM)
{


  set_cosmology_params(OMM,NORM,NS,W0,WA,OMB,OMNUh2,H0, MGSigma, MGmu, THETA_S);
  if (strcmp(pdeltaparams.runmode,"class")==0||strcmp(pdeltaparams.runmode,"CLASS")==0) {
    int status = 0;
    if (H0> 0 &&(OMB*H0*H0 >= 0.04 || OMB*H0*H0 <= 0.005)){printf("BBN limits exceeded\n"); exit(1);}
    p_class(1.,1.,0,&status);
    if (status){printf("CLASS error\n"); exit(1);}
  }
  set_nuisance_shear_calib(M);
  set_nuisance_shear_photoz(SP);
  set_nuisance_ia_mpp(A1,A2,B_TA);

}



void set_params_from_python(input_cosmo_params_y3 ic, input_nuisance_params_y3 in)
{
  double NORM;
  if (ic.A_s > 0. && ic.A_s < 1.e-5){NORM = ic.A_s;}
  else{NORM = ic.sigma_8;}
  if (NORM <= 0){
    printf("log_like_wrapper called with A_s = %e, sigma_8 =%e\nEXIT\n",ic.A_s,ic.sigma_8);
    exit(1);
  }
  set_params(ic.omega_m, NORM, ic.n_s, ic.w0, ic.wa, ic.omega_b,ic.omega_nuh2, ic.h0, ic.MGSigma, ic.MGmu,ic.theta_s,
    in.bias,in.bias2, in.b_mag,
    in.A_z, in.A2_z, in.b_ta,
    in.source_z_bias, in.lens_z_bias, in.shear_m, in.pm);
}

double log_like_wrapper(input_cosmo_params_y3 ic, input_nuisance_params_y3 in)
{
  double NORM;
  if (ic.A_s > 0. && ic.A_s < 1.e-5){NORM = ic.A_s;}
  else{NORM = ic.sigma_8;}
  if (NORM <= 0){
    printf("log_like_wrapper called with A_s = %e, sigma_8 =%e\nEXIT\n",ic.A_s,ic.sigma_8);
    exit(1);
  }
  return log_multi_like(ic.omega_m, NORM, ic.n_s, ic.w0, ic.wa, ic.omega_b,ic.omega_nuh2, ic.h0, ic.MGSigma, ic.MGmu,ic.theta_s,
    in.bias,in.bias2, in.b_mag,
    in.A_z, in.A2_z, in.b_ta,
    in.source_z_bias, in.lens_z_bias, in.shear_m, in.pm);
}

void write_datavector_wrapper(char *filename,input_cosmo_params_y3 ic, input_nuisance_params_y3 in)
{
  printf("write_datavector_wrapper: path to test data vector: %s\n",filename);

  double NORM;
  if (ic.A_s > 0. && ic.A_s < 1.e-5){NORM = ic.A_s;}
  else{NORM = ic.sigma_8;}
  if (NORM <= 0){
    printf("write_datavector_wrapper called with A_s = %e, sigma_8 =%e\nEXIT\n",ic.A_s,ic.sigma_8);
    exit(1);
  }
  compute_data_vector(filename, ic.omega_m, NORM, ic.n_s, ic.w0, ic.wa, ic.omega_b,ic.omega_nuh2, ic.h0, ic.MGSigma, ic.MGmu,ic.theta_s,
    in.bias,in.bias2, in.b_mag,
    in.A_z, in.A2_z, in.b_ta,
    in.source_z_bias, in.lens_z_bias, in.shear_m, in.pm);
}
