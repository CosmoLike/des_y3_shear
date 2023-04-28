#include "like_real_y3.c"
void test_Xis()
{
  clock_t begin, end;
  double time_spent;
  double NORM;
  like.IA = 6;
/* sample parameter values */
  nuisance.A_ia = 0.7;
  nuisance.A2_ia = -1.36;
  nuisance.b_ta_z[0] =1.0;

  nuisance.eta_ia = -1.7;
  nuisance.eta_ia_tt = -2.5;
  nuisance.oneplusz0_ia = 1.62;

/* here, do your time-consuming job */
  init_cosmo_runmode("CLASS");
// init_cosmo_runmode("emu");
//  init_bary("owls_AGN");
  cosmology.Omega_m   = 0.3;
  //NORM    = 0.8433;
  NORM  = 2.19e-09;//0.8433;
  cosmology.n_spec    = 0.97;

  cosmology.w0=-1.0;
  cosmology.wa=0.;
  cosmology.omb=0.048;
  cosmology.h0=0.69;
  double Omega_nuh2 = 0.00083;//0.001743331232934258;
  cosmology.Omega_v=1.0-cosmology.Omega_m;
  double b1[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},b_mag[10], b2[10] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},bias_photoz_l[10],sigma_b_photoz_l[10],pm[10];;
  double A_1[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}, A_2[10] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},b_ta[10] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  A_1[0]= nuisance.A_ia; A_1[1] = nuisance.eta_ia;
  A_2[0]= nuisance.A2_ia; A_2[1] = nuisance.eta_ia_tt;
  b_ta[0] = nuisance.b_ta_z[0];
  double mean_m[10]={-0.0058, -0.0104, -0.0255, -0.032,0.0,0.0,0.0,0.0,0.0,0.0};
  double sigma_m[10]={0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03};
  double bias_photoz_s[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double sigma_b_photoz_s[10]={0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08};


/**** current n(z) for comparison with latest cosmosis config ****/
init_source_sample_mpp("./zdistris/nz_source_Y3_unblinded_final.txt",4);

init_lens_sample_mpp("./zdistris/nz_lens_Y3_unblinded_10_26_20.txt",5,b1,b2,0.0);


  init_binning_mpp(20,2.5,250.);
  init_probes_real_mpp("3x2pt");

  set_shear_priors_mpp(mean_m,sigma_m);
  //set_wlphotoz_priors_mpp(bias_photoz_s,sigma_b_photoz_s);
  //set_clphotoz_priors_mpp(bias_photoz_l,sigma_b_photoz_l);

//  set_LF_GAMA();
//  survey.m_lim = 23.0;
//  nuisance.beta_ia = 1.1;

  //gbias.b1_function = &b1_per_bin_pass_evolv;
  sprintf(like.MASK_FILE,"%s","none");
  printf("PATH TO MASK: %s\n",like.MASK_FILE);
  begin = clock();
  char datavfile[200];
//  sprintf(datavfile,"datav/xi_Y3_baseline_v040_no_nu_TT");
//  sprintf(datavfile,"datav/xi_Y3_baseline_v040+bnl+owlsAGN+1h_no_nu");
  sprintf(datavfile,"datav/xi_Y3_shear_model.txt");
  compute_data_vector(datavfile,cosmology.Omega_m,NORM ,cosmology.n_spec,cosmology.w0,cosmology.wa,cosmology.omb,Omega_nuh2,cosmology.h0,0.0,0.0,cosmology.theta_s,
    b1, b2, b_mag,
    A_1,A_2, b_ta,
    bias_photoz_s, //source photo-z bias
    bias_photoz_l, //lens photo-z bias
    mean_m, //shear calibration
    pm //point mass
    );
      end = clock();

  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  for (double lgl = 1; lgl<3.; lgl+=0.1){
    double l = exp(lgl);
//    printf("%e  %e %e\n", l,dC_deflection(l,0,0)/dC_reduced_shear(l,0,0),dC_deflection(l,3,3)/dC_reduced_shear(l,3,3));//Delta_P_reduced_shear_tab(k,a)/Delta_P_deflection_tab(k,a));

  }
  //set_LF_GAMA();//
/*  for (double z = 0.05; z<3.0; z+=0.05){
    printf("%e  %e\n",z, A_IA_Joachimi(1./(1+z)));///C1_TA(1./(1+z),0));
  }*/
/*   printf("time spent %le\n",time_spent);
   for (double l = 10; l <20000; l*=1.2){
     printf("%e  %e %e %e %e\n",l,C_reduced_shear(l,0,0)/C_EE_TATT(l,0,0),C_reduced_shear(l,1,1)/C_EE_TATT(l,1,1),C_reduced_shear(l,2,2)/C_EE_TATT(l,2,2),C_reduced_shear(l,3,3)/C_EE_TATT(l,3,3));}
*/
/*for (int i = 0; i< 20; i++){
  printf("%d   %e %e %e %e\n",i,w_source_clustering(i,0,0),w_source_clustering(i,1,1),w_source_clustering(i,2,2),w_source_clustering(i,3,3));
}*/
//  printf("sigma8 = %.4f\n",cosmology.sigma_8);

}

//zdistr_photoz(1./a-1.,(int)nz)*hoverh0(a)
int main(void)
{
  test_Xis();
  return 0;
}
