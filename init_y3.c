// these routines need to exposed to python
void init_cosmo_runmode(char *runmode);
void init_probes_5x2pt(char *probes);
void init_probes_real_mpp(char *probes);
void init_source_sample_mpp(char *multihisto_file, int Ntomo);
void init_lens_sample_mpp(char *multihisto_file, int Ntomo, double *b1, double *b2, double ggl_cut);
void init_binning_mpp(int Ntheta,double theta_min_arcmin, double theta_max_arcmin);
void init_IA_mpp(int N);
void init_sample_theta_s (void);
void init_data_real(char *COV_FILE, char *MASK_FILE, char *DATA_FILE);
int get_N_tomo_shear(void);
int get_N_tomo_clustering(void);
int get_N_ggl(void);
int get_N_data_masked(void);
double log_L_shear_calib();
double log_L_wlphotoz();
double log_L_clphotoz();
//void init_filenames(char *COV_FILE, char *MASK_FILE);

// routines internal to this file
int count_rows(char* filename,const char delimiter);
double invcov_mask(int READ, int ci, int cj);
double data_read(int READ, int ci);
double mask(int ci);
void init_lens_sample_();
void init_source_sample_();
void set_equal_tomo_bins(int Ntomo);
double log_L_shear_calib();

int get_N_tomo_shear(void){
  return tomo.shear_Nbin;
}
int get_N_tomo_clustering(void){
  return tomo.clustering_Nbin;
}
int get_N_ggl(void){
  return tomo.ggl_Npowerspectra;
}

int get_N_data_masked(void){
  int N = 0;
  for (int i = 0; i < like.Ndata; i++){ N += mask(i);}
  return N;
}

/*void init_filenames(char *COV_FILE, char *MASK_FILE){
  sprintf(like.MASK_FILE,"%s",MASK_FILE);
  printf("PATH TO MASK: %s\n",like.MASK_FILE);
  double init=mask(1);

  sprintf(like.COV_FILE,"%s",COV_FILE);
  printf("PATH TO COV: %s\n",like.COV_FILE);
}*/

void init_sample_theta_s (void){
  like.theta_s = 1;
}
void init_source_sample_()
{
  int i;
  if(strcmp(survey.sourcephotoz,"none")==0) redshift.shear_photoz=0;
  if(strcmp(survey.sourcephotoz,"voigt")==0) redshift.shear_photoz=1;
  if(strcmp(survey.sourcephotoz,"voigt_out")==0) redshift.shear_photoz=2;
  if(strcmp(survey.sourcephotoz,"gaussian")==0) redshift.shear_photoz=3;
  if(strcmp(survey.sourcephotoz,"multihisto")==0) redshift.shear_photoz=4;
  //printf("Source Sample Redshift Errors set to %s: redshift.shear_photoz=%d\n",survey.sourcephotoz,redshift.shear_photoz);

  if (redshift.shear_photoz!=4)  set_equal_tomo_bins(tomo.shear_Nbin);
  for (i=0;i<tomo.shear_Nbin; i++)
  {
    //printf("zmean_source=%f\n",zmean_source(i));
    nuisance.bias_zphot_shear[i]=0.0;
  }
  tomo.shear_Npowerspectra = tomo.shear_Nbin*(tomo.shear_Nbin+1)/2;
}

void init_lens_sample_()
{
  int i,j,n;

  if(strcmp(survey.lensphotoz,"none")==0) redshift.clustering_photoz=0;
  if(strcmp(survey.lensphotoz,"voigt")==0) redshift.clustering_photoz=1;
  if(strcmp(survey.lensphotoz,"voigt_out")==0) redshift.clustering_photoz=2;
  if(strcmp(survey.lensphotoz,"gaussian")==0) redshift.clustering_photoz=3;
  if(strcmp(survey.lensphotoz,"multihisto")==0) redshift.clustering_photoz=4;
  //  printf("Lens Sample Redshift Errors set to %s: redshift.clustering_photoz=%d\n",survey.lensphotoz,redshift.clustering_photoz);
  tomo.clustering_Npowerspectra = tomo.clustering_Nbin;
  n = 0;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    // printf("zmean_lens=%f\n",zmean(i));
    nuisance.bias_zphot_clustering[i]=0.0;
    for(j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
    }
  }
  tomo.ggl_Npowerspectra = n;

  //call test_kmax once to initialize look-up tables at reference cosmology
  test_kmax(1000.,1);
  // printf("end of lens sample init\n");
}

void set_equal_tomo_bins(int Ntomo)
{
  int k,j;
  double frac, zi;
  tomo.shear_Nbin=Ntomo;
  tomo.shear_Npowerspectra=(int) (Ntomo*(Ntomo+1)/2);
  zdistr_histo_1(0.1, NULL);
  int zbins =2000;
  double da = (redshift.shear_zdistrpar_zmax-redshift.shear_zdistrpar_zmin)/(1.0*zbins);
  double *sum;
  sum=create_double_vector(0, zbins);
  
  sum[0] = 0.0;
  for (k = 0, zi = redshift.shear_zdistrpar_zmin; k<zbins; k++,zi+=da){
    sum[k+1] = sum[k]+zdistr_histo_1(zi, NULL);
  }
  
  tomo.shear_zmin[0] = redshift.shear_zdistrpar_zmin;
  tomo.shear_zmax[tomo.shear_Nbin-1] = redshift.shear_zdistrpar_zmax;
  printf("\n");
  printf("Source Sample - Tomographic Bin limits:\n");
  for(k=0;k<tomo.shear_Nbin-1;k++){
    frac=(k+1.)/(1.*Ntomo)*sum[zbins-1];
    j = 0;
    while (sum[j]< frac){
      j++;
    }
    tomo.shear_zmax[k] = redshift.shear_zdistrpar_zmin+j*da;
    tomo.shear_zmin[k+1] = redshift.shear_zdistrpar_zmin+j*da;
    printf("min=%le max=%le\n",tomo.shear_zmin[k],tomo.shear_zmax[k]);
  }
  printf("min=%le max=%le\n",tomo.shear_zmin[tomo.shear_Nbin-1],tomo.shear_zmax[tomo.shear_Nbin-1]);
  printf("redshift.shear_zdistrpar_zmin=%le max=%le\n",redshift.shear_zdistrpar_zmin,redshift.shear_zdistrpar_zmax);
  free_double_vector(sum,0,zbins);
}

void init_cosmo_runmode(char *runmode)
{
  printf("\n");
  printf("-------------------------------------------\n");
  printf("Initializing Standard Runmode/Cosmology\n");
  printf("-------------------------------------------\n");
  set_cosmological_parameters_to_Planck_15_TT_TE_EE_lowP();
  sprintf(pdeltaparams.runmode,"%s",runmode);
  printf("pdeltaparams.runmode =%s\n",pdeltaparams.runmode);
}

void init_probes_5x2pt(char *probes)
{
  printf("\n");
  printf("------------------------------\n");
  printf("Initializing 5x2pt Probes\n");
  printf("------------------------------\n"); 

  sprintf(like.probes,"%s",probes);
  like.shear_shear=0;
  like.shear_pos=0;
  like.pos_pos=0;
  like.gk=0;
  like.ks=0;
  if (strcmp(probes,"all_2pt")==0){
    printf("init_mpp.c: called init_probes_5x2pt with outdated argument probes=%s\n",probes);
    printf("Update your code to use either \"3x2pt\" or \"3x2pt\"\n");
    exit(1);
  }
  char *type1 = "xi";
  if(strcmp(probes,"shear_shear")==0 || strstr(probes, type1) != NULL){
    like.shear_shear=1;
    printf("Shear-Shear computation initialized\n");
  }

  char *type2 = "wtheta";
  if(strcmp(probes,"pos_pos") == 0 || strstr(probes, type2) != NULL){
    like.pos_pos=1;
    printf("Position-Position computation initialized\n");
  }

  char *type3 = "ggl";
  char *type4 = "gammat";
  if(strcmp(probes,"shear_pos")==0 || strstr(probes, type3) != NULL || strstr(probes, type4) != NULL){
    if(tomo.ggl_Npowerspectra == 0){ void init_ggl_tomo();}
    like.shear_pos=1;
    printf("Position-Shear computation initialized\n");
  }
  char *type5 = "skappa";
  if (strstr(probes, type5) != NULL){
    like.ks = 1; 
    printf("Kappa-Shear computation initialized\n");
  }
  char *type6 = "gkappa";
  if (strstr(probes, type6) != NULL){
    like.gk = 1; 
    printf("Kappa-galaxy computation initialized\n");
  }
  if(strcmp(probes,"3x2pt")==0){
    if(tomo.ggl_Npowerspectra == 0){ void init_ggl_tomo();}
    like.shear_shear=1;
    like.shear_pos=1;
    like.pos_pos=1;
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
  } 
  
  if(strcmp(probes,"5x2pt")==0){
    if(tomo.ggl_Npowerspectra == 0){ void init_ggl_tomo();}
    like.shear_shear=1;
    like.shear_pos=1;
    like.pos_pos=1;
    like.ks = 1;
    like.gk = 1;
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
    printf("Kappa-Shear computation initialized\n");
    printf("Kappa-galaxy computation initialized\n");
  } 
  if(strcmp(probes,"ggl_cl")==0){
    if(tomo.ggl_Npowerspectra == 0){ void init_ggl_tomo();}
    like.shear_pos=1;
    like.pos_pos=1;
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
  }
  printf("like.pos_pos = %d, like.shear_pos = %d,like.shear_shear = %d,like.ks = %d,like.gl = %d\n\n",like.pos_pos,like.shear_pos,like.shear_shear,like.ks,like.gk); 
  printf("like.Ntheta = %d\n", like.Ntheta);
  printf("tomo.shear_Npowerspectra = %d\n", tomo.shear_Npowerspectra);
  printf("tomo.ggl_Npowerspectra = %d\n", tomo.ggl_Npowerspectra);
  printf("tomo.clustering_Npowerspectra = %d\n", tomo.clustering_Npowerspectra);

  like.Ndata=like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra + tomo.shear_Nbin+tomo.clustering_Nbin);

  printf("Total number of data points like.Ndata=%d\n",like.Ndata);
}

void init_probes_real_mpp(char *probes)
{
  printf("\n");
  printf("------------------------------\n");
  printf("Initializing 3x2pt Probes\n");
  printf("------------------------------\n"); 

  sprintf(like.probes,"%s",probes);
  like.shear_shear=0;
  like.shear_pos=0;
  like.pos_pos=0;
  like.gk=0;
  like.ks=0;
  if (strcmp(probes,"all_2pt")==0){
    printf("init_mpp.c: called init_probes_real_mpp with outdated argument probes=%s\n",probes);
    printf("Update your code to use \"3x2pt\" instead\n");
    exit(1);
  }

  char *type1 = "xi";
  if(strcmp(probes,"shear_shear")==0 || strstr(probes, type1) != NULL){
    like.shear_shear=1;
    printf("Shear-Shear computation initialized\n");
  }

  char *type2 = "wtheta";
  if(strcmp(probes,"pos_pos") == 0 || strstr(probes, type2) != NULL){
    like.pos_pos=1;
    printf("Position-Position computation initialized\n");
  }

  char *type3 = "ggl";
  char *type4 = "gammat";
  if(strcmp(probes,"shear_pos")==0 || strstr(probes, type3) != NULL || strstr(probes, type4) != NULL){
    if(tomo.ggl_Npowerspectra == 0){ void init_ggl_tomo();}
    like.shear_pos=1;
    printf("Position-Shear computation initialized\n");
  }

  if(strcmp(probes,"all_2pt")==0 || strcmp(probes,"3x2pt")==0){
    if(tomo.ggl_Npowerspectra == 0){ void init_ggl_tomo();}
    like.shear_shear=1;
    like.shear_pos=1;
    like.pos_pos=1;
    printf("Shear-Shear computation initialized\n");
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
  } 
  
  if(strcmp(probes,"ggl_cl")==0){
    if(tomo.ggl_Npowerspectra == 0){ void init_ggl_tomo();}
    like.shear_pos=1;
    like.pos_pos=1;
    printf("Shear-Position computation initialized\n");
    printf("Position-Position computation initialized\n");
  }
  printf("like.pos_pos = %d, like.shear_pos = %d,like.shear_shear = %d\n\n",like.pos_pos,like.shear_pos,like.shear_shear); 
  printf("like.Ntheta = %d\n", like.Ntheta);
  printf("tomo.shear_Npowerspectra = %d\n", tomo.shear_Npowerspectra);
  printf("tomo.ggl_Npowerspectra = %d\n", tomo.ggl_Npowerspectra);
  printf("tomo.clustering_Npowerspectra = %d\n", tomo.clustering_Npowerspectra);

  like.Ndata=like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);

  printf("Total number of data points like.Ndata=%d\n",like.Ndata);
}


void init_survey_mpp(char *surveyname, double area, double sigma_e){
  sprintf(survey.name,"%s",surveyname);
  survey.area = area;
  survey.sigma_e = sigma_e;
}



void init_source_sample_mpp(char *multihisto_file, int Ntomo)
{
  sprintf(redshift.shear_REDSHIFT_FILE,"%s",multihisto_file);
  redshift.shear_photoz=4;
  tomo.shear_Nbin = Ntomo;
  tomo.shear_Npowerspectra = tomo.shear_Nbin*(tomo.shear_Nbin+1)/2;
  printf("Source redshifts: multi-histo file %s, containing %d tomography bins\n",multihisto_file,tomo.shear_Nbin);
  for (int i=0;i<tomo.shear_Nbin; i++)
  {
    printf("bin %d: <z_s>=%f\n",i,zmean_source(i));
    //tomo.n_source[i]= n_source[i];
    nuisance.bias_zphot_shear[i]=0.0;
  }
  printf("init_source_sample_mpp complete\n");
}

void init_ggl_tomo(){
  if (tomo.clustering_Nbin ==0){
    printf("WARNING! init_mpp.c: init_ggl_tomo called while tomo.clustering_Nbin =0\n");
  }
  if (tomo.shear_Nbin ==0){
    printf("WARNING! init_mpp.c: init_ggl_tomo called while tomo.shear_Nbin =0\n");
  }
  int n = 0;
  for (int i = 0; i < tomo.clustering_Nbin; i++){
    for(int j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
      //printf("GGL combinations zl=%d zs=%d accept=%d; <z_l> = %.3f, <z_s> = %.3f\n",i,j,test_zoverlap(i,j), zmean(i),zmean_source(j));
    }
  }
  tomo.ggl_Npowerspectra = n;
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
}

void init_lens_sample_mpp(char *multihisto_file, int Ntomo, double *b1, double *b2, double ggl_cut)
{
  sprintf(redshift.clustering_REDSHIFT_FILE,"%s",multihisto_file);
  redshift.clustering_photoz=4;
  tomo.clustering_Nbin = Ntomo;
  tomo.clustering_Npowerspectra = tomo.clustering_Nbin;
  survey.ggl_overlap_cut = ggl_cut;
  printf("Lens redshifts: multi-histo file %s, containing %d tomography bins\n",multihisto_file,tomo.clustering_Nbin);
  pf_photoz(0.1,0);
  for (int i=0;i<tomo.clustering_Nbin; i++)
  {
    gbias.b1_function = & b1_per_bin;
 //   tomo.n_lens[i]= n_lens[i];
    gbias.b[i] = b1[i];
    gbias.b2[i] = b2[i];
    nuisance.bias_zphot_clustering[i]=0.0;
 //   printf("bin %d: <z_l>=%.3f, b_1=%.3f, b_2=%.3f\n",i,zmean(i),gbias.b[i],gbias.b2[i]);
  }
  init_ggl_tomo();
  printf("init_lens_sample_mpp complete\n");
}


void init_binning_mpp(int Ntheta,double theta_min_arcmin, double theta_max_arcmin){
  like.Ntheta=Ntheta;
  like.vtmin = theta_min_arcmin*constants.arcmin;
  like.vtmax = theta_max_arcmin*constants.arcmin;
  if (tomo.shear_Npowerspectra + tomo.clustering_Nbin + tomo.ggl_Npowerspectra == 0){
    printf("init_des_real.c: init_binning_mpp called with Npowerspectra = 0.\n EXIT!\n");
    exit(1);
  }
  printf("init_binning_mpp complete\n");
}

void init_IA_mpp(int N)
{  
  if(N ==3){
    like.IA = N;
    printf("Set like.IA =3: NLA with per-z-bin amplitude\nSupply one nuisance.A_z[] parameter per source tomo bin\n");
  }
  else if(N ==4){
    like.IA = N;
    printf("Set like.IA =4; NLA with powerlaw z-evolution\nSupply nuisance.A_ia and nuisance.eta_ia\n");
  }
  else if(N ==5){
    like.IA = N;
    printf("Set like.IA =5; TATT with per-z-bin amplitude\nSupply nuisance.A_z[], nuisance.b_ta_z[], nuisance.A2_z[] parameters per source tomo bin\n");
  }
  else if(N ==6){
    like.IA = N;
    printf("Set like.IA =6; TATT with powerlaw z-evolution\nSupply nuisance.A_ia, nuisance.eta_ia, nuisance.b_ta_z[0], nuisance.A2_ia and nuisance.eta_ia_tt\n");
  }
  else{
    printf("like.IA = %d not supported in des_mpp\nEXIT\n", N);
    exit(1);
  }
}


void init_data_real(char *COV_FILE, char *MASK_FILE, char *DATA_FILE)
{
  double init;
  printf("\n");
  printf("---------------------------------------\n");
  printf("Initializing data vector and covariance\n");
  printf("---------------------------------------\n");

  sprintf(like.MASK_FILE,"%s",MASK_FILE);
  printf("PATH TO MASK: %s\n",like.MASK_FILE);
  init=mask(1);

  sprintf(like.COV_FILE,"%s",COV_FILE);
  printf("PATH TO COV: %s\n",like.COV_FILE);
  init=invcov_mask(0,1,1);

  sprintf(like.DATA_FILE,"%s",DATA_FILE);
  printf("PATH TO DATA: %s\n",like.DATA_FILE);
  init=data_read(0,1);
}


double data_read(int READ, int ci)
{
  int i,intspace;
  static double *data = 0;
  
  if(READ==0 || data ==0){
    data  = create_double_vector(0, like.Ndata-1);      
    FILE *F;
    F=fopen(like.DATA_FILE,"r");
    for (i=0;i<like.Ndata; i++){  
      fscanf(F,"%d %le\n",&intspace,&data[i]);
    }
    fclose(F);
    printf("FINISHED READING DATA VECTOR\n");
  }    
  return data[ci];
}

double mask(int ci)
{
  int i,intspace;
  static double *mask =0;
  if(mask ==0){
    int N = 0;
    FILE *F;
    mask  = create_double_vector(0, like.Ndata-1); 
    double *maskc;
    maskc  = create_double_vector(0, like.Ndata-1); 
    if (strcmp(like.MASK_FILE,"none")==0){
      for (i=0;i<like.Ndata; i++){
        mask[i] = 1.0;
        maskc[i] =1.0;
      }
    }
    else{
      F=fopen(like.MASK_FILE,"r");
      if (!F){
        printf("init.c: invcov_mask: like.MASK_FILE = %s not found!\nEXIT!\n",like.MASK_FILE);
        exit(1);
      }
      for (i=0;i<like.Ndata; i++){
        fscanf(F,"%d %le\n",&intspace,&mask[i]);
        maskc[i] = mask[i];
        N += mask[i];
      }
     fclose(F);
     printf("%d bins within angular mask\n",N);
     printf("like.pos_pos = %d, like.shear_pos = %d,like.shear_shear = %d, like.ks = %d, like.gk = %d\n\n",like.pos_pos,like.shear_pos,like.shear_shear,like.ks,like.gk); 
    }
     int N3x2pt, N5x2pt;
     N3x2pt = like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);    
     N5x2pt = like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra+tomo.shear_Nbin+tomo.clustering_Nbin);
    //test whether Ndata assumes 3x2pt or 5x2pt format
    //if so, mask out probes excluded from the analysis
     if (N == N3x2pt || N== N5x2pt){
      if(like.shear_shear==0){
        printf("masking out shear-shear bins\n");
       for (i = 0; i< like.Ntheta*2*tomo.shear_Npowerspectra;i++){mask[i] = 0.;}
      }
      if(like.pos_pos==0){
        N = like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra);
        printf("masking out clustering bins\n");
        for (i = N; i< N+like.Ntheta*tomo.clustering_Npowerspectra;i++){mask[i] = 0.;}
      }
      if(like.shear_pos==0){
        N = like.Ntheta*2*tomo.shear_Npowerspectra;
        printf("masking out ggl bins\n");
        for (i = N; i <N+like.Ntheta*tomo.ggl_Npowerspectra; i++){mask[i] = 0.;}
      }
    }
    //test whether Ndata 5x2pt format
    //if so, mask out probes excluded from the analysis
    if (like.Ndata == N5x2pt){
      if(like.ks==0){
        printf("masking out shear x kappa bins\n");
        N = like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
        for (i = N; i <N+like.Ntheta*tomo.shear_Nbin; i++){mask[i] = 0.;}
      }
      if(like.gk==0){
        printf("masking out galaxies x kappa bins\n");
        N = like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra+tomo.shear_Nbin);
        for (i = N; i < N+like.Ntheta*tomo.clustering_Nbin; i++){mask[i] = 0.;}
      }
    }
    N = 0;
    for (i=0;i<like.Ndata; i++){
      //printf("mask(%d) = %.1f (was %.1f before probe cut)\n",i,mask[i],maskc[i]);
      N +=  mask[i];
    }
    printf("%d data points left after masking probes\n",N);
    if (N == 0){
      printf("init.c: mask: no data points left\nEXIT\n");
      exit(1);
    }
    printf("READ MASK FILE\n");
  }
  return mask[ci];
}

int count_rows(char* filename,const char delimiter){
  FILE *file = fopen (filename, "r" );
  char line[1000];
  if (file != NULL) {
    fgets(line,sizeof line,file);
    fclose(file);
  }
  else{
    printf("count_rows: file %s not found.\nEXIT\n",filename);
    exit(1);
  }
  int count = 1;
  char *p;

  p = line;
  while (*p != '\0')
  {
    if (*p == delimiter){
        while (*p == delimiter){p++;}
        count++;
    }
      p++;
    }
   return count;
}

double invcov_mask(int READ, int ci, int cj)
{
  static double **inv =0;

  if(READ==0 || inv ==0){
    double cov_G,cov_NG,doublespace,m;
    int i,j,intspace;
    gsl_matrix * cov   = gsl_matrix_calloc(like.Ndata, like.Ndata);      
    gsl_matrix * inv_c   = gsl_matrix_calloc(like.Ndata, like.Ndata);
    gsl_matrix_set_zero(cov);
    inv   = create_double_matrix(0, like.Ndata-1, 0, like.Ndata-1);      
      
    FILE *F;
   int n_rows =count_rows(like.COV_FILE,' ');
   F=fopen(like.COV_FILE,"r");
   if (!F){printf("init.c: invcov_mask: like.COV_FILE = %s not found!\nEXIT!\n",like.COV_FILE);exit(1);}
   switch (n_rows){
    case 3: while (fscanf(F,"%d %d %le\n", &i, &j, &cov_G) ==3) {
            m = 1.0;
            if (i < like.Ndata && j < like.Ndata){
            // apply mask to off-diagonal covariance elements
              if (i!=j){m = mask(i)*mask(j);}
            //  printf("%d %d (/%d) %e\n",i,j,like.Ndata,cov_G);
              gsl_matrix_set(cov,i,j,(cov_G)*m);
              gsl_matrix_set(cov,j,i,(cov_G)*m);
            }
          } break;
    case 4: while (fscanf(F,"%d %d %le %le\n", &i, &j, &cov_G, &cov_NG) ==4) {
            m = 1.0;
            if (i < like.Ndata && j < like.Ndata){
            // apply mask to off-diagonal covariance elements
              if (i!=j){m = mask(i)*mask(j);}
              //printf("%d %d (/%d) %e %e\n",i,j,like.Ndata,cov_G,cov_NG);
              gsl_matrix_set(cov,i,j,(cov_G+cov_NG)*m);
              gsl_matrix_set(cov,j,i,(cov_G+cov_NG)*m);
            }
          } break;
    case 10: while (fscanf(F,"%d %d %le %le %d %d %d %d %le %le\n", &i, &j, &doublespace, &doublespace,&intspace,&intspace,&intspace,&intspace,&cov_G,&cov_NG) ==10) {
            m = 1.0;
            if (i < like.Ndata && j < like.Ndata){
            // apply mask to off-diagonal covariance elements
              if (i!=j){m = mask(i)*mask(j);}
              //printf("%d %d (/%d) %e %e\n",i,j,like.Ndata,cov_G,cov_NG);
              gsl_matrix_set(cov,i,j,(cov_G+cov_NG)*m);
              gsl_matrix_set(cov,j,i,(cov_G+cov_NG)*m);
            }
          } break;
    default: printf("init.c:invcov_mask: covariance file %s has %d columns - unsupported format!\nEXIT\n",like.COV_FILE,n_rows);exit(1);
   }
   fclose(F);
   printf("READ COV_FILE\n");
   //SVD_inversion(cov,inv_c,like.Ndata);
  invert_matrix_colesky(cov);
   // printf("cov inverted\n");
   for (i=0;i<like.Ndata; i++){
    for (j=0;j<like.Ndata; j++){
      //apply mask again, to make sure numerical errors in matrix inversion don't cause problems...
      //also, set diagonal elements corresponding to datavector elements outside mask to zero, so that these elements don't contribute to chi2
      inv[i][j] =gsl_matrix_get(cov,i,j)*mask(i)*mask(j);
    }
   }
   gsl_matrix_free(cov);
   gsl_matrix_free(inv_c);
   printf("FINISHED BUILDING INV COV\n");
 }    
 return inv[ci][cj];
}

double log_L_shear_calib()
{
  int i;
  double log_L = 0.; 
  for (i=0;i<tomo.shear_Nbin; i++){
    if (prior.shear_calibration_m[i][1] == 0.){
      printf("external_prior.c: called log_L_clphotoz while prior.shear_calibration_m[%d][1] not set.\nEXIT\n",i);
      exit(1);
    }
    log_L -= pow((nuisance.shear_calibration_m[i] - prior.shear_calibration_m[i][0])/ prior.shear_calibration_m[i][1],2.0);
  }
  return 0.5*log_L;
}

double log_L_wlphotoz()
{
  int i;
  double log_L = 0.;
  for (i=0;i<tomo.shear_Nbin; i++){
    if (prior.bias_zphot_shear[i][1] == 0.){
      printf("external_prior.c: called log_L_wlphotoz while prior.bias_zphot_shear[%d][1] not set.\nEXIT\n",i);
      exit(1);
    }
    log_L -=  pow((nuisance.bias_zphot_shear[i] - prior.bias_zphot_shear[i][0])/ prior.bias_zphot_shear[i][1],2.0);
  }
  return 0.5*log_L;
}

double log_L_clphotoz()
{
  int i;
  double log_L = 0.;
  for (i=0;i<tomo.clustering_Nbin; i++){
      if (prior.bias_zphot_clustering[i][1] == 0.){
      printf("external_prior.c: called log_L_clphotoz while prior.bias_zphot_clustering[%d][1] not set.\nEXIT\n",i);
      exit(1);
    }
    log_L -= pow((nuisance.bias_zphot_clustering[i] - prior.bias_zphot_clustering[i][0])/ prior.bias_zphot_clustering[i][1],2.0);
  }
  return 0.5*log_L;
}