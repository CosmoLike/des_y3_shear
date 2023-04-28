
cfftlog_dir := ../cosmolike_core/cfftlog/

cfftlog := $(cfftlog_dir)cfftlog.c $(cfftlog_dir)utils.c $(cfftlog_dir)utils_complex.c

cfastpt_dir := ../cosmolike_core/cfastpt/

cfastpt := $(cfastpt_dir)cfastpt.c $(cfastpt_dir)utils.c $(cfastpt_dir)utils_complex.c

opt_home := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers -I/usr/local/include -L/usr/local/lib -L../cosmolike_core/class  -g -std=gnu99 -lgsl -lfftw3 -lgslcblas -lclass -lm

opt_home_CLASS29 := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers -I/usr/local/include -L/usr/local/lib -L../cosmolike_core/class_v29  -g -std=gnu99 -lgsl -lfftw3 -lgslcblas -lclass -lm -DCLASS_V29

opt_ocelote := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I/cm/shared/uaapps/gsl/2.1/include -L/cm/shared/uaapps/gsl/2.1/lib \
-lgsl -lfftw3 -lgslcblas -lm -O0 -g -O3 \
-std=gnu99 -ffast-math -funroll-loops -L../cosmolike_core/class -lclass

opt_puma := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers \
-I/opt/ohpc/pub/libs/gnu8/gsl/2.6/include -L/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib \
-lgsl -lfftw3 -lgslcblas -lm -O0 -g -O3 \
-std=gnu99 -ffast-math -funroll-loops -L../cosmolike_core/class -lclass


home:
	make home_lib
	make home_test

home_lib:
	gcc -shared -o like_real_y3.so -fPIC like_real_y3.c $(cfftlog) $(cfastpt) $(opt_home)

home_test:
	gcc  like_test.c -o ./test_desy3 $(cfftlog) $(cfastpt) $(opt_home)

class_v29:
	gcc  like_test.c -o ./test_desy3 $(cfftlog) $(cfastpt) $(opt_home_CLASS29)
cls:
	gcc  write_cl.c -o ./write_cls $(opt_home) $(cfftlog) $(cfastpt)
	./write_cls

FLASK_Cells:
	gcc  FLASK_Cells.c	-o ./FLASK_Cells $(opt_home) $(cfftlog) $(cfastpt)
	./FLASK_Cells


ocelote:
	make ocelote_lib
	make ocelote_test
	make ocelote_cov


ocelote_lib:
	gcc -shared -o like_real_y3.so -fPIC like_real_y3.c $(opt_ocelote)  $(cfftlog) $(cfastpt)

ocelote_test:
	gcc  like_test.c -o ./test_desy3 $(opt_ocelote) $(cfftlog) $(cfastpt)

ocelote_cov:
	gcc  compute_covariances_real_Y3.c -o ./compute_covariances_real_Y3 $(opt_ocelote) $(cfftlog) $(cfastpt)


puma: 
	make puma_lib
	make puma_test
	make puma_cov

puma_lib:
	gcc -shared -o like_real_y3.so -fPIC like_real_y3.c $(opt_puma)  $(cfftlog) $(cfastpt)

puma_test:
	gcc  like_test.c -o ./test_desy3 $(opt_puma) $(cfftlog) $(cfastpt)

puma_cov:
	gcc  compute_covariances_real_Y3.c -o ./compute_covariances_real_Y3 $(opt_puma) $(cfftlog) $(cfastpt)


class:
	cd ../cosmolike_core/class; $(MAKE)


