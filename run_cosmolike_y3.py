from cosmolike_libs_y3 import *

import os
import numpy as np

def run_cosmolike(params, NOMPI, pool=None):
    do_xip = "xip" in params['twoptnames']
    do_xim = "xim" in params['twoptnames']
    do_wtheta = "wtheta" in params['twoptnames']
    do_gammat = "gammat" in params['twoptnames']

    if do_xip or do_xim:
        assert do_xip and do_xim, "Can only run cosmolike with both of xip and xim or neither.  You could modify the mask if you really want to."

    path ="./"
    if "base_dir" in params:
        path = params['base_dir']

    cov_file = path + params['cov_file']
    data_file = path + params['data_file']
    source_nz = path + params['source_nz']
    lens_nz = path + params['lens_nz']
    if "new_mask_file" in params:
        mask_file = params['new_mask_file']
    else:
        mask_file = path + params['mask_file']
    chain_file = path + params['chain_file']

    ntomo_source = params['ntomo_source']
    ntomo_lens = params['ntomo_lens']
    ntheta = params['tbins']
    theta_min = params['tbounds'][0]
    theta_max = params['tbounds'][1]



    runmode ="Halofit"
    if 'run_mode' in params:
        runmode = params['run_mode']
    probes = "".join(params['twoptnames'])
    initcosmo(runmode.encode('utf-8'))
    initsources(source_nz.encode('utf-8'), ntomo_source)
    initlenses(lens_nz.encode('utf-8'), ntomo_lens, Double10(), Double10(),0.0)
    initbins(ntheta, theta_min, theta_max)
    #only for plotting
    #set_pointer_2PCF(ntheta)
    initprobes(probes.encode('utf-8'))

    initdata_real(cov_file.encode('utf-8'), mask_file.encode('utf-8'), data_file.encode('utf-8'))

    (varied_params,
        cosmo_min, cosmo_fid, cosmo_max,
        nuisance_min, nuisance_fid, nuisance_max) = parse_priors_and_ranges(params)

    #test_datavector = chain_file+".test_datavector"
    # if (get_N_data() != params['mask_checksum']):
    #     print("Number of data points computed from yaml file = %d; N_data from maskfile = %d",params['mask_checksum'],get_N_data())
    #     exit(1)
#    cosmo_fid.print_struct()
#    nuisance_fid.print_struct()
    #write_cosmolike_datavector(test_datavector, cosmo_fid, nuisance_fid)
    print ("will sample over", varied_params)

    nthreads = 1
    if 'n_threads' in params:
        nthreads = params['n_threads']
    iterations = 4000
    if 'iterations' in params:
        iterations = params['iterations']
    nwalkers = 560
    if 'nwalkers' in params:
        nwalkers = params['nwalkers']
    try:
    	import mpi4py
    except:
        NOMPI = True
        print("mpi4py not found\n Run emcee with 1 thread for testing!\n")
    NOMPI = False
    if (NOMPI):
        sample_main(varied_params,
                    iterations, #iterations
                    nwalkers, #walkers
                    1, #nthreads
                    chain_file, #output filename
                    cosmo_min, cosmo_fid, cosmo_max,  #cosmo flat priors
                    nuisance_min, nuisance_fid, nuisance_max #nuisance flat priors
                    )
    else:
        from schwimmbad import MPIPool
        sample_main(varied_params,iterations, #iterations
                    nwalkers, #walkers
                    nthreads, #nthreads
                    chain_file, #output filename
                    cosmo_min, cosmo_fid, cosmo_max,  #cosmo flat priors
                    nuisance_min, nuisance_fid, nuisance_max, #nuisance flat priors
                    pool=MPIPool())

def parse_priors_and_ranges(params):
    do_xip = "xip" in params['twoptnames']
    do_xim = "xim" in params['twoptnames']
    do_wtheta = "wtheta" in params['twoptnames']
    do_gammat = "gammat" in params['twoptnames']

    ntomo_source = params['ntomo_source']
    ntomo_lens = params['ntomo_lens']

    cosmo_min = InputCosmologyParams()
    cosmo_fid = InputCosmologyParams().fiducial()
    cosmo_max = InputCosmologyParams()
    cosmo_names = InputCosmologyParams().names()

    #Loop through the cosmological parameters
    #picking
    varied_params = []
    norm, scale = check_cosmo_duplicates(params)
    for p in cosmo_names:
        if p in ["A_s","sigma_8","h0","theta_s"]:
            if p+"_range" in params:
                p_range = params[p+"_range"]
                min_val, fid_val, max_val, is_var = parse_range(p_range)
                setattr(cosmo_fid, p, fid_val)
                if is_var:
                    varied_params.append(p)
                    setattr(cosmo_min, p, min_val)
                    setattr(cosmo_max, p, max_val)
        else:
            p_range = params[p+"_range"]
            min_val, fid_val, max_val, is_var = parse_range(p_range)
            setattr(cosmo_fid, p, fid_val)
            if is_var:
                varied_params.append(p)
                setattr(cosmo_min, p, min_val)
                setattr(cosmo_max, p, max_val)


    nuisance_min = InputNuisanceParams()
    nuisance_fid = InputNuisanceParams().fiducial()
    nuisance_max = InputNuisanceParams()

    gaussian_prior_params = ["source_z_bias", "lens_z_bias", "shear_m"]

    shear_m_mean = Double10()
    shear_m_sigma = Double10()
    source_z_bias_mean = Double10()
    source_z_bias_sigma = Double10()
    lens_z_bias_mean = Double10()
    lens_z_bias_sigma = Double10()


    if do_wtheta or do_gammat:
        parse_nuisance_flat_prior(params, "bias", ntomo_lens, nuisance_min, nuisance_fid, nuisance_max, varied_params)
        parse_nuisance_flat_prior(params, "bias2", ntomo_lens, nuisance_min, nuisance_fid, nuisance_max, varied_params)
        parse_nuisance_flat_prior(params, "b_mag", ntomo_lens, nuisance_min, nuisance_fid, nuisance_max, varied_params)
        is_var = parse_nuisance_gaussian_prior(params, "lens_z_bias", ntomo_lens, nuisance_fid, lens_z_bias_mean, lens_z_bias_sigma, varied_params)
        if is_var:
            setprior_clusteringphotoz(lens_z_bias_mean, lens_z_bias_sigma)

    if do_xip or do_xim or do_gammat:

        is_var = parse_nuisance_gaussian_prior(params, "source_z_bias", ntomo_source, nuisance_fid, source_z_bias_mean, source_z_bias_sigma, varied_params)
        if is_var:
            setprior_wlphotoz(source_z_bias_mean, source_z_bias_sigma)
        is_var = parse_nuisance_gaussian_prior(params, "shear_m", ntomo_source, nuisance_fid, shear_m_mean, shear_m_sigma, varied_params)
        if is_var:
            setprior_m(shear_m_mean,shear_m_sigma)
        #test which IA model
        #power-law redshift parameterization
        TATT_power_law = parse_IA_TATT_power_law_flat_prior(params, nuisance_min, nuisance_fid, nuisance_max, varied_params)
        #if not, try per-bin redshift parameterization
        if TATT_power_law is False:
            is_var_NLA = parse_nuisance_flat_prior(params, "A_z", ntomo_source, nuisance_min, nuisance_fid, nuisance_max, varied_params)
            is_var_TT = parse_nuisance_flat_prior(params, "A2_z", ntomo_source, nuisance_min, nuisance_fid, nuisance_max, varied_params)
            if is_var_NLA:
                initia(3)
                is_var_b_ta = parse_nuisance_flat_prior(params, "b_ta", ntomo_source, nuisance_min, nuisance_fid, nuisance_max, varied_params)
                if is_var_b_ta:
                    initia(5)
            if is_var_TT:
                initia(5)
    if do_gammat:
        parse_nuisance_flat_prior(params, "pm", ntomo_lens, nuisance_min, nuisance_fid, nuisance_max, varied_params)
    return (varied_params,
            cosmo_min, cosmo_fid, cosmo_max,
            nuisance_min, nuisance_fid, nuisance_max)
def check_cosmo_duplicates(params):
### power spectrum amplitude - either A_s or sigma_8
    norm = 0
    p = "A_s"
    if p+"_range" in params:
        norm +=1
    p = "sigma_8"
    if p+"_range" in params:
        norm +=1
    if (norm ==2):
        print ("yaml file specifies both A_s and sigma_8!\nEXIT")
        exit()
    if (norm == 0):
        print ("yaml file specifies neither A_s nor sigma_8!\nEXIT")
        exit()
### distance scale parameter - either h0 or theta_s
    scale = 0
    p = "h0"
    if p+"_range" in params:
        scale +=1
    p = "theta_s"
    if p+"_range" in params:
        scale +=1
    if (scale ==2):
        print ("\nyaml file specifies both h0 and theta_s\nEXIT")
        exit()
    if (scale == 0):
        print ("yaml file specifies neither h0 nor theta_s!\nEXIT")
        exit()
    return norm, scale


def parse_nuisance_gaussian_prior(params, p, nbin, nuisance_fid, mean_ptr, sigma_ptr, varied_params):
    mean_vector = params[p+"_mean"]
    sigma_vector = params.get(p+"_sigma")
    is_var = (sigma_vector != None)

    if len(mean_vector) != nbin:
        raise ValueError("Wrong length for parameter {}_mean - should be {}".format(p, nbin))
    if is_var and len(sigma_vector) != nbin:
        raise ValueError("Wrong length for parameter {}_sigma - should be {}".format(p, nbin))

    for i in range(nbin):
        getattr(nuisance_fid, p)[i] = mean_vector[i]

    if is_var:
        for i in range(nbin):
            varied_params.append("{}_{}".format(p,i))
            mean_ptr[i] = mean_vector[i]
            sigma_ptr[i] = sigma_vector[i]
    return is_var


def parse_nuisance_flat_prior(params, p, nbin, nuisance_min, nuisance_fid, nuisance_max, varied_params):
    set = 0
    if p+"_range" in params:
        set = 1
        p_range = params[p+"_range"]
        min_val, fid_val, max_val, is_var = parse_range(p_range)
        for i in range(nbin):
            getattr(nuisance_fid, p)[i] = fid_val
            if is_var:
                varied_params.append("{}_{}".format(p,i))
                getattr(nuisance_min, p)[i] = min_val
                getattr(nuisance_max, p)[i] = max_val
    if p+"_fiducial" in params:
        set = 1
        values = params[p+"_fiducial"]
        for i in range(nbin):
            getattr(nuisance_fid, p)[i] = values[i]
    if (set == 0):
        print ("run_cosmolike_mpp.py: %s not found in yaml file, use cosmolike_libs_y3 default value" %(p))
        for i in range(nbin):
            print ("cosmolike_libs_real_mpp.py: %s[%d] =%e" %(p,i,getattr(nuisance_fid,p)[i]))
        is_var = 0
    return is_var

def parse_IA_TATT_power_law_flat_prior(params, nuisance_min, nuisance_fid, nuisance_max, varied_params):
    power_law = False
    if (("A_ia_range" in params) & ("A2_z_range" in params)):
        print("run_cosmolike_y3.py:parse_IA_TATT_power_law_flat_prior: mix of power-law NLA-IA and per-bin TT-IA z-dependence not supported")
        exit(1)
    if (("A_z_range" in params) & ("A2_ia_range" in params)):
        print("run_cosmolike_y3.py:parse_IA_TATT_power_law_flat_prior: mix of power-law TT-IA and per-bin NLA-IA z-dependence not supported")
        exit(1)
    var = 0
#### NLA/TA parameters
    p ="A_z"
    if "A_ia_range" in params:
        power_law = True
        print("Power-law NLA/TA-IA parameterization")
        initia(6)
        p_range = params["A_ia_range"]
        min_val, fid_val, max_val, is_var = parse_range(p_range)
        i = 0
        getattr(nuisance_min, p)[i] = min_val
        getattr(nuisance_fid, p)[i] = fid_val
        getattr(nuisance_max, p)[i] = max_val
        if is_var:
            varied_params.append("{}_{}".format(p,i))
            var = 1
    ### only check for eta_ia if A_ia is set, as default A_ia = 0.
        if "eta_ia_range" in params:
            p_range = params["eta_ia_range"]
            min_val, fid_val, max_val, is_var = parse_range(p_range)
            i = 1
            getattr(nuisance_min, p)[i] = min_val
            getattr(nuisance_fid, p)[i] = fid_val
            getattr(nuisance_max, p)[i] = max_val
            if is_var:
                varied_params.append("{}_{}".format(p,i))
                var +=1
        else:
            print ("run_cosmolike_y3.py: eta_ia not found in yaml file, use cosmolike_libs_y3 defauly value")
    else:
        if "A_z_range" not in params:
            print ("run_cosmolike_y3.py: NLA-IA amplitude not found in yaml file, use cosmolike_libs_y3 default value")
### TT parameters
    p ="A2_z"
    if "A2_ia_range" in params:
        power_law = True
        print("Power-law TT-IA parameterization")
        initia(6)
        p_range = params["A2_ia_range"]
        min_val, fid_val, max_val, is_var = parse_range(p_range)
        i = 0
        getattr(nuisance_min, p)[i] = min_val
        getattr(nuisance_fid, p)[i] = fid_val
        getattr(nuisance_max, p)[i] = max_val
        if is_var:
            varied_params.append("{}_{}".format(p,i))
            var = 1
    ### only check for eta_ia if A_ia is set, as default A_ia = 0.
        if "eta_ia_tt_range" in params:
            p_range = params["eta_ia_tt_range"]
            min_val, fid_val, max_val, is_var = parse_range(p_range)
            i = 1
            getattr(nuisance_min, p)[i] = min_val
            getattr(nuisance_fid, p)[i] = fid_val
            getattr(nuisance_max, p)[i] = max_val
            if is_var:
                varied_params.append("{}_{}".format(p,i))
                var +=1
        else:
            print ("run_cosmolike_y3.py: eta_ia_tt not found in yaml file, use cosmolike_libs_y3 defauly value")

#### TA bias parameter
    p ="b_ta"
    if "b_ta_noz_range" in params:
            p_range = params["b_ta_noz_range"]
            min_val, fid_val, max_val, is_var = parse_range(p_range)
            i = 0
            getattr(nuisance_min, p)[i] = min_val
            getattr(nuisance_fid, p)[i] = fid_val
            getattr(nuisance_max, p)[i] = max_val
            if is_var:
                varied_params.append("{}_{}".format(p,i))
                var +=1
            else:
                print ("run_cosmolike_y3.py: b_ta not found in yaml file, use cosmolike_libs_y3 defauly value")

    return power_law

def parse_range(p_range):
    "return min, fid, max, is_varied"
    if np.isscalar(p_range):
        min_val = p_range
        fid_val = p_range
        max_val = p_range
        is_var  = False
    elif len(p_range)==1:
        min_val = p_range[0]
        fid_val = p_range[0]
        max_val = p_range[0]
        is_var  = False
    else:
        if len(p_range)!=3:
            raise ValueError("Must specify 1 or 3 elements in param ranges")
        min_val = p_range[0]
        fid_val = p_range[1]
        max_val = p_range[2]
        is_var  = True

    return min_val, fid_val, max_val, is_var
