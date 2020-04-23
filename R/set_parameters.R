set_param <- function(enm_model = "ming_wall", ...) {
  if (enm_model == "ming_wall") {
    param <- set_param_ming_wall(...)
  } else if (enm_model == "anm") {
    param <- set_param_anm(...)
  } else {
    print("Warning, undefined option in set_param")
    param <- NA
  }
  param
}


# set model parameter lists
set_param_ming_wall <- function(nmut_per_site,
                                v0,
                                mut_sd_min, mut_dl_sigma,
                                fit_dg_thr, beta,
                                fix_model, fix_mut_rate, fix_n_eff) {
  param <- list(
    trj = list(nmut_per_site = nmut_per_site),
    enm = list(model = "ming_wall", v0 = v0,  d_max = 10.5),
    mut = list(model = "lfenm",  sd_min = mut_sd_min, dl_sigma = mut_dl_sigma),
    fit = list(model = "asm", dg_thr = fit_dg_thr, beta = beta),
    fix = list(model = fix_model, mut_rate = fix_mut_rate,  n_eff = fix_n_eff)
  )
  param
}

# set model parameter lists
set_param_anm <- function(nmut_per_site,
                          v0,
                          mut_sd_min, mut_dl_sigma,
                          fit_dg_thr, beta,
                          fix_model, fix_mut_rate, fix_n_eff) {
  param <- list(
    trj = list(nmut_per_site = nmut_per_site),
    enm = list(model = "anm", v0 = v0,  d_max = 12.5),
    mut = list(model = "lfenm",  sd_min = mut_sd_min, dl_sigma = mut_dl_sigma),
    fit = list(model = "asm", dg_thr = fit_dg_thr, beta = beta),
    fix = list(model = fix_model, mut_rate = fix_mut_rate,  n_eff = fix_n_eff)
  )
  param
}


beta_boltzmann <- function(R = 1.986e-3, T = 298) 1/(R*T)

