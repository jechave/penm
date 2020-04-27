delta_v_min <- function(prot1, prot2)
  get_v_min(prot2) - get_v_min(prot1)

delta_g_entropy <- function(prot1, prot2)
  get_g_entropy(prot2) - get_g_entropy(prot1)

delta_v_stress <- function(prot1, prot2)
  get_v_stress(prot2) - get_v_stress(prot1)

delta_v_activation <- function(prot1, prot2)
  get_dv_activation(prot2) - get_dv_activation(prot1)

delta_g_entropy_activation <- function(prot1, prot2)
  get_g_entropy_activation(prot2) - get_g_entropy_activation(prot1)
