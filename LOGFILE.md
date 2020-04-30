penm development log
================

## 30 April 2020

### Eliminate need to use `ideal` in enm and lfenm calculations

`ideal` is the ideal protein structure, needed to calculate `v_stress`,
`dv_activation`, and `g_entropy_activation`. To eliminate it, I can just
avoid calling these functions to set up and mutate `wt`

1.  Moved `v_stress`, `dv_activation`, and `g_entropy_activation` to
    *activation.R*
2.  Deleted calls to \`enm\_energy everywhere (thus current protein
    setup does not add energies to prot object).

### Eliminate active-site dependent info from enm and lfenm calculations

pdb\_active\_site indices are used to calculate cmat\_activation,
kmat\_activation, and activation energies. I do not need this in the
prot object. If needed, I can add it later.
