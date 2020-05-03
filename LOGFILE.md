penm development log
================

## 30 April 2020

### Eliminate need to use `ideal` in enm and lfenm calculations

`ideal` is the ideal protein structure, needed to calculate `v_stress`,
`dv_activation`, and `g_entropy_activation`. To eliminate it, I can just
avoid calling these functions to set up and mutate `wt`

1.  Moved `v_stress`, `dv_activation`, and `g_entropy_activation` to
    *activation.R*
2.  Split old enm\_energy into two functions: enm\_energy\_activation
    (it calls dv\_activation and g\_entropy\_activation) and enm\_energy
    that calculates just v\_min and g\_entropy

### Eliminate active-site dependent info from enm and lfenm calculations

pdb\_active\_site indices are used to calculate cmat\_activation,
kmat\_activation, and activation energies. I do not need this in the
prot object. If needed, I can add it later.

1.  I moved everything related to active site into file activation.R
    (need to test it, current tests don’t use it).
2.  Everything else is independent of either “pdb\_site\_active” or
    “ideal”. Therefore, it should run for tasks that do not depend on
    defining active sites.
3.  Modified tests accordingly (eliminating dummy active site and
    ideal).
4.  Commited changes to git and github

### Make a one-button setup of prot object and test

right now setting up the protein object is done by first reading a pdb
file into a bio3d pdb object, then calling prot\_sc or prot\_ca to
initialize the prot object, then calling init\_prot to cmplete it. Join
prot\_sc/prot\_ca with init\_prot into a single function `set_prot(pdb)`
Consider adding the enm parameters to prot to be used when passing to
other functions(prot)

  - Replaced previous “pdb\_sc” followed by “init\_prot” by a single
    set\_prot(pdb,…) function that sets up the enm, performs nma etc.
  - Added enm\_param to prot object
  - Set up a prot\_test.R test in tests/testthat directory, which is run
    when devtools::test()
  - Added objects pdb\_2acy\_A and prot\_2acy\_A to data folder for
    testing Tested and commited all changes.

## 1 May 2020

Refactor enm

### beta and energy terms

beta = 1/(k\_boltzman T), depends on temperature. Therefore it is not a
property of the protein but, rather, something that depends on the
protein’s properties and the environment’s temperature. For this reason,
it would be better to calculate energies when needed, rather than
attaching energy terms to the prot object.

Thus, eliminate energy from protein object (and change everything
accordingly)

  - deleted enm\_energy and enm\_energy\_activation functions
  - deleted all calls to enm\_energy and recursively…
  - deleted query functions get\_v\_min, get\_g\_entropy,
    get\_v\_min\_activaiton, get\_g\_entropy\_activation, get\_v\_stress
  - tested
  - commited to git and github

### make plot\_enm functions good enough to move to package

  - added some more plots to test\_enm.Rmd
  - fixed an issue: graph setting (in set\_enm\_xyz) missed some i-(i+1)
    contacts for which dij \> d\_max (not a problem for CA models, but
    it’s a problem for SC models).
  - Changed kij\_anm and kij\_ming\_wall so that they don’t set these
    i-(i+1) kij to 0.
  - Tested and commited.

## 2 May 2020

### Deleted `add_site_indexes()`

  - added `nsites` and `site` to the result returned by `prot_sc()` and
    `prot_ca()`
  - tested and commited

### Deleted `add_enm()`

  - renamed `enm_set_xyz` to `enm_from_xyz`
  - made `enm_from_xyz` call `nma(kmat)` and return also `mode, evalue,
    umat, cmat`
  - made set\_prot call `enm_from_xyz` directly, rather than `add_enm`
  - deleted `add_enm`
  - tested and commited

### Moved non-binary data to `./data_raw`

## 3 May 2020

### merged set\_prot.R and add\_prot.R into single file

### renamed various files in package R directory
