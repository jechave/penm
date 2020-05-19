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

## 4 May 2020

### Changed order of `eval` and `umat` columns in `enm_nma()`

### Changed all enm plotting functions

### Added several functions to *enm\_analysis.R*

### Changed *test\_enm.R*

Now it calculates prot then calls plot functions in the order they are
in the *plot\_enm.R* file.

## 5 May 2020

### Learned useful tools

  - Learnt to use `knitr::purl()` to translate .Rmd into .R
  - Learnt to use `mvbutils::foodweb()` to plot dependencies of
    functions

### Refactored prot getters, calculators, and plotters

  - Changed function name: `cmat(prot)` to `get_reduced_cmat(prot)`
  - Changed function name: `kmat(prot)` to `get_reduced_cmat(prot)`
  - Changed function name: `umat2(prot)` to `get_umat2(prot)`
  - Changed function name: `umat2_matrix(prot)` to
    `get_umat2_matrix(prot)`
  - Changed definition of all matrix plot functions so that they do not
    need to transform to tibble before plotting (they call `plot_matrix`
    instead)
  - Deleted `msf_site_mode`
  - Renamed `msf_site_mode_matrix` to `get_msf_site_mode`
  - Deleted `get_umat2`
  - Renamed `get_umat2_matrix` to `get_umat2`
  - Renamed `rho_matrix`to `get_rho_matrix` and `plot_rho` to
    `plot_rho_matrix`
  - Moved general functions `matrix_to_tibble`and `plot_matrix` to
    package `jefuns`
  - Moved all plot functions of enm into package directory as
    `enm_plot.R`
  - Tested and comitted to git and github

### Other

  - Eliminated `d_max` from `get_cn` and related

## 6 May 2020

### Restructured prot object

  - Added enm parameters to `prot$enm`
  - Wrote getters for prot object
  - Changed `set_prot` to a single `set_enm(pdb, ...)` function that
    sets up the prot object
  - Changed structure of `prot` object (and it’s class is `prot`)
  - Fixed getters according to new structure of prot
  - Changed queries by getters in enm module (but not in penm.R and
    activation.R): e.g. replaced prot\(enm\)umat by get\_umat(prot)

## 7 May 2020

### Rename functions called by `set_enm`

  - Extract function `set_enm_param`
  - `set_nodes` to `set_enm_nodes`
  - `enm_graph_xyz` to `set_enm_graph`
  - `eij_edge` to `set_enm_eij`
  - `kmat_graph` to `set_enm_kmat`
  - `enm_nma` to `set_enm_nma`

### Fix calls to prot object in penm.R and related

The prot object of the enm module was restructured. Therefore, I need to
change everywhere where prot objects are called in penm.R and related.

  - Did all the necessary renaming in penm.R and penm\_analysis.R
  - Checked that get\_penm\_mutant works
  - Wrote automatic tests test\_enm.R test\_penm.R for further
    refactoring
  - Tested automatically set\_enm() and get\_mutant\_site(): they work.
  - Commited to git and
github

### **WARNING: v\_min changes between update\_enm = F and update\_enm = T**

### Refactor penm.R

Get enm parameters from prot object rather than pass it independently in
all penm.R funcitons

  - Eliminate model, d\_max, and frustrated from arguments of
    get\_mutant\_site and functions called from there

## 8 May 2020

### Refactor set\_enm

I refactored set\_enm by adding functions that depend as much as
possible only on prot objects, so that parameters are passed through
prot.

  - new set\_enm\_ functions depend mostly on `prot` rather than
    explicitly on its components
  - added `set_enm_nma()` to file *enm.R* and deleted *enm\_nma.R*
  - old `set_enm_` are now `calculate_enm_`

### Refactor `get_mutant_site()` in *penm.R*

  - replaced all `calculate_enm_` functions by `set_enm_(prot)`
    functions in *penm.R*
  - Changed `get_force`
  - Changed dlij from lij(mut) - dij(wt) to lij(mut) - lij(wt)
  - Made `get_mutant_site` a bit shorter by adding get\_dlij
  - Removed `wt0` (I wasn’t using it, just confusing).
  - Tidied up penm.R file a little bit more

## 9 May 2020

### Put update\_enm on stand by

  - Changed `update_enm` to `mut_model`, that can now be `lfenm` (K
    doesn’t change) or `sclfenm` (the update\_enm = T version previous).
  - Wrote a note regarding my worries about the sclfenm version
  - Separated more clearly the options “lfenm” and “sclfenm” in *penm.R*
    functions
  - Made current “sclfenm” option stop if called because I need to
    revise it.
  - Made “sclefnm” tests skip the test.
  - Tested
  - Commited

## 12 May 2020

  - Run and revisions of test\_mutate\_structure.Rmd
  - created new perturbation\_response\_scanning.Rmd that produces
    slides
  - commited
everything

### Separated data creation and analysis of perturbation response scanning

  - Put all creation in “Rmd/prs\_data\_create.Rmd”
  - Put all analysis in “Rmd/prs\_data\_analyse.Rmd”
  - Kept in these files only “prs”, not pair comparison mut vs. wt
  - Created s pair-comparison file “mut\_vs\_wt.Rmd”, to be completed
    later

## 13 May 2020

### Complete penm\_data\_analyse.Rmd

  - Created some new notes on PRS
  - Completed first full version of penm\_data\_analyse.Rmd
  - Deleted TODO.Rmd (Using Trello and Things… already too much)
  - Tested and commited

### Clean up

  - moved `prs` functions into package: file *R/prs.R*
  - created *./to\_do* directory to hold planned features
  - moved relevant files related to planned features to */to\_do*
  - run all Rmd reports
  - run all tests
  - Commit

## 14 May 2020

### Add energy response to prs module

  - Added theory doc to *docs*, and notes on energy differences
  - Refactored *prs.R*
  - Added calculation of data.frame dfej in *prs.R*
  - Added calculation of dfej to *prs\_create\_data.Rmd*
  - Added analysis of dfej to *prs\_analyse\_data.Rmd*
  - Tested
  - Commited

### Add globality of respnse and influence to prs\_analyse\_data.Rmd

  - Added a few slides analysing response and influence globaility in
    site and mode spaces.
  - Test
  - Commit

## Add v\_stress to prs

  - Developed some more theory, regarding energies and separating them
    into site contributions
  - Added `enm_v_stress()` and `enm_delta_v_stress` to `penm` module
  - Calculated delta\_v\_stress in *prs.R* (and removed `delta_u` and
    `delta_a`)
  - Recalculated response data using *prs\_create\_data.Rmd*
  - Revised and rerun *prs\_analyse\_data.Rmd*
  - Tested
  - Commited

## 18 May 2020

### Superfast calculation of response matrices

  - Wrote the key formulae in *docs/notes*
  - Wrote *prs\_fast.R* that contains functions for fast calcualation of
    site responses: dr2ij, de2ij, and df2ij
  - Wrote *prs\_fast.Rmd* that compares cpu times and results between
    “fast” and “slow” methods
  - Tested
  - Commited

## 19 May 2020

### Better cpu-time comparison of prs\_fast vs. prs

  - Improved prs by making the de2\_site require kmat\_sqrt as input
    that is calculated outside in the calling function only once for the
    whole scan.
  - Made a more relevant cpu-time comparison
  - Verified that simulated prs responses converge to analytical values
    as number of mutations per site increases
  - Tested
  - Commited

### Added \`fast\_delta\_structure\_mode() to *prs\_fast.R*

  - Added fast\_delta\_structure\_mode and needed functions called by
    it.
  - Made `delta_structure_mode()` in *prs.R* faster
  - Tested cpu-time and convergence in *prs\_superfast\_mode.Rmd*
  - Commited
