# penm 0.5.1

## Bug fixes

* `delta_energy_dvs()` gave wrong values when `ideal` has springs that the
  protein lacks, for example a genm `ideal` against an lfenm protein: it took
  `ideal`'s spring lengths by position instead of by spring, with only a
  warning. Proteins with the same springs, such as a wild type and its mutants,
  give the same values as before.

* A coordinate vector whose length is not a multiple of 3 stops with an error
  that says so, instead of printing the explanation and stopping with "Stopping
  execution".

## Documentation

* Formulas in the help pages have plain-text forms, so text help shows them
  readably instead of as LaTeX.
* The `delta_*` help pages say which arguments need normal modes, and their
  examples show the motion and entropy measures under `mut_model = "genm"`.
* `get_enm_param()` is documented as returning the network's and the
  mutational model's parameters.

## Other changes

* `R CMD check` no longer notes undefined global variables.
* Continuous integration uses `actions/checkout` v7.

# penm 0.5.0

## Breaking changes

* `set_enm()`'s `n_alleles` is renamed `alphabet_size`, and so is the entry in
  a protein's `param`: it is the size of the model's sequence alphabet, the
  number of letters a site can carry. Results do not change. A `prot` built by
  0.4.0 is rejected by `get_mutant_site()`: rebuild it with `set_enm()`.

## Documentation

* `?set_enm` has a Details section on sequences and mutations: what a letter
  does to the spring lengths, and that letter 0, the pdb's residue, is the
  reference. A mutation from it changes a spring's length with standard
  deviation `mut_dl_sigma`; a mutation between two other letters, with standard
  deviation `mut_dl_sigma` times the square root of 2.
* `mut_dl_sigma` is defined correctly: the spread of a letter's shift relative
  to the pdb's length, not of any mutation's change.
* `?set_enm` describes the fields of the returned `prot`; `?get_prot_property`
  says that the mode accessors need normal modes, which a genm mutant gets from
  `set_enm_nma()`.
* User-facing text calls a site's value a letter, not an allele.
* The help pages render their markdown. Roxygen's markdown support was never
  turned on, so code, links and emphasis showed as literal backticks, square
  brackets and asterisks in every page written with them.

# penm 0.4.0

## Breaking changes

* Mutants differ from those of 0.3.0. A mutation now gives a site an *allele*,
  from 0 (the pdb's residue) to `n_alleles - 1`, and the equilibrium lengths are
  a function of the alleles: mutating a site back to an earlier allele restores
  its earlier lengths, and a protein's structure no longer depends on the path
  of mutations that led to it. The random draws are new, so the same
  `(ensemble, site_mut, mutation)` gives a different mutant than in 0.3.0;
  results computed with 0.3.0 will not be reproduced exactly. Statistically the
  two agree: the site and mode profiles of a complete single-point scan agree
  as closely as two ensembles of the same version do.

* How a protein mutates is now set when it is built, not when it is mutated.
  `get_mutant_site()` loses its `mut_model`, `mut_dl_sigma`, `mut_sd_min` and
  `ensemble` arguments; `set_enm()` takes them instead, stores them in the
  protein, and mutants inherit them. Code passing them to `get_mutant_site()`
  must pass them to `set_enm()`:

  ```r
  # 0.3.0
  wt <- set_enm(pdb, "ca", "ming_wall", 10.5, frustrated = FALSE)
  mut <- get_mutant_site(wt, 11, 1, mut_sd_min = 1, ensemble = 7)
  # 0.4.0
  wt <- set_enm(pdb, "ca", "ming_wall", 10.5, mut_sd_min = 1, ensemble = 7)
  mut <- get_mutant_site(wt, 11, 1)
  ```

* `mutation` must be an allele, from 0 to `n_alleles - 1` (default 10); it was
  unbounded. The allele a site already has returns the protein unchanged, as
  `mutation = 0` did.

* A `prot` built by an earlier version is rejected by `get_mutant_site()`:
  rebuild it with `set_enm()`. Its layout changed: `eij` moved under
  `internal`, which is not meant to be read; `param` holds all of `set_enm()`'s
  arguments but `pdb`; `nodes` has the `sequence` of alleles; `graph` has the
  pdb lengths `l0ij`. Use the accessors (`?get_prot_property`).

* Networks built with `model = "reach"` change: pairs one, two and three apart
  in sequence now get the model's fixed force constants. Before, the sequence
  separation never reached the force-constant function, and every pair got the
  distance-dependent constant.

* `set_enm()` no longer takes `frustrated`.

* `get_vmin_site()` and its deprecated alias `get_stress()` are removed.

## New features

* A second mutational model, a generalized ENM: `set_enm(..., mut_model =
  "genm")`. The force constants follow the equilibrium lengths, and a mutant's
  structure is the minimum of its energy, so mutations change the normal modes.
  `d_max_graph` sets which pairs of residues get a spring at all.
  `vignette("genm")` uses it for a complete mutational scan and a star tree of
  evolutionary trajectories.

* `set_enm_nma()` is exported: a genm mutant comes back without normal modes,
  which a long trajectory need not compute at every step, and
  `set_enm_nma()` adds them.

* `superpose_prot()` superposes a protein onto a structure, rotating its network
  matrix, normal modes and everything else that depends on orientation with it.

* New models `"anm_smooth"` and `"ming_wall_smooth"` replace the step at
  `d_max` by a smooth switch of width `d_max_width`. `set_enm()` passes
  further named arguments (`...`), such as `d_max_width`, to the
  spring-constant function.

## Other changes

* The mode accessors (`get_evalue()`, `get_umat()`, `get_cmat()`, `get_mode()`,
  `get_nmodes()`) stop with "prot has no normal modes: add them with
  set_enm_nma()" when a protein has none.

* `set_enm()` stops when the network is not a rigid minimum: when the matrix
  has negative eigenvalues, or more or fewer than six null ones.

* An lfenm mutant's coordinates are a vector, as documented, not a one-column
  matrix.

# penm 0.3.0

## Breaking changes

* `delta_structure_dvmi()`, `delta_structure_dvsi()` and
  `delta_structure_dvsi_same_topology()` are removed. All three compare the two
  graphs edge by edge, and none is correct when the contact map changes: `dvmi`
  drops unshared edges, `dvsi` errors, and `dvsi_same_topology` matched rows by
  position and could return wrong numbers silently.

## Other changes

* `get_stress()` is renamed `get_vmin_site()`; `get_stress()` remains as a
  deprecated alias. It is the per-site decomposition of `enm_v_min()`, not the
  stress energy of `delta_energy_dvs()`.
* `get_vmin_site()` always returns a vector of length `nsites`. Previously a
  site with no springs was dropped, shortening and misaligning the result.
* `ddg_tds()` now errors if `wt` and `mut` have a different number of modes,
  instead of subtracting sums over spectra of different length.

# penm 0.2.0

## Breaking changes

* `get_mutant_site()`'s `seed` argument is now `ensemble`, and its default
  changed from `241956` to `1L`. Code passing `seed = ` must be updated.

  The old name was wrong: the value never reached `set.seed()`. It is hashed
  together with `(site_mut, mutation)`, and *that* hash seeds the RNG. Because
  `mutation` is unbounded and the model has no amino acids, the argument names
  **which realization of the mutational process** a mutant belongs to. See
  `?penm_ensemble`.

* **Mutants generated by this version differ from those generated by 0.1.0**,
  even for the same site and mutation index. The internal key dropped a
  redundant component, which changes every hash and so every perturbation drawn.
  Stored results and cached fixtures from 0.1.0 are not comparable with 0.2.0
  output. The *set* of contacts a mutation perturbs is unchanged — that depends
  on the site and `mut_sd_min`, never on the key — only the magnitudes moved.

* A malformed `ensemble` is now an error. Previously `NULL`, `NA`, a string, or
  a vector each produced a valid-looking mutant under a realization nobody
  chose, because the key was built with `paste()`, which stringifies anything.

## Improvements

* `get_mutant_site()` no longer overwrites the caller's `.Random.seed`. Code
  that seeds its own analysis and then generates mutants keeps its RNG stream.
  This changes no value the function itself returns.

* New `?penm_ensemble` help topic explaining what `ensemble` means, why it must
  be held fixed within a scan or an evolutionary trajectory, and when changing
  it is the right thing to do.

# penm 0.1.0

* Initial version, extracted from the package now called `penmscan`, which
  retains the site-by-site scanning layer (the `mrs_*` functions).
