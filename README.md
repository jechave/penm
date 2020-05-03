penm package
================
Julian Echave
5/2/2020

The `penm` package contains functions to build Elastic Network Models
(ENM) of proteins and to perturb them. Thus, `penm` stands for
Perturbing Elastic Network Models.

## Usage

### Set up the ENM for a protein

First, read a pdb file using `bio3d::read.pdb` to generate a `pdb`
object for a protein. Then, create the `prot` object, that contains the
full ENM analysis.

    bio3d::read.pdb("data-raw/2acy_A.pdb") # read a pdb file
    prot <- set_prot(pdb, node = "calpha", model = "anm", d_max = 10.5, frustrated = FALSE)

In this example, network nodes are placed at \(C_\alpha\) coordinates,
the model used is Bahar’s Anisotropic Network Model (“anm”) with a
cut-off distance to define contacts of `d_max = 10.5`. `frustrated`
indicates whether to add frustrations to the model.

## Package architecture

### enm dependencies

  - `set_prot(pdb, node, model, d_max, frustrated)`
      - `nodes(pdb, node)`
          - `prot_ca(pdb)`
          - `prot_sc(pdb)`
              - `residue.coordinates(pdb)`
              - `residue.bfactors(pdb)`
      - `enm_from_xyz(xyz, pdb_site, model, d_max, frustrated)`
          - `enm_graph_xyz(xyz, pdb_site, model, d_max)`
              - `my_as_xyz(xyz)`
              - `kij(dij, d_max)`
              - `dij_edge(xyz, i, j))`
              - `sdij = sdij_edge(pdb_site, i, j))`
          - `eij <- eij_edge(xyz, i, j)`
              - `my_as_xyz(xyz)`
          - `kmat_graph(graph, eij, nsites, frustrated)`

### penm dependencies

  - `get_mutant_site(wt, site_mut, mutation, seed, wt0, mut_sd_min,
    dl_sigma, update_enm, model, d_max, frustrated)`
      - mutates the graph
      - `get_force(wt, mut)`
      - calculates `dxyz`
      - `if (update_enm): * mutate_enm(mut, model, d_max, frustrated)`
          - `mutate_graph(prot, model, d_max)`
              - `enm_graph_xyz(xyz, pdb_site, model, d_max)`
          - `mutate_eij(prot)`
              - `eij_edge(xyz,i, j)`
          - `kmat_graph(graph, eij, nsites, frustrated)`
          - `enm_nma(kmat)` \# recalculate normal modes etc. nma \<-
            enm\_nma(prot\(enm\)kmat) \#returns mode, evalue, cmat, umat
      - `if (!update_enm):`
          - `dij_edge(xyz, i, j)`
