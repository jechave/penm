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
    prot <- set_enm(pdb, node = "calpha", model = "anm", d_max = 10.5, frustrated = FALSE)

In this example, network nodes are placed at \(C_\alpha\) coordinates,
the model used is Bahar’s Anisotropic Network Model (“anm”) with a
cut-off distance to define contacts of `d_max = 10.5`. `frustrated`
indicates whether to add frustrations to the model.

## Package architecture

### enm dependencies

![Dependencies of `set_enm()`.](images/tree_set_enm.png)

### penm dependencies

![Dependencies of `get_mutant_site()`](images/tree_get_mutant_site.png)
