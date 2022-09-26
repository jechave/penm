penm package
================
Julian Echave
26/09/2022

<!-- README.md is generated from README.Rmd. Please edit that file -->

The `penm` package contains functions to build Elastic Network Models
(ENM) of proteins and to perturb them; `penm` stands for Perturbing
Elastic Network Models.

## Overview

In short:

-   `bio3d::read.pdb()`: Set up a **pdb** protein object.
-   `set_enm()`: Set up a `penm` **prot** object, containing protein and
    its ENM info.
-   `amrs()`and `smrs()`: perform single-mutation scans to calculate
    mutation-response matrices.
-   `admrs()`and `sdmrs()`: perform double-mutation scans to calculate
    compensation matrices.

## Installation

## Getting started

### Set up the ENM for a protein

First, read a pdb file using `bio3d::read.pdb` to generate a `pdb`
object for a protein. Then, create the `prot` object, that contains the
full ENM analysis.

    bio3d::read.pdb("data-raw/2acy_A.pdb") # read a pdb file
    prot <- set_enm(pdb, node = "calpha", model = "anm", d_max = 10.5, frustrated = FALSE)

In this example, network nodes are placed at
![C\_\alpha](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C_%5Calpha "C_\alpha")
coordinates, the model used is Bahar’s Anisotropic Network Model (“anm”)
with a cut-off distance to define contacts of `d_max = 10.5`.
`frustrated` indicates whether to add frustrations to the model.
