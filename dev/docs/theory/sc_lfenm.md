---
output:
  pdf_document: default
  html_document: default
---
# Self Consistent Linearly Forced Elastic Network Model

\newcommand{\rvec}{\mathbf{r}}
\newcommand{\wt}{\texttt{wt}}
\newcommand{\mut}{\texttt{mut}}
\newcommand{\rwt}{\rvec^0_\wt}
\newcommand{\rmut}{\rvec^0_\mut}
\newcommand{\kmat}{\mathbf{K}}
\newcommand{\kwt}{\kmat_\wt}
\newcommand{\kmut}{\kmat_\mut}
\newcommand{\kij}{k_{ij}}
\newcommand{\fvec}{\mathbf{f}}
\newcommand{\dij}{d_{ij}}
\newcommand{\lij}{l_{ij}}
\newcommand{\dlij}{\delta l_{ij}}
\newcommand{\half}{\frac{1}{2}}

The [LFENM](lfenm.md) model is defined by:
$$
V_\mut(\rvec) = V_\mut(\rmut)+ \half (\rvec - \rmut)^T\kwt(\rvec - \rmut) \label{vlfenm}
$$
where
$$
V_\mut(\rmut) = \frac{1}{2}\sum \kij \dlij^2 - \half (\rmut - \rwt)^T \kwt(\rmut - \rwt).
$$
The second term of Eq.$\ref{vlfenm}$ is $0$ at the minimum, because any possible frustration has been included in the minimum-energy first term.

A major weakness of the LFENM model is that: 
$$
\kmat_\wt \neq \kmat(\rmut)
$$

Since ENMs are defined in such a way that $\mathbf{K}$ is derived from the minimum-energy conformation $\mathbf{r}^0$, this is inconsistent with the basic assumption of ENM models. 

The previous issue is overcome in the SC-LFENM model that first calculates the mutant's structure using:


$$
\rmut = \rwt - \kwt^{-1}\fvec
$$
Then, we recalculate $\kmat$:
$$
\kmut = \kmat(\rmut)
$$
Thus the full SC-LFENM model is specified by:
$$
V_\mut(\rvec) = V_\mut(\rmut)+ \half (\rvec - \rmut)^T\kmut(\rvec - \rmut) \label{sclfenm}
$$
whith,
$$
V_\mut(\rmut) = \frac{1}{2}\sum \kij \dlij^2 - \half (\rmut - \rwt)^T \kmut(\rmut - \rwt).
$$
Here, I assume that everything is done in the potential energy surface of the mutant: i.e. firs term sums "stress" energies only over the eges of the mutant's network (new edges have $\dlij = 0$, thus they don't contribute, delted edges don't contribute befause their $\kij^\mut = 0$); the second term is relaxation of the mutant's from $\rwt$ to $\rmut$. The second term of Eq.$\ref{sclfenm}$ is the energy of deformation around the minimum $\rmut$. 

We can turn the wild-type into the mutant by a vertical transition at $\rwt$ in which we replace $\kij^\wt$ by $\kij^\mut$ followed by the relaxation of the mutant from $\rwt$ to $\rmut$. 



## Problem: the model is not reversible

If I start with the mutant and do the reverse mutation, I don't get the wild-type's potential... 

#todo
