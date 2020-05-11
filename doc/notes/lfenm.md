# LFENM: linearly forced elastic network model

$$
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
$$

The LFENM keeps only the linear terms of the perturbed potential, discarding the quadratic terms:
$$
V_\mut(\rvec) = V_\mut(\rmut)+ \half (\rvec - \rmut)^T\kwt(\rvec - \rmut)
$$
where the mutant's minimum energy is given by:
$$
V_\mut(\rmut) = \frac{1}{2}\sum \kij \dlij^2 - \half (\rmut - \rwt)^T \kwt(\rmut - \rwt).
$$
Note that $\mathbf{K}$ for the mutant is identical to that of the wild-type (because quadratic terms of the potantial are discarded.) There's a *lateral* and *vertical* shift of the potential energy surface, but no *rotation* (same normal modes) or *deformation* (same eigenvalues).

