# Wild-type and mutant energies

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

## Wild-type energy

$$
V_\wt(\rvec) = V_\wt(\rvec_\wt)+ \half (\rvec - \rwt)^T\kwt(\rvec - \wt)\\
V_\wt(\rvec_\wt) = 0
$$

For this potential, the thermodynamic energies (excluding the kinetic energy term that doesn't change if there're no indels) are:
$$
U = <V> = V_\wt(\rvec_\wt)+ (3N -6)\frac{1}{2\beta},\\
A = U - TS = V_\wt(\rvec_\wt) +\frac{1}{2\beta} \sum_n{\frac{\beta \lambda_n}{2 \pi}}\\
TS = \frac{1}{2\beta} \sum_n({\frac{2 \pi}{\beta \lambda_n} + 1})
$$



## Mutant energy

Let:

$$
V_\mut(\rvec) = V_\mut(\rmut)+ \half (\rvec - \rmut)^T\kmut(\rvec - \rmut)
$$
with,
$$
V_\mut(\rmut) = \frac{1}{2}\sum \kij^\mut \dlij^2 - \half (\rmut - \rwt)^T \kmut(\rmut - \rwt).
$$

Therefore, the thermodynamic energies are:
$$
U = <V> = V_\mut(\rvec_\wt)+ (3N -6)\frac{1}{2\beta},\\
A = U - TS = V_\mut(\rvec_\wt) +\frac{1}{2\beta} \sum_n{\frac{\beta \lambda^\mut_n}{2 \pi}}\\
TS = \frac{1}{2\beta} \sum_n({\frac{2 \pi}{\beta \lambda^\mut_n} + 1})
$$


## Two ways of calculating mutant's minimum

$$
V_\mut(\rmut) = \frac{1}{2}\sum \kij^\mut \dlij^2 - \half (\rmut - \rwt)^T \kmut(\rmut - \rwt)
$$

$$
V_\mut(\rmut) = \half\sum_{\texttt{edge}}k^\mut_{ij}[\dij^\mut- (\lij + \dlij))]^2
$$

