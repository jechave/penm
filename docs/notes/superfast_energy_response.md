# Superfast energy response

$
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
\newcommand{\vstress}{V^\texttt{stress}}
\newcommand{\dvstress}{\delta \vstress}
$

####  $\delta V_\texttt{stress}$

From [prs theory](prs_theory.md), we know that when the mutant adopts the wild-type conformation, it's energy is given by the "stress" energy:
$$
V_\mut(\wt) = V_0 + \half \sum_{i \sim j} V^\mut_{0, ij} + \half \sum \kij^\mut \dlij^2
$$
In PRS, $V_0$ is the minimum energy of the wild-type. Assuming the mutatonal scan does not affect $V^0_{ij}$ (as I have assumed in the stress-model and activity-stability models), then
$$
\delta V_\texttt{stress} =\half \sum \kij^\mut \dlij^2
$$
We can split this into site-stress contributions:
$$
\delta V_\texttt{stress} =\half \sum_i (\half \sum_{j\sim i} \kij^\mut \dlij^2)
$$
where the second $\half$ is introduced to avoid counting twice each edge. We can rewrite:
$$
\dvstress = \half \sum_i \dvstress_i \\
\dvstress_i = \half \sum_{j \sim i} k_{ij}\dlij^2
$$

When mutating a single edge $jk$, 
$$
\dvstress_{i,jk} = \half (\delta_{ij}+\delta_{ik})k_{jk}\delta l_{jk}^2
$$
Averaging over mutations:
$$
<\dvstress_{i,jk}> = \half (\delta_{ij}+\delta_{ik})k_{jk}\sigma_\mut^2
$$
When mutating a site $j$, 
$$
<\dvstress_{i,j}> = \sum_{k\sim j}\half(\delta_{ij} + \delta_{ik})k_{jk}\sigma^2
$$
Therefore:
$$
\begin{eqnarray}
<\dvstress_{i,i}> &=& \sum_{j \sim i} \half \kij \sigma^2\\
<\dvstress_{i,j}>&=& \half \kij \sigma^2 \texttt{for } j \sim i \\
<\dvstress_{i, j}> &=& 0 \texttt{ otherwise }
\end{eqnarray}
$$
These formulas allow the calculation of a "superfast" stress-energy response matrix. The stress energy response to a mutation in $j$ is half the sum of column $j$ of this matrix. The mean over mutations is the profile of site-dependent responses. 
$$
<\dvstress_{.,j}> = \half\sum_{i}<\dvstress_{i,j}> =\sum_{j \sim i} \half \kij \sigma^2
$$
Note that the matrix is symmetric: the response of $i$ to mutating $j$ is identical to the response of $j$ to mutating $i$: just the stretching af a single spring, or zero.





#$\delta V_{min}$

The mutant's minimum energy can be calculated from:
$$
\delta V_\texttt{min} = \delta V_\texttt{stress}-\half\sum\delta e_i^2
$$
where $\delta \evec$ is the energy response vector defined elsewhere.

Therefore the $\delta V_\texttt{min}$ site-by-site response matrix is:
$$
<\delta V^\texttt{min}_{i,j}> = <\dvstress_{i, j}>-\delta e^2_{i,j}
$$


