# Superfast energy response

$
\newcommand{\rvec}{\mathbf{r}}
\newcommand{\evec}{\mathbf{e}}
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
\newcommand{\vmin}{V^\texttt{min}}
\newcommand{\dvmin}{\delta \vmin}
$

##  $\delta V_\texttt{stress}$

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

## $\delta V_{min}$

Rigorously:
$$
\delta V_\texttt{min} = \half\sum\kij^\mut (\dij^\mut - \lij^\mut)^2-\half\sum\kij^\wt (\dij^\wt - \lij^\wt)^2
$$
Assuming that $\lij^\wt = \dij^\wt$ (i.e. the wild-type's springs are all perfectly relaxed at equilibrium):
$$
\delta V_\texttt{min} = \half\sum\kij^\mut (\dij^\mut - \lij^\mut)^2
$$
To use this formula, we we need to calculate the mutant's structure. It is not clear how to average over mutants to obtain an analytical foormula, at this point. However, we can approximate:
$$
\delta V_\texttt{min} = \delta V_\texttt{stress}-\half\sum\delta e_i^2
$$
where $\delta \evec$ is the energy response vector defined elsewhere.

Therefore the $\delta V_\texttt{min}$ site-by-site response matrix is:
$$
\label{dvmin}
<\delta V^\texttt{min}_{i,j}> = <\dvstress_{i, j}>-\delta e^2_{i,j}
$$

## Warning!

 #warning, #todo

 As I defined it, $\delta e^2_i >= 0$. However, $\dvstress_{i} - \dvmin_i$ may be negative. Physically, the local stress at one site may increase to decrease the local stress at other sites when the protein relaxes to minimize the total energy. 

==For consistency, for the fast calculation, I'm using Eq.$\ref{dvmin}$ too, rather than the "exact" calculation of the mutant's $v_{min}$, which I can do numerically but, so far, not analytically (except via. 15).==

## Partial problems to solve

### Where's the minimum? 

#todo

Where's the minimum of a general potential:
$$
V = 1/2 \sum_{ij}\kij (\dij - \lij)^2
$$
when there're no restrictions on $\lij$ (springs not necessarily relaxed at the minimum) and $\dij =\sqrt{||\rvec_j - \rvec_i||^2}$? Note that by definition $\dij$ are not independent. 



### Fast calculation of minimum energy perturbation?

Solving the previous issue, I could solve the issue of finding an analytical expression for average of minimum-energy perturbations over mutations, perhaps. 

### Negative deformation energies?

#todo

When the network deforms, if all springs are completely relaxed at equilibrium, any deformation will make all springs, and thus all nodes more stressed.



However, if the there's frustration, i.e. if at equilibrium some springs are stressed and others relaxed, relieving the stress of some springs may increase that of others. Therefore, even if the total deformation energy is positive (moving away from the minimum increases energy), the contributions may be positive or negative.



The approximation $\half (\rvec - \rvec_0)^T \kmat (\rvec - \rvec_0)$ doesn't seem to account for increases and decreases of stress... 



The problem is, then:
$$
\color{red}{
\half\sum\kij (\dij - \lij)^2-\half\sum\kij (\dij^0 - \lij)^2 \approx \half (\rvec - \rvec_0)^T \kmat (\rvec - \rvec_0)?}
$$








