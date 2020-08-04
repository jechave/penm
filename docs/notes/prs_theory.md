$
\newcommand{\rvec}{\mathbf{r}}
\newcommand{\fmat}{\mathbf{F}}
\newcommand{\svec}{\mathbf{s}}
\newcommand{\smat}{\mathbf{S}}
\newcommand{\wt}{\texttt{wt}}
\newcommand{\mut}{\texttt{mut}}
\newcommand{\rwt}{\rvec^0_\wt}
\newcommand{\rmut}{\rvec^0_\mut}
\newcommand{\kmat}{\mathbf{K}}
\newcommand{\kwt}{\kmat_\wt}
\newcommand{\kmut}{\kmat_\mut}
\newcommand{\amat}{\mathbf{A}}
\newcommand{\kij}{k_{ij}}
\newcommand{\fvec}{\mathbf{f}}
\newcommand{\dij}{d_{ij}}
\newcommand{\lij}{l_{ij}}
\newcommand{\dlij}{\delta l_{ij}}
\newcommand{\half}{\frac{1}{2}}
$

# Structural response theory

Definitions:

* $f_{ij} = k_{ij} \delta l_{ij}$ is the force in the direction of contact $i-j$. $\mathbf{df}$ is the force-fector. 
* $d\mathbf{r} = \mathbf{C df}$ is the change in minimum-energy conformation.
* $\mathbf{de} = K^{1/2}dr = C^{-1/2}df$ is a vector whose norm is the deformation energy: $||\mathbf{de}||^2 = d\mathbf{r}^T \mathbf{K} d\mathbf{r}$.

From the previous definitions, it follwos that in the normal-mode respresentation: 
$$
de_{n} = \sigma_n df_{n} \\ dr_n = \sigma_n^2 df_n
$$
By analogy, in a site representation, it is approximately true that (#check): 
$$
de^2_{i} = \sigma_i^2 df^2_{i} \\ dr^2_i = \sigma_i^4 df^2_i
$$
Also, we have shown previously (Echave & Fernández, 2010)  that, on average: 
$$
df^2_n \propto \frac{1}{\sigma_n^2}
$$
By analogy, it might be that (#check): 
$$
df^2_i \propto 1 / \sigma_i^2
$$

In summary, we expect that, regardless of representation:
$$
df^2 \propto \frac{1}{\sigma^2} \\
de^2 \propto 1 \\ 
dr^2 \propto \sigma^2
$$

## (Scalar) Energy response

A specific mutation $m$ at site $j$ will mutate `wt` to `mut`. As a result of the mutation, there're the following energy differences:
$$
\begin{eqnarray}
\delta V_0 &=& V_\mut(\rmut) - V_\wt(\rwt) \\
\delta U &=& U_\mut - U_\wt \\
\delta A &=& A_\mut - A_\wt 
\end{eqnarray}
$$
Let's generally call any of these scalar responses $\delta X(j, mut)$. The average over mutations will be a vector $\delta X(j)$, a profile representing the influence of site j on property $X$. 



### Turning energy response vectors into matrices

Several energy contributions can be turned into vectors.

#### General sum-over-edges to sum-over-sites

In general, any sum over edges can be express as a double sum over sites:
$$
\sum_{i \sim j} X_{ij} = \half \sum_{i}\sum_{j}\delta(i \sim j) X_{ij} = \sum_{i}X_i
$$
where $\delta(i \sim j)$ is $1$ for connected nodes and $0$ otherwise.
$$
X_i = \sum_{j}\delta(i \sim j) X_{ij}
$$
This is: $X_i$ is the sum over edges of $i$ of $X_{ij}$.

#### Vectorizing the minimum energy
In general, the network's minimum energy can be written:
$$
V_\mut(\rmut) = V_0 + \half \sum_{i \sim j} V^\mut_{0, ij} + \half \sum_{i\sim j}\kij^\mut(\dij^\mut - \lij^\mut)^2
$$


This is a constant, $V_0$ that can be arbitrarily assumed to be zero, plus two terms that can be separated into site contributions (the contribution of a site, is the sum over its contacts). Note that in the present version of the penm models, the second term is zero too.

#### Stress at $\wt$ conformation

When the mutant adopts the wild-type conformation, it's energy is given by:
$$
V_\mut(\wt) = V_0 + \half \sum_{i \sim j} V^\mut_{0, ij} + \half \sum \kij^\mut \dlij^2
$$
Again, this can be separated as a sum of node contributions that represent the local stress of the structure when it adopts the $\wt$ conformation. 

#### Relaxation from $\rwt$ to $\rmut$ 

The difference between the two previous energies is the energy needed to deform the mutant so that it adopts the wild-type conformation:
$$
V_\mut(\wt) - V_\mut(\rmut) \sim \half(\rmut - \rwt)^T\kmut (\rmut - \rwt)
$$
Assuming that $\kmut \sim \kwt$, this is precisely the norm of the vector $\delta \mathbf{e}$ defined previously:
$$
V_\mut(\wt) - V_\mut(\rmut) \sim ||\mathbf{\delta e||^2}
$$


#### Sum-over-sites entropy

#todo

Rigourusly, entropy can be separated only in the normal-mode representation, not in a site representation. 

However, intuitively, we understand that a change in the profile of mean-square fluctuations, $\sigma_i$ entails a change in entropies. It would be possible to derive approximate entropy formulas assuming for instance that sites oscillate independently, kind of a mean-field representation. To use $\sigma_i$ we could just assume spheric oscillators and express entropy as a sum of contributions of these spheric oscillators. 

When I want to pursue this further, I need to:

1. Derive a formula for site entropy, relating it to the 3D distribution of a given site.
2. Implement this formula in the `penm` package.
3. Assess the approximation: Check whether this sum-over-sites entropy is similar to the sum-over-modes entropy.
4. If it works, I entropy can be separated as a sum of site-contributons and entropy changes studied as changes of $\sigma_i$ profiles.





## Tricks for faster response calculation

### Several perturbations at once

Several response vectors are of the form:
$$
\svec = \amat \fvec
$$

where $\fvec$ is *one* perturbation (e.g. one contact, or one site). If we have several perturbations, (several sites, several mutations at same sites, whatever), we can write:
$$
\smat = \amat \fmat
$$

where $\fmat = (\fvec_1 \fvec_2 \ldots)$ is the matrix where each column is a forcing perturbation and $\smat = (\svec_1 \svec_2 \ldots)$ is the matrix of responses.



### Adding perturbations before matrix multiplication

Clearly:
$$
\svec_1 + \svec_2 = \amat (\fvec_1 + \fvec_2)
$$
Thus, for instance, the response to mutating a node is equal to the sum of the responses to mutating each of its contacts, which is what I use in the simulation.

This should also be used to make 2-mutation scans faster.

### Uses for *PRS*

The response to several mutations at a given site $j$ can be written:
$$
\smat^j = \amat \fmat^j
$$
The mean response vector would be:

```R
rowMeans(Sj^2)
```

The response to one mutation at each site (a single mutational scan), I could define a matrix `fmat(3 * nsites, nsites)` where column `j`contains the mutating foce at site `j`. Then, 

```R
smat <- amat %*% fmat
```

has in each column the response to that site's forcing.