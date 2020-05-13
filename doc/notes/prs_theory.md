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