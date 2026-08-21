# Daily notes

## 20 May 2020

In *prs.R* and *prs_fast.R* I'm repeating calculations in `de2_site`, `dr2_site`,  `df2_site`,  and the mode counterparts. I could make the calculation by avoiding this repetition.



## 5 August 2020

### Faster prs using force matrices (rather than vectors)

Attempting to optimize `prs.new`. First, I used a matrix of size `(3 * nsites, nmut)` for each site, and calculated the response by multiplying this matrix by, e.g. the covariance matrix all at once (rather than each mutation at a time, as in the "simulation" approach). 



Second, I generalized this to setting up a force matrix of size `(3 * nsites, nsites * nmut)` representing the whole mutational scan at once, then calculated the responses by multiplying by, e.g. the covariance matrix. 



Third, the "total-scan" force matrix is sparse, so I tried to use this for faster matrix multiplication.



#### Similar cpu times for all approaches

The two sparse versions are the fastest (i.e. options 1 and 3).

Also, probably behind these similarities, the bottleneck is not matrix multiplication (not anymore), but getting the forces that simulate the mutations.





