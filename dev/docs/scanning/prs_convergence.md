# Convergence of simulation-based prs

Do the simulation-based prs response matrices and profile converge to the values predicted by the analytical formulae of [prs theory](superfast_response_theory.md)?

Check *prs_superfast.Rmd* and *prs_superfast.html*: they do. The convergence is rather slow, the inflection point is around 30 mutations per site. The convergence of the matrix is bad (a lot of mutations to reach approx 0.5 relative error), but the profiles converge better (to better than 0.1 relative error with about 20 mutations per site). 

Of course, this just validates the analytical method: it is **much faster** and **much more accurate**!

It is the **superfast and superaccurate prs**