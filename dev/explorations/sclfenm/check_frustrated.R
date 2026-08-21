## keepnet states are STRAINED, so their true Hessian has the transverse term.
## Every keepnet result in the report used frustrated=FALSE. How much does it matter?
source("applications.R")
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max=10.5); N <- length(r0)/3
st <- wt; set.seed(9)
for (i in 1:30) st <- mutate(st, sample(N,1), sigma=.3, model="keepnet")$state
rq <- state_req(st); st$r <- rq$r
strain <- dists(st$r,st$pr) - st$l
cat(sprintf("30-step keepnet state: strain energy %.3f, max|d-l| %.3f\n\n",
    0.5*sum(st$k*strain^2), max(abs(strain))))
Kf <- hessian(st$r, st$pr, st$k, st$l, frustrated=FALSE)
Kt <- hessian(st$r, st$pr, st$k, st$l, frustrated=TRUE)
cat(sprintf("||K_false - K_true||_max = %.4f   (scale max|K| = %.2f, %.1f%%)\n",
    max(abs(Kf-Kt)), max(abs(Kt)), 100*max(abs(Kf-Kt))/max(abs(Kt))))
af <- nma(Kf); at <- nma(Kt)
cat(sprintf("zero modes: false %d, true %d\n", af$nzero, at$nzero))
cat(sprintf("min eigenvalue: false %.4e, true %.4e  %s\n", min(af$evalue), min(at$evalue),
    if (min(at$evalue) < 0) "<-- SADDLE under the true Hessian!" else ""))
n <- min(length(af$evalue), length(at$evalue))
cat(sprintf("softest 10 eigenvalues, relative difference: %.3f%%\n",
    100*mean(abs(sort(at$evalue)[1:10]-sort(af$evalue)[1:10])/sort(at$evalue)[1:10])))
## does it change dV?
set.seed(4); dvF <- dvT <- numeric(150)
stT <- st; stT$frustrated <- TRUE
for (i in 1:150) {
  s <- sample(N,1)
  set.seed(100+i); dvF[i] <- mutate(st,  s, sigma=.3, model="keepnet")$dV
  set.seed(100+i); dvT[i] <- mutate(stT, s, sigma=.3, model="keepnet")$dV
}
cat(sprintf("\ndV over 150 mutations:  frustrated=FALSE mean %.4f  frac<0 %.3f\n",
    mean(dvF), mean(dvF<0)))
cat(sprintf("                        frustrated=TRUE  mean %.4f  frac<0 %.3f\n",
    mean(dvT), mean(dvT<0)))
cat(sprintf("  correlation between the two: %.5f\n", cor(dvF,dvT)))
