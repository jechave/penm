# Functions for frustrated ENMs
# in general V = sum (Vij0 + .5* kij * (dij - lij)^2))

calculate_kijMat <- function(xyz.calpha,model,R0) {
    nsites <- length(xyz.calpha)/3
    dim(xyz.calpha) <- c(3,nsites)
    
    # obtain N x N kmat
    if(model=="hnm"){
        kmat <- K_hnm(xyz.calpha)
    } else if (model == "pfgnm") {
        kmat <- K_pfgnm(xyz.calpha)
    } else if (model == "bgnm") {
        kmat <- K_bgnm(xyz.calpha,site,R0)
    } else if (model == "gnm") {
        kmat <- K_gnm(xyz.calpha,R0)
    } else if (model == "hnm0") {
        kmat <- K_hnm0(xyz.calpha,R0)
    } else if (model == "reach") {
        kmat <- K_reach(xyz.calpha,site)
    } else {
        print("error in calculate_kijMat")
        print(paste("model:",model,"is undefined"))
    }
    kmat <- -kmat
    diag(kmat) = 0
    kmat
}

calculate_sijMat <- function(nsites) {
    sijMat = matrix(0,nsites,nsites)
    for(i in site) {for (j in site) sijMat[i,j] = j-i}
    sijMat
}

myRMSD <- function(xyz1,xyz2) {
    nsites = length(xyz1)/3
    dim(xyz1) = c(3,nsites)
    dim(xyz2) = c(3,nsites)
    msd = mean(colSums(xyz1-xyz2)^2)
    sqrt(msd)
}

calculate_nH = function(overlap) {
    nH = exp( rowSums(-overlap^2*log(overlap^2)) )
    nH
}

calculate_nR = function(overlap) {
    nR = 1/rowSums(overlap^4)
    nR
}

calculate_rmsip = function(u1,u2) {
    n1 = ncol(u1)
    n2 = ncol(u2)
    if(n1 != n2) print("warning: n1 != n2 in calculate_rmsip")
    n = n1
    overlap = t(u1) %*% u2
    msip = sum(overlap^2)/n
    rmsip = sqrt(msip)
    rmsip
}

calculate_rwsip = function(ua,ub,sa2,sb2) {
    nmodes = ncol(ua)
    overlap = t(ua) %*% ub
    wmn = matrix(0,nmodes,nmodes)
    for(m in seq(nmodes)) {
        for (n in seq(nmodes)) {
            wmn[m,n] = sa2[m]*sb2[n]
        }
    }
    rwsip = sqrt(sum(wmn*overlap^2)/sum(wmn))
    rwsip
}


calculateForce <- function(r,kijMat,lijMat) {
    # Calculates the force of ENM at r
    # kijMat is the NxN matrix of spring constans
    diag(kijMat) = 0 # make sure diagonals are 0
    nsites = nrow(kijMat)
    sites <- seq(nsites)
    dim(r) = c(3,nsites)
    force = matrix(0,nrow=3,ncol = nsites)
    for (i in sites) {
        sitesInContact <- sites[kijMat[i,] > 0]
        for (j in sitesInContact) {
            rij = r[,j] - r[,i]
            dij = sqrt(sum(rij^2))
            eij = rij/dij
            force[,i] = force[,i] + kijMat[i,j]*(dij - lijMat[i,j])*eij
        }
    }
    dim(force) = length(force)
    force
}



venm_frustrated <- function(r,vijMat,kijMat,lijMat) {
    require(bio3d)
    r <- as.vector(r)
    nsites <- nrow(kijMat)
    dijMat <- bio3d::dm.xyz(r,mask.lower = F)
    vij <- vijMat + .5*kijMat*(dijMat - lijMat)^2
    venm <- .5*sum(vij)
    venm
}

calculate_hessian <- function(r,kijMat,lijMat)  {
    nsites <- nrow(kijMat)
    dim(r) <- c(3,nsites)
    hessian <- matrix(0,3*nsites,3*nsites)
    sites <- seq(nsites)
    ind.sites <- seq(3*nsites)
    dim(ind.sites) <- c(3,nsites)
    idMat <- diag(1,3,3)
    for (i in sites) {
        inds.i <- as.vector(ind.sites[,i])
        jSites <- sites[kijMat[i,] > 0 ]
        for (j in jSites) {
            if(j > i) {
            inds.j <- as.vector(ind.sites[,j])
            kij <- kijMat[i,j]
            rij <-r[,j] - r[,i]
            dij <- sqrt(sum(rij^2))
            eij <- as.vector(rij)/dij
            eijMat <- eij %*% t(eij)
            hessian[inds.i,inds.j] <- kij*(eijMat + (1-lijMat[i,j]/dij)*(eijMat - idMat))
            hessian[inds.j,inds.i] <- hessian[inds.i,inds.j]
            hessian[inds.i,inds.i] <- hessian[inds.i,inds.i] - hessian[inds.i,inds.j]
            hessian[inds.j,inds.j] <- hessian[inds.j,inds.j] - hessian[inds.i,inds.j] 
            }
        }
    }
    -hessian
}
            


calculate_hessian_k <- function(r,kijMat,lijMat)  {
    nsites <- nrow(kijMat)
    dim(r) <- c(3,nsites)
    hessian <- matrix(0,3*nsites,3*nsites)
    sites <- seq(nsites)
    ind.sites <- seq(3*nsites)
    dim(ind.sites) <- c(3,nsites)
    idMat <- diag(1,3,3)
    for (i in sites) {
        inds.i <- as.vector(ind.sites[,i])
        jSites <- sites[kijMat[i,] > 0 ]
        for (j in jSites) {
            if(j > i) {
            inds.j <- as.vector(ind.sites[,j])
            kij <- kijMat[i,j]
            rij <-r[,j] - r[,i]
            dij <- sqrt(sum(rij^2))
            eij <- as.vector(rij)/dij
            eijMat <- eij %*% t(eij)
            hessian[inds.i,inds.j] <- kij*eijMat 
            hessian[inds.j,inds.i] <- hessian[inds.i,inds.j]
            hessian[inds.i,inds.i] <- hessian[inds.i,inds.i] - hessian[inds.i,inds.j]
            hessian[inds.j,inds.j] <- hessian[inds.j,inds.j] - hessian[inds.i,inds.j] 
            }
        }
    }
    
    
    -hessian
}

calculate_rmin <- function(r0,kijMat,lijMat,maxit=10,TOL=1.e-5) {
    r0 = as.vector(r0) #initial conformation
    nsites = length(r0)/3
    nmodes = length(r0)
    for (n in seq(maxit)) {
        f = calculateForce(r0,kijMat,lijMat)
        H = calculate_hessian_k(r0,kijMat,lijMat)
        H.eigen = eigen(H,symmetric = T)
        U = H.eigen$vectors[,1:(nmodes-6)]
        sigma = diag(1/H.eigen$values[1:(nmodes-6)])
        C = U %*% sigma %*% t(U)
        dr = C %*% f
        r0 = r0 + dr
        print(c("rmsd:",sqrt(3*mean(dr^2)),"rmsf:",sqrt(3*mean(f^2))))
        if (sqrt(sum(f^2)) < TOL) break
        if (sqrt(mean(f^2)) < TOL) break
    }
    r0
}
        
        
         
    
    


    
    
            
        
        
    
    