# Study the effect of frustration on ENMs
# packages
rm(list = ls()) # clear workspace

library(tidyverse) #dplyr, ggplot2, purrr, etc.
library(bio3d)
library(corpcor) # package where pseudoinverse(m,tol) is

library(matrixStats)

# functions
source("scripts/ganm_v2/readCA.R")
source("scripts/ganm_v2/ganm-functions.R")
source("scripts/ganm_v2/reduce_matrix.R")
source("scripts/mutation-model.R")
source("scripts/enm-frustrated-functions.R")




# parameters
TOLERANCE <- 1.e-10
deltaijMax = .5

# files and directories
# input file names
pdbdir <- "data/repaired_pdbs"
models.fname <- "data/models.csv"
dataset.fname <- "data/monomers_dataset.csv"


# input
models <- read.csv(file=models.fname,header=T) # data.frame(modelid,model,R0)
dataset <- read.csv(file=dataset.fname)
dataBenJack <- read.csv(file="data/monomers_data_table.csv.gz")
nmodels <- nrow(models)
nprot <- nrow(dataset)

# parameters

# list of pdbs to calculate
pdbList = dataset$pdb

#debug: to check use only the first protein

prot = 1

    pdbid <- as.character(dataset$pdb[prot])
    chain <- as.character(dataset$chain[prot])
    print(paste("protein:",pdbid,"chain:",chain))
    
    dataProt <- dataBenJack %>%
        filter(pdb == pdbid) %>%
        droplevels()
    nsites <- nrow(dataProt)
    

    # read structure and select CA's
    pdb.fname <- paste("RepairPDB",pdbid,chain,sep="_")
    pdb.fname <- paste(pdb.fname,".pdb",sep="")
    pdb.fname <- file.path(pdbdir,pdb.fname) # structure
    pdb <- readCA(pdb.fname)
    xyz.calpha <- pdb$xyz.calpha
    site <- pdb$site
    if(nsites != pdb$nsites) print("Error: nsites != pdb$nsites")
    
    m <- 1
    mm <- models$modelid[m]
    print(mm)
    R0 = 10
    enm <- ganm(xyz.calpha,site,model=models$model[m],R0=models$R0[m],TOL=TOLERANCE)
    kijMat <- K_gnm(xyz.calpha,R0)
    kijMat <- enm$kijMat
   
    lijMat <- bio3d::dm.xyz(xyz=as.vector(xyz.calpha),mask.lower = FALSE)
    lijMat[kijMat == 0] <- 0 
    vijMat <- matrix(0,nsites,nsites)
    
    r <- as.vector(xyz.calpha)
    venm_frustrated(r,vijMat,kijMat,lijMat)
    
    f <- calculateForce(as.vector(xyz.calpha),kijMat,lijMat)
    
    deltaij <- matrix(rnorm(nsites*nsites,0,.5),nsites,nsites)
    deltaij[kijMat == 0] = 0
    deltaij <- .5*(deltaij + t(deltaij))
    
    
    lijMatMutant <- lijMat + deltaij
    f <- calculateForce(as.vector(xyz.calpha),kijMat,lijMatMutant)
    dr <-  enm$cmat %*% as.vector(f)
    r <- as.vector(xyz.calpha) + dr
    f.new <- calculateForce(as.vector(r),kijMat,lijMatMutant)
    dr <- enm$cmat %*% as.vector(f.new)
    r.new <- as.vector(r) + dr
    f.2<- calculateForce(as.vector(r.new),kijMat,lijMatMutant)
    dr <- enm$cmat %*% as.vector(f.2)
    r.2 <- as.vector(r.new) + dr
    
    venm_frustrated(as.vector(xyz.calpha),vijMat,kijMat,lijMatMutant)
    venm_frustrated(r,vijMat,kijMat,lijMatMutant)
    venm_frustrated(r.new,vijMat,kijMat,lijMatMutant)
    venm_frustrated(r.2,vijMat,kijMat,lijMatMutant)
    
    e <- as.vector(r) - as.vector(xyz.calpha)
    dat = data.frame()
    for (s in seq(from=-1,to=2,by=.1)) {
        xyz = as.vector(xyz.calpha) + s*e
        vwt <- venm_frustrated(xyz,vijMat,kijMat,lijMat)
        vMutant <- venm_frustrated(xyz,vijMat,kijMat,lijMatMutant)
        row = data.frame(s,vwt,vMutant)
        dat = rbind(dat,row)
    }
    dat2 <- dat %>% gather(key = "case",value="v",vwt,vMutant)
    ggplot(dat2,aes(x=s,y=v,col=case)) + geom_path()
    
   # kmat.mutant <- calculate_k(r.2,kijMat,lijMatMutant)
    xyz <- r.2
    dim(xyz) = c(3,nsites)
    enm.mutant <- ganm(xyz,site,model=models$model[m],R0=models$R0[m],TOL=TOLERANCE)
    kmat.mutant <- enm.mutant$kmat
    hmat.mutant <- calculate_hessian(r.2,kijMat,lijMatMutant)
    
    kmat.eigen <- eigen(kmat.mutant,symmetric = T)
    hmat.eigen <- eigen(hmat.mutant,symmetric = T)
    plot(hmat.eigen$values,kmat.eigen$values)
    sk <- t(kmat.eigen$vectors) %*% kmat.eigen$vectors 
    sh <- t(hmat.eigen$vectors) %*% hmat.eigen$vectors 
    shk <- t(hmat.eigen$vectors) %*% kmat.eigen$vectors 
    
    plot(rank(kmat.eigen$values),rowMaxs(shk^2))
    nh = exp(rowSums(-shk^2*log(shk^2)))
    plot(rank(kmat.eigen$values),nh)
    
    kmat.inv = pseudoinverse(kmat.mutant,tol=1.e-3)
    hmat.inv = pseudoinverse(hmat.mutant,tol=1.e-3)
    plot(diag(kmat.inv))
    plot(diag(hmat.inv))
    plot(diag(kmat.inv),diag(hmat.inv))
