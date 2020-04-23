# calculate lpd profiles and join with rsa, bfactors, and empirical rate profiles
# write lpd profiles in output file prof.out.fname (see below)
require(bio3d)

source("sc_lpd_functions.r")

# parameters
enm.model = "wcn" # alternatives: cn and wcn (corresponds to anm and pfanm)
R0 = 150  # use large value for wcn, (e.g. 150 to include all sites, smaller to sppeed up calculation, e.g. R0=15)

# input 
#data = read.csv("input/dataset_monomers_1.csv")
data = read.csv("input/dataset_monomers_209.csv") # read dataset file
prof.in = read.csv("input/prof_monomers_209.csv") # read empirical rates

nprot = nrow(data)
prof.out = data.frame()
for (prot in seq(nprot)) {
    pdb.id = as.character(data$pdb[prot])
    chain  = as.character(data$chain[prot])
    prot.inds = prof.in$pdb == pdb.id
    nsites = sum(prot.inds)
    print(c(prot,"prot",pdb.id,"sites:",nsites))
# keep this protein's profiles
    prof.prot = prof.in[prot.inds,]
# read pdb file
    pdb.file = paste("input/repaired_pdbs/RepairPDB_",pdb.id,"_",chain,".pdb",sep="")
    pdb = read.pdb(pdb.file)
    if (!inspect.connectivity(pdb)){
        print(paste("Warning: ",pdb.id,"has missing residues"))
    } 
# calculate xyz coordinates of ca, cb, com, qm (qb by Micheletti), ql (qb by Levitt)
    r = residue.coordinates(pdb) 
# calculate bfactors
    df.bfactors = residue.bfactors(pdb) # returns bfactos for alpha carbon, beta carbon, and (approximate for) center of mass of side chain
# calculate site-specific LPD profiles
    if (enm.model == "wcn") {
        lpd.prot = calculate.wcn(r,R0) # calculates wcn terms (aa,ab, etc) and returns them in a data.frame
        prof.out.fname = paste("prof_out_",enm.model,".csv",sep="")
    } else if (enm.model == "cn") {
        lpd.prot = calculate.cn(r,R0) # calculates wcn terms (aa,ab, etc) and returns them in a data.frame
        prof.out.fname = paste("prof_out_",enm.model,"_R0_",as.character(R0),".csv",sep="")
    } else {
        print(c("unknown enm.model:",enm.model))
    }
# bind all profiles together
    prof.prot = cbind(prof.prot,df.bfactors,lpd.prot)
# add to profiles of previous proteins
    prof.out = rbind(prof.out,prof.prot)
}

# output all profiles
row.names(prof.out) = NULL

write.csv(prof.out,prof.out.fname,row.names=F)
