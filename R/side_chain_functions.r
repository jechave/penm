residue.coordinates = function(pdb,d=1.5) {
# returns a list with the coordinates:
# alpha carbon, beta carbon, side-chain center of mass
# and a "quasi-beta" atom placed in the direction of Ca(i-1)->Ca(i) + Ca(i+1)->Ca(i) at a distnace l=1.5A
    ca.inds = atom.select(pdb,"calpha",verbose=FALSE)
    cb.inds = atom.select(pdb,"protein", elety = "CB",verbose=FALSE)
    backbone.inds = atom.select(pdb,"backbone",verbose=FALSE)

    site = pdb$atom$resno[ca.inds$atom] #resno of all sites (with CA)
    cb.resno = pdb$atom$resno[cb.inds$atom] #resno of all sites wih CB

    aa = aa321(pdb$atom$resid[ca.inds$atom]) #aa of AA with CA
    nsites = length(aa)

    ca.xyz = c()
    cb.xyz = c()
    com.xyz = c()

    for (resno in site) {
        res.inds = atom.select(pdb,resno=resno,verbose=FALSE) # residue resno
# CA coordinates
        res.ca.inds = combine.select(res.inds,ca.inds,operator = "and", verbose = FALSE) # CA of resno
        res.ca.xyz = pdb$xyz[res.ca.inds$xyz]
        ca.xyz = c(ca.xyz,res.ca.xyz)
# CB coordinates (when there's a CB)
        res.cb.inds = combine.select(res.inds,cb.inds,operator="and",verbose=FALSE) # CB of resno
        res.cb.xyz = c(NA,NA,NA)
        if (length(res.cb.inds$atom != 0))  res.cb.xyz = pdb$xyz[res.cb.inds$xyz]
        cb.xyz = c(cb.xyz,res.cb.xyz)
# side chain center of mass coordinates (when there's a side chain)
        res.sc.inds = combine.select(res.inds,backbone.inds,operator = "not",verbose=FALSE) # Side chain of resno
        res.com.xyz = c(NA,NA,NA)
        if (length(res.sc.inds$atom != 0)) res.com.xyz= com(pdb,inds=res.sc.inds,use.mass=FALSE)
        com.xyz = c(com.xyz,res.com.xyz)
    }

   m.xyz = qb.micheletti(ca.xyz)
   l.xyz = qb.levitt(ca.xyz)

   r = list("site"=site,"aa"=aa,"ca.xyz"=ca.xyz,"cb.xyz"=cb.xyz,"com.xyz"=com.xyz,"m.xyz"=m.xyz,"l.xyz"=l.xyz)
   r
}


residue.bfactors = function(pdb) {
# returns a data.frame with alpha, beta, and com bfactors
    ca.inds = atom.select(pdb,"calpha",verbose=FALSE)
    cb.inds = atom.select(pdb,"protein", elety = "CB",verbose=FALSE)
    backbone.inds = atom.select(pdb,"backbone",verbose=FALSE)

    site = pdb$atom$resno[ca.inds$atom] #resno of all sites (with CA)
    cb.resno = pdb$atom$resno[cb.inds$atom] #resno of all sites wih CB

    ca.bfactor = c()
    cb.bfactor = c()
    com.bfactor = c()
    l.bfactor = c()
    m.bfactor = c()

    for (resno in site) {
        res.inds = atom.select(pdb,resno=resno,verbose=FALSE) # residue resno
# CA bfactor
        res.ca.inds = combine.select(res.inds,ca.inds,operator ="and",verbose=FALSE) # CA of resno
        res.ca.bfactor = pdb$atom$b[res.ca.inds$atom]
        ca.bfactor = c(ca.bfactor,res.ca.bfactor)

# CB bfactor (when there's a CB)
        res.cb.inds = combine.select(res.inds,cb.inds,operator ="and",verbose=FALSE) # CB of resno
        res.cb.bfactor = NA
        if (length(res.cb.inds$atom != 0))  res.cb.bfactor = pdb$atom$b[res.cb.inds$atom]
        cb.bfactor = c(cb.bfactor,res.cb.bfactor)
# side chain center of mass bfactor (when there's a side chain)
        res.sc.inds = combine.select(res.inds,backbone.inds,operator ="not",verbose=FALSE) # Side chain of resno
        res.com.bfactor = NA
        if (length(res.sc.inds$atom != 0)) {
            # res.com.bfactor = mean(sqrt(pdb$atom$b[res.sc.inds$atom]))^2
            res.com.bfactor = mean(pdb$atom$b[res.sc.inds$atom])
        }
        com.bfactor = c(com.bfactor,res.com.bfactor)
    }
# levitt and micheletti bfactors (they are identical, but keep them named differently for now)
    l.bfactor = bfactor.levitt(ca.bfactor)
    m.bfactor = bfactor.micheletti(ca.bfactor)

    df.bfactor = data.frame("b.a"=ca.bfactor,"b.b"=cb.bfactor,"b.c"=com.bfactor,"b.l"=l.bfactor,"b.m"=m.bfactor)

   df.bfactor
}


vnorm = function(v)  sqrt(sum(v^2))

theta = function(v1,v2) {
  th = acos(sum(v1*v2)/sqrt(sum(v1^2)*sum(v2^2)))
  th
}

distance = function(r1,r2) {
  r12 = r2-r1
  d12 = sqrt(sum(r12^2))
  d12
}

qb.levitt = function(ca.xyz,l=3,theta=37.6){
  require(pracma)
  # quasi-beta (centroid) coordinates according to Levitt 1996
  nsites = length(ca.xyz)/3
  dim(ca.xyz) = c(3,nsites)
  qb.xyz = c()
  for (i in seq(nsites)) {
    if(i == 1 || i == nsites) {
      #            q.i = ca.xyz[,i]
      q.i = c(NA,NA,NA)
    } else {
      rm = unlist(ca.xyz[,i] - ca.xyz[,i-1])
      rp = unlist(ca.xyz[,i] - ca.xyz[,i+1])
      x = rp+rm
      x = x/vnorm(x)
      y = cross(rp,rm)
      y = y/vnorm(y)
      q.i = ca.xyz[,i]+x*l*cos(theta)+y*l*sin(theta)
    }
    qb.xyz = c(qb.xyz,q.i)
  }
  qb.xyz
}

qb.micheletti = function(ca.xyz,l=3){
  require(pracma)
  # quasi-beta (centroid) coordinates according to Micheletti 2004 (also Jernigan and HSCa)
  nsites = length(ca.xyz)/3
  dim(ca.xyz) = c(3,nsites)
  qb.xyz = c()
  for (i in seq(nsites)) {
    if(i == 1 || i == nsites) {
      #            q.i = ca.xyz[,i]
      q.i = c(NA,NA,NA)
    } else {
      rm = unlist(ca.xyz[,i] - ca.xyz[,i-1])
      rp = unlist(ca.xyz[,i] - ca.xyz[,i+1])
      x = rp+rm
      x = x/vnorm(x)
      q.i = ca.xyz[,i]+x*l
    }
    qb.xyz = c(qb.xyz,q.i)
  }
  qb.xyz
}


bfactor.levitt = function(ca.bfactor){
  nsites = length(ca.bfactor)
  qb.bfactor = c()
  for (i in seq(nsites)) {
    if(i == 1 || i == nsites) {
      qb.bfactor.i = NA
    } else {
      x = c(ca.bfactor[i-1],ca.bfactor[i],ca.bfactor[i+1])
      qb.bfactor.i = mean(x)
    }
    qb.bfactor = c(qb.bfactor,qb.bfactor.i)
  }
  qb.bfactor
}

bfactor.micheletti = function(ca.bfactor){
  nsites = length(ca.bfactor)
  qb.bfactor = c()
  for (i in seq(nsites)) {
    if(i == 1 || i == nsites) {
      qb.bfactor.i = NA
    } else {
      x = c(ca.bfactor[i-1],ca.bfactor[i],ca.bfactor[i+1])
      qb.bfactor.i = mean(x)
    }
    qb.bfactor = c(qb.bfactor,qb.bfactor.i)
  }
  qb.bfactor
}
