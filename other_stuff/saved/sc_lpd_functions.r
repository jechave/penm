mylm = function(y,x){
# a linear fit that deals with non-finite values
# it returns a list with elements fitted.values, r.squared, and sigma
    # turn +-Inf into NA
    x[ which( !is.finite(x) ) ] = NA
    y[ which( !is.finite(y) ) ] = NA
    fit = lm(y~x,na.action=na.omit)
    fitted.values = predict(fit,as.data.frame(x))
    sfit = summary(fit)
    r.squared = sfit$r.squared
    sigma = sfit$sigma
    results = list("fitted.values"=fitted.values,"r.squared"=r.squared,"sigma"=sigma)
    results
}

mylm2 = function(y,x1,x2){
# a linear fit that deals with non-finite values
# it returns a list with elements fitted.values, r.squared, and sigma
    # turn +-Inf into NA
    x1[ which( !is.finite(x1) ) ] = NA
    x2[ which( !is.finite(x2) ) ] = NA
    y[ which( !is.finite(y) ) ] = NA
    fit = lm(y~x1+x2,na.action=na.omit)
    fitted.values = predict(fit,data.frame(x1,x2))
    sfit = summary(fit)
    r.squared = sfit$r.squared
    sigma = sfit$sigma
    results = list("fitted.values"=fitted.values,"r.squared"=r.squared,"sigma"=sigma)
    results
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

            
residue.coordinates = function(pdb,d=1.5) {
# returns a data.frame with one row/residue with the following residue coordinates:
# alpha carbon, beta carbon, side-chain center of mass
# and a "quasi-beta" atom placed in the direction of Ca(i-1)->Ca(i) + Ca(i+1)->Ca(i) at a distnace l=1.5A
    ca.inds = atom.select(pdb,"calpha",verbose=FALSE)
    cb.inds = atom.select(pdb,"//////CB/",verbose=FALSE)
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
        res.ca.inds = combine.sel(res.inds,ca.inds,op="and",verbose=FALSE) # CA of resno
        res.ca.xyz = pdb$xyz[res.ca.inds$xyz]
        ca.xyz = c(ca.xyz,res.ca.xyz)
# CB coordinates (when there's a CB)
        res.cb.inds = combine.sel(res.inds,cb.inds,op="and",verbose=FALSE) # CB of resno
        res.cb.xyz = c(NA,NA,NA)
        if (length(res.cb.inds$atom != 0))  res.cb.xyz = pdb$xyz[res.cb.inds$xyz]
        cb.xyz = c(cb.xyz,res.cb.xyz)
# side chain center of mass coordinates (when there's a side chain)
        res.sc.inds = combine.sel(res.inds,backbone.inds,op="not",verbose=FALSE) # Side chain of resno
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
    cb.inds = atom.select(pdb,"//////CB/",verbose=FALSE)
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
        res.ca.inds = combine.sel(res.inds,ca.inds,op="and",verbose=FALSE) # CA of resno
        res.ca.bfactor = pdb$atom$b[res.ca.inds$atom]
        ca.bfactor = c(ca.bfactor,res.ca.bfactor)
# CB bfactor (when there's a CB)
        res.cb.inds = combine.sel(res.inds,cb.inds,op="and",verbose=FALSE) # CB of resno
        res.cb.bfactor = NA
        if (length(res.cb.inds$atom != 0))  res.cb.bfactor = pdb$atom$b[res.cb.inds$atom]
        cb.bfactor = c(cb.bfactor,res.cb.bfactor)
# side chain center of mass bfactor (when there's a side chain)
        res.sc.inds = combine.sel(res.inds,backbone.inds,op="not",verbose=FALSE) # Side chain of resno
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


# Calculation of CN

cnij = function(dij,R0=7.5){
# cnij is a step function
# Note: idf dij=NA it returns NA
    if (is.na(dij)) {
        wij = NA
    } else if (dij <= R0) {
        wij = 1
    } else {
        wij = 0
    }
    wij
}

wcnij = function(dij){
# wcnij.max corresopnds to observed maximum value for i-i+1 Calphas
# it prevents artifacts when using quasi-centroids using levitt's or micheletti's approximations
# Note: idf dij=NA it returns NA
    #wcnij.max = 0.07
    wcnij.max = 0.1
    wij = min(1/dij^2,wcnij.max)
    wij

}

sum.lpd = function(wij){
# wij is a vector of wij values for a given i
# if the ith center is undefined, all values are NA: return NA
# if the jth center is undefined, the jth value is NA: don't consider this value, do consider the others
    if(sum(!is.na(wij)) == 0) {
        wi = NA
    } else {
        wi = sum(na.exclude(wij))
    }
    wi
}
    
calculate.wcn = function(r,R0=150) {
# calculates different wcn components, returns them in a data.frame
# R0 is a cut-off distance used to make calculation faster for very large proteins

    ca.xyz = r$ca.xyz
    cb.xyz = r$cb.xyz
    com.xyz = r$com.xyz
    m.xyz = r$m.xyz
    l.xyz = r$l.xyz
    nsites = length(r$ca.xyz)/3
    dim(ca.xyz) = c(3,nsites)
    dim(cb.xyz) = c(3,nsites)
    dim(com.xyz) = c(3,nsites)
    dim(m.xyz) = c(3,nsites)
    dim(l.xyz) = c(3,nsites)

    df = data.frame()

    for (i in seq(nsites)) {
        #print(c("site",i))
        resno.i = r$site[i]
        aa.i = r$aa[i]
# coordinates of network nodes associated to i
        ca.i = ca.xyz[,i]
        cb.i = cb.xyz[,i]
        com.i = com.xyz[,i]
        m.i = m.xyz[,i]
        l.i = l.xyz[,i]

# distances
# of a and b with own site
        d.ai.bi = vnorm(cb.i - ca.i)
        d.ai.ci = vnorm(com.i - ca.i)
        d.ai.mi = vnorm(m.i - ca.i)
        d.ai.li = vnorm(l.i - ca.i)

# distance vectors of i with j
        d.ai.aj = c()
        d.ai.bj = c()
        d.ai.cj = c()
        d.ai.lj = c()
        d.ai.mj = c()
        d.bi.aj = c()
        d.bi.bj = c()
        d.ci.aj = c()
        d.ci.cj = c()
        d.li.aj = c()
        d.li.lj = c()
        d.mi.aj = c()
        d.mi.mj = c()
        for (j in seq(nsites)) {
            if (i != j) {
# coordinates of network nodes associated to j
                ca.j = ca.xyz[,j]
                daa = vnorm(ca.j-ca.i)
                if(daa <= R0) { # don't consider pairs with too large daa
                    cb.j = cb.xyz[,j]
                    com.j = com.xyz[,j]
                    m.j = m.xyz[,j]
                    l.j = l.xyz[,j]

# distances between nodes of i and nodes of j
                    d.ai.aj = c(d.ai.aj,vnorm(ca.j-ca.i))
                    d.ai.bj = c(d.ai.bj,vnorm(cb.j-ca.i))
                    d.ai.cj = c(d.ai.cj,vnorm(com.j-ca.i))
                    d.ai.mj = c(d.ai.mj,vnorm(m.j-ca.i))
                    d.ai.lj = c(d.ai.lj,vnorm(l.j-ca.i))
                    d.bi.aj = c(d.bi.aj,vnorm(ca.j-cb.i))
                    d.bi.bj = c(d.bi.bj,vnorm(cb.j-cb.i))
                    d.ci.aj = c(d.ci.aj,vnorm(ca.j-com.i))
                    d.ci.cj = c(d.ci.cj,vnorm(com.j-com.i))
                    d.li.aj = c(d.li.aj,vnorm(ca.j-l.i))
                    d.li.lj = c(d.li.lj,vnorm(l.j-l.i))
                    d.mi.aj = c(d.mi.aj,vnorm(ca.j-m.i))
                    d.mi.mj = c(d.mi.mj,vnorm(m.j-m.i))

                }
            }
        }
        wcn.ai.bi = 1/d.ai.bi^2
        wcn.ai.ci = 1/d.ai.ci^2
        wcn.ai.li = 1/d.ai.li^2
        wcn.ai.mi = 1/d.ai.mi^2

        wcn.ai.aj = sum.lpd(sapply(d.ai.aj,FUN=wcnij))
        wcn.ai.bj = sum.lpd(sapply(d.ai.bj,FUN=wcnij))
        wcn.ai.cj = sum.lpd(sapply(d.ai.cj,FUN=wcnij))
        wcn.ai.mj = sum.lpd(sapply(d.ai.mj,FUN=wcnij))
        wcn.ai.lj = sum.lpd(sapply(d.ai.lj,FUN=wcnij))
        wcn.bi.aj = sum.lpd(sapply(d.bi.aj,FUN=wcnij))
        wcn.bi.bj = sum.lpd(sapply(d.bi.bj,FUN=wcnij))
        wcn.ci.aj = sum.lpd(sapply(d.ci.aj,FUN=wcnij))
        wcn.ci.cj = sum.lpd(sapply(d.ci.cj,FUN=wcnij))
        wcn.li.aj = sum.lpd(sapply(d.li.aj,FUN=wcnij))
        wcn.li.lj = sum.lpd(sapply(d.li.lj,FUN=wcnij))
        wcn.mi.aj = sum.lpd(sapply(d.mi.aj,FUN=wcnij))
        wcn.mi.mj = sum.lpd(sapply(d.mi.mj,FUN=wcnij))
        df = rbind(df,data.frame(
                            i,resno.i,aa.i,wcn.ai.bi,wcn.ai.ci,wcn.ai.li,wcn.ai.mi,
                            wcn.ai.aj,wcn.ai.bj,wcn.ai.cj,wcn.ai.lj,wcn.ai.mj,
                            wcn.bi.aj,wcn.bi.bj,
                            wcn.ci.aj,wcn.ci.cj,
                            wcn.li.aj,wcn.li.lj,
                            wcn.mi.aj,wcn.mi.mj))
    }
    df
}

calculate.cn = function(r,R0=7.5) {
# calculates different wcn components, returns them in a data.frame
# R0 is a cut-off distance used to make calculation faster for very large proteins

    ca.xyz = r$ca.xyz
    cb.xyz = r$cb.xyz
    com.xyz = r$com.xyz
    m.xyz = r$m.xyz
    l.xyz = r$l.xyz
    nsites = length(r$ca.xyz)/3
    dim(ca.xyz) = c(3,nsites)
    dim(cb.xyz) = c(3,nsites)
    dim(com.xyz) = c(3,nsites)
    dim(m.xyz) = c(3,nsites)
    dim(l.xyz) = c(3,nsites)

    df = data.frame()

    for (i in seq(nsites)) {
        #print(c("site",i))
        resno.i = r$site[i]
        aa.i = r$aa[i]
# coordinates of network nodes associated to i
        ca.i = ca.xyz[,i]
        cb.i = cb.xyz[,i]
        com.i = com.xyz[,i]
        m.i = m.xyz[,i]
        l.i = l.xyz[,i]

# distances
# of a and b with own site
        d.ai.bi = vnorm(cb.i - ca.i)
        d.ai.ci = vnorm(com.i - ca.i)
        d.ai.mi = vnorm(m.i - ca.i)
        d.ai.li = vnorm(l.i - ca.i)

# distance vectors of i with j
        d.ai.aj = c()
        d.ai.bj = c()
        d.ai.cj = c()
        d.ai.lj = c()
        d.ai.mj = c()
        d.bi.aj = c()
        d.bi.bj = c()
        d.ci.aj = c()
        d.ci.cj = c()
        d.li.aj = c()
        d.li.lj = c()
        d.mi.aj = c()
        d.mi.mj = c()
        for (j in seq(nsites)) {
            if (i != j) {
# coordinates of network nodes associated to j
                ca.j = ca.xyz[,j]
                daa = vnorm(ca.j-ca.i)
                if(daa <= 3*R0) { # don't consider pairs with too large daa
                    cb.j = cb.xyz[,j]
                    com.j = com.xyz[,j]
                    m.j = m.xyz[,j]
                    l.j = l.xyz[,j]

# distances between nodes of i and nodes of j
                    d.ai.aj = c(d.ai.aj,vnorm(ca.j-ca.i))
                    d.ai.bj = c(d.ai.bj,vnorm(cb.j-ca.i))
                    d.ai.cj = c(d.ai.cj,vnorm(com.j-ca.i))
                    d.ai.mj = c(d.ai.mj,vnorm(m.j-ca.i))
                    d.ai.lj = c(d.ai.lj,vnorm(l.j-ca.i))
                    d.bi.aj = c(d.bi.aj,vnorm(ca.j-cb.i))
                    d.bi.bj = c(d.bi.bj,vnorm(cb.j-cb.i))
                    d.ci.aj = c(d.ci.aj,vnorm(ca.j-com.i))
                    d.ci.cj = c(d.ci.cj,vnorm(com.j-com.i))
                    d.li.aj = c(d.li.aj,vnorm(ca.j-l.i))
                    d.li.lj = c(d.li.lj,vnorm(l.j-l.i))
                    d.mi.aj = c(d.mi.aj,vnorm(ca.j-m.i))
                    d.mi.mj = c(d.mi.mj,vnorm(m.j-m.i))

                }
            }
        }
        cn.ai.bi = cnij(d.ai.bi,R0)
        cn.ai.ci = cnij(d.ai.ci,R0)
        cn.ai.li = cnij(d.ai.li,R0)
        cn.ai.mi = cnij(d.ai.mi,R0)

        cn.ai.aj = sum.lpd(sapply(d.ai.aj,FUN=cnij,R0=R0))
        cn.ai.bj = sum.lpd(sapply(d.ai.bj,FUN=cnij,R0=R0))
        cn.ai.cj = sum.lpd(sapply(d.ai.cj,FUN=cnij,R0=R0))
        cn.ai.mj = sum.lpd(sapply(d.ai.mj,FUN=cnij,R0=R0))
        cn.ai.lj = sum.lpd(sapply(d.ai.lj,FUN=cnij,R0=R0))
        cn.bi.aj = sum.lpd(sapply(d.bi.aj,FUN=cnij,R0=R0))
        cn.bi.bj = sum.lpd(sapply(d.bi.bj,FUN=cnij,R0=R0))
        cn.ci.aj = sum.lpd(sapply(d.ci.aj,FUN=cnij,R0=R0))
        cn.ci.cj = sum.lpd(sapply(d.ci.cj,FUN=cnij,R0=R0))
        cn.li.aj = sum.lpd(sapply(d.li.aj,FUN=cnij,R0=R0))
        cn.li.lj = sum.lpd(sapply(d.li.lj,FUN=cnij,R0=R0))
        cn.mi.aj = sum.lpd(sapply(d.mi.aj,FUN=cnij,R0=R0))
        cn.mi.mj = sum.lpd(sapply(d.mi.mj,FUN=cnij,R0=R0))
        df = rbind(df,data.frame(
                            i,resno.i,aa.i,cn.ai.bi,cn.ai.ci,cn.ai.li,cn.ai.mi,
                            cn.ai.aj,cn.ai.bj,cn.ai.cj,cn.ai.lj,cn.ai.mj,
                            cn.bi.aj,cn.bi.bj,
                            cn.ci.aj,cn.ci.cj,
                            cn.li.aj,cn.li.lj,
                            cn.mi.aj,cn.mi.mj))
    }
    df
}


distance.matrices = function(r) {
# calculates distances for a,b,c,l,and m
# R0 is a cut-off distance used to make calculation faster for very large proteins
    site=r$site
    ca.xyz = r$ca.xyz
    cb.xyz = r$cb.xyz
    com.xyz = r$com.xyz
    m.xyz = r$m.xyz
    l.xyz = r$l.xyz
    nsites = length(r$ca.xyz)/3
    dim(ca.xyz) = c(3,nsites)
    dim(cb.xyz) = c(3,nsites)
    dim(com.xyz) = c(3,nsites)
    dim(m.xyz) = c(3,nsites)
    dim(l.xyz) = c(3,nsites)

    dmat.aa = matrix(NA,nsites,nsites)
    dmat.ca = matrix(NA,nsites,nsites)
    dmat.cc = matrix(NA,nsites,nsites)
    dmat.ba = matrix(NA,nsites,nsites)
    dmat.bb = matrix(NA,nsites,nsites)
    dmat.la = matrix(NA,nsites,nsites)
    dmat.ll = matrix(NA,nsites,nsites)
    dmat.ma = matrix(NA,nsites,nsites)
    dmat.mm = matrix(NA,nsites,nsites)

    for (i in seq(nsites)) {
        #print(c("site",i))
        resno.i = r$site[i]
        aa.i = r$aa[i]
# coordinates of network nodes associated to i
        ca.i = ca.xyz[,i]
        cb.i = cb.xyz[,i]
        com.i = com.xyz[,i]
        m.i = m.xyz[,i]
        l.i = l.xyz[,i]

# distances
        for (j in seq(nsites)) {
            if (i != j) {
# coordinates of network nodes associated to j
                ca.j = ca.xyz[,j]
                com.j = com.xyz[,j]
                cb.j = cb.xyz[,j]
                l.j = l.xyz[,j]
                m.j = m.xyz[,j]
# distances between nodes of i and nodes of j
                dmat.aa[i,j] = vnorm(ca.j-ca.i)
                dmat.ba[i,j] = vnorm(ca.j-cb.i)
                dmat.bb[i,j] = vnorm(cb.j-cb.i)
                dmat.ca[i,j] = vnorm(ca.j-com.i)
                dmat.cc[i,j] = vnorm(com.j-com.i)
                dmat.la[i,j] = vnorm(ca.j-l.i)
                dmat.ll[i,j] = vnorm(l.j-l.i)
                dmat.ma[i,j] = vnorm(ca.j-m.i)
                dmat.mm[i,j] = vnorm(m.j-m.i)

            } # close if
        } # close j loop
    } # close i loop
    result = list("site"=site,"daa"=dmat.aa,
                              "dba"=dmat.ba,"dbb"=dmat.bb,
                              "dca"=dmat.ca,"dcc"=dmat.cc,
                              "dla"=dmat.la,"dll"=dmat.ll,
                              "dma"=dmat.ma,"dmm"=dmat.mm)

    result
}
calc.wmat = function(xyz){
    n.nodes = length(xyz)/3
    dim(xyz)=c(3,n.nodes)
    wmat = matrix(0,n.nodes,n.nodes)
    for (i in seq(n.nodes-1)){
        ri = xyz[,i]
        for (j in c((i+1):(n.nodes))){
            rj = xyz[,j]
            dij = vnorm(rj-ri)
            wmat[i,j] = 1/dij^2
            wmat[j,i] = wmat[i,j]
        }
    }
    wmat
}

