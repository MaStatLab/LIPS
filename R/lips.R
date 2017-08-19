lips = function( y.vec, x.mat, T = 10, kstep = 0, alpha = NA, method = "g", nparticle = 100, rho0=0.5, rho0.method ="const",rho.prop=0.1,p.mix=1,resample.param=1,n.top.vars=5){
  ## p.mix is the mixing probability for boostrap filter
  ## p.mix = 1 --> always boostrap filter
  ## p.mix = 0 --> always k-step filter
  ## p.mix = p in [0,1] --> p proportion boostrap filter, 1 - p proportion k-step filteró
  ## p.mix = 2 --> hybrid: prior stopping + k-step selection
  ## p.mix = 3 --> hybrid: k-step stopping + prior selection
  ## p.mix = 4 --> hybrid: prespecified small stopping probability + k-step selection
  ## p.mix = 5 --> SSS style proposal: selection probabilities for the top MAXVAR variables evenly distributed
  ## p.mix = 6 --> k-step filter (same as 0) implemented using reordering in the selection step

  nobs = length(y.vec)
  p = ncol(x.mat)
  kstep = min(c(kstep,T), na.rm=TRUE)
  T = min(p,T);



  if (method == "g") method.int = as.integer(0)
  if (method == "hyper-g") method.int = as.integer(1)

  if (is.na(alpha)) {
    if (method == "g") alpha = nobs ## unit information g-prior
    if (method == "hyper-g") alpha = 3
  }

  if (rho0.method == "const") rho0.method.int = as.integer(0)
  if (rho0.method == "exp") rho0.method.int = as.integer(1)
  if (rho0.method == "beta-binomial") rho0.method.int = as.integer(2)

  colnames(x.mat) = as.character(1:ncol(x.mat))

  x.mat.centered <- apply(x.mat, 2, scale, scale=FALSE, center=TRUE)
  y.vec.centered = y.vec - mean(y.vec)


  ans = .Call('modelsel_pf_C', PACKAGE='LIPS',y.vec.centered, x.mat.centered, p, as.integer(T), as.integer(kstep), as.numeric(alpha), method.int, as.integer(nparticle), as.numeric(rho0), rho0.method.int, as.numeric(rho.prop), as.numeric(p.mix),as.numeric(resample.param),as.integer(min(n.top.vars,T)))

  names(ans)=c("models","logweights")

  return(ans)
}
