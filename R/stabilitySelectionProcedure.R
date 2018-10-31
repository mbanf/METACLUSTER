stabilitySelectionProcedure <- function(x,y,nbootstrap=500,nsteps=5,alpha=0.2,plotme=FALSE){
  
  dimx <- dim(x)
  n <- dimx[1]
  p <- dimx[2]
  halfsize <- as.integer(n/2)
  freq <- matrix(0,nsteps+1,p)
  
  for (i in seq(nbootstrap)) {
    
    # Randomly reweight each variable
    xs <- t(t(x)*runif(p,alpha,1))
    
    # Ramdomly split the sample in two sets
    perm <- sample(dimx[1])
    i1 <- perm[1:halfsize]
    i2 <- perm[(halfsize+1):n]
    
    # run the randomized lasso on each sample and check which variables are selected
    r <- lars(xs[i1,],y[i1], max.steps=nsteps,normalize=FALSE, use.Gram=FALSE) #, eps = .Machine$double.eps)
    freq <- freq + abs(sign(coef.lars(r)))
    r <- lars(xs[i2,],y[i2],max.steps=nsteps,normalize=FALSE, use.Gram=FALSE)
    freq <- freq + abs(sign(coef.lars(r)))
    
  }
  
  # normalize frequence in [0,1]
  freq <- freq/(2*nbootstrap)
  
  if (plotme) {
    matplot(freq,type='l',xlab="LARS iteration",ylab="Frequency")
  }
  
  # the final stability score is the maximum frequency over the steps
  result <- apply(freq,2,max)
}