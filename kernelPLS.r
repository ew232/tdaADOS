
kernelPLS = function (kernelG, responseY) {
  # responseY has to already be scaled to mean 0
  output=NULL
  nSubjects = nrow(kernelG)
  u = matrix(rnorm(nSubjects))
  for (iter in 1:5) {
    t = kernelG%*%u
    t = t/norm(t, 'f')
    c = t(responseY)%*%t
    u = responseY%*%c
    u = u/norm(u, 'f')
  }
  B = u%*%solve(t(t)%*%kernelG%*%u)%*%t(t)%*%responseY
  output = list(t=t, c=c, u=u, B=B)
  output
}

kernelPLS.loocv = function (kern, Y) {
  # run kernel PLS over leave one out cross validation
  require(pracma)
  output = NULL
  for (i in 1:nrow(kern)) {
    kern.temp = kern[-i,-i,drop=F] # get the column mean without the testing row
    cMu = t(colMeans(kern.temp))
    kern.temp2 = kern[,-i,drop=F] # include the testing row back in
    kern.temp2 = kern.temp2 - repmat(cMu, nrow(kern.temp2), 1)
    rMu = matrix(rowMeans(kern.temp2))
    kern.temp2 = kern.temp2 - repmat(rMu, 1, ncol(kern.temp2))
    
    Ymean = mean(Y[-i])
    Y.centered = Y-Ymean
    kern.train = kern.temp2[-i,,drop=FALSE]
    kern.test = kern.temp2[i,,drop=FALSE]
    kplsOut = kernelPLS(kern.train, Y.centered[-i])
    fit.train = kern.train%*%kplsOut$B
    fit.test = kern.test%*%kplsOut$B
    output[i] = fit.test + Ymean
  }
  output
}

kernelPLS.test = function (test.idx, kern, Y) {
  # testing new data
  require(pracma)
  kern.temp = kern[-test.idx, -test.idx, drop=F]
  cMu = t(colMeans(kern.temp))
  kern.temp2 = kern[,-test.idx, drop=F]
  kern.temp2 = kern.temp2 - repmat(cMu, nrow(kern.temp2), 1)
  rMu = matrix(rowMeans(kern.temp2))
  kern.temp2 = kern.temp2 - repmat(rMu,1,ncol(kern.temp2))
  Ymean = mean(Y[-test.idx])
  Y.centered = Y-Ymean
  kern.train = kern.temp2[-test.idx,, drop=FALSE]
  kern.test = kern.temp2[test.idx,, drop=FALSE]
  kplsOut = kernelPLS(kern.train, Y.centered[-test.idx])
  fit.train = kern.train%*%kplsOut$B
  fit.test = kern.test%*%kplsOut$B
  output = fit.test + Ymean
  output
}

permTest = function (Y, kern) {
  # runs permutation test on kernel PLS
  output = NULL
  kern = centerKernel(kern)
  Y = scale(Y, center=TRUE, scale=F)
  as.numeric(Sys.time())-> t; set.seed((t - floor(t)) * 1e8 -> seed);
  kernCovPermTest = NULL
  nIters = 1000
  for (iter in 1:nIters) {
    Y.temp = Y
    if (iter == 1) {	  
      kernelPLSOut = kernelPLS(kern, Y.temp)  
      kernCovPermTest[1] = cov(kernelPLSOut$t, kernelPLSOut$u)
    } else {
      shuff = sample(nrow(Y.temp))
      Y.temp = Y.temp[shuff,]
      kernelPLSOut = kernelPLS(kern, Y.temp) 
      kernCovPermTest[iter] = cov(kernelPLSOut$t, kernelPLSOut$u)
    }
  }
  trueCov = kernCovPermTest[1]
  permPval = sum(kernCovPermTest>kernCovPermTest[1])/nIters
  output = c(trueCov, permPval)
  output
}

buildBaselineKernel = function (corMats) {
  # linear kernel of correlation matrices
  output = matrix(0, nrow(corMats), nrow(corMats))
  for (i in 1:nrow(corMats)) {
    c.i = corMats[i,]
    for (j in i:nrow(corMats)) {
      c.j = corMats[j,]
      output[i,j] = output[j,i] = c.i%*%c.j
    }
  }
  output
}

centerKernel = function (K) {
  # centers a kernel
  ell = nrow(K);
  D = colSums(K) / ell;
  E = sum(D) / ell;
  J = matrix(rep(1,ell)) %*% t(D);
  K = K - J - t(J) + E * matrix(1, ell, ell);
  output = K
  output
}

centerY = function (Y) {
	Y = scale(Y, center=T, scale=F)
	Y
}

