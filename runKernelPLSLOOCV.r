require(AtmRay)
require(foreach)
require(doSNOW)
source('kernelPLS.r') # this loads convID

runKernelPLSLOOCV = function (Y, corMats) {
    ### Y: ados score
    ### corMats: triangularized subject correlation matrices
    
    dropSubs = which(is.na(Y))
    Y = as.matrix(Y[-dropSubs])
    scaleY = FALSE
    Y = scale(Y, center=FALSE, scale=scaleY)
    corMats = corMats[-dropSubs,]
    corMats = scale(corMats, center=FALSE, scale=scaleY)
    corMatKernel = buildBaselineKernel(corMats)
    bsLOOCVOut = kernelPLS.loocv (corMatKernel, Y)
    
    # do crossvalidation over kernel sizes 0,1, weights 0,1,baseline
    weight0 = seq(0, 1, by=0.05)
    weight1 = seq(0, 1, by=0.05)
    weightmg = meshgrid(weight0, weight1)
    weight0 = as.numeric(weightmg$x)
    weight1 = as.numeric(weightmg$y)
    weightbaseline = 1 - (weight0 + weight1)
    weights = cbind(weight0, weight1, weightbaseline)
    weights = weights[(weights[,3]>=0), ] # restrict the weights to sum to 1
    # weights = weights[-1,] # remove the baseline row
    kernelDir = "newKernels/"
    sigSize = (seq(-8, 6, by=0.2))
    sigSizemg = meshgrid(sigSize, sigSize)
    sigSize0 = as.numeric(sigSizemg$x)
    sigSize1 = as.numeric(sigSizemg$y)
    sigSizes = cbind(sigSize0, sigSize1)
    print(bsLOOCVOut)
    
    baselineR2 = NULL
    baselinePredCor = NULL
    cl = makeCluster(1)
    registerDoSNOW(cl)
    foreach(i=1:nrow(sigSizes)) %dopar% {
      loocvOut.i = NULL
      mat0 = as.matrix(read.csv(paste0(kernelDir,'dim_0_1e', round(sigSizes[i,1]/0.1)*0.1, '.csv'), header=FALSE))
      mat1 = as.matrix(read.csv(paste0(kernelDir,'dim_1_1e', round(sigSizes[i,2]/0.1)*0.1, '.csv'), header=FALSE))
      
      mat0Med = median(abs(mat0))
      mat1Med = median(abs(mat1))
      corMatMed = median(abs(corMatKernel))
      mat0 = mat0 / mat0Med # normalizing the scale of the values by its median
      mat1 = mat1 / mat1Med
      corMatKernel = corMatKernel / corMatMed
      for (j in 1:nrow(weights)) {
        w0 = weights[j,1]
        w1 = weights[j,2]
        wbs = weights[j,3]
        w0w1 = w0*mat0 + w1*mat1
        w0w1 = w0w1[-dropSubs, -dropSubs]
        kernCombo = w0w1 + wbs*corMatKernel
        loocvOut = kernelPLS.loocv (kernCombo, Y)
        loocvOut.i = rbind(loocvOut.i, c(sigSizes[i,],weights[j,],loocvOut))
      }
      
      # print table of results
      write.table(loocvOut.i, file=paste0('Results/ss_',i,'.csv'),sep=',', row.names=FALSE, col.names=FALSE)
      print(i)
    }
    stopCluster(cl)
}
