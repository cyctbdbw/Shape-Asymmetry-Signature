args = commandArgs(trailingOnly=TRUE)

require(OpenMx)   
require(MASS) 
require(R.matlab) 

library(R.matlab)
library(MASS)
library(OpenMx)

set.seed(200)

# load the data
setwd("....") # Change to the path with twin_SAS.mat   
data_covar=readMat('twinCovariates.mat')
data=readMat('twin_SAS.mat') # the file can be changed to individual hemisphere

# remove second additional sibling from the data
data$Output.MZ <- data$Output.MZ[,-4,];
data$Output.DZ <- data$Output.DZ[,-4,];

data_covar$MZ.ID <- data_covar$MZ.ID[,-4];
data_covar$MZ.age <- data_covar$MZ.age[,-4];
data_covar$MZ.sex <- data_covar$MZ.sex[,-4];

data_covar$DZ.ID <- data_covar$DZ.ID[,-4];
data_covar$DZ.age <- data_covar$DZ.age[,-4];
data_covar$DZ.sex <- data_covar$DZ.sex[,-4];

# replace years with centuries (recommended by openMx administrators)
data_covar$DZ.age = data_covar$DZ.age/100
data_covar$MZ.age = data_covar$MZ.age/100

# replace missing covariates with a very different number setting them to "pseudo-missing" values as openMx doesn't allow missing covariate values (recommended by openMx administrators)
data_covar$MZ.ID[is.na(data_covar$MZ.ID)] <- -999
data_covar$MZ.age[is.na(data_covar$MZ.age)] <- -999
data_covar$MZ.sex[is.na(data_covar$MZ.sex)] <- -999

data_covar$DZ.ID[is.na(data_covar$DZ.ID)] <- -999
data_covar$DZ.age[is.na(data_covar$DZ.age)] <- -999
data_covar$DZ.sex[is.na(data_covar$DZ.sex)] <- -999


# all edges
numEdges = dim(data$Output.DZ)[3]
heritabilityA <- numeric(numEdges)
heritabilityC <- numeric(numEdges)
heritabilityT <- numeric(numEdges)
heritabilityE <- numeric(numEdges)
heritabilityS <- numeric(numEdges)
heritabilitySp <- numeric(numEdges)
best_model <- numeric(numEdges)
pValue <- numeric(numEdges)
pValueModelC <- numeric(numEdges)
saturated <- numeric(numEdges)

# (recommended by openMx administrators)
mxOption(NULL,"Default optimizer","SLSQP")

for (edge in c(1:numEdges)){
  #for (edge in c(1:20)){
  # each connection (phenotype of interest) is denoted as edge
  myDataMZ<-data.frame(data$Output.MZ[,,edge], data_covar$MZ.age[,1], data_covar$MZ.age[,2],data_covar$MZ.age[,3],data_covar$MZ.sex[,1],data_covar$MZ.sex[,2],data_covar$MZ.sex[,3])
  myDataDZ<-data.frame(data$Output.DZ[,,edge], data_covar$DZ.age[,1], data_covar$DZ.age[,2],data_covar$DZ.age[,3],data_covar$DZ.sex[,1],data_covar$DZ.sex[,2],data_covar$DZ.sex[,3])
  
  myDataMZ_measure<-data$Output.MZ[,,edge]
  myDataDZ_measure<-data$Output.DZ[,,edge]
  # replaces outlier with NA
  outliersValue<- boxplot.stats(myDataMZ_measure)$out
  myDataMZ_measure[myDataMZ_measure %in% outliersValue]=NA;
  
  outliersValue<- boxplot.stats(myDataDZ_measure)$out
  myDataDZ_measure[myDataDZ_measure %in% outliersValue]=NA;
  
  colnames(myDataMZ) <- c('twin1', 'twin2', 'sib','ageT1MZ', 'ageT2MZ', 'ageSIBMZ', 'sexT1MZ', 'sexT2MZ', 'sexSIBMZ')
  colnames(myDataDZ) <- c('twin1', 'twin2', 'sib','ageT1DZ', 'ageT2DZ', 'ageSIBDZ', 'sexT1DZ', 'sexT2DZ', 'sexSIBDZ')
  selVars <- c('twin1','twin2', 'sib')
  
  # "complete" specifies use only cases with data in all columns, otherwise subjects would be excluded if they're missing an additional sibling.
  CovMZ = cov(data$Output.MZ[,,edge],use="complete")
  CovDZ = cov(data$Output.DZ[,,edge],use="complete")
  
  # mean across MZ twins
  MeanMZ = colMeans(data$Output.MZ[,1:2,edge],na.rm=TRUE)
  MeanMZ = mean(MeanMZ)
  # mean across DZ twins
  MeanDZ = colMeans(data$Output.DZ[,1:2,edge],na.rm=TRUE)
  MeanDZ = mean(MeanDZ)
  # mean across siblings
  MeanSIBMZ = mean(data$Output.MZ[,3,edge],na.rm=TRUE)
  MeanSIBDZ = mean(data$Output.DZ[,3,edge],na.rm=TRUE)
  
  
  mylabels=c("twin1","twin2","sib")
  MZsat <- mxModel("MZsat",
                   
                   mxMatrix( type = "Full", nrow=1, ncol=3, free=T, c(MeanMZ,MeanMZ,MeanSIBMZ), labels =c("b0_mz1","b0_mz2","b0_mzsib"), name="Intercepts" ),
                   mxMatrix( type = "Full", nrow=1, ncol=2, free=TRUE, values= 0, labels=c("betaAge","betaSex"), name="beta"),
                   mxMatrix( type = "Full", nrow=2, ncol=3, free=F, labels=c("data.ageT1MZ","data.sexT1MZ","data.ageT2MZ","data.sexT2MZ","data.ageSIBMZ","data.sexSIBMZ"), name="MZDefVars"),
                   mxAlgebra( expression=Intercepts + beta %*% MZDefVars, name="expMeanMZ"),
                   mxMatrix( type = "Lower", nrow=3, ncol=3, free=T, 0.5, name="CholMZ" ),
                   mxAlgebra( CholMZ %*% t(CholMZ), name="expCovMZ"),
                   mxData( myDataMZ, type="raw"), mxFitFunctionML(),
                   mxExpectationNormal( "expCovMZ", "expMeanMZ", mylabels))
  
  DZsat <- mxModel("DZsat",
                   mxMatrix( type = "Full", nrow=1, ncol=3, free=T, c(MeanDZ,MeanDZ,MeanSIBDZ), labels=c("b0_dz1","b0_dz2","b0_dzsib"), name="Intercepts" ),
                   mxMatrix( type = "Full", nrow=1, ncol=2, free=TRUE, values= 0, labels=c("betaAge","betaSex"), name="beta"),
                   mxMatrix( type = "Full", nrow=2, ncol=3, free=F, labels=c("data.ageT1DZ","data.sexT1DZ","data.ageT2DZ","data.sexT2DZ","data.ageSIBDZ","data.sexSIBDZ"), name="DZDefVars"),
                   mxAlgebra( expression=Intercepts + beta %*% DZDefVars, name="expMeanDZ"),
                   mxMatrix( type = "Lower", nrow=3, ncol=3, free=T, 0.5, name="CholDZ" ),
                   mxAlgebra( CholDZ %*% t(CholDZ),name="expCovDZ"),
                   mxData( myDataDZ, type="raw"), mxFitFunctionML(),
                   mxExpectationNormal( "expCovDZ", "expMeanDZ", mylabels))
  
  # Model specification starts here
  SatModel <- mxModel("twinSat", MZsat, DZsat, mxFitFunctionMultigroup(c('MZsat', 'DZsat')))
  #---------------------------------------------------------------------------------------------------------------------------
  SatModelFit <- mxTryHard(SatModel)
  summary(SatModelFit)
  
  # -----------------------------------------------------------------------
  #Fit ACE Model with RawData and Matrices Input
  # -----------------------------------------------------------------------
  twinACTE <- mxModel("twinACTE",
                      # Matrices X, Y, and Z to store a, c, and e path coefficients
                      mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=sqrt((CovDZ[2,2]/3)), label="a", name="X" ),
                      mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=sqrt((CovDZ[2,2]/3)), label="c", name="Y" ),
                      mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=sqrt((CovDZ[2,2]/3)), label="e", name="Z" ),
                      mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=sqrt((CovDZ[2,2]/3)), label="t", name="W" ),
                      # Matrices A, C, and E compute variance components
                      mxAlgebra( expression=X %*% t(X), name="A" ),
                      mxAlgebra( expression=Y %*% t(Y), name="C" ),
                      mxAlgebra( expression=Z %*% t(Z), name="E" ),
                      mxAlgebra( expression=W %*% t(W), name="T" ),
                      
                      # Algebra for expected variance/covariance matrix in MZ+sibling
                      mxAlgebra(
                        expression= rbind  (cbind(A+C+T+E , A+C+T, 0.5%x%A+C),
                                            cbind(A+C+T   , A+C+T+E, 0.5%x%A+C),
                                            cbind(0.5%x%A+C, 0.5%x%A+C, A+C+E)),
                        name="expCovMZ"),
                      
                      # Algebra for expected variance/covariance matrix in DZ+siblig
                      mxAlgebra(
                        expression= rbind  (cbind(A+C+T+E , 0.5%x%A+C+T, 0.5%x%A+C),
                                            cbind(0.5%x%A+C+T , A+C+T+E, 0.5%x%A+C),
                                            cbind(0.5%x%A+C , 0.5%x%A+C, A+C+E)),
                        name="expCovDZ"),
                      
                      mxModel("MZ", mxData( observed=myDataMZ, type="raw" ),
                              # Algebra for making the means a function of the definition variables age and sex
                              mxMatrix( type="Full", nrow=1, ncol=3, free=T, c(MeanMZ,MeanMZ,MeanSIBMZ), labels =c("b0_mz1","b0_mz2","b0_mzsib"), name="Intercepts"),
                              mxMatrix( type="Full", nrow=1, ncol=2, free=T, values= 0, labels=c("betaAge","betaSex"), name="beta"),
                              mxMatrix( type="Full", nrow=2, ncol=3, free=F, labels=c("data.ageT1MZ","data.sexT1MZ","data.ageT2MZ","data.sexT2MZ","data.ageSIBMZ","data.sexSIBMZ"), name="MZDefVars"),
                              mxAlgebra( expression=Intercepts + beta %*% MZDefVars, name="expMeanMZ"),
                              mxExpectationNormal( covariance="twinACTE.expCovMZ", means="expMeanMZ", dimnames=selVars ),
                              mxFitFunctionML()),
                      
                      mxModel("DZ", mxData( observed=myDataDZ, type="raw" ),
                              mxMatrix( type="Full", nrow=1, ncol=3, free=T, c(MeanDZ,MeanDZ,MeanSIBDZ), labels=c("b0_dz1","b0_dz2","b0_dzsib"), name="Intercepts"),
                              mxMatrix( type="Full", nrow=1, ncol=2, free=T, values= 0, labels=c("betaAge","betaSex"), name="beta"),
                              mxMatrix( type="Full", nrow=2, ncol=3, free=F, labels=c("data.ageT1DZ","data.sexT1DZ","data.ageT2DZ","data.sexT2DZ","data.ageSIBDZ","data.sexSIBDZ"), name="DZDefVars"),
                              mxAlgebra( expression=Intercepts + beta %*% DZDefVars, name="expMeanDZ"),
                              mxExpectationNormal( covariance="twinACTE.expCovDZ", means="expMeanDZ", dimnames=selVars ),
                              mxFitFunctionML()),
                      
                      mxFitFunctionMultigroup( c("MZ.fitfunction",  "DZ.fitfunction"))
  )
  twinACTEFit<-mxTryHard(twinACTE)
  
  #Run ACE model
  # -----------------------------------------------------------------------
  #estCovMZ  <- mxEval(twinACE.expCovMZ, twinACEFit)      # expected covariance matrix for MZ's
  #estCovDZ  <- mxEval(twinACE.expCovDZ, twinACEFit)      # expected covariance matrix for DZ's
  estVA1     <- mxEval(a*a, twinACTEFit)              # additive genetic variance, a^2
  estVC1     <- mxEval(c*c, twinACTEFit)              # shared enviromnemtal variance, c^2
  estVT1     <- mxEval(t*t, twinACTEFit)
  estVE1     <- mxEval(e*e, twinACTEFit)              # unique environmental variance, e^2
  estVP1     <- (estVA1+estVC1+estVT1+estVE1)                  # total variance
  
  heritmodels <- matrix(nrow = 8, ncol = 5)
  
  heritmodels[1,1] <- estVA1/estVP1                          # standardized additive genetic variance
  heritmodels[1,2] <- estVC1/estVP1
  heritmodels[1,3] <- estVT1/estVP1
  heritmodels[1,4] <- estVE1/estVP1
  heritmodels[1,5] <- twinACTEFit@output$status$code
  
  # Generate ACE Model - T=0
  twinACE      <- twinACTE
  twinACE      <- mxRename(twinACE, "twinACE")
  twinACE      <- omxSetParameters(twinACE, labels="t", free=FALSE, values=0 )
  twinACEFit   <- mxTryHard(twinACE)
  ACESumm      <- summary(twinACEFit)
  
  estVA2    <- mxEval(a*a, twinACEFit)             
  estVC2    <- mxEval(c*c, twinACEFit)    
  estVT2    <- mxEval(t*t, twinACEFit) 
  estVE2    <- mxEval(e*e, twinACEFit)             
  estVP2    <- (estVA2+estVC2+estVT2+estVE2)                
  
  heritmodels[2,1] <- estVA2/estVP2                          # standardized additive genetic variance
  heritmodels[2,2] <- estVC2/estVP2
  heritmodels[2,3] <- estVT2/estVP2
  heritmodels[2,4] <- estVE2/estVP2
  heritmodels[2,5] <- twinACEFit@output$status$code
  
  # Generate AE Model - C=0
  twinAE      <- twinACE
  twinAE      <- mxRename(twinAE, "twinAE")
  twinAE      <- omxSetParameters(twinAE, labels="c", free=FALSE, values=0 )
  twinAEFit   <- mxTryHard(twinAE)
  AESumm      <- summary(twinAEFit)
  
  estVA3    <- mxEval(a*a, twinAEFit)              # additive genetic variance, a^2
  estVC3    <- mxEval(c*c, twinAEFit)              # shared enviromnemtal variance, c^2
  estVT3    <- mxEval(t*t, twinAEFit)
  estVE3    <- mxEval(e*e, twinAEFit)              # unique environmental variance, e^2
  estVP3    <- (estVA3+estVC3+estVT3+estVE3)                  # total variance
  
  heritmodels[3,1] <- estVA3/estVP3                          # standardized additive genetic variance
  heritmodels[3,2] <- estVC3/estVP3
  heritmodels[3,3] <- estVT3/estVP3
  heritmodels[3,4] <- estVE3/estVP3
  heritmodels[3,5] <- twinAEFit@output$status$code
  
  #Generate CE Model - A=0
  twinCE      <- twinACE
  twinCE      <- mxRename(twinCE, "twinCE")
  twinCE      <- omxSetParameters(twinCE, labels="a", free=FALSE, values=0 )
  twinCEFit   <- mxTryHard(twinCE)
  
  estVA4    <- mxEval(a*a, twinCEFit)              # additive genetic variance, a^2
  estVC4    <- mxEval(c*c, twinCEFit)              # shared enviromnemtal variance, c^2
  estVT4    <- mxEval(t*t, twinCEFit) 
  estVE4    <- mxEval(e*e, twinCEFit)              # unique environmental variance, e^2
  estVP4    <- (estVA4+estVC4+estVT4+estVE4)              # total variance
  
  heritmodels[4,1] <- estVA4/estVP4                          # standardized additive genetic variance
  heritmodels[4,2] <- estVC4/estVP4
  heritmodels[4,3] <- estVT4/estVP4
  heritmodels[4,4] <- estVE4/estVP4
  heritmodels[4,5] <- twinCEFit@output$status$code
  
  #Generate E Model - A & C=0
  twinE      <- twinCE
  twinE      <- mxRename(twinE, "twinE")
  twinE      <- omxSetParameters(twinE, labels="c", free=FALSE, values=0 )
  twinEFit   <- mxTryHard(twinE)
  
  estVA5    <- mxEval(a*a, twinEFit)              # additive genetic variance, a^2
  estVC5    <- mxEval(c*c, twinEFit)   
  estVT5    <- mxEval(t*t, twinEFit)# shared enviromnemtal variance, c^2
  estVE5    <- mxEval(e*e, twinEFit)              # unique environmental variance, e^2
  estVP5    <- (estVA5+estVC5+estVT5+estVE5)              # total variance
  
  heritmodels[5,1] <- estVA5/estVP5                          # standardized additive genetic variance
  heritmodels[5,2] <- estVC5/estVP5
  heritmodels[5,3] <- estVT5/estVP5
  heritmodels[5,4] <- estVE5/estVP5
  heritmodels[5,5] <- twinEFit@output$status$code
  
  # Generate CTE Model - A=0
  twinCTE      <- twinACTE
  twinCTE      <- mxRename(twinCTE, "twinCTE")
  twinCTE      <- omxSetParameters(twinCTE, labels="a", free=FALSE, values=0 )
  twinCTEFit   <- mxTryHard(twinCTE)
  
  estVA6    <- mxEval(a*a, twinCTEFit)             
  estVC6    <- mxEval(c*c, twinCTEFit)    
  estVT6    <- mxEval(t*t, twinCTEFit) 
  estVE6    <- mxEval(e*e, twinCTEFit)             
  estVP6    <- (estVA6+estVC6+estVT6+estVE6)                
  
  heritmodels[6,1] <- estVA6/estVP6                          # standardized additive genetic variance
  heritmodels[6,2] <- estVC6/estVP6
  heritmodels[6,3] <- estVT6/estVP6
  heritmodels[6,4] <- estVE6/estVP6
  heritmodels[6,5] <- twinCTEFit@output$status$code
  
  # Generate ATE Model - C=0
  twinATE      <- twinACTE
  twinATE      <- mxRename(twinATE, "twinATE")
  twinATE      <- omxSetParameters(twinATE, labels="c", free=FALSE, values=0 )
  twinATEFit   <- mxTryHard(twinATE)
  
  estVA7    <- mxEval(a*a, twinATEFit)             
  estVC7    <- mxEval(c*c, twinATEFit)    
  estVT7    <- mxEval(t*t, twinATEFit) 
  estVE7    <- mxEval(e*e, twinATEFit)             
  estVP7    <- (estVA7+estVC7+estVT7+estVE7)                
  
  heritmodels[7,1] <- estVA7/estVP7                          # standardized additive genetic variance
  heritmodels[7,2] <- estVC7/estVP7
  heritmodels[7,3] <- estVT7/estVP7
  heritmodels[7,4] <- estVE7/estVP7
  heritmodels[7,5] <- twinATEFit@output$status$code
  
  # Generate TE Model - A=0
  twinTE      <- twinATE
  twinTE      <- mxRename(twinTE, "twinTE")
  twinTE      <- omxSetParameters(twinTE, labels="a", free=FALSE, values=0 )
  twinTEFit   <- mxTryHard(twinTE)
  
  estVA8    <- mxEval(a*a, twinTEFit)             
  estVC8    <- mxEval(c*c, twinTEFit)    
  estVT8    <- mxEval(t*t, twinTEFit) 
  estVE8    <- mxEval(e*e, twinTEFit)             
  estVP8    <- (estVA8+estVC8+estVT8+estVE8)                
  
  heritmodels[8,1] <- estVA8/estVP8                          # standardized additive genetic variance
  heritmodels[8,2] <- estVC8/estVP8
  heritmodels[8,3] <- estVT8/estVP8
  heritmodels[8,4] <- estVE8/estVP8
  heritmodels[8,5] <- twinTEFit@output$status$code
  
  # model comparison
  options('digits' = 8)
  # compare saturated to other models
  compValuesACE = mxCompare(SatModelFit, c(twinACTEFit, twinACEFit,twinAEFit,twinCEFit,twinEFit,twinCTEFit,twinATEFit,twinTEFit))
  #compValuesACE = mxCompare(twinACEFit, c(twinAEFit,twinCEFit,twinEFit))
  
  # find model with the lowest AIC
  INDmin = which.min(compValuesACE$AIC[1:9])
  AICmin = compValuesACE$AIC[INDmin]
  
  # check if this model is not significantly different from saturated
  pmin = compValuesACE$p[INDmin]
  
  if (INDmin==1) {
    INDmin2 = which.min(compValuesACE$AIC[2:9])
    heritabilityA[edge] <- heritmodels[INDmin2,1]
    heritabilityC[edge] <- heritmodels[INDmin2,2]
    heritabilityT[edge] <- heritmodels[INDmin2,3]
    heritabilityE[edge] <- heritmodels[INDmin2,4]
    heritabilityS[edge] <- 1; # 1 if saturated model gave lower AIC value
    heritabilitySp[edge] <- compValuesACE$p[INDmin2+1] # selecting from a table where first row is excluded
    saturated[edge] <- 1  
  } else {
    heritabilityA[edge] <- heritmodels[INDmin-1,1]
    heritabilityC[edge] <- heritmodels[INDmin-1,2]
    heritabilityT[edge] <- heritmodels[INDmin-1,3]
    heritabilityE[edge] <- heritmodels[INDmin-1,4]
    heritabilityS[edge] <- heritmodels[INDmin-1,5]
    heritabilitySp[edge] <- compValuesACE$p[INDmin] # selecting from a table where values in first row are NaNs
    saturated[edge] <- 0
  }
  
  INDmin3 = which.min(compValuesACE$AIC[2:9])
  best_model[edge] <- INDmin3
  if (INDmin3==1) {
    compforP_A = mxCompare(twinACTEFit, twinCTEFit)
    pValue[edge] <- compforP_A$p[2]
    compforP_C = mxCompare(twinACTEFit, twinATEFit)
    pValueModelC[edge] <- compforP_C$p[2]
  } else if (INDmin3==2) {
    compforP_A = mxCompare(twinACEFit, twinCEFit)
    pValue[edge] <- compforP_A$p[2]
    compforP_C = mxCompare(twinACEFit, twinAEFit)
    pValueModelC[edge] <- compforP_C$p[2]
  } else if (INDmin3==3) {
    compforP_A = mxCompare(twinAEFit, twinEFit)
    pValue[edge] <- compforP_A$p[2]
    compforP_C = mxCompare(twinACEFit, twinAEFit)
    pValueModelC[edge] <- compforP_C$p[2]
  } else if (INDmin3==4) {
    compforP_A = mxCompare(twinACEFit, twinCEFit)
    pValue[edge] <- compforP_A$p[2]
    compforP_C = mxCompare(twinCEFit, twinEFit)
    pValueModelC[edge] <- compforP_C$p[2]
  } else if (INDmin3==5) {
    compforP_A = mxCompare(twinAEFit, twinEFit)
    pValue[edge] <- compforP_A$p[2]
    compforP_C = mxCompare(twinCEFit, twinEFit)
    pValueModelC[edge] <- compforP_C$p[2]
  } else if (INDmin3==6) {
    compforP_A = mxCompare(twinACTEFit, twinCTEFit)
    pValue[edge] <- compforP_A$p[2]
    compforP_C = mxCompare(twinCTEFit, twinTEFit)
    pValueModelC[edge] <- compforP_C$p[2]
  } else if (INDmin3==7) {
    compforP_A = mxCompare(twinATEFit, twinTEFit)
    pValue[edge] <- compforP_A$p[2]
    compforP_C = mxCompare(twinACTEFit, twinATEFit)
    pValueModelC[edge] <- compforP_C$p[2]
  } else if (INDmin3==8) {
    compforP_A = mxCompare(twinATEFit, twinTEFit)
    pValue[edge] <- compforP_A$p[2]
    compforP_C = mxCompare(twinCTEFit, twinTEFit)
    pValueModelC[edge] <- compforP_C$p[2]
  }
  
   
  #if (pmin>0.05) {
  # heritabilityA[edge] <- heritmodels[INDmin,1]
  # heritabilityC[edge] <- heritmodels[INDmin,2]
  #  heritabilityT[edge] <- heritmodels[INDmin,3]
  #  heritabilityE[edge] <- heritmodels[INDmin,4]
  #  heritabilityS[edge] <- heritmodels[INDmin,5]
  #} else {
  #  heritabilityA[edge] <- NA; 
  #  heritabilityC[edge] <- NA; 
  #  heritabilityT[edge] <- NA;
  #  heritabilityE[edge] <- NA; 
  #  heritabilityS[edge] <- NA; 
  #}
}

heritabilityACE <- data.frame(heritabilityA,heritabilityC,heritabilityT,heritabilityE,heritabilityS,heritabilitySp,best_model,pValue,pValueModelC,saturated)
setwd("...")# Change to the output path
write.csv(heritabilityACE,sprintf("heritability_SAS%s.txt",args[2]),row.names=FALSE)
