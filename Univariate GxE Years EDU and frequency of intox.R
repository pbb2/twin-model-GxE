# ------------------------------------------------------------------------------
# Program: Univariate G x E Intearction
#  Author: Peter Barr
#    Date: 10 06 2015 

# -------|---------|---------|---------|---------|---------|---------|---------|

setwd("~/Desktop/R/FT12/data/")

# Load Library -------------------------------------------------------------------
require(OpenMx)
require(psych)


#NPSOL, SLQSP, CSOLNP
mxOption( NULL, "Default optimizer", "NPSOL" )

# PREPARE DATA -------------------------------------------------------------------

# Restructuring data on family instead of twins
dat <-read.csv('twin_data_FT12_age22.csv', header = T)
describe(dat)

# Standardizing frequency of intoxication
#dat$freq_intox <-scale(dat$zfreq_intox)

twin1 <-subset(dat, (firstborn==1 & zyg!=6) | (firstborn==0 & zyg ==6))
twin2 <-subset(dat, firstborn==0 & zyg !=6 | zyg==6 & firstborn==1)
fam <-merge(twin1, twin2, by=c('fam_nb','zyg'), suffixes=c('_1','_2'), all=T)
dim(fam);head(fam)
remove(twin1,twin2,dat)

fam$age <- fam$age_1
fam$age_2 <- NULL
twinData <-fam
describe(twinData, skew=F)
remove(fam)


# Select Variables for Analysis
Vars      <- 'freq_intox'
Defs      <- c('yrsed')
nv        <- 1       # number of variables
ntv       <- nv*2    # number of total variables
selVars   <- paste(Vars,c(rep(1,nv),rep(2,nv)),sep="_")
defVars   <- paste(Defs,c(rep(1,nv),rep(2,nv)),sep="_")
covVars   <- 'age'

# Select Data for Analysis
twinData[,'yrsed_1'] <- (twinData[,'yrsed_1'])
twinData[,'yrsed_2'] <- (twinData[,'yrsed_2'])
mzData    <- subset(twinData, zyg==1 | zyg==2, c(selVars,defVars,covVars))
dzData    <- subset(twinData, zyg==3 | zyg==4, c(selVars,defVars,covVars))
describe(twinData)
describe(mzData)

# Eliminate cases with missingness on the definition variables
mzData$missingDef <- ifelse( is.na(mzData$yrsed_1) | is.na(mzData$yrsed_2), 1, 0)                    
mzData      <- subset(mzData, missingDef==0 )
dzData$missingDef <- ifelse( is.na(dzData$yrsed_1) | is.na(dzData$yrsed_2), 1, 0)                    
dzData      <- subset(dzData, missingDef==0 )
#mzData    <- mzData[!is.na(mzData[,selVars)&!is.na(mzData[,defVars),]


# Generate Descriptive Statistics
colMeans(mzData,na.rm=TRUE)
colMeans(dzData,na.rm=TRUE)
cov(mzData,use="complete")
cov(dzData,use="complete")


# Set Starting Values
svMe      <-  3.2                  # start value for means
svPa      <- .6                   # start value for path coefficients (sqrt(variance/#ofpaths))



# ------------------------------------------------------------------------------
# PREPARE MODEL

# ACE Model
# Matrices declared to store a, c, and e Path Coefficients
pathA     <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values= svPa, label="a11", lbound=.001, name="a" ) 
pathC     <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values= svPa, label="c11", lbound=.001, name="c" )
pathE     <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values= svPa, label="e11", lbound=.001, name="e" )



# Matrices generated to hold A, C, and E computed Variance Components
covA      <- mxAlgebra( expression=a %*% t(a), name="A" )
covC      <- mxAlgebra( expression=c %*% t(c), name="C" ) 
covE      <- mxAlgebra( expression=e %*% t(e), name="E" )

# Algebra to compute total variances and standard deviations (diagonal only)
covP      <- mxAlgebra( expression=A+C+E, name="V" )
matI      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I")
invSD     <- mxAlgebra( expression=solve(sqrt(I*V)), name="iSD")

# Algebras generated to create summary Table of Derived Variance Components
rowVars   <- rep('vars',nv)
colVars   <- rep(c('A','C','E','SA','SC','SE'),each=nv)
estVars   <- mxAlgebra( expression=cbind(A,C,E,A/V,C/V,E/V), name="Vars", dimnames=list(rowVars,colVars) )

# Matrix for expected Mean 
meanG     <- mxMatrix( type="Full", nrow=1, ncol=nv, free=TRUE, values= svMe, label="mean", name="Mean" )

# Matrix for moderating/interacting variable
defEdu1    <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.yrsed_1"), name="Edu1")
defEdu2    <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.yrsed_2"), name="Edu2")
defAge    <- mxMatrix( type="Full", nrow=1, ncol=1, free=FALSE, labels=c("data.age"), name="Age" )


# Matrices declared to store moderated a, c, and e Path Coefficients
pathAI    <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.01, label="aI11", name="aI" ) 
pathCI    <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.01, label="cI11", name="cI" )
pathEI    <- mxMatrix( type="Full", nrow=nv, ncol=nv, free=TRUE, values=.01, label="eI11", name="eI" )

# Matrices generated to hold moderated A, C, and E computed Variance Components
covAI1     <- mxAlgebra( expression=(a+ Edu1%*%aI) %*% t(a+ Edu1%*%aI), name="AI1" )
covCI1     <- mxAlgebra( expression=(c+ Edu1%*%cI) %*% t(c+ Edu1%*%cI), name="CI1" )
covEI1     <- mxAlgebra( expression=(e+ Edu1%*%eI) %*% t(e+ Edu1%*%eI), name="EI1" )
covAI2     <- mxAlgebra( expression=(a+ Edu2%*%aI) %*% t(a+ Edu2%*%aI), name="AI2" )
covCI2     <- mxAlgebra( expression=(c+ Edu2%*%cI) %*% t(c+ Edu2%*%cI), name="CI2" )
covEI2     <- mxAlgebra( expression=(e+ Edu2%*%eI) %*% t(e+ Edu2%*%eI), name="EI2" )
covAI12    <- mxAlgebra( expression=(a+ Edu1%*%aI) %*% t(a+ Edu2%*%aI), name="AI12" )
covCI12    <- mxAlgebra( expression=(c+ Edu1%*%cI) %*% t(c+ Edu2%*%cI), name="CI12" )
covAI21    <- mxAlgebra( expression=(a+ Edu2%*%aI) %*% t(a+ Edu1%*%aI), name="AI21" )
covCI21    <- mxAlgebra( expression=(c+ Edu2%*%cI) %*% t(c+ Edu1%*%cI), name="CI21" )

# Algebra to compute total variances and standard deviations (diagonal only)
covPI1      <- mxAlgebra( expression=AI1+CI1+EI1, name="VI1" )
covPI2      <- mxAlgebra( expression=AI2+CI2+EI2, name="VI2" )

# Matrices declared to store linear and quadratic Regression Coefficients for covariate
pathB     <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=.01, label=c("l11"), name="b" )
pathAge     <- mxMatrix( type="Full", nrow=1, ncol=1, free=TRUE, values=.01, label="b11", name="bage" )

# Algebra for expected Mean and Variance/Covariance Matrices in MZ & DZ twins
meanEdu1   <- mxAlgebra( expression= b%*%Edu1, name="EduR1")
meanEdu2   <- mxAlgebra( expression= b%*%Edu2, name="EduR2")
meanAge    <- mxAlgebra( expression= bage%*%Age, name="AgeR1")
meanGI    <- mxAlgebra( expression= cbind((Mean + EduR1 + AgeR1),(Mean + EduR2 + AgeR1)), name="expMean")
covMZI    <- mxAlgebra( expression= rbind( cbind(VI1, AI12+CI12), cbind(AI21+CI21, VI2)), name="expCovMZ" )
covDZI    <- mxAlgebra( expression= rbind( cbind(VI1, 0.5%x%AI12+CI12), cbind(0.5%x%AI21+CI21, VI2)), name="expCovDZ" )

# Data objects for Multiple Groups
dataMZ    <- mxData( observed=mzData, type="raw" )
dataDZ    <- mxData( observed=dzData, type="raw" )

# Objective objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="expMean", dimnames=selVars )
expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="expMean", dimnames=selVars )
funML     <- mxFitFunctionML()

# Matrices & Algebra to plot Means and Variance Components by Wt (not required for model fitting)
Wts       <- mxMatrix( type="Full", nrow=5, ncol=1, values=c(9,11,13,15,17), name="Wts")
unit      <- mxMatrix( type="Unit", nrow=5, ncol=1, name="unit")
meanI     <- mxAlgebra( expression=unit%x%Mean+ t(b%*%t(Wts)), name="Mi")
compAI    <- mxAlgebra( expression=diag2vec((unit%x%a+ Wts%x%aI) %*% t(unit%x%a+ Wts%x%aI)), name="Ai" )
compCI    <- mxAlgebra( expression=diag2vec((unit%x%c+ Wts%x%cI) %*% t(unit%x%c+ Wts%*%cI)), name="Ci" )
compEI    <- mxAlgebra( expression=diag2vec((unit%x%e+ Wts%x%eI) %*% t(unit%x%e+ Wts%*%eI)), name="Ei" )
compPI    <- mxAlgebra( expression=Ai+Ci+Ei, name="Vi" )

# Algebras generated to create summary Table of Derived Variance Components by Wt
rowsWt   <- c('9','11','13','15','17')
colsWt   <- c("MI","AI","CI","EI","VI")
estWt    <- mxAlgebra( expression=cbind(Mi,Ai,Ci,Ei,Vi), name="byWt", dimnames=list(rowsWt,colsWt) )
estPrWt  <- mxAlgebra( expression=cbind(Mi,Ai/Vi,Ci/Vi,Ei/Vi,Vi), name="byPrWt", dimnames=list(rowsWt,colsWt) )
propWt   <- list( Wts, unit, meanI, compAI, compCI, compEI, compPI, estWt, estPrWt)

# Combine Groups
pars      <- list( pathAge, pathA, pathC, pathE, pathAI, pathCI, pathEI, pathB, meanG, covA, covC, covE, covP, matI, invSD, estVars )
defs      <- list( defAge, defEdu1, defEdu2, meanEdu1, meanEdu2, meanAge,
                   covAI1, covCI1, covEI1, covAI2, covCI2, covEI2, covAI12, covCI12, covAI21, covCI21, covPI1, covPI2 )
ciVars    <- mxCI("Vars")
ciEdu     <-mxCI("byPrWt")
modelMZ   <- mxModel( pars, defs, meanGI, covMZI, dataMZ, expMZ, funML, name="MZ" )
modelDZ   <- mxModel( pars, defs, meanGI, covDZI, dataDZ, expDZ, funML, name="DZ" )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )
ModAceLModel  <- mxModel( "ModACEl", pars, modelMZ, modelDZ, multi, propWt, ciVars, ciEdu )

# ------------------------------------------------------------------------------
# RUN MODELS

# Run Moderated ACE model + Linear Moderated Means
ModAceLFit    <- mxRun(ModAceLModel, intervals = T)
ModAceLSum    <- summary(ModAceLFit)
ModAceLSum$pa
round(ModAceLFit@output$estimate,4)
round(ModAceLFit$Vars@result,4)
round(ModAceLFit$byWt@result,4)
mxCompare(ModAceLFit)

# ------------------------------------------------------------------------------
# FIT SUBMODELS

# Run non-Moderated ACE model
OneAceLModel  <- mxModel( ModAceLFit, name="OneACEL" )
OneAceLModel  <- omxSetParameters( OneAceLModel, labels=c("aI11","cI11","eI11"), free=FALSE, values=0 )
OneAceLFit    <- mxRun(OneAceLModel)
round(OneAceLFit@output$estimate,4)
round(OneAceLFit$Vars@result,4)
round(OneAceLFit$byWt@result,4)
mxCompare(ModAceLFit, OneAceLFit)

# Run No Moderation on A model
NoAAceLModel  <- mxModel( ModAceLFit, name="NoAACEL" )
NoAAceLModel  <- omxSetParameters(NoAAceLModel, labels="aI11", free=FALSE, values=0)
NoAAceLFit    <- mxRun(NoAAceLModel)
round(NoAAceLFit@output$estimate,4)
round(NoAAceLFit$Vars@result,4)
round(NoAAceLFit$byWt@result,4)
mxCompare(ModAceLFit, NoAAceLFit)

# Run No Moderation on C model
NoCAceLModel  <- mxModel( ModAceLFit, name="NoCACEL" )
NoCAceLModel  <- omxSetParameters(NoCAceLModel, labels="cI11", free=FALSE, values=0)
NoCAceLFit    <- mxRun(NoCAceLModel)
round(NoCAceLFit@output$estimate,4)
round(NoCAceLFit$Vars@result,4)
round(NoCAceLFit$byWt@result,4)
mxCompare(ModAceLFit, NoCAceLFit)

# Run No Moderation on E model
NoEAceLModel  <- mxModel( ModAceLFit, name="NoEACEL" )
NoEAceLModel  <- omxSetParameters(NoEAceLModel, labels="eI11", free=FALSE, values=0)
NoEAceLFit    <- mxRun(NoEAceLModel)
round(NoEAceLFit@output$estimate,4)
round(NoEAceLFit$Vars@result,4)
round(NoEAceLFit$byWt@result,4)
mxCompare(ModAceLFit, NoEAceLFit)

# Fit Moderated ACE model + no Moderated Means
ModAceModel    <- mxModel( OneAceLFit, name="ModACE" )
ModAceModel    <- omxSetParameters( ModAceModel, labels="l11", free=FALSE, values=0 )
ModAceFit      <- mxRun(ModAceModel)
round(ModAceFit@output$estimate,4)
round(ModAceFit$Vars@result,4)
round(ModAceFit$byWt@result,4)

# ------------------------------------------------------------------------------

# Print Comparative Fit Statistics
AceNested <- list(NoAAceLFit, NoCAceLFit, NoEAceLFit, OneAceLFit, ModAceFit)
mxCompare(ModAceLFit,AceNested)






#-------------------------------------------------------------------------------

# Graphing Variance components from the GxE model


# R graphics
results <-as.data.frame(round(ModAceLFit$byWt@result,4))

# Variance Estimates
par(mfrow=c(1,2), oma = c(0, 0, 2, 0), xpd=NA)

plot(results$AI, xlab="Years Of Education", ylab="Variance Estimates",ylim=c(0,3), axes = FALSE, type = "l",lwd=4 ,col="red")
axis(1,1:5,labels=rowsWt)
axis(2,0:3)
lines(results$CI,lwd=4,col="green")
lines(results$EI,lwd=4,col="blue")
legend(x=1, y=4.25,legend=c("Additive Genetic", "Shared Environment","Unique Environment"),lty=c(1,1,1),lwd=c(3,3,3),col=c("red","green","blue"))
# Changes in mean  
plot(results$MI ,xlab="Years Of Education", ylab="Mean Estimates",ylim=c(2,7), axes = FALSE, type = "l", lwd=4, col="purple" )
axis(1,1:5,labels=rowsWt)
axis(2,2:7)
#
legend(x=3, y=8.25,legend="Mean", lty=1, lwd=3, col="purple")
#
mtext("Moderation of Frequency of Intoxication by Education", outer = TRUE, cex = 1.5)
#




