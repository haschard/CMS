
###############################################################################################
#																							  #
# CMS - 2017.04.19                                                                            #
# Hugues Aschard                                                                              #
#																							  #
# This software is free for academic use. 													  #
# Please contact MTAnalytics (david.skurnik@mt-analytics.com) if you are interested in using  #
# the software for commercial purposes. This software must not be distributed or modified 	  #
# without prior permission of the author.											 		  #
#																							  #
###############################################################################################
 

MCcomp <- function(DATAMAT,myOUTCOME,myPREDICTOR,myFIXCOV,listCOVARIATES,myTEST,optionStand,optionImpute,minTotSample,sigmaMax,dweight,Tmul,verbose){

		version 			= 1.0
		# INPUT ========================================
		# DATAMAT	     : matrix data
		# myOUTCOME	     : outcome label
		# myPredictor    : predictor label
		# myFIXCOV	     : table of fixed covariates names (e.g. confounding factors)
		# listCOVARIATES : table of candidate covariates names (that would be assessed by CMS)
		# myTEST         : NA (for future dev)
		# optionStand    : 1 to standardize of all variable (necessary if variable are not standardized already), 0 otherwise
		# optionImpute   : 1 to impute missing values in candidate covariate, o otherwise
		# verbose        : 1 to display additional output, 0 otherwise

		#------parameters----------------------------#
		# minTotSample   : minimum sample size with complete outcome and predictor data
		# sigmaMax       : recommended value = 2
		# dweight        : recommended value = 8
		# Tmul           : recommended value = 0.05
		
		minCovSample		= 50
		minSmallCat			= 1
	
		#------result tables-------------------------#
		resMA						= matrix(NA,1,3)
		resMC						= matrix(NA,1,8)
		names(resMA)		= c("beta.MA","sd.MA","pval.MA")
		names(resMC)		= c("beta.MC","sd.MC","pval.MC","Yrsq","Ncov","Lcov","pMUL","PcutoffBias")

		#------check covariate missingness-----------#
		reduceSample		<- complete.cases(DATAMAT[,c(myPREDICTOR,myOUTCOME)])
		
		if(verbose == 1){ writeLines(paste("reduceSample = ",sum(reduceSample)))} 
		
		if(sum(reduceSample) > minTotSample){
		DATAMAT					<- DATAMAT[reduceSample,]
		missingCOV			<- apply(as.matrix(DATAMAT[,listCOVARIATES]),2,function(x) sum(!is.na(x)))
		varianceCOV			<- apply(as.matrix(DATAMAT[,listCOVARIATES]),2,function(x) var(x, na.rm=TRUE))
		listCOVARIATES	<- listCOVARIATES[missingCOV > minCovSample & varianceCOV > 0]
		missingFIX			<- apply(as.matrix(DATAMAT[,myFIXCOV]),2,function(x) sum(!is.na(x)))
		varianceFIX			<- apply(as.matrix(DATAMAT[,myFIXCOV]),2,function(x) var(x, na.rm=TRUE))
		myFIXCOV				<- myFIXCOV[missingFIX > minCovSample & varianceFIX > 0]
		VARpred					<- NROW(unique(sort(rank(DATAMAT[,myPREDICTOR]))))
		VARout					<- NROW(unique(sort(rank(DATAMAT[,myOUTCOME]))))
		smallSampPred		<- (NROW(DATAMAT) - max(table(DATAMAT[,myPREDICTOR])))
		smallSampOut		<- (NROW(DATAMAT) - max(table(DATAMAT[,myOUTCOME])))
		Nc							<- NROW(listCOVARIATES)
		N								<- NROW(DATAMAT)

		if(verbose == 1){ 
				writeLines(paste("VARpred       = ",VARpred))
				writeLines(paste("VARout        = ",VARout))
				writeLines(paste("smallSampPred = ",smallSampPred))
				writeLines(paste("Nc            = ",Nc))
				print(cbind(missingCOV,varianceCOV))
		} 


		#------enough data to run--------------------#
		if(!(myPREDICTOR %in% myOUTCOME | myPREDICTOR %in% listCOVARIATES) & VARpred > 1 & VARout > 1 & smallSampPred > minSmallCat & smallSampOut > minSmallCat){
		if(Nc > 1){
			#####################################################
			#==OPTION1==========================================#
			#standardized all variables
			if(optionStand ==1){
					if(verbose == 1){ writeLines("Standardized all variables")} 
					TMP	= apply(DATAMAT,2,function(x) (x-mean(x, na.rm=TRUE))/sqrt(var(x, na.rm=TRUE)))
					DATAMAT	= TMP
			}
			
			#####################################################
			#==REDEFINITION=====================================#
			COVAR			= as.matrix(DATAMAT[,listCOVARIATES])
			OUTCO			= DATAMAT[,myOUTCOME]
			PREDI			= DATAMAT[,myPREDICTOR]
			if(sum(!is.na(myFIXCOV))>0){
			FIXCOV		= as.matrix(DATAMAT[,myFIXCOV])
			}
			
			#####################################################
			#==OPTION2==========================================#
			#impute missing covariates
			if(optionImpute ==1){
				if(verbose == 1){ writeLines("Imputed missing data for covariates")}
				for(x in  1:NCOL(COVAR)){		COVAR[is.na(COVAR[,x]),x]		= mean(COVAR[,x],na.rm=TRUE) }
				if(sum(!is.na(myFIXCOV))>0){
				for(x in  1:NCOL(FIXCOV)){	FIXCOV[is.na(FIXCOV[,x]),x] = mean(FIXCOV[,x],na.rm=TRUE) }
				}
			}
			#===================================================#
			#####################################################
		
		
			#####################################################
			#==STEP 1===========================================#
			#marginal effect of the predictor on the outcome
			fitMarg 		= lm(OUTCO~PREDI)
			betaPY			= summary(fitMarg)$coef[2,1]
			sdPY				= summary(fitMarg)$coef[2,2]
			pvalPY			= summary(fitMarg)$coef[2,4]		
			chi2ASSO		= (betaPY/sdPY)^2
			
			if(sum(!is.na(myFIXCOV))>0){
					fitAdj 			= lm(OUTCO ~ PREDI + as.matrix(FIXCOV))
					resMA[1:3]	= summary(fitAdj)$coef[2,c(1,2,4)]
			}else{
					resMA[1:3]	= c(betaPY,sdPY,pvalPY)
			}
			
			#Marginal effect of the predictor on the covariate
			betaPC = pvalPC = 0
			for(i in 1:Nc){
		 			fit <- lm(COVAR[,i]~PREDI)
		  		betaPC[i] = summary(fit)$coef[2,1];
		  		pvalPC[i] = summary(fit)$coef[2,4];
			}
		
			#Marginal effect of the covariate on the outcome
			gammaCY = pvalCY = 0
			for(i in 1:Nc){
		  		fit 				<- lm(OUTCO~COVAR[,i])
		  		gammaCY[i]	<- summary(fit)$coef[2,1];
		  		pvalCY[i] 	<- summary(fit)$coef[2,4];
			}
			
			gammaFX = 0
			if(sum(!is.na(myFIXCOV))>0){
			for(i in 1:NROW(myFIXCOV)){
					fit 				<- lm(OUTCO~FIXCOV[,i])
					gammaFX[i]	<- summary(fit)$coef[2,1];
	 		}}
	 		
			selGCY = which(abs(gammaCY) >= 1)
			gammaCY[selGCY] = 0.99*sign(gammaCY[selGCY])
			#===================================================#
			#####################################################
	
	
			#####################################################
			#==STEP 2: filter multivariate======================#
			pMUL 				= 0
			ss					= 1
			myCOV				= seq(1:Nc)
			PcutoffBias = 0.001
			while(pMUL< Tmul & NROW(myCOV)>1){
				 	fitMAN <- manova(COVAR[,myCOV] ~ PREDI)
					if(inherits(try(pMUL <- summary(fitMAN, test="Wilks")$stats[1,6],silent=TRUE), "try-error")){ pMUL = 0.000001 }
					if(pMUL< Tmul){
							ss=ss+1
							myCOV = which(rank(pvalPC) >= ss)
					}
			}
			if(ss==1){		PcutoffBias	= 0.001
			}else if(ss>1){	PcutoffBias	= pvalPC[rank(pvalPC,ties.method = "first")==ss] }
			
			if(verbose == 1){
					print(paste("ss=",ss))
			}
	
			#===================================================#
			#####################################################
	
	
			#####################################################
			#==STEP 3: filter univariate========================#
			w					<- 0.1 * pMUL * (1-gammaCY^2)
			
			key1=key2=chi2ASSO/dweight
			if(key2 <= sigmaMax){ key2 = key2}
			if(key2 >  sigmaMax){ key2 = max(2*sigmaMax - key2,0)}
			sdThres1	<- min(key1,sigmaMax)  		
			sdThres2	<- min(key2,sigmaMax)  		
			expPC   	<- betaPY*gammaCY - betaPC
			expSD			<- sqrt((1-gammaCY^2)/N)
			
			#joint estimates
			if(sum(!is.na(myFIXCOV))>0){
				corMat		<- cor(cbind(COVAR,FIXCOV), use="pairwise")
				CC				<- corMat*N
				BB				<- (t(cbind(COVAR,FIXCOV)) %*% OUTCO)
				locFX			<- (Nc+1):(Nc+NCOL(FIXCOV))
			}else{
				corMat		<- cor(COVAR, use="pairwise")
				CC				<- corMat*N
				BB				<- (t(COVAR) %*% OUTCO)
			}
			Nsel			<- NROW(myCOV)
			NselP			<- Nsel+1
			while((Nsel > 0) & (Nsel < NselP )){
					NselP						<- NROW(myCOV)
					if(sum(!is.na(myFIXCOV))>0){
							jointBetaYC 		<- solve(CC[myCOV,myCOV]) %*% BB[myCOV]
					 		rsqY				<- t(c(gammaCY[myCOV],gammaFX)) %*% solve(corMat[c(myCOV,locFX),c(myCOV,locFX)]) %*% c(gammaCY[myCOV],gammaFX)
				 	}else{
							jointBetaYC 		<- solve(CC[myCOV,myCOV]) %*% BB[myCOV]
					 		rsqY				<- t(gammaCY[myCOV]) %*% solve(corMat[myCOV,myCOV]) %*% gammaCY[myCOV]
				 	}
					weightJoint					<- rep(0,Nc)
					weightJoint[myCOV] 			<- abs(jointBetaYC)
				 	# interval around 0 ====================================
				 	bound1					= as.numeric(lapply(w * (1-rsqY) / weightJoint^2, function(x) min(sdThres1,x)))
					W1maxThres      		= as.numeric(rep(sqrt(1/N),Nc) *  bound1)
					W1minThres      		= as.numeric(rep(sqrt(1/N),Nc) * -bound1)
				 	# interval around the expected =========================
					bound2					= as.numeric(lapply(w * (1-rsqY) / weightJoint^2, function(x) min(sdThres2,x)))
					W2maxThres      		= as.numeric(+ expSD * bound2)
					W2minThres     			= as.numeric(- expSD * bound2)
					# final list of covariate ==============================
					myCOV					= which(((betaPC > W1minThres & betaPC < W1maxThres) |  (expPC > W2minThres & expPC < W2maxThres)) & pvalPC > PcutoffBias)
					Nsel					= NROW(myCOV)
			}
			#===================================================#
			#####################################################
			
		
			#####################################################
			#==STEP 4: Derive final test========================#
			if(NROW(myCOV)==0){ 
			  	myCOV				= NULL
			  	resMC[1:3]	= resMA[1:3]
			  	resMC[4]		= 0
			  	resMC[5]		= 0
			  	resMC[6]		= 0
			  	resMC[7]		= 0
			  	resMC[8]		= 0
			}else{
			 		if(sum(!is.na(myFIXCOV))>0){
							fitF 				<- lm(OUTCO ~ PREDI + COVAR[,myCOV] + as.matrix(FIXCOV))
							fitC 				<- lm(OUTCO ~ COVAR[,myCOV] + as.matrix(FIXCOV))
							rsF 				<- summary(lm(OUTCO ~ as.matrix(FIXCOV)))$adj
					}else{
							fitF 				<- lm(OUTCO ~ PREDI + COVAR[,myCOV])
							fitC 				<- lm(OUTCO ~ COVAR[,myCOV])
							rsF 				<- 0
							Yres 				<- fitC$res
							#fitF 				<- lm(Yres ~ PREDI)
					}
			  	resMC[1:3]	= summary(fitF)$coef[2,c(1,2,4)]
			  	resMC[4]		= summary(fitC)$r.squ - rsF
			  	resMC[5]		= NROW(myCOV)
			  	resMC[6]		= paste(myCOV, collapse = ",")
			  	resMC[7]		= pMUL
			  	resMC[8]		= PcutoffBias
			}
 ############# if no covariates ############################################################
	}else{
				OUTCO			= DATAMAT[,myOUTCOME]
				PREDI			= DATAMAT[,myPREDICTOR]
				if(sum(!is.na(myFIXCOV))>0){
					FIXCOV		<- as.matrix(DATAMAT[,myFIXCOV])
					fitF 			<- lm(OUTCO ~ PREDI + as.matrix(FIXCOV))
				}else{
					fitF 			<- lm(OUTCO ~ PREDI)
				}
				resMC[1:3] = resMA[1:3] = summary(fitF)$coef[2,c(1,2,4)]
	}
	}}
	
	if(verbose == 1){
		#selC		<- (((betaPC > W1minThres & betaPC < W1maxThres) |  (betaPC > W2minThres & betaPC < W2maxThres)) & (pvalPC > PcutoffBias))
		#selC		<- ((pvalPC > PcutoffBias))
		#detRes	<- cbind(selC,betaPC,pvalPC,W1minThres,W1maxThres,W2minThres,W2maxThres,jointBetaYC,gammaCY)
		#colnames(detRes) = c("sel","betaPC","pvalPC","minT1","maxT1","minT2","maxT2","jointBetaYC","gammaCY")
		#print(detRes)
		#print(PcutoffBias)
		#print(cbind(selC,betaPC,pvalPC))
	}

	return(c(resMA,resMC))

}

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
 



