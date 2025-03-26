setClass("PCA2STEP",
	representation = representation(
						strEqcCommand		=	"character",
						rcdStep1Crit		=	"character",
						colStep2Pval		=	"character",
						numPcaThrs			= 	"numeric",
						filePcaBed			= 	"character",
						filePcaSample		= 	"character",
						numR2PosSize		= 	"numeric",
						numR2VarSize		= 	"numeric",
						colInChr			=	"character",
						colInMarker		 	= 	"character",
						colOutPval			= 	"character",
						blnParal			= 	"logical",
						pathLibLoc 			=	"character",
						blnCrossChr 		= 	"logical",
						strTag				=	"character"
						),
	prototype = prototype(
						strEqcCommand		=	"",
						rcdStep1Crit		=	"",
						colStep2Pval		=	"",
						numPcaThrs			= 	0.995, 
						filePcaBed			= 	"",
						filePcaSample		= 	"",
						numR2PosSize		= 	-9,
						numR2VarSize		= 	-9,
						colInChr			=	"",
						colInMarker			=	"",
						colOutPval			= 	"",
						blnParal			= 	FALSE,
						pathLibLoc 			=	"",
						blnCrossChr 		= 	FALSE,
						strTag				=	"PCA2STEP"
						)
	#contains = c("EcfReader")
)

setGeneric("setPCA2STEP", function(object) standardGeneric("setPCA2STEP"))
setMethod("setPCA2STEP", signature = (object = "PCA2STEP"), function(object) {
	
	aEqcSlotNamesIn = c("rcdStep1Crit","colStep2Pval","numPcaThrs","filePcaBed","filePcaSample","numR2PosSize","numR2VarSize","colInChr","colInMarker","colOutPval","blnParal","pathLibLoc","blnCrossChr","strTag")
	
	objEqcReader <- EqcReader(object@strEqcCommand,aEqcSlotNamesIn)
	
	if(length(objEqcReader@lsEqcSlotsOut) > 0) {
		for(i in 1:length(objEqcReader@lsEqcSlotsOut)) {
			tmpSlot <- names(objEqcReader@lsEqcSlotsOut)[i]
			tmpSlotVal <- objEqcReader@lsEqcSlotsOut[[i]]
			
			if(all(!is.na(tmpSlotVal))) slot(object, tmpSlot) <- tmpSlotVal
		}
	}
	return(object)
})

#############################################################################################################################
validPCA2STEP <- function(objPCA2STEP) {
	
	### Valid with GWADATA?
	
	# if(objPCA2STEP@rcdStep1Crit == "") 
		# cat(paste(" EASY WARNING:PCA2STEP\n No criterion rcdStep1Crit defined. All data will be used for independentisation.", sep=""))
		
	if(objPCA2STEP@colStep2Pval == "") 
		stop(paste(" EASY ERROR:PCA2STEP\n No column colStep2Pval defined. Please set colStep2Pval.", sep=""))
	if(objPCA2STEP@colInChr == "") 
		stop(paste(" EASY ERROR:PCA2STEP\n No column colInChr defined. Please set colInChr.", sep=""))

	return(TRUE)
}
PCA2STEP.GWADATA.valid <- function(objPCA2STEP, objGWA) {

	isAv <- objPCA2STEP@colStep2Pval %in% objGWA@aHeader
	if(!isAv)
		stop(paste(" EASY ERROR:PCA2STEP\n Defined column colStep2Pval \n",objPCA2STEP@colStep2Pval, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	isAv <- objPCA2STEP@colInChr %in% objGWA@aHeader
	if(!isAv)
		stop(paste(" EASY ERROR:PCA2STEP\n Defined column colInChr \n",objPCA2STEP@colInChr, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	if(grepl("<CHR>",objPCA2STEP@filePcaBed,fixed=T)) {
		filePcaBed = c()
		for(chr in 1:22) filePcaBed = c(filePcaBed, gsub("<CHR>",chr,objPCA2STEP@filePcaBed,fixed=T))
		objPCA2STEP@filePcaBed = filePcaBed
	} else {
		objPCA2STEP@filePcaBed = rep(objPCA2STEP@filePcaBed,22)
	}
	
	for(filePcaBedi in objPCA2STEP@filePcaBed) {
		if(!file.exists(paste0(filePcaBedi,".bed"))) {
			stop(paste("EASY ERROR:PCA2STEP\n File filePcaBed\n ",filePcaBedi,"\n does not exist.\n Please check path or remove --filePcaBed.", sep=""))
		}
	}

	if(!objPCA2STEP@filePcaSample=="") {
		if(!file.exists(objPCA2STEP@filePcaSample)) {
				stop(paste("EASY ERROR:PCA2STEP\n File filePcaSample\n ",objPCA2STEP@filePcaSample,"\n does not exist.\n Please check path or remove --filePcaSample.", sep=""))
		}
	}
	
	if(!objPCA2STEP@pathLibLoc=="") {
		if(!file.exists(objPCA2STEP@pathLibLoc)) {
				stop(paste("EASY ERROR:PCA2STEP\n Path pathLibLoc\n ",objPCA2STEP@pathLibLoc,"\n does not exist.\n Please check path or remove --pathLibLoc.", sep=""))
		}
	}
	

	return(objPCA2STEP)
	
}

#############################################################################################################################
PCA2STEP.run <- function(objPCA2STEP, objGWA, objREPORT, isValidScript) {
	
	rcdStep1Crit 	<- objPCA2STEP@rcdStep1Crit
	colStep2Pval	<- objPCA2STEP@colStep2Pval
	colInChr		<- objPCA2STEP@colInChr
	colInMarker 	<- objPCA2STEP@colInMarker
	strTag			<- objPCA2STEP@strTag
	numPcaThrs  <- objPCA2STEP@numPcaThrs
	filePcaBed 	<- objPCA2STEP@filePcaBed
	filePcaSample <- objPCA2STEP@filePcaSample
	numR2PosSize <- objPCA2STEP@numR2PosSize
	numR2VarSize <- objPCA2STEP@numR2VarSize
	blnParal	<- objPCA2STEP@blnParal
	pathLibLoc	<- objPCA2STEP@pathLibLoc
	blnCrossChr	<- objPCA2STEP@blnCrossChr
	colOutPval	<- objPCA2STEP@colOutPval
	
	#### 
	
	if(nchar(strTag)>0) strTag = paste(strTag,".",sep="") 

	objRCD 	<- RCD(rcdStep1Crit)
	isStep1 	<- RCD.eval(objRCD, objGWA)
	
	isStep1 = isStep1 & !is.na(GWADATA.getcol(objGWA, colStep2Pval))
	
	numVarStep1Crit = length(which(isStep1))
	
	objGWA.pca <- GWADATA.getrows(objGWA, which(isStep1))
	
	objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numVarStep1Crit",sep=""),NA)
	objREPORT <- REPORT.addval(objREPORT,paste(strTag,"Meff",sep=""),NA)
	objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numVarMissing",sep=""),NA)
	objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numVarStep2Signif",sep=""),NA)
	
	if(numVarStep1Crit == 0 | !isValidScript) {
		
		objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numVarStep1Crit",sep=""),0)
		objREPORT <- REPORT.setval(objREPORT,paste(strTag,"Meff",sep=""),0)
		objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numVarMissing",sep=""),0)
		objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numVarStep2Signif",sep=""),0)
		
		strNewColName <- ifelse(colOutPval == "", paste(colStep2Pval,".pca",sep=""), colOutPval)
		objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), strNewColName)
		
		return(list(objGWA,objGWA.pca,objGWA.pca,objREPORT))
	}
	
	objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numVarStep1Crit",sep=""),numVarStep1Crit)
	
	##############
	### PCA analysis to derive number of independent tests 

	asnps = objGWA.pca@tblGWA[[colInMarker]]
	achr = objGWA.pca@tblGWA[[colInChr]]
	
	achruni = sort(as.integer(unique(achr)))
	nchruni = length(achruni)
	
	### read sample file
	tfam1 = fread(paste0(filePcaBed[1],".fam"),header=F)
	if(filePcaSample!="") {
		asampleused = scan(filePcaSample,sep="\n",what="character",quiet=TRUE)
		aidxsample = which(tfam1$V2%in%asampleused)
	} else {
		aidxsample = seq(1,nrow(tfam1))
	}
	
	if(blnCrossChr) {
		
		fnCor <- function(fileBedi, asnpsi, aidxsample, numR2PosSize) {
			
			tbimi = fread(paste0(fileBedi,".bim"),header=F)
			
			aidx = match(asnpsi,tbimi$V2)
			
			## there should be no missings 
			
			asnpsmiss = asnpsi[is.na(aidx)]
			aidx = aidx[!is.na(aidx)]
			
			aidx = sort(aidx)
			
			if(length(aidx)==0) return(list(NULL,asnpsmiss))
			if(length(aidx)==1) return(list(as(matrix(1,1,1),"dgCMatrix"),asnpsmiss))
			
			# snp_readBed2(paste0(fileBedi,".bed"), ind.col = aidx, ind.row = aidxsample)	
			# bin <- snp_attach(paste0(fileBedi,".rds"))
			
			# matcor = bed_cor(bed(paste0(filePcaBedi,".bed")), ind.col = aidx, ind.row = aidxsample, size=length(aidx))
			# matcor = bed_cor(bed(paste0(fileBedi,".bed")), ind.col = aidx, ind.row = aidxsample, size=500)
			
			if(numR2PosSize!=-9) {
				matcor = bed_cor(bed(paste0(fileBedi,".bed")), ind.col = aidx, ind.row = aidxsample, infos.pos = tbimi$V4[aidx], size = floor(numR2PosSize/1000))
			} else if(numR2VarSize!=-9) {
				matcor = bed_cor(bed(paste0(fileBedi,".bed")), ind.col = aidx, ind.row = aidxsample, size = numR2VarSize)
			} else {
				matcor = bed_cor(bed(paste0(fileBedi,".bed")), ind.col = aidx, ind.row = aidxsample, size=length(aidx))
			}
			
			
			# matcor = bed_cor(bed(paste0(fileBedi,".bed")), ind.col = aidx, ind.row = aidxsample, size=length(aidx))
			## size = 500 is default
			
			# matcor[upper.tri(matcor)] = 0
			matcor[is.na(matcor)] = 0
			matcor = as(matcor,"dgCMatrix")
			
			return(list(matcor,asnpsmiss))
		}
		
		
		if(blnParal) {
			cl<-makeCluster(nchruni)
			registerDoParallel(cl) 
			
			if(pathLibLoc == "") {
				clusterCall(cl, function(x) .libPaths(x), .libPaths())
			} else {
				clusterCall(cl, function(x) .libPaths(x), c(.libPaths(),pathLibLoc))
			}
			
			lscor = foreach(i=1:nchruni,.packages=c("bigsnpr","data.table")) %dopar% {
				fnCor(filePcaBed[as.integer(achruni[i])], asnps[which(achr==achruni[i])], aidxsample, numR2PosSize)
			}
			stopCluster(cl)
			
		} else {
			
			lscor = list()
			for(i in 1:nchruni) {
				print(i)
				lscor[[i]] = fnCor(filePcaBed[as.integer(achruni[i])], asnps[which(achr==achruni[i])], aidxsample, numR2PosSize)
			}
			
		}
		
		asnpmiss = unlist(lapply(lscor,function(x) x[[2]]))
		nmiss = length(asnpmiss)
		
		ansnp = unlist(lapply(lscor,function(x) ifelse(is.null(x[[1]]),0,nrow(x[[1]]))))
		ansnpcum = cumsum(ansnp)
		
		###################################
		### v1
		
		# start_time1 <- Sys.time()
		
		# nsnps = sum(asnpsmat)
		# matcorx2 = as(emptySparse(nrow=nsnps, ncol=nsnps),"dgCMatrix")
		# # for(i in 1:nchruni) {
		# for(i in 1:4) {
			# print(i)
			# if(i==1) {
				# matcorx2[(1:ansnpcum[1]),(1:ansnpcum[1])] = lscor[[i]][[1]]
			# } else { 
				# matcorx2[((ansnpcum[i-1]+1):ansnpcum[i]), ((ansnpcum[i-1]+1):ansnpcum[i])] = lscor[[i]][[1]]
			# }
		# }
		# stop_time1 <- Sys.time()
		# stop_time1 - start_time1
		
		### v2
		
		# start_time2 <- Sys.time()
		isStart = TRUE
		for(i in 1:nchruni) {
		# for(i in 1:4) {
			print(i)
			matcori = lscor[[i]][[1]]
			
			if(!is.null(matcori)) {
				# if(i==1) {
				if(isStart) {
					matcorx = matcori
					isStart = FALSE
				} else {
					
					## using MatrixExtra:emptySparse()
					# matcori = cbind(as(emptySparse(nrow=ansnp[i], ncol=ansnpcum[i-1]),"dgCMatrix"), matcori)
					# matcorx = cbind(matcorx,as(emptySparse(nrow=ansnpcum[i-1], ncol=ansnp[i]),"dgCMatrix"))
					# matcorx = rbind(matcorx,matcori)
					
					
					ncoli = ncol(matcori)
					nrowi = nrow(matcori)				
					ncolx = ncol(matcorx)
					nrowx = nrow(matcorx)				
					
					matcori = cbind(as(matrix(0,nrowi,ncolx),"dgCMatrix"), matcori)
					matcorx = cbind(matcorx,as(matrix(0,nrowx,ncoli),"dgCMatrix"))
					matcorx = rbind(matcorx,matcori)

				}
			}
		}
		# stop_time2 <- Sys.time()
		# stop_time2 - start_time2
		
		
		#ei = eigen(matcorx,symmetric=TRUE,only.values=TRUE)
		
		
		ei = eigen(matcorx,symmetric=TRUE,only.values=TRUE)
		# aev = sort(ei$values,decreasing=TRUE)/sum(ei$values)
		aeigval = (ei$values)
		aeigval[aeigval<0] = 0
		aexplvar = sort(aeigval,decreasing=TRUE)/sum(aeigval)
		acumexplvar = cumsum(aexplvar)
		meff = length(which(acumexplvar<numPcaThrs))+1
		
		# ## https://stats.stackexchange.com/questions/254592/calculating-pca-variance-explained
		# pca <- prcomp(USArrests, scale = TRUE)
		
		# pca = prcomp(matcori)
		# aexplvar <- pca$sdev^2/sum(pca$sdev^2)
		# acexplvar = cumsum(aexplvar)
		# meff = length(which(acexplvar<numPcaThrs))+1
		
		# t1=Sys.time()
		# ei = eigen(matcorx,symmetric=TRUE,only.values=TRUE)
		# aev = sort(ei$values,decreasing=TRUE)/sum(ei$values)
		# acev = cumsum(aev)
		# meff = length(which(acev<numPcaThrs))+1
		# t2 = Sys.time()
		# t2-t1
		
		
		# t1=Sys.time()
		# h = as.spam.dgCMatrix(matcorx)
		# # ei = eigen.spam(h,symmetric=TRUE,only.values=TRUE)
		# ei = eigen.spam(h,nev=ceiling(9077/2),symmetric=TRUE,only.values=TRUE,control=list(ncv=(ceiling(9077/2)+1)))
		# aev = sort(ei$values,decreasing=TRUE)/sum(ei$values,na.rm=T)
		# acev = cumsum(aev)
		# meff = length(which(acev<numPcaThrs))+1
		# t2 = Sys.time()
		# t2-t1
		
	} else {
		
		fnPca <- function(fileBedi, asnpsi, aidxsample, numPcaThrs, numR2PosSize) {
			
			tbimi = fread(paste0(fileBedi,".bim"),header=F)
			
			aidx = match(asnpsi,tbimi$V2)
			
			## there should be no missings 
			
			asnpsmiss = asnpsi[is.na(aidx)]
			aidx = aidx[!is.na(aidx)]
			
			aidx = sort(aidx)
			
			if(length(aidx)<=1) return(list(length(aidx),asnpsmiss))
			
			# matcor = bed_cor(bed(paste0(filePcaBedi,".bed")), ind.col = aidx, ind.row = aidxsample, size=length(aidx))
			# matcor = bed_cor(bed(paste0(fileBedi,".bed")), ind.col = aidx, ind.row = aidxsample, size=500)
			# matcor = bed_cor(bed(paste0(fileBedi,".bed")), ind.col = aidx, ind.row = aidxsample, size=length(aidx))
			if(numR2PosSize!=-9) {
				matcor = bed_cor(bed(paste0(fileBedi,".bed")), ind.col = aidx, ind.row = aidxsample, infos.pos = tbimi$V4[aidx], size = floor(numR2PosSize/1000))
			} else if(numR2VarSize!=-9){
				matcor = bed_cor(bed(paste0(fileBedi,".bed")), ind.col = aidx, ind.row = aidxsample, size = numR2VarSize)
			} else {
				matcor = bed_cor(bed(paste0(fileBedi,".bed")), ind.col = aidx, ind.row = aidxsample, size = length(aidx))
				
			}

			# set missings to 0 cor 
			matcor[is.na(matcor)] = 0
			
			ei = eigen(matcor,symmetric=TRUE,only.values=TRUE)
			aeigval = (ei$values)
			aeigval[aeigval<0] = 0
			aexplvar = sort(aeigval,decreasing=TRUE)/sum(aeigval)
			acumexplvar = cumsum(aexplvar)
			em = length(which(acumexplvar<numPcaThrs))+1
			## +1 to make sure explained var is indeed >numPcaThrs !
		
			
			return(list(em,asnpsmiss))
		}
		
		if(blnParal) {
			cl<-makeCluster(nchruni)
			registerDoParallel(cl) 
			# clusterCall(cl, function(x) .libPaths(x), .libPaths())
			
			if(pathLibLoc == "") {
				clusterCall(cl, function(x) .libPaths(x), .libPaths())
			} else {
				clusterCall(cl, function(x) .libPaths(x), c(.libPaths(),pathLibLoc))
			}
			
			lspca = foreach(i=1:nchruni,.packages=c("bigsnpr","data.table")) %dopar% {
				fnPca(filePcaBed[as.integer(achruni[i])], asnps[which(achr==achruni[i])], aidxsample, numPcaThrs, numR2PosSize)
			}
			stopCluster(cl)
			
		} else {
			
			lspca = list()
			for(i in 1:nchruni) {
				lspca[[i]] = fnPca(filePcaBed[as.integer(achruni[i])], asnps[which(achr==achruni[i])], aidxsample, numPcaThrs, numR2PosSize)
			}
			
		}
		
		ameff = unlist(lapply(lspca,function(x) x[[1]]))
		asnpmiss = unlist(lapply(lspca,function(x) x[[2]]))
		
		meff = sum(ameff)
		nmiss = length(asnpmiss)
	
	}
	
	### add on missings add individual variants
	# meff = meff + nmiss
	
	###
	objREPORT <- REPORT.setval(objREPORT,paste(strTag,"Meff",sep=""),meff)
	objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numVarMissing",sep=""),nmiss)
	
	asnpsin = GWADATA.getcol(objGWA, colInMarker)
	isMiss = asnpsin%in%asnpmiss
	
	apvalout = GWADATA.getcol(objGWA, colStep2Pval)
	apvalout[!isStep1 | is.na(isStep1) | isMiss] = NA
	apvalout = apvalout*meff
	apvalout[!is.na(apvalout)&apvalout>1]<-1
	
	strNewColName <- ifelse(colOutPval == "", paste(colStep2Pval,".pca",sep=""), colOutPval)
	objGWA <- GWADATA.cbind(objGWA, apvalout, strNewColName)
	
	isStep2Signif = apvalout<0.05
	objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numVarStep2Signif",sep=""),length(which(isStep2Signif)))
	
	objGWA.signif <- GWADATA.getrows(objGWA, which(isStep2Signif))
	objGWA.missbed <- GWADATA.getrows(objGWA, which(isMiss))
	
	return(list(objGWA,objGWA.signif,objGWA.missbed,objREPORT))

}

PCA2STEP <- function(strEqcCommand){ 
	## Wrapper for class definition
	PCA2STEPout <- setPCA2STEP(new("PCA2STEP", strEqcCommand = strEqcCommand))
	validPCA2STEP(PCA2STEPout)
	#PCA2STEPout.valid <- validPCA2STEP(PCA2STEPout)
	return(PCA2STEPout)
}

