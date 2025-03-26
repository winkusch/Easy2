setClass("INDEP2",
	representation = representation(
						strEqcCommand		=	"character",
						rcdCriterion		=	"character",
						acolPval			=	"character",
						astrPvalTag			=	"character",
						anumPvalLim			=	"numeric",
						colInChr			=	"character",
						colInPos			=	"character",
						numPosLim			=	"numeric",
						numPosRegionExtension =	"numeric",
						fileRecombRate		=	"character",
						colRecombRateChr	=	"character",
						colRecombRatePos	=	"character",
						colRecombRate		=	"character",
						numRecombRateLim 	= 	"numeric",
						fileGene			=   "character",
						colGeneChr			=	"character",
						colGenePosStart		=	"character",
						colGenePosStop		=	"character",
						colGeneName			=	"character",
						fileAnnot			=	"character",
						strAnnotTag			=	"character",
						colAnnotTag			=	"character",
						colAnnotChr			=	"character",
						colAnnotPos			=	"character",
						colAnnotCoord		=	"character",
						numAnnotPosLim		=	"numeric",
						blnAddIndepInfo 	= 	"logical",
						blnRunByAnalysis 	= 	"logical",
						strTag				=	"character"
						),
	prototype = prototype(
						strEqcCommand		=	"",
						rcdCriterion		=	"",
						acolPval			=	"",
						astrPvalTag			=	"",
						anumPvalLim			=	5e-8,
						colInChr			=	"",
						colInPos			=	"",
						numPosLim			=	500000,
						numPosRegionExtension =	-1, # by default will be set to numPosLim/2
						fileRecombRate		=	"",
						colRecombRateChr	=	"chr",
						colRecombRatePos	=	"position",
						colRecombRate		=	"COMBINED_rate.cM.Mb.",
						numRecombRateLim 	= 	20,
						fileGene			= 	"",
						colGeneChr			=	"Chr",
						colGenePosStart		=	"Pos1",
						colGenePosStop		=	"Pos2",
						colGeneName			=	"Gene",
						fileAnnot			=	"",
						strAnnotTag			=	"Annot",
						colAnnotTag			=	"",
						colAnnotChr			=	"Chr",
						colAnnotPos			=	"Pos",
						colAnnotCoord		=	"",
						numAnnotPosLim		=	-1,
						blnAddIndepInfo 	= 	FALSE,
						blnRunByAnalysis 	= 	TRUE,
						strTag				=	""
						)
	#contains = c("EcfReader")
)

setGeneric("setINDEP2", function(object) standardGeneric("setINDEP2"))
setMethod("setINDEP2", signature = (object = "INDEP2"), function(object) {
	
	aEqcSlotNamesIn = c("rcdCriterion","acolPval","astrPvalTag","anumPvalLim","colInChr","colInPos","numPosLim","numPosRegionExtension","fileRecombRate","colRecombRateChr","colRecombRatePos","colRecombRate","numRecombRateLim",
						"fileGene","colGeneChr","colGenePosStart","colGenePosStop","colGeneName",
						"fileAnnot","strAnnotTag","colAnnotTag","colAnnotChr","colAnnotPos","colAnnotCoord","numAnnotPosLim",
						"blnAddIndepInfo","blnRunByAnalysis","strTag")
	
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
validINDEP2 <- function(objINDEP2) {
	
	### Valid with GWADATA?
	
	# if(objINDEP2@rcdCriterion == "") 
		# cat(paste(" EASY WARNING:INDEP2\n No criterion rcdCriterion defined. All data will be used for independentisation.", sep=""))
		
	if(objINDEP2@acolPval[1] == "") 
		stop(paste(" EASY ERROR:INDEP2\n No column acolPval defined. Please set acolPval.", sep=""))
	if(objINDEP2@colInChr == "") 
		stop(paste(" EASY ERROR:INDEP2\n No column colInChr defined. Please set colInChr.", sep=""))
	if(objINDEP2@colInPos == "") 
		stop(paste(" EASY ERROR:INDEP2\n No column colInPos defined. Please set colInPos.", sep=""))


	# ## check fileRecombRate
	# if(objINDEP2@fileRecombRate=="") 
		# stop(paste("EASY ERROR:INDEP2\n File fileRecombRate undefined.\n Please set --fileRecombRate or reset --blnRecombRate 0.", sep=""))
	
	if(!objINDEP2@fileRecombRate=="") {
	
		if(!file.exists(objINDEP2@fileRecombRate))
			stop(paste("EASY ERROR:INDEP2\n File fileRecombRate\n ",objINDEP2@fileRecombRate,"\n does not exist.\n Please check path or remove --fileRecombRate or reset --blnRecombRate 0.", sep=""))
		
		tr = read.table(objINDEP2@fileRecombRate, sep="\t", header=T, stringsAsFactors=F, nrows=1)
			
		isAv = objINDEP2@colRecombRateChr %in% names(tr)
		if(!isAv)
			stop(paste(" EASY ERROR:INDEP2\n Defined column --colRecombRateChr \n",objINDEP2@colRecombRateChr, "\n is not available in recombination rate file \n",objINDEP2@fileGene,"\n PLease specify correct column name.", sep=""))
		
		isAv = objINDEP2@colRecombRatePos %in% names(tr)
		if(!isAv)
			stop(paste(" EASY ERROR:INDEP2\n Defined column --colRecombRatePos \n",objINDEP2@colRecombRatePos, "\n is not available in recombination rate file \n",objINDEP2@fileGene,"\n PLease specify correct column name.", sep=""))

		isAv = objINDEP2@colRecombRate %in% names(tr)
		if(!isAv)
			stop(paste(" EASY ERROR:INDEP2\n Defined column --colRecombRate \n",objINDEP2@colRecombRate, "\n is not available in recombination rate file \n",objINDEP2@fileGene,"\n PLease specify correct column name.", sep=""))
	}
	
	## check fileGene
	if(objINDEP2@fileGene != "") {
		if(!file.exists(objINDEP2@fileGene))
			stop(paste("EASY ERROR:INDEP2\n File fileGene\n ",objINDEP2@fileGene,"\n does not exist.\n Please check path or remove --fileGene.", sep=""))
			
		tg = read.table(objINDEP2@fileGene, sep="\t", header=T, stringsAsFactors=F, nrows=1)
		
		isAv = objINDEP2@colGeneChr %in% names(tg)
		if(!isAv)
			stop(paste(" EASY ERROR:INDEP2\n Defined column --colGeneChr \n",objINDEP2@colGeneChr, "\n is not available in gene file \n",objINDEP2@fileGene,"\n PLease specify correct column name.", sep=""))
		
		isAv = objINDEP2@colGenePosStart %in% names(tg)
		if(!isAv)
			stop(paste(" EASY ERROR:INDEP2\n Defined column --colGenePosStart \n",objINDEP2@colGenePosStart, "\n is not available in gene file \n",objINDEP2@fileGene,"\n PLease specify correct column name.", sep=""))

		isAv = objINDEP2@colGenePosStop %in% names(tg)
		if(!isAv)
			stop(paste(" EASY ERROR:INDEP2\n Defined column --colGenePosStop \n",objINDEP2@colGenePosStop, "\n is not available in gene file \n",objINDEP2@fileGene,"\n PLease specify correct column name.", sep=""))
	
		isAv = objINDEP2@colGeneName %in% names(tg)
		if(!isAv)
			stop(paste(" EASY ERROR:INDEP2\n Defined column --colGeneName \n",objINDEP2@colGeneName, "\n is not available in gene file \n",objINDEP2@fileGene,"\n PLease specify correct column name.", sep=""))
		
	}
	
	
	## check fileAnnot
	if(objINDEP2@fileAnnot != "") {
		if(!file.exists(objINDEP2@fileAnnot))
			stop(paste("EASY ERROR:INDEP\n File --fileAnnotGene\n ",objINDEP2@fileAnnot,"\n does not exist.\n Please check path or remove --fileAnnot.", sep=""))
			
		tk = read.table(objINDEP2@fileAnnot, sep="\t", header=T, stringsAsFactors=F, nrows=1)
		
		isUseCoord = objINDEP2@colAnnotCoord != ""
		
		if(isUseCoord) {
		
			isAv = objINDEP2@colAnnotCoord %in% names(tk)
				if(!isAv)
					stop(paste(" EASY ERROR:INDEP2\n Defined column --colAnnotCoord \n",objINDEP2@colAnnotCoord, "\n is not available in known loci file \n",objINDEP2@fileAnnot,"\n PLease specify correct column name.", sep=""))
		
		} else {
		
			isAv = objINDEP2@colAnnotChr %in% names(tk)
			if(!isAv)
				stop(paste(" EASY ERROR:INDEP2\n Defined column --colAnnotChr \n",objINDEP2@colAnnotChr, "\n is not available in known loci file \n",objINDEP2@fileAnnot,"\n PLease specify correct column name.", sep=""))
			
			isAv = objINDEP2@colAnnotPos %in% names(tk)
			if(!isAv)
				stop(paste(" EASY ERROR:INDEP2\n Defined column --colAnnotPos \n",objINDEP2@colAnnotPos, "\n is not available in known loci file \n",objINDEP2@fileAnnot,"\n PLease specify correct column name.", sep=""))
			
		}
		
		if(objINDEP2@colAnnotTag != "") {
			isAv = objINDEP2@colAnnotTag %in% names(tk)
			if(!isAv)
				stop(paste(" EASY ERROR:INDEP2\n Defined column --colAnnotTag \n",objINDEP2@colAnnotTag, "\n is not available in known loci file \n",objINDEP2@fileAnnot,"\n PLease specify correct column name.", sep=""))
		}
	}
	
	return(TRUE)
}
INDEP2.GWADATA.valid <- function(objINDEP2, objGWA) {
	
	aisMatch = objINDEP2@acolPval %in% objGWA@aHeader
	if(any(!aisMatch))
		stop(paste("EASY ERROR:INDEP2\n Column \n",objINDEP2@acolPval[which(!aisMatch)[1]]," does not exist in file\n",objGWA@fileInShortName,"\n !!!" ,sep=""))
	
	isAv <- objINDEP2@colInChr %in% objGWA@aHeader
	if(!isAv)
		stop(paste(" EASY ERROR:INDEP2\n Defined column colInChr \n",objINDEP2@colInChr, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	isAv <- objINDEP2@colInPos %in% objGWA@aHeader
	if(!isAv)
		stop(paste(" EASY ERROR:INDEP2\n Defined column colInPos \n",objINDEP2@colInPos, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	if(length(objINDEP2@anumPvalLim)==1 & length(objINDEP2@acolPval)>1) objINDEP2@anumPvalLim = rep(objINDEP2@anumPvalLim[1], length(objINDEP2@acolPval))
	
	if(all(objINDEP2@astrPvalTag=="")) objINDEP2@astrPvalTag = objINDEP2@acolPval
	
	if(objINDEP2@numPosRegionExtension == -1) objINDEP2@numPosRegionExtension = objINDEP2@numPosLim/2
	
	if(length(objINDEP2@acolPval)==1) objINDEP2@blnRunByAnalysis = FALSE
	
	return(objINDEP2)
	
}

INDEP2.read <- function(objINDEP2, blnReadAll) {
	
	fR = objINDEP2@fileRecombRate
	
	if(blnReadAll) {
		tblR = read.table(fR, header=T, stringsAsFactors=FALSE, strip.white = TRUE, comment.char = "")
	} else {
		tblR = read.table(fR, header=T, stringsAsFactors=FALSE, strip.white = TRUE, comment.char = "", nrows = 10)
	}
	
	tblRShort = tblR[,c(objINDEP2@colRecombRateChr,objINDEP2@colRecombRatePos,objINDEP2@colRecombRate)]
	names(tblRShort) <- c("chr","pos","recomb_rate_cM_Mb")
	
	###sort by chr and pos
	tblRShortOrdered = tblRShort[order(as.integer(tblRShort$chr),as.integer(tblRShort$pos)),]
	
	return(tblRShortOrdered)
}

#############################################################################################################################
INDEP2.run <- function(objINDEP2, objGWA, objREPORT, tblRR, isValidScript) {
	
	# objINDEP2, objGWA, objREPORT, tblRR, isValidScript
	
	
	
	rcdCriterion 	<- objINDEP2@rcdCriterion
	acolPval		<- objINDEP2@acolPval
	astrPvalTag 	<- objINDEP2@astrPvalTag
	anumPvalLim		<- objINDEP2@anumPvalLim
	colInChr		<- objINDEP2@colInChr
	colInPos		<- objINDEP2@colInPos
	numPosLim		<- objINDEP2@numPosLim
	numPosRegionExtension <- objINDEP2@numPosRegionExtension
	
	fileRecombRate <- objINDEP2@fileRecombRate
	numRecombRateLim <- objINDEP2@numRecombRateLim
	
	fileGene 	<- objINDEP2@fileGene
	colGeneChr 	<- objINDEP2@colGeneChr
	colGenePosStart <- objINDEP2@colGenePosStart
	colGenePosStop 	<- objINDEP2@colGenePosStop
	colGeneName 	<- objINDEP2@colGeneName
	
	fileAnnot 	<- objINDEP2@fileAnnot
	strAnnotTag 	<- objINDEP2@strAnnotTag
	colAnnotTag 	<- objINDEP2@colAnnotTag
	colAnnotChr 	<- objINDEP2@colAnnotChr
	colAnnotPos 	<- objINDEP2@colAnnotPos
	colAnnotCoord  	<- objINDEP2@colAnnotCoord
	numAnnotPosLim 	<- objINDEP2@numAnnotPosLim
	
	blnAddIndepInfo	<- objINDEP2@blnAddIndepInfo
	
	blnRunByAnalysis <- objINDEP2@blnRunByAnalysis
	
	strTag			<- objINDEP2@strTag
	
	#### 
	
	if(nchar(strTag)>0) strTag = paste(strTag,".",sep="") 
	
	blnRecombRate = !(fileRecombRate=="")
	blnGene = !(fileGene=="")
	blnAnnot = !(fileAnnot=="")
	if(numAnnotPosLim == -1) numAnnotPosLim = numPosLim
	
	### 
	
	if(rcdCriterion != "") {
		objRCD 	<- RCD(rcdCriterion)
		out 	<- RCD.eval(objRCD, objGWA)
		numIndepCrit = length(which(out))
		objGWA.indep <- GWADATA.getrows(objGWA, which(out))
	} else {
		## use anumPvalLim's
		
		isUsed = rep(FALSE, nrow(objGWA@tblGWA))
		for(i in 1:length(anumPvalLim)) {
			
			isUsedi = GWADATA.getcol(objGWA, acolPval[i]) < anumPvalLim[i]
			isUsed = isUsed | isUsedi
		}
		
		numIndepCrit = length(which(isUsed))
		objGWA.indep <- GWADATA.getrows(objGWA, which(isUsed))
		
		# numIndepCrit <- nrow(objGWA@tblGWA)
		# objGWA.indep <- objGWA
	}
	
	objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numIndepIn",sep=""),NA)
	objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numPosMissing",sep=""),NA)
	for(i in 1:length(acolPval)) objREPORT <- REPORT.addval(objREPORT,paste(strTag,astrPvalTag[i],".numVarSignif",sep=""),NA)	
	objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numRegion",sep=""),NA)
	for(i in 1:length(acolPval)) objREPORT <- REPORT.addval(objREPORT,paste(strTag,astrPvalTag[i],".numRegion",sep=""),NA)
	if(blnRecombRate) objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numSignal",sep=""),NA)

	if(numIndepCrit == 0) {
	
		objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numIndepIn",sep=""),0)
		objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numPosMissing",sep=""),0)
		for(i in 1:length(acolPval)) objREPORT <- REPORT.setval(objREPORT,paste(strTag,astrPvalTag[i],".numVarSignif",sep=""),0)	
		objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numRegion",sep=""),0)
		for(i in 1:length(acolPval)) objREPORT <- REPORT.setval(objREPORT,paste(strTag,astrPvalTag[i],".numRegion",sep=""),0)
		if(blnRecombRate) objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numSignal",sep=""),0)
		
		
		if(blnAddIndepInfo) {
			
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"pMin",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionId",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionLead",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionTag",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionGwsVariants",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionGwsCoordinates",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionGwsSize",sep=""))
			# objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionExtCoordinates",sep=""))
			# objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionExtSize",sep=""))
			
			if(blnRecombRate) {
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalId",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalLead",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalTag",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalGwsVariants",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalCoordinates",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalSize",sep=""))
			}
			
			if(blnAnnot) {
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionAnnot",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"annotDistance",sep=""))
				if(blnRecombRate) objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalAnnot",sep=""))
			}
			
			if(blnGene) {
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"NearestGene",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"NearestGeneDistance",sep=""))
			}

		}
		
		return(list(objGWA,objGWA.indep,objGWA.indep,objGWA.indep,objREPORT))
	}
	
	objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numIndepIn",sep=""),numIndepCrit)
	
	
	##############
	### remove missing positions
	
	isNaPos = is.na(GWADATA.getcol(objGWA.indep, colInChr)) | is.na(GWADATA.getcol(objGWA.indep, colInPos))
	numPosMissing = length(which(isNaPos))

	objGWA.indep@tblGWA <- 	objGWA.indep@tblGWA[!isNaPos,]

	objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numPosMissing",sep=""),numPosMissing)


	##############
	### reduce to significant variants and count significant variants by tag
	
	isSignif = rep(FALSE,nrow(objGWA.indep@tblGWA))
	
	for(i in 1:length(acolPval)) { 
		
		apvalTmp = GWADATA.getcol(objGWA.indep, acolPval[i])
		
		isSignifTmp = apvalTmp < anumPvalLim[i] & !is.na(apvalTmp)
		
		objREPORT <- REPORT.setval(objREPORT,paste(strTag,astrPvalTag[i],".numVarSignif",sep=""),length(which(isSignifTmp)))
		
		isSignif = isSignif | isSignifTmp
	}
	
	if(all(!isSignif)) {
		
		objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numRegion",sep=""),0)
		for(i in 1:length(acolPval)) objREPORT <- REPORT.setval(objREPORT,paste(strTag,astrPvalTag[i],".numRegion",sep=""),0)
		if(blnRecombRate) objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numSignal",sep=""),0)
		
		
		if(blnAddIndepInfo) {
			
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"pMin",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionId",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionLead",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionTag",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionGwsVariants",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionGwsCoordinates",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionGwsSize",sep=""))
			# objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionExtCoordinates",sep=""))
			# objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionExtSize",sep=""))
			
			if(blnRecombRate) {
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalId",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalLead",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalTag",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalGwsVariants",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalCoordinates",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalSize",sep=""))
			}
			
			if(blnAnnot) {
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionAnnot",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"annotDistance",sep=""))
				if(blnRecombRate) objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalAnnot",sep=""))
			}
			
			if(blnGene) {
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"NearestGene",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"NearestGeneDistance",sep=""))
			}

		}
		
		return(list(objGWA,objGWA.indep,objGWA.indep,objGWA.indep,objREPORT))
		
	}
	### return if no variants are left
	
	objGWA.indep <- GWADATA.getrows(objGWA.indep, which(isSignif))
	
	##############
	### start independetization
		
	tIn = as.data.frame(objGWA.indep@tblGWA)
	
	if(length(acolPval)>1) {
		aPmin = apply(tIn[,acolPval],1,min,na.rm=T)
	} else {
		aPmin = tIn[,acolPval]
	}
	aidxsort = order(aPmin)
	
	tInSort = tIn[aidxsort,]
	aPminSort = aPmin[aidxsort]
	
	aChr = tInSort[,colInChr]
	aPos = tInSort[,colInPos]
	
	### work with aChr, aPos, aPminSort and tInSort[,acolPval] to obtain : 
	
	####################################
	### 1. Obtain regionId, regionLead across analyses
	
	regionId <- regionLead <- regionGwsCoordinates <- regionGwsSize <- regionExtCoordinates <- regionExtSize <- rep(NA, nrow(tInSort))
	
	regioncount = 1
	
	aPminSortBackup = aPminSort
	
	while(any(is.na(regionId))) {
			
		# print(regioncount)
		
		iTmpExtr = which(aPminSort == min(aPminSort))[1]
				
		chrExtr =  aChr[iTmpExtr]
		posExtr =  aPos[iTmpExtr]
		
		### use position mapping to obtain regions
		
		isCurRegion = NA
		isCurRegionNew = aChr==chrExtr & aPos>=(posExtr-numPosLim) & aPos<=(posExtr+numPosLim)
		
		while(!identical(isCurRegion,isCurRegionNew)) {
			
			isCurRegion = isCurRegionNew
			
			regionPos1 = min(aPos[isCurRegion])
			regionPos2 = max(aPos[isCurRegion])
			
			isCurRegionNew = aChr==chrExtr & aPos>=(regionPos1-numPosLim) & aPos<=(regionPos2+numPosLim)
			
		}
		
		regionId[isCurRegionNew] = regioncount
		regionLead[iTmpExtr] = regioncount
		# regionGwsCoordinates[isCurRegionNew] = paste(chrExtr,":",regionPos1-numPosLim,"_",regionPos2+numPosLim,sep="")
		# INDEP210 update: It is important to increase the regions by numPosLim/2 because otherwise regions will overlap !
		regionGwsCoordinates[isCurRegionNew] = paste(chrExtr,":",regionPos1,"_",regionPos2,sep="")
		regionGwsSize[isCurRegionNew] = regionPos2 - regionPos1
		regionExtCoordinates[isCurRegionNew] = paste(chrExtr,":",regionPos1-numPosRegionExtension,"_",regionPos2+numPosRegionExtension,sep="")
		regionExtSize[isCurRegionNew] = (regionPos2+numPosRegionExtension) - (regionPos1-numPosRegionExtension)
		
		regioncount = regioncount + 1
		
		aPminSort[isCurRegionNew] = Inf
		
	}
	
	aPminSort = aPminSortBackup
	
	
	####################################
	### 2. Obtain signalId, signalLead, signalCoordinates across analyses
	
	if(blnRecombRate) {
	
		signalId <- signalLead <- signalCoordinates <- signalSize <- rep(NA, nrow(tInSort))
		
		tblRR = tblRR[order(tblRR$chr,tblRR$pos),]
		
		for(rid in unique(regionId)) {
			
			# if(rid==82) stop()
			
			isRegion = regionId == rid
			
			signalcount = 1
			
			aPminSortBackup = aPminSort
			
			while(any(is.na(signalId[isRegion]))) {
				
				# print(paste(rid,signalcount,sep="."))
				
				iTmpExtr = which(aPminSort == min(aPminSort[isRegion]) & isRegion)[1]
				
				chrExtr =  aChr[iTmpExtr]
				posExtr =  aPos[iTmpExtr]
				
				isChr = tblRR$chr==chrExtr
				isPosLow = tblRR$pos<=posExtr
				isPosHigh = tblRR$pos>=posExtr
				
				if(any(isChr & isPosLow)) {
					pos1 = tblRR$pos[max(which(isChr & isPosLow))]
				} else {
					pos1 = 0
				}
				
				if(any(isChr & isPosHigh)) {
					pos2 = tblRR$pos[min(which(isChr & isPosHigh))]
				} else {
					pos2 = Inf
				}
				
				isCurSignal = aChr == chrExtr & aPos >= pos1 & aPos <= pos2 & is.na(signalId) & isRegion
				
				# signalPos1 = min(aPos[isCurSignal])
				# signalPos2 = min(aPos[isCurSignal])
				
				signalId[isCurSignal] = signalcount
				signalLead[iTmpExtr] = signalcount
				signalCoordinates[isCurSignal] = paste(chrExtr,":",pos1,"_",pos2,sep="")
				if(pos2==Inf) {
					signalSize[isCurSignal] = max(aPos[isCurSignal])-pos1
				} else {
					signalSize[isCurSignal] = pos2-pos1
				}
				
				signalcount = signalcount + 1
				
				aPminSort[isCurSignal] = Inf
				
			}
			
			aPminSort = aPminSortBackup
			
		}
		
		signalCoordinates = gsub("Inf","EoCHR",signalCoordinates)
		
		signalId = paste(regionId,signalId,sep=".")
		signalLead = ifelse(!is.na(signalLead), signalId, NA)
	}
	
	####################################
	### 3. Obtain signalTag and regionTag 

	regionTag <- regionGwsVariants <- regionNumSignals <- rep(NA, nrow(tInSort))
	
	for(rid in unique(regionId)) {
		
		# print(rid)
		
		isRegion = regionId == rid
		
		numvar = length(which(isRegion))
		if(blnRecombRate) numsig = length(unique(signalId[isRegion]))
		
		artag = c()
		arnumvar = c()
		arnumsig = c()
		
		for(i in 1:length(acolPval)) { 
			
			colp = acolPval[i]
			tag = astrPvalTag[i]
			plim = anumPvalLim[i]
			
			isRidSignif = isRegion & tInSort[,colp]<plim
			isRidSignif[is.na(isRidSignif)] = FALSE
			
			if(any(isRidSignif)) { 
				artag = c(artag, astrPvalTag[i])
				arnumvar = c(arnumvar, length(which(isRidSignif)))
				if(blnRecombRate) arnumsig = c(arnumsig, length(unique(signalId[isRidSignif])))
			} 
		}
		
		strrtag = ifelse(length(artag)==0,NA,ifelse(length(artag)==1, artag, paste(artag,collapse=";")))
		regionTag[isRegion] = strrtag
		
		strNumVar = ifelse(length(arnumvar)==0,NA,ifelse(length(arnumvar)==1, as.character(arnumvar), paste(numvar,"(",paste(arnumvar,collapse=";"),")",sep="")))
		regionGwsVariants[isRegion] = strNumVar
		
		if(blnRecombRate) {
			strNumSig = ifelse(length(arnumsig)==0,NA,ifelse(length(arnumsig)==1, as.character(arnumsig), paste(numsig,"(",paste(arnumsig,collapse=";"),")",sep="")))
			regionNumSignals[isRegion] = strNumSig
		}
	}
	
	if(blnRecombRate) {
		signalTag <- signalGwsVariants <- rep(NA, nrow(tInSort))
		
		for(sid in unique(signalId)) {
			
			isSignal = signalId == sid
			
			numvar = length(which(isSignal))
			
			astag = c()
			asnumvar = c()
			
			for(i in 1:length(acolPval)) { 
				
				colp = acolPval[i]
				tag = astrPvalTag[i]
				plim = anumPvalLim[i]
				
				isSidSignif = isSignal & tInSort[,colp]<plim
				isSidSignif[is.na(isSidSignif)] = FALSE
				
				if(any(isSidSignif)) {
					astag = c(astag, astrPvalTag[i])
					asnumvar = c(asnumvar, length(which(isSidSignif)))
				}
			}
			
			strstag = ifelse(length(astag)==0,NA,ifelse(length(astag)==1, astag, paste(astag,collapse=";")))
			signalTag[isSignal] = strstag
			
			strNumVar = ifelse(length(asnumvar)==0,NA,ifelse(length(asnumvar)==1, as.character(asnumvar), paste(numvar,"(",paste(asnumvar,collapse=";"),")",sep="")))
			signalGwsVariants[isSignal] = strNumVar
		}
	}
	
	####################################
	### 4. Optional Annotation
	
	aNearestGene <- aNearestGeneDistance <- rep(NA,length(aPos))
	
	if(blnGene) {
		
		cat("\n   -> Starting gene annotation ... ")
		
		tG = read.table(fileGene,sep="\t",stringsAsFactors=F,header=T)
				
		gchr = tG[,colGeneChr]
		gpos1 = as.integer(tG[,colGenePosStart])
		gpos2 = as.integer(tG[,colGenePosStop])
		gene = tG[,colGeneName]
		
		tG = data.frame(gchr,gpos1,gpos2,gene,stringsAsFactors=F)
		
		## existant from above
		# aChr 		= GWADATA.getcol(objGWA.indep, colInChr)
		# aPos 		= GWADATA.getcol(objGWA.indep, colInPos)
		
		for(i in 1:length(aChr)) {
			
			chr=aChr[i]
			pos=aPos[i]
			
			isWithinGene = (chr == tG$gchr) & (pos >= tG$gpos1) & (pos <= tG$gpos2)
			
			if(any(isWithinGene)) {
				aNearestGene[i] = paste(tG$gene[isWithinGene],collapse=",")
				aNearestGeneDistance[i] = 0
			} else {
				isChr = 
				tGChr = tG[tG$gchr==chr,]
				adismin = pmin(abs(pos-tGChr$gpos1),abs(pos-tGChr$gpos2))
				idxMin = which.min(adismin)
				aNearestGene[i] = tGChr$gene[idxMin[1]] # 1 if variant has equal dis to two genes
				aNearestGeneDistance[i] = min(adismin)
				## duplicates in glist -> take union !!!
				# if(pos-tGChr$gpos1[idxMin[1]]<0) aNearestGeneDistance[i] = - aNearestGeneDistance[i]
				if(pos-tGChr$gpos1[idxMin[1]]>0) aNearestGeneDistance[i] = - aNearestGeneDistance[i]
			}
		}
	}
	
	regionAnnot <- annotDistance <- signalAnnot <- rep(NA,length(aPos))
	
	if(blnAnnot) {
	
		cat("\n   -> Starting annotation ... ")
		
		tAn<-read.table(fileAnnot, header=T, sep="\t", stringsAsFactors=FALSE)
		
		# colAnnotCoord
		if(colAnnotCoord!="") {
			## extract values from coord 1:123_466
			acoord = tAn[,colAnnotCoord]
			anchr = unlist(lapply(strsplit(acoord,":"),function(x) x[1]))
			strpos = unlist(lapply(strsplit(acoord,":"),function(x) x[2]))
			anpos1 = as.integer(unlist(lapply(strsplit(strpos,"_"),function(x) x[1])))
			anpos2 = as.integer(unlist(lapply(strsplit(strpos,"_"),function(x) x[2])))
		} else {
			anchr = tAn[,colAnnotChr]
			anpos = as.integer(tAn[,colAnnotPos])
			anpos1 = anpos
			anpos2 = anpos
		}
		
		if(colAnnotTag != "") {
			antag = tAn[, colAnnotTag]
		} else {
			antag = rep(strAnnotTag,length(anpos1))
		}
		
		## if known lead variant is < numAnnotPosLim distant from region border (outer gws variant!; not (outer gws variant + numPosLim) """")
		
		for(rid in unique(regionId)) {
			
			# if(rid==135) stop()
			
			isRegion = regionId == rid
			rchr = unique(aChr[isRegion])
			rpos1 = min(aPos[isRegion],na.rm=T)
			rpos2 = max(aPos[isRegion],na.rm=T)
			
			## extend by annot distance
			rpos1 = rpos1-numAnnotPosLim
			rpos2 = rpos2+numAnnotPosLim
			
			## do not use regionGwsCoordinates because those were extended by numPoslim
			
			isKnown = anchr == rchr & ((rpos1>=anpos1 & rpos1<=anpos2) | (rpos2>=anpos1 & rpos2<=anpos2) | (rpos1<=anpos1 & rpos2>=anpos2))
			
			if(any(isKnown)) {
				
				atagknown = antag[isKnown]
				unitag = unique(atagknown)
				
				aanpos1tmp = anpos1[isKnown]
				aanpos2tmp = anpos2[isKnown]
				
				countannots = 1
				
				for(ik in 1:length(which(isKnown))) {
					if(countannots == 1) {
						is1 = abs(aanpos1tmp[ik]-aPos[isRegion])<abs(aanpos2tmp[ik]-aPos[isRegion])
						adistmp = ifelse(is1, aanpos1tmp[ik] - aPos[isRegion], aanpos2tmp[ik] - aPos[isRegion])
					} else {
						is1 = abs(aanpos1tmp[ik]-aPos[isRegion])<abs(aanpos2tmp[ik]-aPos[isRegion])
						adistmp = c(adistmp, ifelse(is1, aanpos1tmp[ik] - aPos[isRegion], aanpos2tmp[ik] - aPos[isRegion]))
					}
					countannots = countannots + 1
				}
				
				if(length(unitag)>0) {
					aidxsorttag = order(unitag)
					regionAnnot[isRegion] = paste(unitag[aidxsorttag],collapse=";")	
					annotDistance[isRegion] = paste(adistmp[aidxsorttag],collapse=";")	
				}
				
			}
		}
		
		if(blnRecombRate) {
			for(sid in unique(signalId)) {
				
				isSignal = signalId == sid
				schr = unique(aChr[isSignal])
				spos = unique(unlist(lapply(strsplit(signalCoordinates[isSignal],":"), function(x) x[2])))
				spos = gsub("EoCHR","Inf",spos)
				spos1 = as.numeric(unlist(lapply(strsplit(spos,"_"), function(x) x[1])))
				spos2 = as.numeric(unlist(lapply(strsplit(spos,"_"), function(x) x[2])))
				
				isKnown = anchr == schr & ((anpos1 >= (spos1-numAnnotPosLim) & anpos1 <= (spos2+numAnnotPosLim)) | (anpos2 >= (spos1-numAnnotPosLim) & anpos2 <= (spos2+numAnnotPosLim)))
				
				isKnown = anchr == schr & ((spos1>=anpos1 & spos1<=anpos2) | (spos2>=anpos1 & spos2<=anpos2) | (spos1<=anpos1 & spos2>=anpos2))
				
				if(any(isKnown)) {
					unitag = unique(antag[isKnown])
					unitag = unitag[!is.na(unitag)]	
					if(length(unitag)>0) signalAnnot[isSignal] = paste(sort(unitag),collapse=";")
				}
			}
		}

	}
	
	####################################
	### 5. Add variables

	## reset objGWA to sorted table
	objGWA.indep@tblGWA <- as.data.table(tInSort)
		
	objGWA.indep = GWADATA.cbind(objGWA.indep, regionId, paste(strTag,"regionId",sep=""))
	objGWA.indep = GWADATA.cbind(objGWA.indep, regionLead, paste(strTag,"regionLead",sep=""))
	objGWA.indep = GWADATA.cbind(objGWA.indep, regionTag, paste(strTag,"regionTag",sep=""))
	
	if(blnRecombRate) objGWA.indep = GWADATA.cbind(objGWA.indep, regionNumSignals, paste(strTag,"regionNumSignals",sep=""))
	
	objGWA.indep = GWADATA.cbind(objGWA.indep, regionGwsVariants, paste(strTag,"regionGwsVariants",sep=""))
	objGWA.indep = GWADATA.cbind(objGWA.indep, regionGwsCoordinates, paste(strTag,"regionGwsCoordinates",sep=""))
	objGWA.indep = GWADATA.cbind(objGWA.indep, regionGwsSize, paste(strTag,"regionGwsSize",sep=""))
	# objGWA.indep = GWADATA.cbind(objGWA.indep, regionExtCoordinates, paste(strTag,"regionExtCoordinates",sep=""))
	# objGWA.indep = GWADATA.cbind(objGWA.indep, regionExtSize, paste(strTag,"regionExtSize",sep=""))
	
	
	if(blnAnnot) {
		objGWA.indep = GWADATA.cbind(objGWA.indep, regionAnnot, paste(strTag,"regionAnnot",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, annotDistance, paste(strTag,"annotDistance",sep=""))
	}	
	
	if(blnRecombRate) {
		objGWA.indep = GWADATA.cbind(objGWA.indep, signalId, paste(strTag,"signalId",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, signalLead, paste(strTag,"signalLead",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, signalTag, paste(strTag,"signalTag",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, signalGwsVariants, paste(strTag,"signalGwsVariants",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, signalCoordinates, paste(strTag,"signalCoordinates",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, signalSize, paste(strTag,"signalSize",sep=""))
		
		if(blnAnnot) objGWA.indep = GWADATA.cbind(objGWA.indep, signalAnnot, paste(strTag,"signalAnnot",sep=""))
	}
	
	if(blnGene) {
		objGWA.indep = GWADATA.cbind(objGWA.indep, aNearestGene, paste(strTag,"NearestGene",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, aNearestGeneDistance, paste(strTag,"NearestGeneDistance",sep=""))
	}
	
	
	
	########################################################################################################################################################################################################################
	########################################################################################################################################################################################################################
	########################################################################################################################################################################################################################
	########################################################################################################################################################################################################################
	########################################################################################################################################################################################################################
	########################################################################################################################################################################################################################
	########################################################################################################################################################################################################################
	
	###########################################################################################
	### 6. do the indep by analysis separately
	
	if(blnRunByAnalysis) {
	
		mergeid = seq(1,nrow(tInSort))
		
		objGWA.indep = GWADATA.cbind(objGWA.indep, mergeid, "mergeid")
		
		for(i in 1:length(acolPval)) { 
			
			colPvali = acolPval[i]
			tagi = astrPvalTag[i]
			
			objGWA.indepi = GWADATA.getcols(objGWA.indep, c("mergeid", colInChr, colInPos, colPvali))
			tIni = as.data.frame(objGWA.indepi@tblGWA)
						
			isSignif = tIni[,colPvali] < anumPvalLim[i] & !is.na(tIni[,colPvali])
			
			if(any(isSignif)) {
				tIni = tIni[which(isSignif),]
				
				aPi = tIni[,colPvali]
				aidxsort = order(aPi)
				
				tIniSort = tIni[aidxsort,]
				aPiSort = aPi[aidxsort]
				
				aChr = tIniSort[,colInChr]
				aPos = tIniSort[,colInPos]
				
				regionIdi <- regionLeadi <- regionGwsCoordinatesi <- regionGwsVariantsi <- rep(NA, nrow(tIniSort))
			
				regioncount = 1
				
				while(any(is.na(regionIdi))) {
						
					# print(regioncount)
					
					iTmpExtr = which(aPiSort == min(aPiSort))[1]
							
					chrExtr =  aChr[iTmpExtr]
					posExtr =  aPos[iTmpExtr]
																	
					### use position mapping to obtain regions
					
					isCurRegion = NA
					isCurRegionNew = aChr==chrExtr & aPos>=(posExtr-numPosLim) & aPos<=(posExtr+numPosLim)
					
					while(!identical(isCurRegion,isCurRegionNew)) {
						
						isCurRegion = isCurRegionNew
						
						regionPos1 = min(aPos[isCurRegion])
						regionPos2 = max(aPos[isCurRegion])
						
						isCurRegionNew = aChr==chrExtr & aPos>=(regionPos1-numPosLim) & aPos<=(regionPos2+numPosLim)
						
					}
					
					regionIdi[isCurRegionNew] = regioncount
					regionLeadi[iTmpExtr] = regioncount

					regionGwsCoordinatesi[isCurRegionNew] = paste(chrExtr,":",regionPos1,"_",regionPos2,sep="")
					regionGwsVariantsi[isCurRegionNew] = length(which(isCurRegionNew))
					
					regioncount = regioncount + 1
					
					aPiSort[isCurRegionNew] = Inf
					
				}
				
				objGWA.indepi@tblGWA = as.data.table(tIniSort)
				
				objGWA.indepi = GWADATA.cbind(objGWA.indepi, regionIdi, paste(strTag,tagi,".regionId",sep=""))
				objGWA.indepi = GWADATA.cbind(objGWA.indepi, regionLeadi, paste(strTag,tagi,".regionLead",sep=""))
				objGWA.indepi = GWADATA.cbind(objGWA.indepi, regionGwsVariantsi, paste(strTag,tagi,".regionGwsVariants",sep=""))
				objGWA.indepi = GWADATA.cbind(objGWA.indepi, regionGwsCoordinatesi, paste(strTag,tagi,".regionGwsCoordinates",sep=""))
				
				
				objGWA.indepi = GWADATA.removecols(objGWA.indepi, c(colInChr, colInPos, colPvali))
				## mergeid and regioni cols
				
				objGWA.indep = GWADATA.merge(objGWA.indep, objGWA.indepi, blnAll.In = TRUE, blnAll.Add = TRUE, strBy.In = "mergeid", strBy.Add = "mergeid", strSuffix.In = "", strSuffix.Add = "")
				
			} else {
				objGWA.indep = GWADATA.cbind(objGWA.indep, rep(NA,nrow(objGWA.indep@tblGWA)), paste(strTag,tagi,".regionId",sep=""))
				objGWA.indep = GWADATA.cbind(objGWA.indep, rep(NA,nrow(objGWA.indep@tblGWA)), paste(strTag,tagi,".regionLead",sep=""))
				objGWA.indep = GWADATA.cbind(objGWA.indep, rep(NA,nrow(objGWA.indep@tblGWA)), paste(strTag,tagi,".regionGwsCoordinates",sep=""))
				objGWA.indep = GWADATA.cbind(objGWA.indep, rep(NA,nrow(objGWA.indep@tblGWA)), paste(strTag,tagi,".regionGwsVariantsi",sep=""))
			}
		}
		
		objGWA.indep = GWADATA.removecols(objGWA.indep, "mergeid")
	
	}
	########################################################################################################################################################################################################################
	########################################################################################################################################################################################################################
	########################################################################################################################################################################################################################
	########################################################################################################################################################################################################################
	########################################################################################################################################################################################################################
	########################################################################################################################################################################################################################
	########################################################################################################################################################################################################################
	
	
	####################################
	### 7. Resort and prepare output

	## no resort: 
	objGWA.indep@tblGWA <- objGWA.indep@tblGWA
	
	### obtain file with any region lead (sometimes different leads by analysis)
	objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numRegion",sep=""),length(unique(regionId)))
	
	isRegionLead = !is.na(GWADATA.getcol(objGWA.indep, paste(strTag,"regionLead",sep="")))
	
	if(blnRunByAnalysis) {
		for(i in 1:length(astrPvalTag)) {
			tagi = astrPvalTag[i]
			aregionleadi = GWADATA.getcol(objGWA.indep, paste(strTag,tagi,".regionLead",sep=""))
			objREPORT <- REPORT.setval(objREPORT,paste(strTag,tagi,".numRegion",sep=""),length(unique(aregionleadi)))
			
			isRegionLead = isRegionLead | !is.na(aregionleadi)
		}
	}
	
	objGWA.indep.regionLeads <- GWADATA.getrows(objGWA.indep, isRegionLead)
	
	# ## regionTag
	# artleads = regionTag[which(!is.na(regionLead))]
	# artleads = artleads[grepl(";",artleads,fixed=T)]
	# if(length(artleads)>0) {
		# for(rt in sort(unique(artleads))) {
			# nrt = length(which(artleads==rt))
			# objREPORT <- REPORT.addval(objREPORT,paste(strTag,gsub(";","_AND_",rt,fixed=T),".numRegion",sep=""),nrt)
		# }
	# }
	#### THIS CAUSED REPORT PROBLEMS WHEN 0 VAR were signif in one data
	
	if(blnRecombRate) {
		objGWA.indep.signalLeads <- GWADATA.getrows(objGWA.indep, which(!is.na(signalLead)))
		objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numSignal",sep=""),length(unique(signalId)))
		## signalTag
		# astleads = signalTag[which(!is.na(signalLead))]
		# for(st in sort(unique(astleads))) {
			# nst = length(which(astleads==st))
			# objREPORT <- REPORT.addval(objREPORT,paste(strTag,gsub(";","_AND_",st,fixed=T),".numSignal",sep=""),nst)
		# }
		#### THIS CAUSED REPORT PROBLEMS WHEN 0 VAR were signif in one data
	} else {
		objGWA.indep.signalLeads <- objGWA.indep.regionLeads
	}
	
	if(blnAddIndepInfo) {
		### merges by intersect(names)		
		objGWA <- GWADATA.merge(objGWA,objGWA.indep, strSuffix.In = "", strSuffix.Add = "", blnAll.In = TRUE, blnAll.Add = TRUE, strBy.In = NA, strBy.Add = NA)
	}	

	return(list(objGWA,objGWA.indep,objGWA.indep.regionLeads,objGWA.indep.signalLeads,objREPORT))

}

INDEP2 <- function(strEqcCommand){ 
	## Wrapper for class definition
	INDEP2out <- setINDEP2(new("INDEP2", strEqcCommand = strEqcCommand))
	validINDEP2(INDEP2out)
	#INDEP2out.valid <- validINDEP2(INDEP2out)
	return(INDEP2out)
}

