setClass("INDEP",
	representation = representation(
						strEqcCommand		=	"character",
						rcdCriterion		=	"character",
						arcdCriterion		=	"character",
						acolPval			=	"character",
						astrTag				=	"character", # used for multiple acolPval
						colTag				=	"character", # used for one acolPval with different tags in the column
						anumPvalLim			=	"numeric",
						acolIndep			=	"character",
						anumIndepLim		=	"numeric",
						strIndepDir			=	"character",
						colInChr			=	"character",
						colInPos			=	"character",
						numPosLim			=	"numeric",
						numPosRegionExtension =	"numeric",
						colInGPos			=	"character",
						numGPosLim			=	"numeric",
						fileRecombRate		=	"character",
						colRecombRateChr	=	"character",
						colRecombRatePos	=	"character",
						colRecombRate		=	"character",
						numRecombRateLim 	= 	"numeric",
						fileClumpBed		= 	"character",
						fileClumpSample		= 	"character",
						numR2Thrs			= 	"numeric",
						blnRetainRegionLead = 	"logical",
						blnParal			= 	"logical",
						pathLibLoc			= 	"character",
						numR2PosSize		=	"numeric",
						numR2VarSize		=	"numeric",
						blnClumpInSignal	= 	"logical",
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
						#blnAddIndepInfoExtCoord 	= 	"logical",
						colInMarker		 	= 	"character",
						strTag				=	"character"
						),
	prototype = prototype(
						strEqcCommand		=	"",
						rcdCriterion		=	"",
						arcdCriterion		=	"",
						acolPval			=	"",
						astrTag				=	"",
						colTag				=	"",
						anumPvalLim			=	1,
						acolIndep			=	"",
						anumIndepLim		=	1,
						strIndepDir			=	"min",
						colInChr			=	"",
						colInPos			=	"",
						numPosLim			=	500000,
						numPosRegionExtension =	-1, # by default will be set to numPosLim/2
						colInGPos			=	"",
						numGPosLim			=	0.1,
						fileRecombRate		=	"",
						colRecombRateChr	=	"chr",
						colRecombRatePos	=	"position",
						colRecombRate		=	"COMBINED_rate.cM.Mb.",
						numRecombRateLim 	= 	20,
						fileClumpBed		= 	"",
						fileClumpSample		= 	"",
						numR2Thrs			= 	0.2,
						numR2PosSize		=	-9,
						numR2VarSize		=	-9,
						blnClumpInSignal	= 	TRUE,
						blnParal			= 	FALSE,
						pathLibLoc			= 	"",
						blnRetainRegionLead = 	FALSE,
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
						#blnAddIndepInfoExtCoord 	= 	FALSE,
						colInMarker			=	"",
						strTag				=	"INDEP"
						)
	#contains = c("EcfReader")
)

setGeneric("setINDEP", function(object) standardGeneric("setINDEP"))
setMethod("setINDEP", signature = (object = "INDEP"), function(object) {
	
	aEqcSlotNamesIn = c("rcdCriterion","arcdCriterion","acolPval","astrTag","colTag","anumPvalLim","colInChr","colInPos","numPosLim","numPosRegionExtension",
						"acolIndep", "anumIndepLim","strIndepDir",
						"fileRecombRate","colRecombRateChr","colRecombRatePos","colRecombRate","numRecombRateLim",
						"fileClumpBed","fileClumpSample","numR2Thrs","blnParal","pathLibLoc","numR2PosSize","numR2VarSize","blnClumpInSignal","blnRetainRegionLead",
						"fileGene","colGeneChr","colGenePosStart","colGenePosStop","colGeneName",
						"fileAnnot","strAnnotTag","colAnnotTag","colAnnotChr","colAnnotPos","colAnnotCoord","numAnnotPosLim",
						"blnAddIndepInfo","colInMarker","strTag",
						"colInGPos","numGPosLim")
	
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
validINDEP <- function(objINDEP) {
	
	### Valid with GWADATA?
	
	# if(objINDEP@rcdCriterion == "") 
		# cat(paste(" EASY WARNING:INDEP\n No criterion rcdCriterion defined. All data will be used for independentisation.", sep=""))
		
	if(objINDEP@acolPval[1] == "" & objINDEP@acolIndep[1] == "") 
		stop(paste(" EASY ERROR:INDEP\n No column acolPval or acolIndep defined. Please set acolPval or acolIndep.", sep=""))
	if(objINDEP@colInChr == "") 
		stop(paste(" EASY ERROR:INDEP\n No column colInChr defined. Please set colInChr.", sep=""))
	if(objINDEP@colInPos == "") 
		stop(paste(" EASY ERROR:INDEP\n No column colInPos defined. Please set colInPos.", sep=""))
	
	# ## check fileRecombRate
	# if(objINDEP@fileRecombRate=="") 
		# stop(paste("EASY ERROR:INDEP\n File fileRecombRate undefined.\n Please set --fileRecombRate or reset --blnRecombRate 0.", sep=""))
	
	if(!objINDEP@fileRecombRate=="") {
		
		if(objINDEP@colInGPos!="") stop("EASY ERROR:INDEP\n You cannot set --colInGPos AND --fileRecombRate. Choose either of the two. ")
		
		if(!file.exists(objINDEP@fileRecombRate))
			stop(paste("EASY ERROR:INDEP\n File fileRecombRate\n ",objINDEP@fileRecombRate,"\n does not exist.\n Please check path or remove --fileRecombRate or reset --blnRecombRate 0.", sep=""))
		
		tr = read.table(objINDEP@fileRecombRate, sep="\t", header=T, stringsAsFactors=F, nrows=1)
			
		isAv = objINDEP@colRecombRateChr %in% names(tr)
		if(!isAv)
			stop(paste(" EASY ERROR:INDEP\n Defined column --colRecombRateChr \n",objINDEP@colRecombRateChr, "\n is not available in recombination rate file \n",objINDEP@fileGene,"\n PLease specify correct column name.", sep=""))
		
		isAv = objINDEP@colRecombRatePos %in% names(tr)
		if(!isAv)
			stop(paste(" EASY ERROR:INDEP\n Defined column --colRecombRatePos \n",objINDEP@colRecombRatePos, "\n is not available in recombination rate file \n",objINDEP@fileGene,"\n PLease specify correct column name.", sep=""))

		isAv = objINDEP@colRecombRate %in% names(tr)
		if(!isAv)
			stop(paste(" EASY ERROR:INDEP\n Defined column --colRecombRate \n",objINDEP@colRecombRate, "\n is not available in recombination rate file \n",objINDEP@fileGene,"\n PLease specify correct column name.", sep=""))
	}
	
	## check fileGene
	if(objINDEP@fileGene != "") {
		if(!file.exists(objINDEP@fileGene))
			stop(paste("EASY ERROR:INDEP\n File fileGene\n ",objINDEP@fileGene,"\n does not exist.\n Please check path or remove --fileGene.", sep=""))
			
		tg = read.table(objINDEP@fileGene, sep="\t", header=T, stringsAsFactors=F, nrows=1)
		
		isAv = objINDEP@colGeneChr %in% names(tg)
		if(!isAv)
			stop(paste(" EASY ERROR:INDEP\n Defined column --colGeneChr \n",objINDEP@colGeneChr, "\n is not available in gene file \n",objINDEP@fileGene,"\n PLease specify correct column name.", sep=""))
		
		isAv = objINDEP@colGenePosStart %in% names(tg)
		if(!isAv)
			stop(paste(" EASY ERROR:INDEP\n Defined column --colGenePosStart \n",objINDEP@colGenePosStart, "\n is not available in gene file \n",objINDEP@fileGene,"\n PLease specify correct column name.", sep=""))

		isAv = objINDEP@colGenePosStop %in% names(tg)
		if(!isAv)
			stop(paste(" EASY ERROR:INDEP\n Defined column --colGenePosStop \n",objINDEP@colGenePosStop, "\n is not available in gene file \n",objINDEP@fileGene,"\n PLease specify correct column name.", sep=""))
	
		isAv = objINDEP@colGeneName %in% names(tg)
		if(!isAv)
			stop(paste(" EASY ERROR:INDEP\n Defined column --colGeneName \n",objINDEP@colGeneName, "\n is not available in gene file \n",objINDEP@fileGene,"\n PLease specify correct column name.", sep=""))
		
	}
	
	
	## check fileAnnot
	if(objINDEP@fileAnnot != "") {
		if(!file.exists(objINDEP@fileAnnot))
			stop(paste("EASY ERROR:INDEP\n File --fileAnnotGene\n ",objINDEP@fileAnnot,"\n does not exist.\n Please check path or remove --fileAnnot.", sep=""))
			
		tk = read.table(objINDEP@fileAnnot, sep="\t", header=T, stringsAsFactors=F, nrows=1)
		
		isUseCoord = objINDEP@colAnnotCoord != ""
		
		if(isUseCoord) {
		
			isAv = objINDEP@colAnnotCoord %in% names(tk)
				if(!isAv)
					stop(paste(" EASY ERROR:INDEP\n Defined column --colAnnotCoord \n",objINDEP@colAnnotCoord, "\n is not available in known loci file \n",objINDEP@fileAnnot,"\n PLease specify correct column name.", sep=""))
		
		} else {
		
			isAv = objINDEP@colAnnotChr %in% names(tk)
			if(!isAv)
				stop(paste(" EASY ERROR:INDEP\n Defined column --colAnnotChr \n",objINDEP@colAnnotChr, "\n is not available in known loci file \n",objINDEP@fileAnnot,"\n PLease specify correct column name.", sep=""))
			
			isAv = objINDEP@colAnnotPos %in% names(tk)
			if(!isAv)
				stop(paste(" EASY ERROR:INDEP\n Defined column --colAnnotPos \n",objINDEP@colAnnotPos, "\n is not available in known loci file \n",objINDEP@fileAnnot,"\n PLease specify correct column name.", sep=""))
			
		}
		
		if(objINDEP@colAnnotTag != "") {
			isAv = objINDEP@colAnnotTag %in% names(tk)
			if(!isAv)
				stop(paste(" EASY ERROR:INDEP\n Defined column --colAnnotTag \n",objINDEP@colAnnotTag, "\n is not available in known loci file \n",objINDEP@fileAnnot,"\n PLease specify correct column name.", sep=""))
		}
	}
	
	if(!objINDEP@pathLibLoc=="") {
		if(!file.exists(objINDEP@pathLibLoc)) {
				stop(paste("EASY ERROR:PCA2STEP\n Path pathLibLoc\n ",objINDEP@pathLibLoc,"\n does not exist.\n Please check path or remove --pathLibLoc.", sep=""))
		}
	}
	
	return(TRUE)
}
INDEP.GWADATA.valid <- function(objINDEP, objGWA) {
	

	if(objINDEP@acolIndep[1] == "") {
		aisMatch = objINDEP@acolPval %in% objGWA@aHeader
		if(any(!aisMatch))
			stop(paste("EASY ERROR:INDEP\n Column \n",objINDEP@acolPval[which(!aisMatch)[1]]," does not exist in file\n",objGWA@fileInShortName,"\n !!!" ,sep=""))
		
		if(length(objINDEP@anumPvalLim)==1 & length(objINDEP@acolPval)>1) objINDEP@anumPvalLim = rep(objINDEP@anumPvalLim[1], length(objINDEP@acolPval))
		
		if(all(objINDEP@astrTag=="")) objINDEP@astrTag = objINDEP@acolPval
		
		if(any(duplicated(objINDEP@astrTag))) {
			stop(paste("EASY ERROR:INDEP\n Created duplicated astrTag \n ",paste(objINDEP@astrTag,collapse=","),"\n Please use --astrTag to set unique tags.", sep=""))
		}

		# if(all(objINDEP@astrTag=="")) objINDEP@astrTag = objINDEP@astrPvalTag
		
		if(length(objINDEP@arcdCriterion)>1&objINDEP@arcdCriterion[1]!="") {
			if(length(objINDEP@acolPval)!=length(objINDEP@arcdCriterion)) {
				stop(paste("EASY ERROR:INDEP\n Length of --acolPval must be identical to length of --arcdCriterion. \n Please correct.", sep=""))
			}
		}
		
	} else {
	
		aisMatch = objINDEP@acolIndep %in% objGWA@aHeader
		if(any(!aisMatch))
			stop(paste("EASY ERROR:INDEP\n Column \n",objINDEP@acolIndep[which(!aisMatch)[1]]," does not exist in file\n",objGWA@fileInShortName,"\n !!!" ,sep=""))	
		
		if(objINDEP@strIndepDir=="max" & length(objINDEP@anumIndepLim)==1 & objINDEP@anumIndepLim[1]==1) {
			objINDEP@anumIndepLim[1]==0
		}		
		
		if(length(objINDEP@anumIndepLim)==1 & length(objINDEP@acolIndep)>1) objINDEP@anumIndepLim = rep(objINDEP@anumIndepLim[1], length(objINDEP@acolIndep))
		
		if(all(objINDEP@astrTag=="")) objINDEP@astrTag = objINDEP@acolIndep
		
		if(any(duplicated(objINDEP@astrTag))) {
			stop(paste("EASY ERROR:INDEP\n Created duplicated astrTag \n ",paste(objINDEP@astrTag,collapse=","),"\n Please use --astrTag to set unique tags.", sep=""))
		}
		
		#if(all(objINDEP@astrTag=="")) objINDEP@astrTag = objINDEP@astrIndepTag
		
		if(length(objINDEP@arcdCriterion)>1&objINDEP@arcdCriterion[1]!="") {
			if(length(objINDEP@acolIndep)!=length(objINDEP@arcdCriterion)) {
				stop(paste("EASY ERROR:INDEP\n Length of --acolIndep must be identical to length of --arcdCriterion. \n Please correct.", sep=""))
			}
		}
		
	}
	
	isAv <- objINDEP@colInChr %in% objGWA@aHeader
	if(!isAv)
		stop(paste(" EASY ERROR:INDEP\n Defined column colInChr \n",objINDEP@colInChr, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	isAv <- objINDEP@colInPos %in% objGWA@aHeader
	if(!isAv)
		stop(paste(" EASY ERROR:INDEP\n Defined column colInPos \n",objINDEP@colInPos, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	if(objINDEP@colInGPos!="") {
		isAv <- objINDEP@colInGPos %in% objGWA@aHeader
		if(!isAv)
			stop(paste(" EASY ERROR:INDEP\n Defined column colInGPos \n",objINDEP@colInGPos, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))	
	}
	
	isAv <- objINDEP@colTag %in% objGWA@aHeader
	if(objINDEP@colTag!="" & !(objINDEP@colTag %in% objGWA@aHeader))
		stop(paste(" EASY ERROR:INDEP\n Defined column colTag \n",objINDEP@colTag, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name or remove --colTag.", sep=""))
	

	if(any(duplicated(objINDEP@astrTag))) {
		stop(paste("EASY ERROR:INDEP\n Created duplicated astrTag \n ",paste(objINDEP@astrTag,collapse=","),"\n Please use --astrTag to set unique tags.", sep=""))
	}

	
	if(objINDEP@numPosRegionExtension == -1) objINDEP@numPosRegionExtension = objINDEP@numPosLim/2
	
	if(!objINDEP@fileClumpBed=="") {
	
		if(grepl("<CHR>",objINDEP@fileClumpBed,fixed=T)) {
			fileClumpBed = c()
			for(chr in c(1:22,"X")) fileClumpBed = c(fileClumpBed, gsub("<CHR>",chr,objINDEP@fileClumpBed,fixed=T))
			objINDEP@fileClumpBed = fileClumpBed
		} 
		for(k in 1:length(objINDEP@fileClumpBed)) {
		#for(fileClumpBedi in objINDEP@fileClumpBed) {
			fileClumpBedi = objINDEP@fileClumpBed[k]
			if(k<23 & !file.exists(paste0(fileClumpBedi,".bed"))) {
				stop(paste("EASY ERROR:INDEP\n File fileClumpBed\n ",fileClumpBedi,"\n does not exist.\n Please check path or remove --fileClumpBed.", sep=""))
			}
		}
	}
	
	if(!objINDEP@fileClumpSample=="") {
		if(!file.exists(objINDEP@fileClumpSample)) {
				stop(paste("EASY ERROR:INDEP\n File fileClumpSample\n ",objINDEP@fileClumpSample,"\n does not exist.\n Please check path or remove --fileClumpBed.", sep=""))
		}
	}
	
	### set Indep vals from Pvals 
	if(objINDEP@acolIndep[1] == "") {
		objINDEP@acolIndep = objINDEP@acolPval
		#objINDEP@astrIndepTag = objINDEP@astrPvalTag
		objINDEP@anumIndepLim = objINDEP@anumPvalLim
	}
	
	#if(blnAddIndepInfoExtCoord) objINDEP@blnAddIndepInfo <- TRUE
	
	return(objINDEP)
	
}

INDEP.read <- function(objINDEP, blnReadAll) {
	
	fR = objINDEP@fileRecombRate
	
	if(blnReadAll) {
		tblR = read.table(fR, header=T, stringsAsFactors=FALSE, strip.white = TRUE, comment.char = "")
	} else {
		tblR = read.table(fR, header=T, stringsAsFactors=FALSE, strip.white = TRUE, comment.char = "", nrows = 10)
	}
	
	tblRShort = tblR[,c(objINDEP@colRecombRateChr,objINDEP@colRecombRatePos,objINDEP@colRecombRate)]
	names(tblRShort) <- c("chr","pos","recomb_rate_cM_Mb")
	
	###sort by chr and pos
	tblRShortOrdered = tblRShort[order(as.integer(tblRShort$chr),as.integer(tblRShort$pos)),]
	
	return(tblRShortOrdered)
}

#############################################################################################################################
INDEP.run <- function(objINDEP, objGWA, objREPORT, tblRR, isValidScript) {
	#objINDEP, objGWA, objREPORT, tblRRX, isValidScript
	
	rcdCriterion 	<- objINDEP@rcdCriterion
	arcdCriterion 	<- objINDEP@arcdCriterion
	astrTag 		<- objINDEP@astrTag
	colTag 			<- objINDEP@colTag
	## astrTag is tag for arcdCriterion
	
	# deprecated: 
	# acolPval		<- objINDEP@acolPval
	# astrPvalTag 	<- objINDEP@astrPvalTag
	# anumPvalLim		<- objINDEP@anumPvalLim
	
	acolIndep 		<- objINDEP@acolIndep
	#astrIndepTag 	<- objINDEP@astrIndepTag
	anumIndepLim 	<- objINDEP@anumIndepLim
	
	strIndepDir		<- objINDEP@strIndepDir
	
	colInChr		<- objINDEP@colInChr
	colInPos		<- objINDEP@colInPos
	numPosLim		<- objINDEP@numPosLim
	numPosRegionExtension <- objINDEP@numPosRegionExtension
	
	## signaling by genetic position:
	colInGPos		<- objINDEP@colInGPos
	numGPosLim		<- objINDEP@numGPosLim
	
	## signaling by recomb rate file: 
	fileRecombRate <- objINDEP@fileRecombRate
	numRecombRateLim <- objINDEP@numRecombRateLim
	
	## clumping to loci by r2:
	fileClumpBed <- objINDEP@fileClumpBed
	fileClumpSample <- objINDEP@fileClumpSample
	numR2Thrs 	 <- objINDEP@numR2Thrs
	blnParal <- objINDEP@blnParal
	
	pathLibLoc	<- objINDEP@pathLibLoc
	
	numR2PosSize <- objINDEP@numR2PosSize
	numR2VarSize <- objINDEP@numR2VarSize
	blnClumpInSignal <- objINDEP@blnClumpInSignal
	## default is to clump within a cM signal
	blnRetainRegionLead <- objINDEP@blnRetainRegionLead
	
	fileGene 	<- objINDEP@fileGene
	colGeneChr 	<- objINDEP@colGeneChr
	colGenePosStart <- objINDEP@colGenePosStart
	colGenePosStop 	<- objINDEP@colGenePosStop
	colGeneName 	<- objINDEP@colGeneName
	
	fileAnnot 	<- objINDEP@fileAnnot
	strAnnotTag 	<- objINDEP@strAnnotTag
	colAnnotTag 	<- objINDEP@colAnnotTag
	colAnnotChr 	<- objINDEP@colAnnotChr
	colAnnotPos 	<- objINDEP@colAnnotPos
	colAnnotCoord  	<- objINDEP@colAnnotCoord
	numAnnotPosLim 	<- objINDEP@numAnnotPosLim
	
	blnAddIndepInfo	<- objINDEP@blnAddIndepInfo
	#blnAddIndepInfoExtCoord	<- objINDEP@blnAddIndepInfoExtCoord
	colInMarker 	<- objINDEP@colInMarker
	strTag			<- objINDEP@strTag
	
	#### 
	
	if(nchar(strTag)>0) strTag = paste(strTag,".",sep="") 
	
	ncrit = 1
	if(arcdCriterion[1]!="") {
		ncrit = length(arcdCriterion)
		if(ncrit == 1) {
			rcdCriterion = arcdCriterion
		}
	}
	
	blnRecombRate = !(fileRecombRate=="")
	blnGene = !(fileGene=="")
	blnAnnot = !(fileAnnot=="")
	blnClump = any(!(fileClumpBed==""))
	blnGpos = colInGPos!=""
	blnIndepMin = strIndepDir == "min"
	
	if(!blnRecombRate & !blnGpos) blnClumpInSignal = FALSE
	
	if(numAnnotPosLim == -1) numAnnotPosLim = numPosLim
	### 
	
	if(rcdCriterion != "") {
		objRCD 	<- RCD(rcdCriterion)
		out 	<- RCD.eval(objRCD, objGWA)
		out[is.na(out)] <- FALSE
		numIndepCrit = length(which(out))
		objGWA.indep <- GWADATA.getrows(objGWA, which(out))
	} else if(arcdCriterion[1]!="") {
		tCrit = data.frame()
		out = rep(FALSE,nrow(objGWA@tblGWA))
		for(i in 1:length(arcdCriterion)) {
			objRCDi 	<- RCD(arcdCriterion[i])
			outi 	<- RCD.eval(objRCDi, objGWA)
			outi[is.na(outi)] <- FALSE
			if(i==1) {
				tCrit = data.frame(outi)
			} else {
				tCrit = cbind(tCrit, outi)
			}
			names(tCrit)[ncol(tCrit)] = astrTag[i]
			out = out | outi
		}
		numIndepCrit = length(which(out))
		objGWA.indep <- GWADATA.getrows(objGWA, which(out))
		tCrit = tCrit[which(out),]
		
	} else {
		numIndepCrit <- nrow(objGWA@tblGWA)
		objGWA.indep <- objGWA
	}
	
	objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numIndepIn",sep=""),NA)
	objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numPosMissing",sep=""),NA)
	for(i in 1:length(acolIndep)) objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numVarSignif.",astrTag[i],sep=""),NA)
	objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numRegion",sep=""),NA)
	if(blnRecombRate | blnGpos) objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numSignal",sep=""),NA)
	if(blnClump) {
		objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numLoci",sep=""),NA)
		objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numLociMiss",sep=""),NA)
	}
	
	if(numIndepCrit == 0) {
		
		objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numIndepIn",sep=""),0)
		objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numPosMissing",sep=""),0)
		for(i in 1:length(acolIndep)) objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numVarSignif.",astrTag[i],sep=""),0)
		objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numRegion",sep=""),0)
		if(blnRecombRate | blnGpos) objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numSignal",sep=""),0)
		if(blnClump) {
			objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numLoci",sep=""),0)
			objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numLociMiss",sep=""),0)
		}
		
		if(blnAddIndepInfo) {
			
			# objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"pMin",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionId",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionLead",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionTag",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionNumVariants",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionCoordinates",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionSize",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionExtCoordinates",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionExtSize",sep=""))
			## for all (also non-signif variants): TRUE/FALSE value
			# objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"inregionExt",sep=""))
			
			if(blnRecombRate | blnGpos) {
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionNumSignals",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalId",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalLead",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalTag",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalNumVariants",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalCoordinates",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalSize",sep=""))
			}
			if(blnClump) {
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionNumLoci",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusId",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusLead",sep=""))
				# objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusMiss",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusTag",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusNumVariants",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusCoordinates",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusSize",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusR2",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusNumMiss",sep=""))
			}
			
			if(blnAnnot) {
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionAnnot",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"annotDistance",sep=""))
				if(blnRecombRate | blnGpos) objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalAnnot",sep=""))
			}
			
			if(blnGene) {
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"NearestGene",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"NearestGeneDistance",sep=""))
			}
			
		}
		
		return(list(objGWA,objGWA.indep,objGWA.indep,objGWA.indep,objGWA.indep,objGWA.indep,objREPORT))
	}
	
	objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numIndepIn",sep=""),numIndepCrit)
	
	
	##############
	### remove missing positions
	
	isNaPos = is.na(GWADATA.getcol(objGWA.indep, colInChr)) | is.na(GWADATA.getcol(objGWA.indep, colInPos))
	numPosMissing = length(which(isNaPos))

	objGWA.indep@tblGWA <- 	objGWA.indep@tblGWA[!isNaPos,]


	objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numPosMissing",sep=""),numPosMissing)

	
	##############
	### reduce to significant variants and count significant variants by tag / criterion
	
	isSignif = rep(FALSE,nrow(objGWA.indep@tblGWA))
	
	if(rcdCriterion != "") {
		tCrit = data.frame()
		for(i in 1:length(acolIndep)) { 
			avalTmp = GWADATA.getcol(objGWA.indep, acolIndep[i])
			if(blnIndepMin) {
				isSignifTmp = avalTmp < anumIndepLim[i] & !is.na(avalTmp)
			} else {
				isSignifTmp = avalTmp > anumIndepLim[i] & !is.na(avalTmp)
			}
			if(i == 1) {
				tCrit = data.frame(isSignifTmp)
			} else {
				tCrit = cbind(tCrit, isSignifTmp)
			}
			names(tCrit)[ncol(tCrit)] = astrTag[i]
			objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numVarSignif.",astrTag[i],sep=""),length(which(isSignifTmp)))
			isSignif = isSignif | isSignifTmp
		}
	} else {
		tCrit = tCrit[!isNaPos,]
		for(i in 1:length(arcdCriterion)) { 
	
			isSignifTmp = tCrit[,i]
			objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numVarSignif.",astrTag[i],sep=""),length(which(isSignifTmp)))
			isSignif = isSignif | isSignifTmp
		}
	}
	
	if(all(!isSignif)) {
		
		objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numRegion",sep=""),0)
		if(blnRecombRate | blnGpos) objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numSignal",sep=""),0)
		if(blnClump) {
			objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numLoci",sep=""),0)
			objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numLociMiss",sep=""),0)
		}
		
		if(blnAddIndepInfo) {
			
			# objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"pMin",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionId",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionLead",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionTag",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionNumVariants",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionCoordinates",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionSize",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionExtCoordinates",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionExtSize",sep=""))
			
			if(blnRecombRate | blnGpos) {
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionNumSignals",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalId",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalLead",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalTag",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalNumVariants",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalCoordinates",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalSize",sep=""))
			}
			if(blnClump) {
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionNumLoci",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusId",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusLead",sep=""))
				# objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusMiss",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusTag",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusNumVariants",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusCoordinates",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusSize",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusR2",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"locusNumMiss",sep=""))
				
			}
			
			if(blnAnnot) {
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"regionAnnot",sep=""))
				objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"annotDistance",sep=""))
				if(blnRecombRate | blnGpos) objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"signalAnnot",sep=""))
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
	
	if(length(acolIndep)>1) {
		if(blnIndepMin) {
			aValext = apply(tIn[,acolIndep],1,min,na.rm=T)
			aidxsort = order(aValext)
		} else {
			aValext = apply(tIn[,acolIndep],1,max,na.rm=T)
			aidxsort = order(aValext,decreasing=TRUE)
		}
	} else {
		aValext = tIn[,acolIndep]
		aidxsort = order(aValext)
	}
	
	tInSort = tIn[aidxsort,]
	# aPminSort = aPmin[aidxsort]
	aValextSort = aValext[aidxsort]
	tCritSort = as.data.frame(tCrit[aidxsort,])
	names(tCritSort) = names(tCrit)
	
	aChr = tInSort[,colInChr]
	aPos = tInSort[,colInPos]
	
	if(blnGpos) {
		aGPos = tInSort[,colInGPos]
	} else {
		aGPos = rep(NA,nrow(tInSort))
	}
	
	
	### work with aChr, aPos, aPminSort and tInSort[,acolIndep] to obtain : 
	
	####################################
	### 1. Obtain regionId, regionLead
	
	regionId <- regionLead <- regionCoordinates <- regionSize <- regionExtCoordinates <- regionExtSize <- rep(NA, nrow(tInSort))
	
	regioncount = 1
	
	# aPminSortBackup = aPminSort
	aValextSortBackup = aValextSort
	
	while(any(is.na(regionId))) {
			
		# print(regioncount)
		
		# iTmpExtr = which(aPminSort == min(aPminSort))[1]
		
		if(blnIndepMin) {
			iTmpExtr = which(aValextSort == min(aValextSort))[1]
		} else {
			iTmpExtr = which(aValextSort == max(aValextSort))[1]
		}
		
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
		# regionCoordinates[isCurRegionNew] = paste(chrExtr,":",regionPos1-numPosLim,"_",regionPos2+numPosLim,sep="")
		# INDEP10 update: It is important to increase the regions by numPosLim/2 because otherwise regions will overlap !
		regionCoordinates[isCurRegionNew] = paste(chrExtr,":",regionPos1,"_",regionPos2,sep="")
		regionSize[isCurRegionNew] = regionPos2 - regionPos1 + 1
		regionExtCoordinates[isCurRegionNew] = paste(chrExtr,":",regionPos1-numPosRegionExtension,"_",regionPos2+numPosRegionExtension,sep="")
		regionExtSize[isCurRegionNew] = (regionPos2+numPosRegionExtension) - (regionPos1-numPosRegionExtension)
		
		regioncount = regioncount + 1
		
		if(blnIndepMin) {
			# aPminSort[isCurRegionNew] = Inf
			aValextSort[isCurRegionNew] = Inf
		} else {
			aValextSort[isCurRegionNew] = -Inf
		}
		
	}
	
	# aPminSort = aPminSortBackup
	aValextSort = aValextSortBackup
	
	####################################
	### 2. Obtain signalId, signalLead, signalCoordinates
	if(blnGpos) {
		## repeat on genetic signal positions within regions
		signalId <- signalLead <- signalCoordinates <- signalSize <- rep(NA, nrow(tInSort))
		
		for(rid in unique(regionId)) {
		
			isRegion = regionId == rid
		
			signalcount = 1
		
			# aPminSortBackup = aPminSort
			aValextSortBackup = aValextSort
			
			while(any(is.na(signalId[isRegion]))) {
					
				# print(signalcount)
				
				# iTmpExtr = which(aPminSort == min(aPminSort))[1]
				
				#iTmpExtr = which(aPminSort == min(aPminSort[isRegion]) & isRegion)[1]
				
				if(blnIndepMin) {
					iTmpExtr = which(aValextSort == min(aValextSort[isRegion]) & isRegion)[1]
				} else {
					iTmpExtr = which(aValextSort == max(aValextSort[isRegion]) & isRegion)[1]
				}
				
				
				chrExtr =  aChr[iTmpExtr]
				posExtr =  aGPos[iTmpExtr]
				
				### use position mapping to obtain signals
				
				isCurSignal = NA
				isCurSignalNew = aChr==chrExtr & aGPos>=(posExtr-numGPosLim) & aGPos<=(posExtr+numGPosLim) & isRegion
				
				while(!identical(isCurSignal,isCurSignalNew)) {
					
					isCurSignal = isCurSignalNew
					
					signalGPos1 = min(aGPos[isCurSignal])
					signalGPos2 = max(aGPos[isCurSignal])
					
					signalPos1 = min(aPos[isCurSignal])
					signalPos2 = max(aPos[isCurSignal])
					
					isCurSignalNew = aChr==chrExtr & aGPos>=(signalGPos1-numGPosLim) & aGPos<=(signalGPos2+numGPosLim) & isRegion
					
				}
				
				signalId[isCurSignalNew] = signalcount
				signalLead[iTmpExtr] = signalcount

				signalCoordinates[isCurSignalNew] = paste(chrExtr,":",signalPos1,"_",signalPos2,sep="")
				signalSize[isCurSignalNew] = signalPos2 - signalPos1 + 1

				signalcount = signalcount + 1
				
				#aPminSort[isCurSignalNew] = Inf
				
				if(blnIndepMin) {
					aValextSort[isCurSignalNew] = Inf
				} else {
					aValextSort[isCurSignalNew] = -Inf
				}
				
			}
			
			# aPminSort = aPminSortBackup
			aValextSort = aValextSortBackup
		}
		
		signalId = paste(regionId,signalId,sep=".")
		signalLead = ifelse(!is.na(signalLead), signalId, NA)
		
	}
	
	
	if(blnRecombRate) {
	
		signalId <- signalLead <- signalCoordinates <- signalSize <- rep(NA, nrow(tInSort))
		
		tblRR = tblRR[order(tblRR$chr,tblRR$pos),]
		
		for(rid in unique(regionId)) {
			
			# if(rid==82) stop()
			
			isRegion = regionId == rid
			
			signalcount = 1
			
			# aPminSortBackup = aPminSort
			aValextSortBackup = aValextSort
			
			while(any(is.na(signalId[isRegion]))) {
				
				# print(paste(rid,signalcount,sep="."))
				
				# iTmpExtr = which(aPminSort == min(aPminSort[isRegion]) & isRegion)[1]
				
				if(blnIndepMin) {
					iTmpExtr = which(aValextSort == min(aValextSort[isRegion]) & isRegion)[1]
				} else {
					iTmpExtr = which(aValextSort == max(aValextSort[isRegion]) & isRegion)[1]
				}
				
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
					signalSize[isCurSignal] = max(aPos[isCurSignal])-pos1+1
				} else {
					signalSize[isCurSignal] = pos2-pos1+1
				}
				
				signalcount = signalcount + 1
				
				# aPminSort[isCurSignal] = Inf
				
				if(blnIndepMin) {
					aValextSort[isCurSignal] = Inf
				} else {
					aValextSort[isCurSignal] = -Inf
				}
				
			}
			
			# aPminSort = aPminSortBackup
			aValextSort = aValextSortBackup
			
		}
		
		signalCoordinates = gsub("Inf","EoCHR",signalCoordinates)
		
		signalId = paste(regionId,signalId,sep=".")
		signalLead = ifelse(!is.na(signalLead), signalId, NA)
	}
	
	if(blnClump) {
		
		achruni = unique(aChr)
		
		lsidxsample = list()
		
		for(k in 1:length(achruni)) {
			chrk = achruni[k]
			
			tfamk = fread(paste0(fileClumpBed[as.integer(ifelse(chrk=="X",23,chrk))],".fam"),header=F)
			
			# tfam1 = fread(paste0(fileClumpBed[1],".fam"),header=F)
			
			if(fileClumpSample!="") {
				asampleused = scan(fileClumpSample,sep="\n",what="character",quiet=TRUE)
				#aidxsample = which(tfam1$V2%in%asampleused)
				lsidxsample[[k]] = which(tfamk$V2%in%asampleused)
			} else {
				#aidxsample = seq(1,nrow(tfam1))
				lsidxsample[[k]] = seq(1,nrow(tfamk))
			}
		}
		# locusMiss <- rep(0, nrow(tInSort))
		
		lsClump = list()
		
		## lsClump[[k]] = fnClumpChr(chrk, aChr, aPos, tInSort, colInMarker, areaId, aValextSort, fileClumpBed, aidxsample, numR2PosSize, numR2VarSize, numR2Thrs, blnIndepMin)
		
		# fnClumpChr <- function(chrk, aChr, aPos, tInSort, colInMarker, regionId, aValextSort, fileClumpBed, aidxsample, numR2PosSize, numR2VarSize, numR2Thrs, blnIndepMin, blnRetainRegionLead) {
		fnClumpChr <- function(chrk, aChr, aPos, tInSort, colInMarker, tArea, aValextSort, fileClumpBed, aidxsample, numR2PosSize, numR2VarSize, numR2Thrs, blnIndepMin, blnRetainRegionLead) {
			
			idxk = which(aChr==chrk)
			
			aChrk = aChr[idxk]
			aPosk = aPos[idxk]
			tInSortk = tInSort[idxk,]
			aSnpsk = tInSortk[[colInMarker]]
			
			regionIdk = tArea[[1]][idxk]
			regionLeadk = tArea[[2]][idxk]
			regionCoordinatesk = tArea[[3]][idxk]
			regionSizek = tArea[[4]][idxk]
			
			# aPminSortk = aPminSort[idxk]
			aValextSortk = aValextSort[idxk]
			
			locusIdk <- locusLeadk <- locusCoordinatesk <- locusSizek <- locusR2k <- locusNumMissk <- rep(NA, nrow(tInSortk))
			
			if(length(fileClumpBed)>1) {
				tbimk = fread(paste0(fileClumpBed[as.integer(ifelse(chrk=="X",23,chrk))],".bim"),header=F)
			} else {
				tbimk = fread(paste0(fileClumpBed,".bim"),header=F)
			}
			
			
			for(m in 1:length(unique(regionIdk))) {
			# for(rid in unique(regionIdk)) {
				rid = unique(regionIdk)[m]
				
				#if(rid==6) stop()
				
				isRegion = regionIdk == rid
				
				print(paste("Clumping chromsome",chrk,"regions:",m,"of",length(unique(regionIdk))))
				
				# aPminSortBackup = aPminSortk
				aValextSortBackup = aValextSortk
				
				chrregion = unique(aChrk[isRegion])
				
				if(length(chrregion)>1) stop("Region with multiple chrs.")
				
				# if(length(fileClumpBed)>1) {
					# tbimchr = lsbim[[as.integer(chrregion)]]
				# } else {
					# tbimchr = lsbim[[1]]
				# }
				tbimchr = tbimk 
				
				aSnpsRegion = aSnpsk[isRegion]
				aidxbimused = which(tbimchr$V2%in%aSnpsRegion)
				
				if(length(aidxbimused)==0) {
					## no region SNP is in reference panel
					if(blnRetainRegionLead) {
						locusNumMissk[isRegion] = length(aSnpsRegion)
						locusIdk[isRegion] = 0
						locusLeadk[isRegion & !is.na(regionLeadk)] = 0
						locusCoordinatesk[isRegion] = regionCoordinatesk[isRegion]
						locusSizek[isRegion] = regionSizek[isRegion]
						
					} else {
						locusNumMissk[isRegion] = length(aSnpsRegion)
						# locusId[isRegion] = (-1)*seq(1,sum(isRegion))
						locusIdk[isRegion] = -9
						# locusMiss[isRegion] = 1
					}
					
				} else if(length(aidxbimused)==1) {
					
					locusNumMissk[isRegion] = length(unique(aSnpsRegion)) - 1
					isOneSnp = which(aSnpsk==tbimchr$V2[aidxbimused])
					locusIdk[isOneSnp] = 1
					locusLeadk[isOneSnp[1]] = 1
					locusCoordinatesk[isOneSnp] = paste(chrregion,":",aPosk[isOneSnp],"_",aPosk[isOneSnp],sep="")
					locusSizek[isOneSnp] = 1
					locusR2k[isOneSnp] = 1
					
				} else {
					
					if(length(fileClumpBed)>1) { 
						fileClumpBedused = fileClumpBed[as.integer(ifelse(chrregion=="X",23,chrregion))]
					} else {
						fileClumpBedused = fileClumpBed
					}
					
					# matcorused = bed_cor(bed(paste0(fileClumpBedused,".bed")), ind.col=aidxbimused, ind.row = aidxsample, size=length(aidxbimused))
					# matcorused = bed_cor(bed(paste0(fileClumpBedused,".bed")), ind.col = aidxbimused, ind.row = aidxsample, size = length(aidxbimused))
					
					if(numR2PosSize!=-9) {
						matcorused = bed_cor(bed(paste0(fileClumpBedused,".bed")), ind.col = aidxbimused, ind.row = aidxsample, infos.pos = tbimchr$V4[aidxbimused], size = floor(numR2PosSize/1000))
					} else if(numR2VarSize!=-9) {
						matcorused = bed_cor(bed(paste0(fileClumpBedused,".bed")), ind.col = aidxbimused, ind.row = aidxsample, size = numR2VarSize)
					} else {
						matcorused = bed_cor(bed(paste0(fileClumpBedused,".bed")), ind.col = aidxbimused, ind.row = aidxsample, size = length(aidxbimused))
					}
					
					# t1 = Sys.time()
					# matcorused = bed_cor(bed(paste0(fileClumpBedused,".bed")), ind.col = aidxbimused, ind.row = aidxsample, size = 100)
					# t2 = Sys.time()
					# t2-t1
					
					# t1 = Sys.time()
					# matcorused = bed_cor(bed(paste0(fileClumpBedused,".bed")), ind.col = aidxbimused, ind.row = aidxsample, size = 200)
					# t2 = Sys.time()
					# t2-t1
					
					# t1 = Sys.time()
					# matcorused = bed_cor(bed(paste0(fileClumpBedused,".bed")), ind.col = aidxbimused, ind.row = aidxsample, size = 400)
					# t2 = Sys.time()
					# t2-t1
					
					# t1 = Sys.time()
					# matcorused = bed_cor(bed(paste0(fileClumpBedused,".bed")), ind.col = aidxbimused, ind.row = aidxsample, size = 800)
					# t2 = Sys.time()
					# t2-t1
					
					# apos = tbimchr$V4[aidxbimused]
					
					# an500=rep(NA,length(apos))
					# for(i in 1:length(apos)) an500[i]=sum(abs(apos[i]-apos)<250000)
					
					
					# numR2PosSize
					# infos.pos = NULL,

					aisMissMatrix = apply(matcorused,1,function(x) all(is.na(x)|x==1)&(sum(x==1,na.rm=T)==1))
					if(any(aisMissMatrix)) {
						aidxbimused = aidxbimused[which(!aisMissMatrix)]
						matcorused = matcorused[which(!aisMissMatrix),which(!aisMissMatrix)]
					}
					
					tbimregion = tbimchr[aidxbimused,]
					
					## set NAs to -1 in locusId
					aisMiss = !aSnpsRegion%in%tbimregion$V2
					if(any(aisMiss)) {
						aSnpsRegionMiss = aSnpsRegion[aisMiss]
						isSetMinus = aSnpsk%in%aSnpsRegionMiss
						locusIdk[isSetMinus] = -9
						# locusId[isSetMinus] = (-1)*seq(1,sum(isSetMinus))  ## negative locusId indicates missing in bed 
						locusNumMissk[isRegion] = sum(isSetMinus)
						# locusMiss[isSetMinus] = 1
						
						# aPminSortk[isSetMinus] = Inf
						if(blnIndepMin) {
							aValextSortk[isSetMinus] = Inf
						} else {
							aValextSortk[isSetMinus] = -Inf
						}
						
					} else {
						# isSetMinus = rep(FALSE,length(aPminSortk))
						isSetMinus = rep(FALSE,length(aValextSortk))
						locusNumMissk[isRegion] = 0
					}
					
					##
					## "16:20383049:C_G"
					
					locuscount = 1
					
					matcorused = matcorused^2
					
					while(any(is.na(locusIdk[isRegion]))) {
						
						# print(paste(rid,locuscount,sep="."))
						
						# iTmpExtr = which(aPminSortk == min(aPminSortk[isRegion & !isSetMinus]) & isRegion & !isSetMinus)[1]
						
						if(blnIndepMin) {
							iTmpExtr = which(aValextSortk == min(aValextSortk[isRegion & !isSetMinus]) & isRegion & !isSetMinus)[1]	
						} else {
							iTmpExtr = which(aValextSortk == max(aValextSortk[isRegion & !isSetMinus]) & isRegion & !isSetMinus)[1]	
						}
						
						snpExtr = aSnpsk[iTmpExtr]
						
						idxbimleadsnp = which(tbimregion$V2 == snpExtr)
						
						isbimlocus = matcorused[idxbimleadsnp,] >= numR2Thrs
						
						asnpLocus = tbimregion$V2[which(isbimlocus)]
						
						isCurLocus = aSnpsk%in%asnpLocus & is.na(locusIdk) & isRegion
						
						locusPos1 = min(aPosk[which(isCurLocus)])
						locusPos2 = max(aPosk[which(isCurLocus)])
						locusCoordinatesk[which(isCurLocus)] = paste(chrregion,":",locusPos1,"_",locusPos2,sep="")
						locusSizek[which(isCurLocus)] = locusPos2-locusPos1+1
						
						locusIdk[which(isCurLocus)] = locuscount
						locusLeadk[iTmpExtr] = locuscount
						
						aidxmatch2bim = match(aSnpsk[which(isCurLocus)],tbimregion$V2)
						locusR2k[which(isCurLocus)] = matcorused[idxbimleadsnp,aidxmatch2bim]
						
						locuscount = locuscount + 1
						
						# aPminSortk[which(isCurLocus)] = Inf
						
						if(blnIndepMin) {
							aValextSortk[which(isCurLocus)] = Inf
						} else {
							aValextSortk[which(isCurLocus)] = -Inf
						}
						
					}
					
					#aPminSortk = aPminSortBackup
					aValextSortk = aValextSortBackup
				}
			}
			
			return(data.frame(idxk, locusIdk, locusLeadk, locusCoordinatesk, locusSizek, locusR2k, locusNumMissk))
		}
		
		nchruni = length(achruni)
		if(blnClumpInSignal) {
			# areaId = signalId
			tarea = data.frame(signalId, signalLead, signalCoordinates, signalSize, stringsAsFactors=F)
		} else {
			#areaId = regionId
			
			tarea = data.frame(regionId, regionLead, regionCoordinates, regionSize, stringsAsFactors=F)
		}
		
		if(blnParal) {
			cl<-makeCluster(nchruni)
			registerDoParallel(cl) 
			#clusterCall(cl, function(x) .libPaths(x), .libPaths())
			
			if(pathLibLoc == "") {
				clusterCall(cl, function(x) .libPaths(x), .libPaths())
			} else {
				clusterCall(cl, function(x) .libPaths(x), c(.libPaths(),pathLibLoc))
			}
			
			lsClump = foreach(k=1:nchruni,.packages=c("bigsnpr","data.table")) %dopar% {
#				fnCor(filePcaBed[as.integer(achruni[i])], asnps[which(achr==achruni[i])], aidxsample, numR2PosSize)
				chrk = achruni[k]
				aidxsample = lsidxsample[[k]]
				# fnClumpChr(chrk, aChr, aPos, tInSort, colInMarker, areaId, aPminSort, fileClumpBed, aidxsample, numR2PosSize, numR2VarSize, numR2Thrs, blnIndepMin)
				# fnClumpChr(chrk, aChr, aPos, tInSort, colInMarker, areaId, aValextSort, fileClumpBed, aidxsample, numR2PosSize, numR2VarSize, numR2Thrs, blnIndepMin, blnRetainRegionLead)
				fnClumpChr(chrk, aChr, aPos, tInSort, colInMarker, tarea, aValextSort, fileClumpBed, aidxsample, numR2PosSize, numR2VarSize, numR2Thrs, blnIndepMin, blnRetainRegionLead)
			}
			stopCluster(cl)
		} else {
			for(k in 1:length(achruni)) {
				chrk = achruni[k]
				aidxsample = lsidxsample[[k]]
				# lsClump[[k]] = fnClumpChr(chrk, aChr, aPos, tInSort, colInMarker, areaId, aPminSort, fileClumpBed, aidxsample, numR2PosSize, numR2VarSize, numR2Thrs, blnIndepMin)
				# lsClump[[k]] = fnClumpChr(chrk, aChr, aPos, tInSort, colInMarker, areaId, aValextSort, fileClumpBed, aidxsample, numR2PosSize, numR2VarSize, numR2Thrs, blnIndepMin, blnRetainRegionLead)
				lsClump[[k]] = fnClumpChr(chrk, aChr, aPos, tInSort, colInMarker, tarea, aValextSort, fileClumpBed, aidxsample, numR2PosSize, numR2VarSize, numR2Thrs, blnIndepMin, blnRetainRegionLead)
			}			
		}
		
		
		locusId <- locusLead <- locusCoordinates <- locusSize <- locusR2 <- locusNumMiss <- rep(NA, nrow(tInSort))
		
		for(k in 1:length(achruni)) {
			tclump = lsClump[[k]]
			locusId[tclump$idxk] = tclump$locusId
			locusLead[tclump$idxk] = tclump$locusLead
			locusCoordinates[tclump$idxk] = tclump$locusCoordinates
			locusSize[tclump$idxk] = tclump$locusSize
			locusR2[tclump$idxk] = tclump$locusR2
			locusNumMiss[tclump$idxk] = tclump$locusNumMiss
		}
		
		locusId[which(locusId==-9)] = NA
			
		# locusId = paste(areaId,locusId,sep=".") ## 1.NA, for variants missing in bed
		locusId = paste(tarea$regionId,locusId,sep=".") ## 1.NA, for variants missing in bed
		
		#locusLead = ifelse(is.na(locusLead), NA, paste(areaId,locusLead,sep="."))  # only set for real locusLead
		locusLead = ifelse(is.na(locusLead), NA, paste(tarea$regionId,locusLead,sep="."))  # only set for real locusLead
		
	}
	
	####################################
	### 3. Obtain signalTag and regionTag 

	regionTag <- regionNumVariants <- regionNumSignals <- regionNumLoci <- rep(NA, nrow(tInSort))
	
	for(rid in unique(regionId)) {
		
		# print(rid)
		
		isRegion = regionId == rid
		
		#numvar = length(which(isRegion))
		numvaruni = length(unique(tInSort[[colInMarker]][isRegion]))
		
		if(blnRecombRate | blnGpos) numsig = length(unique(signalId[isRegion]))
		if(blnClump) numloci = length(unique(locusId[isRegion & !grepl("NA",locusId)]))
		
		artag = c()
		arnumvar = c()
		arnumsig = c()
		arnumloc = c()
		
		if(colTag == "") {
		
			for(i in 1:length(acolIndep)) { 
				
				# colp = acolIndep[i]
				#plim = anumIndepLim[i]
				#tag = astrTag[i]
				
				# isRidSignif = isRegion & tInSort[,colp]<plim
				isRidSignif = isRegion & tCritSort[,i]
				
				if(any(isRidSignif)) { 
					artag = c(artag, astrTag[i])
					arnumvar = c(arnumvar, length(which(isRidSignif)))
					if(blnRecombRate | blnGpos) arnumsig = c(arnumsig, length(unique(signalId[isRidSignif])))
					if(blnClump & !blnClumpInSignal) arnumloc = c(arnumloc, length(unique(locusId[isRidSignif & !grepl("NA",locusId)])))
				} 
			}
		} else {
			## use colTag for ONE acolIndep
			artag = sort(unique(tInSort[[colTag]][isRegion]))
			for(rtag in artag) {
				isRtagRegion = tInSort[[colTag]]==rtag & isRegion
				arnumvar = c(arnumvar, length(which(isRtagRegion)))
				if(blnRecombRate | blnGpos) arnumsig = c(arnumsig, length(unique(signalId[isRtagRegion])))
				if(blnClump & !blnClumpInSignal) arnumloc = c(arnumloc, length(unique(locusId[isRtagRegion & !grepl("NA",locusId)])))
			}
			
		}
		
		strrtag = ifelse(length(artag)==0,NA,ifelse(length(artag)==1, artag, paste(artag,collapse=";")))
		regionTag[isRegion] = strrtag
		
		#strNumVar = ifelse(length(arnumvar)==0,NA,ifelse(length(arnumvar)==1, as.character(arnumvar), paste(numvar,"(",paste(arnumvar,collapse=";"),")",sep="")))
		strNumVar = ifelse(length(arnumvar)==0,NA,ifelse(length(arnumvar)==1, as.character(arnumvar), paste(numvaruni,"(",paste(arnumvar,collapse=";"),")",sep="")))
		regionNumVariants[isRegion] = strNumVar
		
		if(blnRecombRate | blnGpos) {
			strNumSig = ifelse(length(arnumsig)==0,NA,ifelse(length(arnumsig)==1, as.character(arnumsig), paste(numsig,"(",paste(arnumsig,collapse=";"),")",sep="")))
			regionNumSignals[isRegion] = strNumSig
		}
		if(blnClump & !blnClumpInSignal) {
			strNumLoc = ifelse(length(arnumloc)==0,NA,ifelse(length(arnumloc)==1, as.character(arnumloc), paste(numloci,"(",paste(arnumloc,collapse=";"),")",sep="")))
			regionNumLoci[isRegion] = strNumLoc
		}		
		
	}
	
	if(blnRecombRate | blnGpos) {
		signalTag <- signalNumVariants <- signalNumLoci <- rep(NA, nrow(tInSort))
		
		for(sid in unique(signalId)) {
			
			isSignal = signalId == sid
			
			#numvar = length(which(isSignal))
			numvaruni = length(unique(tInSort[[colInMarker]][isSignal]))
			
			if(blnClump) numloci = length(unique(locusId[isSignal & !grepl("NA",locusId)]))
			
			astag = c()
			asnumvar = c()
			asnumloc = c()
			
			if(colTag == "") {
			
				for(i in 1:length(acolIndep)) { 
					
					# colp = acolIndep[i]
					# tag = astrTag[i]
					# plim = anumIndepLim[i]
					
					# isSidSignif = isSignal & tInSort[,colp]<plim
					isSidSignif = isSignal & tCritSort[,i]
					
					if(any(isSidSignif)) {
						astag = c(astag, astrTag[i])
						asnumvar = c(asnumvar, length(which(isSidSignif)))
						
						if(blnClump & blnClumpInSignal) asnumloc = c(asnumloc, length(unique(locusId[isSidSignif & !grepl("NA",locusId)])))
					}
				}
			
			} else {
			
				## use colTag for ONE acolIndep
				astag = sort(unique(tInSort[[colTag]][isSignal]))
				for(stag in astag) {
					isStagSignal = tInSort[[colTag]]==stag & isSignal
					asnumvar = c(asnumvar, length(which(isStagSignal)))
					if(blnClump & !blnClumpInSignal) asnumloc = c(asnumloc, length(unique(locusId[isStagSignal & !grepl("NA",locusId)])))
				}
			
			
			
			}
			
			
			strstag = ifelse(length(astag)==0,NA,ifelse(length(astag)==1, astag, paste(astag,collapse=";")))
			signalTag[isSignal] = strstag
			
			#strNumVar = ifelse(length(asnumvar)==0,NA,ifelse(length(asnumvar)==1, as.character(asnumvar), paste(numvar,"(",paste(asnumvar,collapse=";"),")",sep="")))
			strNumVar = ifelse(length(asnumvar)==0,NA,ifelse(length(asnumvar)==1, as.character(asnumvar), paste(numvaruni,"(",paste(asnumvar,collapse=";"),")",sep="")))
			signalNumVariants[isSignal] = strNumVar
			
			if(blnClump & blnClumpInSignal) {
				strNumLoc = ifelse(length(asnumloc)==0,NA,ifelse(length(asnumloc)==1, as.character(asnumloc), paste(numloci,"(",paste(asnumloc,collapse=";"),")",sep="")))
				signalNumLoci[isSignal] = strNumLoc
			}	
		}
	}
	
	if(blnClump) {
		locusTag <- locusNumVariants <- rep(NA, nrow(tInSort))
		
		liduni = unique(locusId)
		liduni = liduni[!grepl("NA",liduni)]
		
		for(lid in liduni) {
			
			isLocus = locusId == lid
			
			#numvar = length(which(isLocus))
			numvaruni = length(unique(tInSort[[colInMarker]][isLocus]))
			
			altag = c()
			alnumvar = c()
			
			if(colTag == "") {
			
				for(i in 1:length(acolIndep)) { 
					
					# colp = acolIndep[i]
					# tag = astrTag[i]
					# plim = anumIndepLim[i]
					
					# isLidSignif = isLocus & tInSort[,colp]<plim
					isLidSignif = isLocus & tCritSort[,i]
					
					if(any(isLidSignif)) {
						altag = c(altag, astrTag[i])
						alnumvar = c(alnumvar, length(which(isLidSignif)))
					}
				}
			} else {
			
				altag = sort(unique(tInSort[[colTag]][isLocus]))
				for(ltag in altag) {
					isLtagLocus = tInSort[[colTag]]==ltag & isLocus
					alnumvar = c(alnumvar, length(which(isLtagLocus)))
				}
			
			}
			
			
			strltag = ifelse(length(altag)==0,NA,ifelse(length(altag)==1, altag, paste(altag,collapse=";")))
			locusTag[isLocus] = strltag
			
			#strNumVar = ifelse(length(alnumvar)==0,NA,ifelse(length(alnumvar)==1, as.character(alnumvar), paste(numvar,"(",paste(alnumvar,collapse=";"),")",sep="")))
			strNumVar = ifelse(length(alnumvar)==0,NA,ifelse(length(alnumvar)==1, as.character(alnumvar), paste(numvaruni,"(",paste(alnumvar,collapse=";"),")",sep="")))
			locusNumVariants[isLocus] = strNumVar
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
		
		if(blnRecombRate | blnGpos) {
			for(sid in unique(signalId)) {
				
				isSignal = signalId == sid
				schr = unique(aChr[isSignal])
				spos = unique(unlist(lapply(strsplit(signalCoordinates[isSignal],":"), function(x) x[2])))
				spos = gsub("EoCHR","Inf",spos)
				spos1 = as.numeric(unlist(lapply(strsplit(spos,"_"), function(x) x[1])))
				spos2 = as.numeric(unlist(lapply(strsplit(spos,"_"), function(x) x[2])))
				
				#isKnown = anchr == schr & ((anpos1 >= (spos1-numAnnotPosLim) & anpos1 <= (spos2+numAnnotPosLim)) | (anpos2 >= (spos1-numAnnotPosLim) & anpos2 <= (spos2+numAnnotPosLim)))
				
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
	objGWA.indep = GWADATA.cbind(objGWA.indep, regionNumVariants, paste(strTag,"regionNumVariants",sep=""))
	objGWA.indep = GWADATA.cbind(objGWA.indep, regionCoordinates, paste(strTag,"regionCoordinates",sep=""))
	objGWA.indep = GWADATA.cbind(objGWA.indep, regionSize, paste(strTag,"regionSize",sep=""))
	objGWA.indep = GWADATA.cbind(objGWA.indep, regionExtCoordinates, paste(strTag,"regionExtCoordinates",sep=""))
	objGWA.indep = GWADATA.cbind(objGWA.indep, regionExtSize, paste(strTag,"regionExtSize",sep=""))
	
	
	if(blnRecombRate | blnGpos) {
		objGWA.indep = GWADATA.cbind(objGWA.indep, regionNumSignals, paste(strTag,"regionNumSignals",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, signalId, paste(strTag,"signalId",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, signalLead, paste(strTag,"signalLead",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, signalTag, paste(strTag,"signalTag",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, signalNumVariants, paste(strTag,"signalNumVariants",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, signalCoordinates, paste(strTag,"signalCoordinates",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, signalSize, paste(strTag,"signalSize",sep=""))
	}
	
	if(blnClump) {
		if(!blnClumpInSignal) {
			objGWA.indep = GWADATA.cbind(objGWA.indep, regionNumLoci, paste(strTag,"regionNumLoci",sep=""))
		} else {
			objGWA.indep = GWADATA.cbind(objGWA.indep, signalNumLoci, paste(strTag,"signalNumLoci",sep=""))
		}
		objGWA.indep = GWADATA.cbind(objGWA.indep, locusId, paste(strTag,"locusId",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, locusLead, paste(strTag,"locusLead",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, locusTag, paste(strTag,"locusTag",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, locusNumVariants, paste(strTag,"locusNumVariants",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, locusCoordinates, paste(strTag,"locusCoordinates",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, locusSize, paste(strTag,"locusSize",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, locusR2, paste(strTag,"locusR2",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, locusNumMiss, paste(strTag,"locusNumMiss",sep=""))
	}
	
	if(blnAnnot) {
		objGWA.indep = GWADATA.cbind(objGWA.indep, regionAnnot, paste(strTag,"regionAnnot",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, annotDistance, paste(strTag,"annotDistance",sep=""))
	}	
	
	if(blnGene) {
		objGWA.indep = GWADATA.cbind(objGWA.indep, aNearestGene, paste(strTag,"NearestGene",sep=""))
		objGWA.indep = GWADATA.cbind(objGWA.indep, aNearestGeneDistance, paste(strTag,"NearestGeneDistance",sep=""))
	}
	
	####################################
	### 6. Resort and prepare output

	## no resort: 
	objGWA.indep@tblGWA <- objGWA.indep@tblGWA
	objGWA.indep.regionLeads <- GWADATA.getrows(objGWA.indep, which(!is.na(regionLead)))
	
	objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numRegion",sep=""),length(unique(regionId)))
	# ## regionTag
	# artleads = regionTag[which(!is.na(regionLead))]
	# for(rt in sort(unique(artleads))) {
		# nrt = length(which(artleads==rt))
		# objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numRegion.",rt,sep=""),nrt)
	# }
	### -> CAUSES ERR IN REPORT 
	
	if(blnRecombRate | blnGpos) {
		objGWA.indep.signalLeads <- GWADATA.getrows(objGWA.indep, which(!is.na(signalLead)))
		objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numSignal",sep=""),length(unique(signalId)))
		## signalTag
		# astleads = signalTag[which(!is.na(signalLead))]
		# for(st in sort(unique(astleads))) {
			# nst = length(which(astleads==st))
			# objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numSignal.",st,sep=""),nst)
		# }
		### -> CAUSES ERR IN REPORT 
		
	} else {
		objGWA.indep.signalLeads <- objGWA.indep.regionLeads
	}
	
	if(blnClump) {
		objGWA.indep.locusLeads <- GWADATA.getrows(objGWA.indep, which(!is.na(locusLead)))
		objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numLoci",sep=""),length(unique(locusId[!grepl("NA",locusId)])))
		objGWA.indep.locusMiss <- GWADATA.getrows(objGWA.indep, which(grepl("NA",locusId)))
		objREPORT <- REPORT.setval(objREPORT,paste(strTag,"numLociMiss",sep=""),length(which(grepl("NA",locusId))))
	} else {
		objGWA.indep.locusLeads <- objGWA.indep.locusMiss <- objGWA.indep.regionLeads
	}
	
	if(blnAddIndepInfo) {
		### merges by intersect(names)		
		
		if(colInMarker!="") {
			aColAdd = names(objGWA.indep@tblGWA)[ ! names(objGWA.indep@tblGWA) %in% names(objGWA@tblGWA)]
			#objGWA.indep = GWADATA.getcols(objGWA.indep, c(colInMarker, aColAdd))
			objGWA.indep_for_merge = GWADATA.getcols(objGWA.indep, c(colInMarker, aColAdd))
			objGWA <- GWADATA.merge(objGWA,objGWA.indep_for_merge, strSuffix.In = "", strSuffix.Add = "", blnAll.In = TRUE, blnAll.Add = TRUE, strBy.In = colInMarker, strBy.Add = colInMarker)
		} else {
			objGWA <- GWADATA.merge(objGWA,objGWA.indep, strSuffix.In = "", strSuffix.Add = "", blnAll.In = TRUE, blnAll.Add = TRUE, strBy.In = NA, strBy.Add = NA)
		}
		
		# if(blnAddIndepInfoExtCoord) {
			# ## extend coord info column to all variants in the Ext Region
			# aChrIn = GWADATA.getcol(objGWA, colInChr)
			# aPosIn = as.integer(GWADATA.getcol(objGWA, colInPos))
			
			# aExtCoord = GWADATA.getcol(objGWA, paste(strTag,"regionExtCoordinates",sep=""))
			# aExtCoordOut = aExtCoord
			# ## need to replace ton signif variants in regions with coord !
			
			# aExtCoordUni = unique(aExtCoord)
			# aExtCoordUni = aExtCoordUni[!is.na(aExtCoordUni)]
			
			# aChrExtUni = unlist(lapply(strsplit(aExtCoordUni,":"),function(x) x[1]))
			# strpos = unlist(lapply(strsplit(aExtCoordUni,":"),function(x) x[2]))
			# aPos1ExtUni = as.integer(unlist(lapply(strsplit(strpos,"_"),function(x) x[1])))
			# aPos2ExtUni = as.integer(unlist(lapply(strsplit(strpos,"_"),function(x) x[2])))
			
			# for(iregion in 1:length(aChrExtUni)) {
				# isin = aChrIn == aChrExtUni[iregion] & aPosIn>=aPos1ExtUni[iregion] & aPosIn<=aPos2ExtUni[iregion]
				# aExtCoordOut[isin] = aExtCoordUni[iregion]
			# }
			
			# objGWA.indep = GWADATA.cbind(objGWA.indep, aExtCoordOut, paste(strTag,"regionExtCoordinates",sep=""),blnOverwrite=TRUE)
		
		# }
		
	}	
	
	
	return(list(objGWA,objGWA.indep,objGWA.indep.regionLeads,objGWA.indep.signalLeads,objGWA.indep.locusLeads,objGWA.indep.locusMiss,objREPORT))

}

INDEP <- function(strEqcCommand){ 
	## Wrapper for class definition
	INDEPout <- setINDEP(new("INDEP", strEqcCommand = strEqcCommand))
	validINDEP(INDEPout)
	#INDEPout.valid <- validINDEP(INDEPout)
	return(INDEPout)
}

