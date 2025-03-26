setClass("ANNOTATE",
	representation = representation(
						strEqcCommand	=	"character",
						colInChr		=	"character",
						colInPos		=	"character",
						colInPval		=	"character",
						colInCoord		=	"character",
						colInSkipAnnot	=	"character",
						numPvalLim		= 	"numeric",
						fileAnnot		=	"character",
						colAnnotChr		=	"character",
						colAnnotPos		=	"character",
						strAnnotTag		=	"character",
						colAnnotTag		=	"character",
						colAnnotCoord	= 	"character",
						colOutAnnot		=	"character",
						colOutCoord		=	"character",
						numAnnotPosLim	=	"numeric",
						blnMerge		=	"logical"
						),
	prototype = prototype(
						strEqcCommand	=	"",
						colInChr		=	"",
						colInPos		=	"",
						colInPval		=	"",
						colInCoord		=	"",
						colInSkipAnnot	=	"",
						numPvalLim		= 	1,
						fileAnnot		=	"",
						colAnnotChr		=	"Chr",
						colAnnotPos		=	"Pos",
						strAnnotTag		=	"",
						colAnnotTag		=	"",
						colAnnotCoord	=   "",
						colOutAnnot		=	"AnnotTag",
						colOutCoord		=	"AnnotCoord",
						numAnnotPosLim	=	500000,
						blnMerge		=	FALSE
						)
	#contains = c("EcfReader")
)

setGeneric("setANNOTATE", function(object) standardGeneric("setANNOTATE"))
setMethod("setANNOTATE", signature = (object = "ANNOTATE"), function(object) {
	
	#aEqcSlotNamesIn = c("colInChr","colInPos","colInPval","numPvalLim","fileAnnot","colAnnotChr","colAnnotPos","strAnnotTag","colOutAnnot","numAnnotPosLim")
	aEqcSlotNamesIn = c("colInChr","colInPos","colInPval","colInSkipAnnot","colInCoord","numPvalLim","fileAnnot","colAnnotChr","colAnnotPos","strAnnotTag","colOutAnnot","colOutCoord","numAnnotPosLim","colAnnotTag","colAnnotCoord","blnMerge")

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
validANNOTATE <- function(objANNOTATE) {
	
	if(objANNOTATE@colInChr == "") 
		stop(paste(" EASY ERROR:ANNOTATE\n No column colInChr defined. Please set colInChr.", sep=""))
	if(objANNOTATE@colInPos == "") 
		stop(paste(" EASY ERROR:ANNOTATE\n No column colInPos defined. Please set colInPos.", sep=""))
	if(objANNOTATE@fileAnnot == "" & objANNOTATE@colInCoord == "") 
		stop(paste(" EASY ERROR:ANNOTATE\n No reference file defined and no colInCoord defined. Please set fileAnnot OR colInCoord to obtain annotation info.", sep=""))
	
	if(objANNOTATE@strAnnotTag == "" & objANNOTATE@colAnnotTag == "") 
		stop(paste(" EASY ERROR:ANNOTATE\n No strAnnotTag or colAnnotTag defined. Please set one of the two.", sep=""))
			
	if(objANNOTATE@strAnnotTag != "" & objANNOTATE@colAnnotTag != "") 		 
		warning(paste(" EASY ERROR:ANNOTATE\n Both strAnnotTag and colAnnotTag were defined. ANNOTATE will use colAnnotTag and ignore the strAnnotTag.", sep=""))
	
	
	### Valid with GWADATA?
	
	if(objANNOTATE@fileAnnot != "") {
		if(!file.exists(objANNOTATE@fileAnnot))
			stop(paste("EASY ERROR:ANNOTATE\n File fileAnnot\n ",objANNOTATE@fileAnnot,"\n does not exist.", sep=""))
		### Cols exist?
		
		tblAnnot<-read.table(objANNOTATE@fileAnnot,header=T, sep="\t",  stringsAsFactors=FALSE)
		
		isUseCoord = objANNOTATE@colAnnotCoord != ""
		
		if(isUseCoord) {
			isAv = objANNOTATE@colAnnotCoord %in% names(tblAnnot)
				if(!isAv)
					stop(paste(" EASY ERROR:ANNOTATE\n Defined column --colAnnotCoord \n",objANNOTATE@colAnnotCoord, "\n is not available in known loci file \n",objANNOTATE@fileAnnot,"\n PLease specify correct column name.", sep=""))
		} else {
			isAv <- objANNOTATE@colAnnotChr %in% names(tblAnnot)
			if(!isAv)
				stop(paste(" EASY ERROR:ANNOTATE\n Defined column colAnnotChr \n",objANNOTATE@colAnnotChr, "\n is not available in fileAnnot. PLease specify correct column name.", sep=""))
			isAv <- objANNOTATE@colAnnotPos %in% names(tblAnnot)	
			if(!isAv)
				stop(paste(" EASY ERROR:ANNOTATE\n Defined column colAnnotPos \n",objANNOTATE@colAnnotPos, "\n is not available in fileAnnot. PLease specify correct column name.", sep=""))
			
		}
		
		if(objANNOTATE@colAnnotTag != "") {
			isAv <- objANNOTATE@colAnnotTag %in% names(tblAnnot)
			if(!isAv)
				stop(paste(" EASY ERROR:ANNOTATE\n Defined column colAnnotTag \n",objANNOTATE@colAnnotTag, "\n is not available in fileAnnot. PLease specify correct column name.", sep=""))
		}
	}
	
	return(TRUE)
}
ANNOTATE.GWADATA.valid <- function(objANNOTATE, objGWA) {
	
	isNotAv <- !(objANNOTATE@colInChr %in% objGWA@aHeader)
	if(isNotAv)
		stop(paste(" EASY ERROR:ANNOTATE\n Defined column colInChr \n",objANNOTATE@colInChr, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	isNotAv <- !(objANNOTATE@colInPos %in% objGWA@aHeader)
	if(isNotAv)
		stop(paste(" EASY ERROR:ANNOTATE\n Defined column colInPos \n",objANNOTATE@colInPos, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	iPos = match(objANNOTATE@colInPos, objGWA@aHeader)
	isPosNumeric <- objGWA@aClasses[iPos] == "numeric" | objGWA@aClasses[iPos] == "integer" | objGWA@aClasses[iPos] == "double"

	if(!isPosNumeric)
		stop(paste(" EASY ERROR:ANNOTATE\n Defined column colInPos \n",objANNOTATE@colInPos, "\n is not numeric for GWA data-set \n",objGWA@fileIn,"\n . Please cast colInPos to numeric, integer or double.", sep=""))
	
	
	isDefinedNotAv <- objANNOTATE@colInPval != "" & !(objANNOTATE@colInPval %in% objGWA@aHeader)
	if(isDefinedNotAv)
		stop(paste(" EASY ERROR:ANNOTATE\n Defined column colInPval \n",objANNOTATE@colInPval, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	if(objANNOTATE@numPvalLim < 1) {
		if(!(objANNOTATE@colInPval %in% objGWA@aHeader))
			stop(paste(" EASY ERROR:ANNOTATE\n To use the Pvalue threshold for the annotation of loci, you have to define colInPval as well.","\n PLease specify colInPval with the ANNOTATE command !!!", sep=""))
		if(objANNOTATE@numPvalLim <=0)
			stop(paste(" EASY ERROR:ANNOTATE\n PLease specify numPvalLim within the range ]0,1[ within the ANNOTATE command !!!", sep=""))
	}
	
	if(objANNOTATE@colInCoord != "") {
		isNotAv <- !(objANNOTATE@colInCoord %in% objGWA@aHeader)
		if(isNotAv)
		stop(paste(" EASY ERROR:ANNOTATE\n Defined column colInCoord \n",objANNOTATE@colInCoord, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	}
	
}
#############################################################################################################################
ANNOTATE.run <- function(objANNOTATE, objGWA) {
	
	colInChr		=	objANNOTATE@colInChr
	colInPos		=	objANNOTATE@colInPos
	colInPval		=	objANNOTATE@colInPval
	colInCoord		=	objANNOTATE@colInCoord
	colInSkipAnnot	=	objANNOTATE@colInSkipAnnot
	numPvalLim		=	objANNOTATE@numPvalLim
	colAnnotChr		=	objANNOTATE@colAnnotChr
	colAnnotPos		=	objANNOTATE@colAnnotPos
	colAnnotTag		=	objANNOTATE@colAnnotTag
	strAnnotTag		=	objANNOTATE@strAnnotTag
	colOutAnnot		= 	objANNOTATE@colOutAnnot
	colOutCoord    	= 	objANNOTATE@colOutCoord
	
	numAnnotPosLim	=	objANNOTATE@numAnnotPosLim
	colAnnotCoord  	= objANNOTATE@colAnnotCoord
	
	blnMerge = objANNOTATE@blnMerge
	
	aInChr = GWADATA.getcol(objGWA, colInChr)
	aInPos = GWADATA.getcol(objGWA, colInPos)
	
	
	if(numPvalLim < 1) {
		aInPval = GWADATA.getcol(objGWA, colInPval)
	} else {
		aInPval = rep(0, dim(objGWA@tblGWA)[1])
	}
	
	if(colInSkipAnnot!="") {
		aInSkip = GWADATA.getcol(objGWA, colInSkipAnnot)
	} else {
		aInSkip = rep(NA,length(aInPval))
	}
	

	if(objANNOTATE@fileAnnot != "") {
		tblLoci<-read.table(objANNOTATE@fileAnnot, header=T, sep="\t", stringsAsFactors=FALSE)
		numLoci = dim(tblLoci)[1]
		
		if(colAnnotCoord!="") {
			## extract values from coord 1:123_466
			acoord = tblLoci[,colAnnotCoord]
			aRefChr = unlist(lapply(strsplit(acoord,":"),function(x) x[1]))
			strpos = unlist(lapply(strsplit(acoord,":"),function(x) x[2]))
			aRefPos1 = as.integer(unlist(lapply(strsplit(strpos,"_"),function(x) x[1])))
			aRefPos2 = as.integer(unlist(lapply(strsplit(strpos,"_"),function(x) x[2])))
		} else {
			aRefChr = tblLoci[,colAnnotChr]
			anpos = as.integer(tblLoci[,colAnnotPos])
			aRefPos1 = anpos
			aRefPos2 = anpos
		}
		
		# aRefChr = tblLoci[, which(names(tblLoci) == colAnnotChr)]
		# aRefPos = tblLoci[, which(names(tblLoci) == colAnnotPos)]
		
		if(colAnnotTag != "") {
			aRefTag = tblLoci[, which(names(tblLoci) == colAnnotTag)]
		} else {
			aRefTag = rep(strAnnotTag,length(aRefPos1))
		}
	} else {
		# use colInCoord to obtain annotation
		acoord = unique(GWADATA.getcol(objGWA, colInCoord))
		acoord = acoord[!is.na(acoord)]
		if(length(acoord)>0) {
			aRefChr = unlist(lapply(strsplit(acoord,":"),function(x) x[1]))
			strpos = unlist(lapply(strsplit(acoord,":"),function(x) x[2]))
			aRefPos1 = as.integer(unlist(lapply(strsplit(strpos,"_"),function(x) x[1])))
			aRefPos2 = as.integer(unlist(lapply(strsplit(strpos,"_"),function(x) x[2])))
			aRefTag = rep(strAnnotTag,length(aRefPos1))
		}
		numLoci = length(acoord)
	}
	
	if(colOutAnnot%in%objGWA@aHeader) {
		aTagOut = GWADATA.getcol(objGWA, colOutAnnot) # is already annotated
		if(colOutCoord%in%objGWA@aHeader) {
			aCoordOut = GWADATA.getcol(objGWA, colOutCoord) # is already annotated
		} else {
			aCoordOut = rep(NA, dim(objGWA@tblGWA)[1]) # not used before
		}
	} else {
		aTagOut = rep(NA, dim(objGWA@tblGWA)[1]) # not used before
		aCoordOut = rep(NA, dim(objGWA@tblGWA)[1]) # not used before
	}
	
	if(numLoci>0) {
		for(i in 1:numLoci) {
			chrTmp = aRefChr[i]
			# posTmp = aRefPos[i]
			pos1Tmp = aRefPos1[i]
			pos2Tmp = aRefPos2[i]
			tagTmp = aRefTag[i]
			
			pos1Tmp = pos1Tmp-numAnnotPosLim
			pos2Tmp = pos2Tmp+numAnnotPosLim
			
			coordTmp = paste0(chrTmp,":",pos1Tmp,"_",pos2Tmp)
			
			#isCurrentLocus = aInChr == chrTmp & abs(aInPos - posTmp) <= numAnnotPosLim & aInPval <= numPvalLim
			# isCurrentLocus = aInChr == chrTmp & abs(aInPos - posTmp) <= numAnnotPosLim
			# isCurrentLocus[is.na(isCurrentLocus)] = FALSE
			
			isCurrentLocus = aInChr == chrTmp & aInPos>=pos1Tmp & aInPos<=pos2Tmp
			isCurrentLocus[is.na(isCurrentLocus)] = FALSE
			
			if(any(isCurrentLocus)) {
				isAnyPvalLow = any(aInPval[isCurrentLocus] <= numPvalLim)
				isSkip = any(!is.na(aInSkip[isCurrentLocus]))
				# if value!=NA in aInSkip then skip annotation !
				
				if(isAnyPvalLow & !isSkip) {
					
					if(any(!is.na(aTagOut[isCurrentLocus]))) {
						# isNa = is.na(aTagOut[isCurrentLocus])
						
						if(!blnMerge) {
							isAvail = unlist(lapply(strsplit(aTagOut[isCurrentLocus],";"),function(x) any(x==tagTmp)))
							## NA, TRUE or FALSE
							aTagOut[isCurrentLocus] = ifelse(is.na(isAvail), tagTmp, ifelse(isAvail, aTagOut[isCurrentLocus], paste(aTagOut[isCurrentLocus],tagTmp, sep=";")))
							aCoordOut[isCurrentLocus] = ifelse(is.na(isAvail), coordTmp, ifelse(isAvail, aCoordOut[isCurrentLocus], paste(aCoordOut[isCurrentLocus],coordTmp, sep=";")))
							
						} else {
							## merge tags and coord
							
							aCoordUniCl = unique(c(aCoordOut[isCurrentLocus],coordTmp))
							aCoordUniCl = aCoordUniCl[!is.na(aCoordUniCl)]
							
							strposUniCl = unlist(lapply(strsplit(aCoordUniCl,":"),function(x) x[2]))
							aPos1UniCl = as.integer(unlist(lapply(strsplit(strposUniCl,"_"),function(x) x[1])))
							aPos2UniCl = as.integer(unlist(lapply(strsplit(strposUniCl,"_"),function(x) x[2])))
							
							pos1TmpCl = min(aPos1UniCl,na.rm=T)
							pos2TmpCl = max(aPos2UniCl,na.rm=T)
							
							coordTmpNew = paste0(chrTmp,":",pos1TmpCl,"_",pos2TmpCl)
							
							isCurrentLocusExt = aInChr == chrTmp & aInPos>=pos1TmpCl & aInPos<=pos2TmpCl
							isCurrentLocusExt[is.na(isCurrentLocusExt)] = FALSE
							
							aCoordOut[isCurrentLocusExt] = coordTmpNew
							
							aTagUniCl = unique(c(unlist(strsplit(aTagOut[isCurrentLocusExt],";")),tagTmp))
							aTagUniCl = aTagUniCl[!is.na(aTagUniCl)]
							aTagOut[isCurrentLocusExt] = paste(aTagUniCl,collapse=";")
						
						}
					} else {
						aTagOut[isCurrentLocus] = tagTmp
						aCoordOut[isCurrentLocus] = paste0(chrTmp,":",pos1Tmp,"_",pos2Tmp)
					}
				}
			}
		}
	}
	
	objGWA <- GWADATA.cbind(objGWA, aTagOut, colOutAnnot, blnOverwrite=TRUE)
	
	objGWA <- GWADATA.cbind(objGWA, aCoordOut, colOutCoord, blnOverwrite=TRUE)
	
	return(objGWA)
}

ANNOTATE <- function(strEqcCommand){ 
	## Wrapper for class definition
	ANNOTATEout <- setANNOTATE(new("ANNOTATE", strEqcCommand = strEqcCommand))
	validANNOTATE(ANNOTATEout)
	#ANNOTATEout.valid <- validANNOTATE(ANNOTATEout)
	return(ANNOTATEout)
}

