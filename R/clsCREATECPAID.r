setClass("CREATECPAID",
	representation = representation(
						strEqcCommand		=	"character",
						colInA1				=	"character",
						colInA2				=	"character",
						colInEA				=	"character",
						colInOA				=	"character",
						colInChr			=	"character",
						colInPos			=	"character",
						acolInFreq			=	"character",
						acolInBeta			=	"character",
						strTag				=	"character"
						),
	prototype = prototype(
						strEqcCommand		=	"",
						colInA1				=	"",
						colInA2				=	"",
						colInEA				=	"",
						colInOA				=	"",
						colInChr			=	"",
						colInPos			=	"",
						acolInFreq			=	"",
						acolInBeta			=	"",
						strTag				=	""
						)
)

CREATECPAID.set <- function(strEqcCommand, objCPA) {

	aEqcSlotNamesIn = c("colInA1","colInA2","colInEA","colInOA","colInChr","colInPos","acolInFreq","acolInBeta","strTag")
	
	## astrPatterns
	## *_CHR_POS
	## chrCHR:POS
	
	### Last 4 are inherited from class GWADATA and can be used with CREATECPAID for reference file!
	
	objEqcReader <- EqcReader(strEqcCommand,aEqcSlotNamesIn)
	
	if(length(objEqcReader@lsEqcSlotsOut) > 0) {
		for(i in 1:length(objEqcReader@lsEqcSlotsOut)) {
			tmpSlot <- names(objEqcReader@lsEqcSlotsOut)[i]
			tmpSlotVal <- objEqcReader@lsEqcSlotsOut[[i]]
			
			if(all(!is.na(tmpSlotVal))) slot(objCPA, tmpSlot) <- tmpSlotVal
		}
	}
	
	return(objCPA)
}

# setGeneric("setCREATECPAID", function(object) standardGeneric("setCREATECPAID"))
# setMethod("setCREATECPAID", signature = (object = "CREATECPAID"), function(object, objGWADATA.default) {
	
	# object@blnUseFastRead <- objGWADATA.default@blnUseFastRead
	
	# aEqcSlotNamesIn = c("colInMarker","colInA1","colInA2","fileMap","colMapMarker","colMapChr","colMapPos","strTag","colInChr","colInPos","blnUseInMarker","blnUseFastRead")
	
	# ## astrPatterns
	# ## *_CHR_POS
	# ## chrCHR:POS
	
	# ### Last 4 are inherited from class GWADATA and can be used with CREATECPAID for reference file!
	
	# objEqcReader <- EqcReader(object@strEqcCommand,aEqcSlotNamesIn)
	
	# if(length(objEqcReader@lsEqcSlotsOut) > 0) {
		# for(i in 1:length(objEqcReader@lsEqcSlotsOut)) {
			# tmpSlot <- names(objEqcReader@lsEqcSlotsOut)[i]
			# tmpSlotVal <- objEqcReader@lsEqcSlotsOut[[i]]
			
			# if(all(!is.na(tmpSlotVal))) slot(object, tmpSlot) <- tmpSlotVal
		# }
	# }
	
	# return(object)
# })

#############################################################################################################################

CREATECPAID.GWADATA.valid <- function(objCPA, objGWA) {
	
	if(objCPA@colInA1!="" & objCPA@colInA2!="") {
	
		isAv <- objCPA@colInA1 %in% objGWA@aHeader
		if(!isAv)
			stop(paste(" EASY ERROR:CREATECPAID\n Defined column colInA1 \n",objCPA@colInA1, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
		
		isAv <- objCPA@colInA2 %in% objGWA@aHeader
		if(!isAv)
			stop(paste(" EASY ERROR:CREATECPAID\n Defined column colInA2 \n",objCPA@colInA2, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
		
		if(any(objCPA@acolInFreq!="")) stop(paste(" EASY ERROR:CREATECPAID\n --acolInFreq cannot be used with --colInA1 and --colInA2. \n PLease specify --colInEA (effect allele) and --colInOA (other allele) instead !", sep=""))
		if(any(objCPA@acolInBeta!="")) stop(paste(" EASY ERROR:CREATECPAID\n --acolInBeta cannot be used with --colInA1 and --colInA2. \n PLease specify --colInEA (effect allele) and --colInOA (other allele) instead !", sep=""))
		
	} else if(objCPA@colInEA!="" & objCPA@colInOA!="") {
	
		isAv <- objCPA@colInEA %in% objGWA@aHeader
		if(!isAv)
			stop(paste(" EASY ERROR:CREATECPAID\n Defined column colInEA \n",objCPA@colInEA, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
		
		isAv <- objCPA@colInOA %in% objGWA@aHeader
		if(!isAv)
			stop(paste(" EASY ERROR:CREATECPAID\n Defined column colInOA \n",objCPA@colInOA, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	} else {
		stop(paste(" EASY ERROR:CREATECPAID\n PLease specify (--colInA1 AND --colInA2) OR (--colInEA AND --colInOA)) with correct column names.", sep=""))
	}
	
	isAv <- objCPA@colInChr %in% objGWA@aHeader
	if(!isAv)
		stop(paste(" EASY ERROR:CREATECPAID\n Defined column colInChr \n",objCPA@colInChr, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct --colInChr OR remove --colInchr and --colInPos.", sep=""))

	isAv <- objCPA@colInPos %in% objGWA@aHeader
	if(!isAv)
		stop(paste(" EASY ERROR:CREATECPAID\n Defined column colInPos \n",objCPA@colInPos, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct --colInPos OR remove --colInchr and --colInPos.", sep=""))
	
	if(any(objCPA@acolInFreq!="")) {
		for(colInFreq in objCPA@acolInFreq) {
			isAv <- colInFreq %in% objGWA@aHeader
			if(!isAv)
				stop(paste(" EASY ERROR:CREATECPAID\n Defined column acolInFreq \n",colInFreq, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct --acolInFreq.", sep=""))
		}
	}

	if(any(objCPA@acolInBeta!="")) {
		for(colInBeta in objCPA@acolInBeta) {
			isAv <- colInBeta %in% objGWA@aHeader
			if(!isAv)
				stop(paste(" EASY ERROR:CREATECPAID\n Defined column acolInBeta \n",colInBeta, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct --acolInBeta.", sep=""))
		}
	}
	
	
}

#############################################################################################################################
CREATECPAID.run <- function(objCPA, objGWA, isValidScript) {
	
	colInChr 	<- objCPA@colInChr
	colInPos 	<- objCPA@colInPos
	colInA1 	<- objCPA@colInA1
	colInA2 	<- objCPA@colInA2
	colInEA 	<- objCPA@colInEA
	colInOA 	<- objCPA@colInOA
	acolInFreq 	<- objCPA@acolInFreq
	acolInBeta 	<- objCPA@acolInBeta
	
	strTag 		<- objCPA@strTag
	
	if(nchar(strTag)>0) strTag = paste(strTag,".",sep="") 
		
	cpaid <- rep(NA,nrow(objGWA@tblGWA))
	
	if(colInA1!="" & colInA2!="") {
		a1 <- as.character(objGWA@tblGWA[[colInA1]])
		a2 <- as.character(objGWA@tblGWA[[colInA2]])
	} else {
		a1 <- as.character(objGWA@tblGWA[[colInEA]])
		a2 <- as.character(objGWA@tblGWA[[colInOA]])
	}
	
	#tcp[[colInChr]] = gsub("^0+","",tcp[[colInChr]])
	#tcp[[colInPos]] = gsub("^0+","",tcp[[colInPos]])
	
	chr <- gsub("^0+","",objGWA@tblGWA[[colInChr]])
	pos <- gsub("^0+","",objGWA@tblGWA[[colInPos]])
	
	a1[is.na(a1)] <- "NA"
	a2[is.na(a2)] <- "NA"
	
	isA2First = a2 < a1 
	## switch dir
	
	cpaid <- ifelse(isA2First, paste(chr,":",pos,":",a2,"_",a1,sep=""), paste(chr,":",pos,":",a1,"_",a2,sep=""))
	
	if(objCPA@colInEA!="" & objCPA@colInOA!="") {
		objGWA@tblGWA[[colInEA]][isA2First] <- a2[isA2First]
		objGWA@tblGWA[[colInOA]][isA2First] <- a1[isA2First]
	
		if(any(acolInFreq!="")) {
			for(colInFreq in acolInFreq) {
				objGWA@tblGWA[[colInFreq]][isA2First] = 1 - objGWA@tblGWA[[colInFreq]][isA2First]
			}
		} else {
			warning(paste(" EASY WARNING:CREATECPAID\n Alleles colInEA =",colInEA, "and colInOA =",colInOA," have been changed to alphabetical order (colInEA first) without changing any allele frequency columns ! \n If required, specific --acolInFreq to allign allele frequencies to --colInEA.", sep=""))
		}
		
		if(any(acolInBeta!="")) {
			for(colInBeta in acolInBeta) {
				objGWA@tblGWA[[colInBeta]][isA2First] = - objGWA@tblGWA[[colInBeta]][isA2First]
			}
		} else {
			warning(paste(" EASY WARNING:CREATECPAID\n Alleles colInEA =",colInEA, "and colInOA =",colInOA," have been changed to alphabetical order (colInEA first) without changing any genetic effect size columns ! \n If required, specific --acolInBeta to allign betas to --colInEA.", sep=""))
		}
	}
		
	objGWA <- GWADATA.cbind(objGWA, cpaid, "cpaid")
	
	return(objGWA)
}
#############################################################################################################################
CREATECPAID <- function(strEqcCommand){ 
	## Wrapper for class definition
	CREATECPAIDout <- new("CREATECPAID")
	CREATECPAIDout <- CREATECPAID.set(strEqcCommand, CREATECPAIDout)
	
	return(CREATECPAIDout)

}
