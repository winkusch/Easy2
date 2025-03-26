setClass("GWASCATALOG",
	representation = representation(
						strEqcCommand	=	"character",
						colInChr		=	"character",
						colInPos		=	"character",
						colInCoord		=	"character",
						colOutAnnot		=	"character",
						numPosLim		=	"numeric",
						fileGwasCatalogue = "character",
						fileMap 		= "character",
						blnSimple = 	"logical",
						blnGwsPval = 	"logical",
						astrTraits = 	"character"
						),
	prototype = prototype(
						strEqcCommand	=	"",
						colInChr		=	"",
						colInPos		=	"",
						colInCoord		=	"",
						colOutAnnot		=	"AnnotTag",
						numPosLim		=	500000,
						fileGwasCatalogue = "",
						fileMap 		= "",
						blnSimple = 	FALSE,
						blnGwsPval = 	TRUE,
						astrTraits =    ""
						)
	#contains = c("EcfReader")
)

setGeneric("setGWASCATALOG", function(object) standardGeneric("setGWASCATALOG"))
setMethod("setGWASCATALOG", signature = (object = "GWASCATALOG"), function(object) {
	
	#aEqcSlotNamesIn = c("colInChr","colInPos","colInPval","numPvalLim","fileAnnot","colRefChr","colRefPos","strAnnotTag","colOutAnnot","numPosLim")
	aEqcSlotNamesIn = c("colInChr","colInPos","colInCoord","colOutAnnot","numPosLim","fileGwasCatalogue","fileMap","blnSimple","blnGwsPval","astrTraits")

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
validGWASCATALOG <- function(objGWASCATALOG) {
	
	if(objGWASCATALOG@colInCoord == "" & (objGWASCATALOG@colInChr == "" | objGWASCATALOG@colInPos == "") ) {
		stop(paste(" EASY ERROR:GWASCATALOG\n No column colInChr, colInPos or colInCoord defined. Please set colInCoord OR colInChr & colInPos.", sep=""))
	}
	
	if(objGWASCATALOG@fileGwasCatalogue == "") 
		stop(paste(" EASY ERROR:GWASCATALOG\n No reference file defined. Please set fileGwasCatalogue.", sep=""))
	if(objGWASCATALOG@fileMap == "") 
		stop(paste(" EASY ERROR:GWASCATALOG\n No reference file defined. Please set fileMap.", sep=""))
		
	### Valid with GWADATA?
	

	if(!file.exists(objGWASCATALOG@fileGwasCatalogue))
		stop(paste("EASY ERROR:GWASCATALOG\n File fileGwasCatalogue\n ",objGWASCATALOG@fileGwasCatalogue,"\n does not exist.", sep=""))
	
	if(!file.exists(objGWASCATALOG@fileMap))
		stop(paste("EASY ERROR:GWASCATALOG\n File fileMap\n ",objGWASCATALOG@fileMap,"\n does not exist.", sep=""))
	
	return(TRUE)
}
GWASCATALOG.GWADATA.valid <- function(objGWASCATALOG, objGWA) {
	
	if(objGWASCATALOG@colInCoord == "") {
		isNotAv <- !(objGWASCATALOG@colInChr %in% objGWA@aHeader)
		if(isNotAv)
			stop(paste(" EASY ERROR:GWASCATALOG\n Defined column colInChr \n",objGWASCATALOG@colInChr, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
		
		isNotAv <- !(objGWASCATALOG@colInPos %in% objGWA@aHeader)
		if(isNotAv)
			stop(paste(" EASY ERROR:GWASCATALOG\n Defined column colInPos \n",objGWASCATALOG@colInPos, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
		
		iPos = match(objGWASCATALOG@colInPos, objGWA@aHeader)
		isPosNumeric <- objGWA@aClasses[iPos] == "numeric" | objGWA@aClasses[iPos] == "integer" | objGWA@aClasses[iPos] == "double"
		
		if(!isPosNumeric)
			stop(paste(" EASY ERROR:GWASCATALOG\n Defined column colInPos \n",objGWASCATALOG@colInPos, "\n is not numeric for GWA data-set \n",objGWA@fileIn,"\n . Please cast colInPos to numeric, integer or double.", sep=""))
		
	} else {
		isNotAv <- !(objGWASCATALOG@colInCoord %in% objGWA@aHeader)
		if(isNotAv)
			stop(paste(" EASY ERROR:GWASCATALOG\n Defined column colInCoord \n",objGWASCATALOG@colInCoord, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	}
	
}
#############################################################################################################################
GWASCATALOG.run <- function(objGWASCATALOG, objGWA) {
	
	colInChr		=	objGWASCATALOG@colInChr
	colInPos		=	objGWASCATALOG@colInPos
	colInCoord		=	objGWASCATALOG@colInCoord
	colOutAnnot		= 	objGWASCATALOG@colOutAnnot
	numPosLim		=	objGWASCATALOG@numPosLim
	fileGwasCatalogue = objGWASCATALOG@fileGwasCatalogue
	fileMap 		= objGWASCATALOG@fileMap
	
	blnSimple = objGWASCATALOG@blnSimple
	blnGwsPval = objGWASCATALOG@blnGwsPval
	astrTraits = objGWASCATALOG@astrTraits
	#### 
	
	tcat<-read.csv(objGWASCATALOG@fileGwasCatalogue, header=T, stringsAsFactors=FALSE)
	
	tcat$rsid = unlist(lapply(strsplit(tcat$Variant.and.risk.allele,"-"),function(x) x[1]))
	
	ptmp = gsub("10-","e-",tcat$P.value)
	ptmp = gsub(" x ",".",ptmp)
	tcat$P.value = as.numeric(ptmp)
	if(blnGwsPval) tcat = tcat[tcat$P.value<=5e-8, ]
	
	if(astrTraits!="") {
		atraits = strsplit(astrTraits,";")[[1]]
		atraits = gsub("'","",atraits,fixed=T)
		isok = tcat$Reported.trait %in% atraits
		tcat = tcat[isok, ]
	}
	
	## add chr, pos to gwas catalogue data-set
	tmap<-fread(objGWASCATALOG@fileMap, header=T, stringsAsFactors=FALSE, data.table=F)
		
	tblLoci = merge(tcat, tmap, by = "rsid", all = F)
	colRefChr = "chr"
	colRefPos = "pos"
	colRefTag1 = "First.Author"
	colRefTag2 = "Reported.trait"
	
	numLoci = dim(tblLoci)[1]
	
	if(colInCoord=="") {
		aInChr = GWADATA.getcol(objGWA, colInChr)
		aInPos1 = GWADATA.getcol(objGWA, colInPos)
		aInPos2 = GWADATA.getcol(objGWA, colInPos)
	} else {
		aInCoord = GWADATA.getcol(objGWA, colInCoord)
		aInChr = unlist(lapply(strsplit(aInCoord,":"), function(x) x[1]))
		aInPosReg = unlist(lapply(strsplit(aInCoord,":"), function(x) x[2]))
		aInPos1 = as.integer(unlist(lapply(strsplit(aInPosReg,"_"), function(x) x[1])))
		aInPos2 = as.integer(unlist(lapply(strsplit(aInPosReg,"_"), function(x) x[2])))
	}
	
	aRefChr = tblLoci[, which(names(tblLoci) == colRefChr)]
	aRefPos = tblLoci[, which(names(tblLoci) == colRefPos)]
	aRefTag = paste(tblLoci[, which(names(tblLoci) == colRefTag1)], "(", tblLoci[, which(names(tblLoci) == colRefTag2)], ")",sep = "")
	
	aTagOut = rep(NA, dim(objGWA@tblGWA)[1]) # not used before
	
	for(i in 1:numLoci) {
		chrTmp = aRefChr[i]
		posTmp = aRefPos[i]
		
		print(paste("GWASCATALOG annotation step ",i,"of",numLoci, " ... "))
		
		# isCurrentLocus = aInChr == chrTmp & abs(aInPos - posTmp) <= numPosLim
		# isCurrentLocus = aInChr == chrTmp & (abs(aInPos1 - posTmp) <= numPosLim | abs(aInPos2 - posTmp) <= numPosLim)
		
		isCurrentLocus = aInChr == chrTmp & posTmp >= (aInPos1 - numPosLim) & posTmp <= (aInPos2 + numPosLim)
		
		
		isCurrentLocus[is.na(isCurrentLocus)] = FALSE
		
		if(blnSimple) {
			tagTmp = "y"
			if(any(isCurrentLocus)) {
				aTagOut[isCurrentLocus] = tagTmp
			}
		} else {
			tagTmp = aRefTag[i]

			if(any(isCurrentLocus)) {
				if(any(!is.na(aTagOut[isCurrentLocus]))) {
					# isNa = is.na(aTagOut[isCurrentLocus])
					isAvail = unlist(lapply(strsplit(aTagOut[isCurrentLocus],";"),function(x) any(x==tagTmp)))
					## NA, TRUE or FALSE
					aTagOut[isCurrentLocus] = ifelse(is.na(isAvail), tagTmp, ifelse(isAvail, aTagOut[isCurrentLocus], paste(aTagOut[isCurrentLocus],tagTmp, sep=";")))
				} else {
					aTagOut[isCurrentLocus] = tagTmp
				}
			}
		}
	}
	
	objGWA <- GWADATA.cbind(objGWA, aTagOut, colOutAnnot, blnOverwrite=TRUE)
	
	return(objGWA)
}

GWASCATALOG <- function(strEqcCommand){ 
	## Wrapper for class definition
	GWASCATALOGout <- setGWASCATALOG(new("GWASCATALOG", strEqcCommand = strEqcCommand))
	validGWASCATALOG(GWASCATALOGout)
	#GWASCATALOGout.valid <- validGWASCATALOG(GWASCATALOGout)
	return(GWASCATALOGout)
}

