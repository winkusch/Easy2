setClass("ADDGENE",
	representation = representation(
						strEqcCommand	= "character",
						rcdCriterion	= "character",
						colInChr		= "character",
						colInPos		= "character",
						fileAnnotGene	= "character",
						colInMarker		= "character",
						strTag			= "character"
						),
	prototype = prototype(
						strEqcCommand	= "",
						rcdCriterion	= "",
						colInChr		= "",
						colInPos		= "",
						fileAnnotGene	= "",
						colInMarker		= "",
						strTag			= ""
						)
	#contains = c("EcfReader")
)

setGeneric("setADDGENE", function(object) standardGeneric("setADDGENE"))
setMethod("setADDGENE", signature = (object = "ADDGENE"), function(object) {
	
	#aEqcSlotNamesIn = c("colInChr","colInPos","colInPval","numPvalLim","fileAnnot","colRefChr","colRefPos","strAnnotTag","colOutAnnot","numPosLim")
	aEqcSlotNamesIn = c("rcdCriterion","colInChr","colInPos","fileAnnotGene","colInMarker","strTag")

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
validADDGENE <- function(objADDGENE) {
	
	if(objADDGENE@rcdCriterion == "") 
		cat(paste(" EASY WARNING:ADDGENE\n No criterion rcdCriterion defined. All data will be used for gene adding.", sep=""))
	
	if(objADDGENE@colInChr == "") 
		stop(paste(" EASY ERROR:ADDGENE\n No column colInChr defined. Please set colInChr.", sep=""))
	if(objADDGENE@colInPos == "") 
		stop(paste(" EASY ERROR:ADDGENE\n No column colInPos defined. Please set colInPos.", sep=""))
	
	if(objADDGENE@fileAnnotGene == "" | !file.exists(objADDGENE@fileAnnotGene)) {
			stop(paste("EASY ERROR:ADDGENE\n File fileAnnotGene\n ",objADDGENE@fileAnnotGene,"\n does not exist.\n Please check definition of --fileAnnotGene.", sep=""))
	}
	
	return(TRUE)
}

ADDGENE.GWADATA.valid <- function(objADDGENE, objGWA) {
	
	isNotAv <- !(objADDGENE@colInChr %in% objGWA@aHeader)
	if(isNotAv)
		stop(paste(" EASY ERROR:ADDGENE\n Defined column colInChr \n",objADDGENE@colInChr, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	isNotAv <- !(objADDGENE@colInPos %in% objGWA@aHeader)
	if(isNotAv)
		stop(paste(" EASY ERROR:ADDGENE\n Defined column colInPos \n",objADDGENE@colInPos, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	return(TRUE)
	
}

#############################################################################################################################
ADDGENE.run <- function(objADDGENE, objGWA) {
	
	rcdCriterion	= objADDGENE@rcdCriterion
	colInChr		= objADDGENE@colInChr
	colInPos		= objADDGENE@colInPos
	fileAnnotGene	= objADDGENE@fileAnnotGene
	colInMarker		= objADDGENE@colInMarker
	strTag			= objADDGENE@strTag
	
	if(nchar(strTag)>0) strTag = paste(strTag,".",sep="") 
	
	if(rcdCriterion == "") {
		objGWA.addgene <- objGWA
	} else {
		objRCD 	= RCD(rcdCriterion)
		out 	= RCD.eval(objRCD, objGWA)
		numCrit = length(which(out))
		if(numCrit == 0) {
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"NearestGene",sep=""))
			objGWA <- GWADATA.cbind(objGWA, rep(NA,nrow(objGWA@tblGWA)), paste(strTag,"NearestGeneDistance",sep=""))
			return(objGWA)
		}
		objGWA.addgene <- GWADATA.getrows(objGWA, which(out))
	}
	
	
	## read gene file 
	tG = read.table(fileAnnotGene,sep="\t",stringsAsFactors=F,header=T,colClasses=c("character","integer","integer","character"))
	
	## t:\M_GENETICS\ReferenceData\glist-hg19.txt
	## 19 58858171 58864865 A1BG
	
	names(tG) <- c("chr","pos1","pos2","gene")
	
	aChr 		= GWADATA.getcol(objGWA.addgene, colInChr)
	aPos 		= GWADATA.getcol(objGWA.addgene, colInPos)
	
	aNearestGene <- aNearestGeneDistance <- rep(NA,length(aPos))
	
	print("Starting annotation ... ")
	for(i in 1:length(aChr)) {
		
		chr=aChr[i]
		pos=aPos[i]
		
		isWithinGene = (chr == tG$chr) & (pos >= tG$pos1) & (pos <= tG$pos2)
		
		if(any(isWithinGene)) {
			aNearestGene[i] = paste(tG$gene[isWithinGene],collapse=",")
			aNearestGeneDistance[i] = 0
		} else {
			if(any(tG$chr==chr)) {
				tGChr = tG[tG$chr==chr,]
				adismin = pmin(abs(pos-tGChr$pos1),abs(pos-tGChr$pos2))
				idxMin = which.min(adismin)
				aNearestGene[i] = tGChr$gene[idxMin[1]] # 1 if variant has equal dis to two genes
				aNearestGeneDistance[i] = min(adismin)
				## duplicates in glist -> take union !!!
				if(pos-tGChr$pos1[idxMin[1]]<0) aNearestGeneDistance[i] = - aNearestGeneDistance[i]
			}
		}
	}
	
	objGWA.addgene = GWADATA.cbind(objGWA.addgene, aNearestGene, paste(strTag,"NearestGene",sep=""))
	objGWA.addgene = GWADATA.cbind(objGWA.addgene, aNearestGeneDistance, paste(strTag,"NearestGeneDistance",sep=""))

	# objGWA.addgene = GWADATA.getcols(objGWA.addgene, c(colInMarker, paste(strTag,"NearestGene",sep=""), paste(strTag,"NearestGeneDistance",sep="")))
	
	aColAdd = names(objGWA.addgene@tblGWA)[ ! names(objGWA.addgene@tblGWA) %in% names(objGWA@tblGWA)]
	objGWA.addgene = GWADATA.getcols(objGWA.addgene, c(colInMarker, aColAdd))
	
	objGWA.addgene <- GWADATA.getrows(objGWA.addgene, which(!duplicated(objGWA.addgene@tblGWA)))
	
	objGWA <- GWADATA.merge(objGWA,objGWA.addgene, strSuffix.In = "", strSuffix.Add = "", blnAll.In = TRUE, blnAll.Add = TRUE, strBy.In = colInMarker, strBy.Add = colInMarker)
	
	return(objGWA)
}

ADDGENE <- function(strEqcCommand){ 
	## Wrapper for class definition
	ADDGENEout <- setADDGENE(new("ADDGENE", strEqcCommand = strEqcCommand))
	validADDGENE(ADDGENEout)
	#ADDGENEout.valid <- validADDGENE(ADDGENEout)
	return(ADDGENEout)
}

