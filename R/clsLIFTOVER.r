setClass("LIFTOVER",
	representation = representation(
						strEqcCommand		=	"character",
						colInChr			=	"character",
						colInPos			=	"character",
						fileLiftOver  		=   "character",
						tblLift				=	"data.table",
						colLiftFromChr		=	"character",
						colLiftFromPos		=	"character",
						colLiftToChr		=	"character",
						colLiftToPos		=	"character",
						strBuildFrom		=	"character",
						strBuildTo			=	"character",
						blnCheckLift		=	"logical",
						blnLift				=	"logical",
						blnRead				=	"logical",
						blnUseFastRead		=	"logical",
						strTag				=	"character"
						),
	prototype = prototype(
						strEqcCommand		=	"",
						colInChr			=	"",
						colInPos			=	"",
						fileLiftOver  		=   "",
						tblLift				=	data.table(),
						colLiftFromChr		=	"chr37",
						colLiftFromPos		=	"pos37",
						colLiftToChr		=	"chr38",
						colLiftToPos		=	"pos38",
						strBuildFrom		=	"37",
						strBuildTo			=	"38",
						blnCheckLift		=	TRUE,
						blnLift				=	TRUE,
						blnRead				=	FALSE,
						blnUseFastRead		=	TRUE,
						strTag				=	""
						)
)

LIFTOVER.set <- function(strEqcCommand, objLO, objGWA.default) {
	
	objLO@blnUseFastRead <- objGWA.default@blnUseFastRead
	
	aEqcSlotNamesIn = c("colInChr","colInPos","fileLiftOver","colLiftFromChr","colLiftFromChr","colLiftFromPos","colLiftToChr","colLiftToPos","blnCheckLift","strBuildFrom","strBuildTo","strTag","blnUseFastRead")
	
	## astrPatterns
	## *_CHR_POS
	## chrCHR:POS
	
	### Last 4 are inherited from class GWADATA and can be used with LIFTOVER for reference file!
	
	objEqcReader <- EqcReader(strEqcCommand,aEqcSlotNamesIn)
	
	if(length(objEqcReader@lsEqcSlotsOut) > 0) {
		for(i in 1:length(objEqcReader@lsEqcSlotsOut)) {
			tmpSlot <- names(objEqcReader@lsEqcSlotsOut)[i]
			tmpSlotVal <- objEqcReader@lsEqcSlotsOut[[i]]
			
			if(all(!is.na(tmpSlotVal))) slot(objLO, tmpSlot) <- tmpSlotVal
		}
	}
	
	return(objLO)
}

#############################################################################################################################

LIFTOVER.GWADATA.valid <- function(objLO, objGWA) {
	
	
	if(!(objLO@colInChr %in% objGWA@aHeader)) stop(paste(" EASY ERROR:LIFTOVER\n Defined column colInChr \n",objLO@colInChr, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct --colInChr OR remove --colInchr and --colInPos.", sep=""))
	if(!(objLO@colInPos %in% objGWA@aHeader)) stop(paste(" EASY ERROR:LIFTOVER\n Defined column colInPos \n",objLO@colInPos, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct --colInPos OR remove --colInchr and --colInPos.", sep=""))
	
	if(!file.exists(objLO@fileLiftOver)) stop(paste("EASY ERROR:LIFTOVER\n Liftover file \n ",objLO@fileLiftOver,"\n does not exist!!!\n", sep=""))
	
	tlift10 = read.table(objLO@fileLiftOver,header=T,nrows=10,stringsAsFactors=FALSE)
	if(!(objLO@colLiftFromChr %in% names(tlift10))) stop(paste(" EASY ERROR:LIFTOVER\n Defined column colLiftFromChr \n",objLO@colLiftFromChr, "\n is not available in Liftover file \n",objLO@fileLiftOver,"\n PLease specify correct --colLiftFromChr.", sep=""))
	if(!(objLO@colLiftFromPos %in% names(tlift10))) stop(paste(" EASY ERROR:LIFTOVER\n Defined column colLiftFromPos \n",objLO@colLiftFromPos, "\n is not available in Liftover file \n",objLO@fileLiftOver,"\n PLease specify correct --colLiftFromPos.", sep=""))
	if(!(objLO@colLiftToChr %in% names(tlift10))) stop(paste(" EASY ERROR:LIFTOVER\n Defined column colLiftToChr \n",objLO@colLiftToChr, "\n is not available in Liftover file \n",objLO@fileLiftOver,"\n PLease specify correct --colLiftToChr.", sep=""))
	if(!(objLO@colLiftToPos %in% names(tlift10))) stop(paste(" EASY ERROR:LIFTOVER\n Defined column colLiftToPos \n",objLO@colLiftToPos, "\n is not available in Liftover file \n",objLO@fileLiftOver,"\n PLease specify correct --colLiftToPos.", sep=""))
	
	
	return(objLO)
	
}
LIFTOVER.check <- function(objLO, isValidScript, objGWA, objREPORT) {
	
	colInChr			<- objLO@colInChr
	colInPos			<- objLO@colInPos
	fileLiftOver  		<- objLO@fileLiftOver
	colLiftFromChr		<- objLO@colLiftFromChr
	colLiftFromPos		<- objLO@colLiftFromPos
	colLiftToChr		<- objLO@colLiftToChr
	colLiftToPos		<- objLO@colLiftToPos
	strBuildFrom 		<- objLO@strBuildFrom
	strBuildTo 			<- objLO@strBuildTo
	blnCheckLift		<- objLO@blnCheckLift
	strTag				<- objLO@strTag
	
	
	blnUseFastRead <- objLO@blnUseFastRead
	
	### check whether data should be lifted 
	
	if(isValidScript) {
			
		if(blnCheckLift) {
			
			### check whether data should be lifted
			objGWAcp = GWADATA.getcols(objGWA, c(colInChr,colInPos))
			tcp = objGWAcp@tblGWA
			tcp[[colInChr]] = as.character(tcp[[colInChr]])
			tcp[[colInPos]] = as.character(tcp[[colInPos]])
			
			tcp[[colInChr]] = gsub("^0+","",tcp[[colInChr]])
			tcp[[colInPos]] = gsub("^0+","",tcp[[colInPos]])
			
			rm(objGWAcp)
			
			tlift1000 = as.data.table(read.table(fileLiftOver, header=T, stringsAsFactors=FALSE, strip.white = TRUE, comment.char = "", colClasses = "character", sep = "\t", nrows = 1000))
			
			t1000testFrom = merge(tlift1000,tcp,by.x=c(colLiftFromChr,colLiftFromPos),by.y=c(colInChr,colInPos),all=F)
			t1000testTo = merge(tlift1000,tcp,by.x=c(colLiftToChr,colLiftToPos),by.y=c(colInChr,colInPos),all=F)
			
			numVarFrom = nrow(t1000testFrom)
			numVarTo = nrow(t1000testTo)
			prVarFrom_1000informative = numVarFrom/1000
			prVarTo_1000informative = numVarTo/1000
			
			objREPORT <- REPORT.addval(objREPORT,paste0("LIFTOVER.",strTag,"prVar",strBuildFrom,"_among_1000informative"),prVarFrom_1000informative)
			objREPORT <- REPORT.addval(objREPORT,paste0("LIFTOVER.",strTag,"prVar",strBuildTo,"_among_1000informative"),prVarTo_1000informative)
			
			if(numVarFrom<numVarTo) {
				objLO@blnLift = FALSE
				warning(paste("LIFTOVER: No variants will be lifted because", prVarTo_1000informative*100, "% of 1000 informative variants are already on build",strBuildTo," comparerd to",prVarFrom_1000informative*100,"% on build",strBuildFrom,"in file",objGWA@fileIn," !!! "))
			}
			
			if(numVarFrom==0 & numVarTo==0) {
				objLO@blnLift = TRUE
				warning(paste("LIFTOVER: No variants available for blnCheckLift. All variants will be lifted to build ",strBuildTo,"."))				
			}
			
		}
		
		# if(objLO@blnLift & !objLO@blnRead) {
			# cat(paste("   + Reading ",fileLiftOver, "... \n"))
			# if(!blnUseFastRead) {
				# objLO@tblLift 	<- as.data.table(read.table(fileLiftOver, header=T, sep = "\t", stringsAsFactors=FALSE, strip.white = TRUE, comment.char = "", colClasses = "character"))
			# } else {
				
				# if(tolower(substring(fileLiftOver,nchar(fileLiftOver)-2,nchar(fileLiftOver)))==".gz") strFread <- paste("zcat ",fileLiftOver,sep="")
				# else strFread <- fileLiftOver
						
				# objLO@tblLift 	<- tryCatch(
					# #fread(strFread, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses = c("character","character","integer")),
					# fread(strFread, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses = c("character","character","character","character")),
					# error = function(err) {
						# strError = err$message
						# val=strsplit(strError,"'",fixed=T)[[1]][length(strsplit(strError,"'",fixed=T)[[1]])]
						# g=scan(file = fileLiftOver, what=character(0), n = -1,sep = "\n",quiet=TRUE)
						# iRow = which(grepl(paste(val,"\t",sep=""),g,fixed=T) | grepl(paste(val,"\n",sep=""),g,fixed=T))[1]
						# stop(paste(strError,"\n EASY ERROR:\n Cannot read '",val,"' from row '",iRow,"' !!!\n ", sep=""))
					# }
				# )
			# }
			# objLO@blnRead <- TRUE
			# ## file only be read once
			
			# if(dim(objLO@tblLift)[1]==0) stop(paste("EASY ERROR:LIFTOVER\n There are no rows available in \n",fileLiftOver,"\n The file is empty!!!\n", sep=""))
			
		# } 
	} 
	# else {
		# objLO@tblLift = as.data.table(read.table(fileLiftOver, header=T, stringsAsFactors=FALSE, strip.white = TRUE, comment.char = "", colClasses = "character", sep = "\t", nrows = 10))
	# }
	
	
	lsOut <- list(objLO, objREPORT)
	return(lsOut)

}
LIFTOVER.read <- function(objLO, isValidScript, objGWA, objREPORT) {
	
	colInChr			<- objLO@colInChr
	colInPos			<- objLO@colInPos
	fileLiftOver  		<- objLO@fileLiftOver
	colLiftFromChr		<- objLO@colLiftFromChr
	colLiftFromPos		<- objLO@colLiftFromPos
	colLiftToChr		<- objLO@colLiftToChr
	colLiftToPos		<- objLO@colLiftToPos
	strBuildFrom 		<- objLO@strBuildFrom
	strBuildTo 			<- objLO@strBuildTo
	blnCheckLift		<- objLO@blnCheckLift
	strTag				<- objLO@strTag
	
	
	blnUseFastRead <- objLO@blnUseFastRead
	
	### check whether data should be lifted 
	
	if(isValidScript) {
			
		
		if(objLO@blnLift & !objLO@blnRead) {
			cat(paste("   + Reading ",fileLiftOver, "... \n"))
			if(!blnUseFastRead) {
				objLO@tblLift 	<- as.data.table(read.table(fileLiftOver, header=T, sep = "\t", stringsAsFactors=FALSE, strip.white = TRUE, comment.char = "", colClasses = "character"))
			} else {
				
				if(tolower(substring(fileLiftOver,nchar(fileLiftOver)-2,nchar(fileLiftOver)))==".gz") strFread <- paste("zcat ",fileLiftOver,sep="")
				else strFread <- fileLiftOver
						
				objLO@tblLift 	<- tryCatch(
					#fread(strFread, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses = c("character","character","integer")),
					fread(strFread, sep="\t", header=TRUE, stringsAsFactors=FALSE, colClasses = c("character","character","character","character")),
					error = function(err) {
						strError = err$message
						val=strsplit(strError,"'",fixed=T)[[1]][length(strsplit(strError,"'",fixed=T)[[1]])]
						g=scan(file = fileLiftOver, what=character(0), n = -1,sep = "\n",quiet=TRUE)
						iRow = which(grepl(paste(val,"\t",sep=""),g,fixed=T) | grepl(paste(val,"\n",sep=""),g,fixed=T))[1]
						stop(paste(strError,"\n EASY ERROR:\n Cannot read '",val,"' from row '",iRow,"' !!!\n ", sep=""))
					}
				)
			}
			objLO@blnRead <- TRUE
			## file only be read once
			
			if(dim(objLO@tblLift)[1]==0) stop(paste("EASY ERROR:LIFTOVER\n There are no rows available in \n",fileLiftOver,"\n The file is empty!!!\n", sep=""))
			
		} 
	} else {
		objLO@tblLift = as.data.table(read.table(fileLiftOver, header=T, stringsAsFactors=FALSE, strip.white = TRUE, comment.char = "", colClasses = "character", sep = "\t", nrows = 10))
	}
	
	
	lsOut <- list(objLO, objREPORT)
	return(lsOut)
	
}

#############################################################################################################################
LIFTOVER.run <- function(objLO, objGWA, objREPORT, isValidScript) {
	
	colInChr			<- objLO@colInChr
	colInPos			<- objLO@colInPos
	fileLiftOver  		<- objLO@fileLiftOver
	colLiftFromChr		<- objLO@colLiftFromChr
	colLiftFromPos		<- objLO@colLiftFromPos
	colLiftToChr		<- objLO@colLiftToChr
	colLiftToPos		<- objLO@colLiftToPos
	strBuildFrom 		<- objLO@strBuildFrom
	strBuildTo 			<- objLO@strBuildTo
	strTag				<- objLO@strTag
	
	blnLift 			<- objLO@blnLift
	
	if(nchar(strTag)>0) strTag = paste(strTag,".",sep="")
	
	if(isValidScript) {
		# cpid = paste(objGWA@tblGWA[[colInChr]],objGWA@tblGWA[[colInPos]],sep=":")
	
		
		if(blnLift) {
			# remap all to hg38
			cat("Lifting positions ... \n")
			
			objGWAcp = GWADATA.getcols(objGWA, c(colInChr,colInPos))
			tcp = objGWAcp@tblGWA
			
			tcp[[colInChr]] = as.character(tcp[[colInChr]])
			tcp[[colInPos]] = as.character(tcp[[colInPos]])
			
			tcp[[colInChr]] = gsub("^0+","",tcp[[colInChr]])
			tcp[[colInPos]] = gsub("^0+","",tcp[[colInPos]])
			rm(objGWAcp)
			
			tcp = merge(tcp,objLO@tblLift,by.x=c(colInChr,colInPos),by.y=c(colLiftFromChr,colLiftFromPos),all.x=T,all.y=F,sort=F)
			
			isLifted = !is.na(tcp[[colLiftToPos]])
						
			if(any(!isLifted)) {
				objGWA.removed = GWADATA.getrows(objGWA, which(!isLifted))
				GWADATA.write(objGWA.removed, strSuffix = paste0(".",strTag,"notlifted"))
				rm(objGWA.removed)
			}
			
			objGWA = GWADATA.getrows(objGWA, which(isLifted))
			
			objGWA@tblGWA[[colInChr]] = tcp[[colLiftToChr]][which(isLifted)]
			objGWA@tblGWA[[colInPos]] = tcp[[colLiftToPos]][which(isLifted)]
			
			objREPORT <- REPORT.addval(objREPORT,paste0("LIFTOVER.",strTag,"numVarLifted_",strBuildFrom,"to",strBuildTo),length(which(isLifted)))
			objREPORT <- REPORT.addval(objREPORT,paste0("LIFTOVER.",strTag,"numVarRemoved"),length(which(!isLifted)))
						
			warning(paste("LIFTOVER:", length(which(isLifted)), "variants from",objGWA@fileInShortName,"were lifted from Hg19 to Hg38.",length(which(!isLifted)),"variants were removed due to unavailabilty in the liftover file."))
			rm(tcp)
			
		} else {
		
			objREPORT <- REPORT.addval(objREPORT,paste0("LIFTOVER.",strTag,"numVarLifted_",strBuildFrom,"to",strBuildTo),0)
			objREPORT <- REPORT.addval(objREPORT,paste0("LIFTOVER.",strTag,"numVarRemoved"),0)
						
		}
		
	} 
		
	lsOut <- list(objGWA, objREPORT)
	
	return(lsOut)
	
	# return(objGWA)
}
#############################################################################################################################
LIFTOVER <- function(strEqcCommand, objGWA.default){ 
	## Wrapper for class definition
	LIFTOVERout <- new("LIFTOVER")
	LIFTOVERout <- LIFTOVER.set(strEqcCommand, LIFTOVERout, objGWA.default)
	
	return(LIFTOVERout)

}
