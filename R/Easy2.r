# rm(list = ls(all = TRUE))
# graphics.off()
# closeAllConnections()
# options(show.error.messages = TRUE)
# options(warn=0)

# setwd("/data/M_GENETICS/ScriptsDBDebug/easy_x/git/250328_addgene/")

# #fileECF	<-	"2_cpaid.ecf"
# # fileECF	<-	"rbind_maf001.ecf"
# fileECF	<-	"5_3_add_genes.ecf"

# pathOut=NA
# fileMerge=NA

# pathClasses <- 	"/data/M_GENETICS/ScriptsDBDebug/easy_x/git/Easy2/R/"
# blnLogMemoryExtern <- FALSE
# blnValidityCheckOnly <- FALSE
# blnReturnGwadata <- FALSE
# blnReturnReport <- FALSE

# aFileIn=c()

# library(Cairo)
# library(plotrix)
# library(data.table)
# library(forestplot)
# library(bigsnpr)
# library(foreach)
# library(doParallel)
# # library(MatrixExtra)

# source(paste(pathClasses,"clsSPLOT.r",sep=""))
# arfun = list.files(pathClasses,pattern=".r")
# arfun = arfun[arfun != "Easy2.r"]
# for (rfun in arfun) source(paste(pathClasses,rfun,sep=""))

# # source("clsCREATECPAID.r")

# if(!file.exists(fileECF))
	# stop(paste("EASY ERROR:\n ECF-file \n ",fileECF,"\n does not exist!!!\n", sep=""))

# fileEcfOut<-paste(fileECF,".out",sep="")
# connection.Logfile <- file(fileEcfOut,open='w')
# sink(connection.Logfile, split=TRUE)


# #fileMemoryLogExtern <- paste(fileECF,".mem",sep="")
# shMemLog<-"/d/GENETICS/ScriptsDBDebug/easy_x/bin/fnLogMem.sh"
# if(blnLogMemoryExtern) 
	# system2(shMemLog, args=c(Sys.getpid(),Sys.getenv("USER"),fileEcfOut), wait=FALSE)
	
# # # ### NEED TO SET STOP CONDITION AFTER EASY2 FUN CALL !!!!!!!!!


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
fnSetNumCol <- function(objGWA) {
	if(objGWA@astrSetNumCol[1]!="") {
		for(strAddCol in objGWA@astrSetNumCol) {
			## strAddCol N=200
			aAddColSplit = strsplit(strAddCol, "=", fixed=TRUE)[[1]]
			newcolname=aAddColSplit[1]
			newcolval=aAddColSplit[2]
			objGWA <- GWADATA.cbind(objGWA, rep(as.numeric(newcolval),nrow(objGWA@tblGWA)), newcolname)
		}
	}
	return(objGWA)
}
fnSetFileInTagCol <- function(objGWA) {

	if(objGWA@fileInTag!="") {
		objGWA <- GWADATA.cbind(objGWA, rep(objGWA@fileInTag,nrow(objGWA@tblGWA)), "fileInTag")
	}
	
	return(objGWA)
}



Easy2.run <- function(fileECF,blnValidityCheckOnly=FALSE,blnReturnGwadata=FALSE,blnReturnReport=FALSE,aFileIn=c(),fileMerge=NA,pathOut=NA){ 
	
	cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
	cat("|       Easy2	   |      1.3.4.b42  |     23/October/2025    |\n")
	cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
	cat("|  (C) 2013 Thomas Winkler, GNU General Public License, v3   |\n")
	cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
	cat("|  For bug-report, please e-mail:                            |\n")
	cat("|  thomas.winkler@ur.de         					          |\n")
	cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
	#### b3 everytime something changes the build gets up
	####   eg bugfix 
	#### 1.X.Y -> Y: new parameter or new function added or bugfix
	#### 1.X.Y -> X: completely newly exported 
	#### 1.1.1.b5 improved LIFTOVER class and Easy2 call
	#### 1.1.2.b5 added isregionExt to INDEPX
	####          added colInCoord to ANNOTATE
	#### 1.1.8.b6 added fixed=T to STRSPLITCOL
	#### 1.1.8.b7 added PCA to INDEPX / deprecated !
	#### 1.1.9 added PCA2STEP
	#### 1.1.12 add +1 to PCASTEP
	#### 1.1.13 debuged CLUMPX
	#### b15: added annot and gene to INDEP
	#### b17: added annot to MHPLOT
	#### b18: added --colInGPos and signaling to INDEP
	#### b19: signalleads output
	#### b20: signalleads debug
	#### b22: updated INDEP
	#### b23: debuged LIFTOVER --blnCheckLift no variants
	#### b24: do not run JOINTTEST in validity check
	#### b25: add colInSkipAnnot to ANNOTATE
	#### b26: add maf filter for Easyrbind
	#### b27: add paramter fileMerge
	#### b28: fix paral issue with missing lib on worker
	#### b29: add tempdir to fread
	#### b30: correct tempdir default value to fread, correct blnHeader EASYMERGE
	#### b31: remove dups from ADDGENE
	#### b32: fix empty forestplot by adding print(forestplot(...))
	#### 1.2.6.b33: add --strMarkerTag to GETVAR
	#### 1.2.7.b34: add X chr to INDEP
	#### 1.2.8.b35: allow "X" in manhattan plot (recode X to 23)
	#### 1.2.10.b37: bugfix chr x check 
	#### 1.3.0.b38: add rbind tag to GWADATA; add one-column tags to INDEP
	#### 1.3.1.b39: mergeeasysin using isMergedSet
	#### 1.3.4.b42: INDEP blnRetainRegionLead to retain regionLead in locusLead when no region variant is in refpanel
	
	
	cat("\n")
	cat("+++++\n")
	.timeStart <- Sys.time()
	cat(paste("Starting Easy2:", .timeStart,"\n"))
	cat(paste("Running script ",fileECF, "...\n"))
	
	idxDeviceOffsetIn <- as.integer(dev.cur()) - 1
	
	fileECFName <- strsplit(fileECF,"/")[[1]][length(strsplit(fileECF,"/")[[1]])]
	fileEcfOut  <- paste(fileECF,".out",sep="")
	##########################
	#print("Reading Easy2 config-data...")
	#aEcfRows	<- scan(file = fileECF, what=character(0), n = -1, comment.char = "#",  strip.white = TRUE,sep = "\n",quiet=TRUE)
	aEcfRows	<- scan(file = fileECF, what=character(0), n = -1, comment.char = "",  strip.white = TRUE,sep = "\n",quiet=TRUE)
	isRemove = substring(aEcfRows,1,1)=="#"
	aEcfRows = aEcfRows[!isRemove]

	#iScript = match(c("START EASYQC", "STOP EASYQC"), aEcfRows)

	iStart = which(aEcfRows == "START EASYX" | aEcfRows == "START EASY2" | aEcfRows == "START EASYQC" | aEcfRows == "START EASYSTRATA" | aEcfRows == "START EASYQC2" | aEcfRows == "START EASYSTRATA2")
	iMergeEasyin = which(substring(aEcfRows,1,11) == "MERGEEASYIN")
	iMergeStrat = which(substring(aEcfRows,1,10) == "MERGESTRAT")
	iRbindTrait = which(substring(aEcfRows,1,11) == "RBINDTRAITS")
	iStop = which(aEcfRows == "STOP EASYX" | aEcfRows == "STOP EASY2" | aEcfRows == "STOP EASYQC" | aEcfRows == "STOP EASYSTRATA" | aEcfRows == "STOP EASYQC2" | aEcfRows == "STOP EASYSTRATA2")

	#if(any(is.na(iScript)) | iScript[2] < iScript[1]) 
	if(!all(length(iStart)==1 & length(iStop)==1 & iStart < iStop))
		stop(paste("EASY ERROR:\n 'START EASYX' and 'STOP EASYX' not properly defined in \n  ",fileECFName,"\n !!!\n", sep=""))
	if(length(iMergeEasyin) > 0) {
		if(!all(length(iMergeEasyin) == 1 & iMergeEasyin > iStart & iMergeEasyin < iStop))
			stop(paste("EASY ERROR:\n 'MERGEEASYIN' not properly defined in \n  ",fileECFName,"\n This command may only be stated once !!!\n", sep=""))
	}
	if(length(iMergeStrat) > 0) {
		if(!all(length(iMergeStrat) == 1 & iMergeStrat > iStart & iMergeStrat < iStop))
			stop(paste("EASY ERROR:\n 'MERGESTRAT' not properly defined in \n  ",fileECFName,"\n This command may only be stated once !!!\n", sep=""))
	}
	if(length(iRbindTrait) > 0) {
		if(!all(length(iRbindTrait) == 1 & iRbindTrait > iStart & iRbindTrait < iStop))
			stop(paste("EASY ERROR:\n 'RBINDTRAITS' not properly defined in \n  ",fileECFName,"\n This command may only be stated once !!!\n", sep=""))
	}
	cat("\n")
	cat("+++++\n")
	cat("Getting list of input files ...\n")
	
	#if(iScript[1] == 1)
	#	stop(paste("EASY ERROR:\n No DEFINE and no EASYIN argument found in configuration part of ECF-file \n ",fileECFName,"\n !!!\n", sep=""))
	
	### Get aEcfConfigCommands
	
	aindConfigRows = c(1,iStart-1)
	aEcfConfigRows <- aEcfRows[aindConfigRows[1]:aindConfigRows[2]]
	strEcfConfigRows=paste(aEcfConfigRows,collapse="\n")
	strEcfConfigRows=gsub("\t"," ",strEcfConfigRows)
	strEcfConfigRows.ByCommand=gsub("\n--"," --",strEcfConfigRows)
	while(grepl("  ", strEcfConfigRows.ByCommand)) strEcfConfigRows.ByCommand <- sub("  "," ",strEcfConfigRows.ByCommand)
	aEcfConfigCommands=strsplit(strEcfConfigRows.ByCommand,"\n")[[1]]
	
	### Get aEqcScriptCommands
	aindScriptRows = c(iStart+1,iStop-1)
	aEqcScriptRows <- aEcfRows[aindScriptRows[1]:aindScriptRows[2]]
	strEqcScriptRows=paste(aEqcScriptRows,collapse="\n")
	strEqcScriptRows=gsub("\t"," ",strEqcScriptRows)
	strEqcScriptRows.ByCommand=gsub("\n--"," --",strEqcScriptRows)
	while(grepl("  ", strEqcScriptRows.ByCommand)) strEqcScriptRows.ByCommand <- sub("  "," ",strEqcScriptRows.ByCommand)
	aEqcScriptCommands=strsplit(strEqcScriptRows.ByCommand,"\n")[[1]]
	
	#### Check correct usage of functions
	astrEasyFuns.used <- unlist(lapply(strsplit(aEqcScriptCommands," "),function(x) x[1]))
	
	astrEasyFuns.allowed <- c("MERGE", "MERGEEASYIN","SPLOT", "ADDCOL","ADJUSTALLELES",
								"CALCULATE","CLEAN","CLEANDUPLICATES","CRITERION","EDITCOL","EVALSTAT","EXTRACTSNPS","FILTER","JOINTTEST","METAANALYSIS","SORT","STRSPLITCOL","WRITE",
								"FLIPSTRAND","GC","GETCOLS","GETNUM","QQPLOT","REMOVECOL","REMOVESNPS","RENAMECOL","RPLOT",
								"RENAMEMARKER","CREATECPTID","HARMONIZEALLELES","AFCHECK","PZPLOT","CREATECPAID","ADJUSTALLELESSIMPLE","RADDCOL","LOADINFO","GETVARTYPE","BRPLOT", 
								"ANNOTATE","BONFERRONI","CALCPDIFF","CALCPHET","CLUMP","FDR","INDEP","MHPLOT","MIAMIPLOT","FANCYPLOT","MULTINDEP","MULTCLUMP","EFFECTPLOT","ENRICHMENT",
								"GETVAR", "FORESTPLOT","WEIGHTEDHYPOTHESIS","ADDGENE","ADDGMAP","INDEPX","GWASCATALOG","INDEP2","CALCPDIFFDIFF","LIFTOVER", "PCA2STEP", "CLUMPX")
	
	isEasyFunMatch = astrEasyFuns.used%in%astrEasyFuns.allowed
	if(any(!isEasyFunMatch)) 
		stop(paste("EASY ERROR:\n Functions \n ",paste(astrEasyFuns.used[which(!isEasyFunMatch)],collapse=","),"\n are not allowed Easy2 functions. Please revise or remove the function !!!\n", sep=""))
	
	
	ls_aEqcScriptCommands = list()
	aEqcScriptCodeType = c()
	iStartRegular = 1

	iMergeCommand = which(substring(aEqcScriptCommands,1,11) == "MERGEEASYIN")
	isMergeCommand = length(iMergeCommand)==1
	nEqcCommands = length(aEqcScriptCommands)
	
	if(isMergeCommand) {

		ls_aEqcScriptCommands = list(aEqcScriptCommands[1:iMergeCommand])
		aEqcScriptCodeType = "GWA"
		if(nEqcCommands > iMergeCommand) {
			ls_aEqcScriptCommands = c(ls_aEqcScriptCommands, list(aEqcScriptCommands[(iMergeCommand+1):nEqcCommands]))
			aEqcScriptCodeType = c(aEqcScriptCodeType, "MERGED")
		}
	
	} else {
			ls_aEqcScriptCommands = list(aEqcScriptCommands[1:nEqcCommands])
			aEqcScriptCodeType = "GWA"
	}
	## EASYMERGE will be applied when GWADATA is read
	
	################################################################################################################################################
	################################################################################################################################################
	
	container.GWADATA <- container.EASYMERGE <- container.EASYRBIND <- list()
	objGWADATA.default 		<- GWADATA()
	objEASYMERGE.default 	<- EASYMERGE()
	objEASYRBIND.default 	<- EASYRBIND()
	isEasyMerge = FALSE
	isEasyRbind = FALSE
	icount.GWADATA <- 0
	
	for(iConfigCommand in 1:length(aEcfConfigCommands)) {

		strConfigCommand = aEcfConfigCommands[iConfigCommand]
		strConfigFun = strsplit(strConfigCommand," ")[[1]][1]
		
		if(strConfigFun == "DEFINE") {
			objGWADATA.default 			<- GWADATA.define(objGWADATA.default, strConfigCommand, pathOut)
			objEASYMERGE.default 		<- EASYMERGE.define(objEASYMERGE.default, strConfigCommand)
			objEASYRBIND.default 		<- EASYRBIND.define(objEASYRBIND.default, strConfigCommand)
		}
		if(strConfigFun == "EASYIN") {
			objGWADATA 	<- GWADATA.easyin(objGWADATA.default, strConfigCommand)	
			
			if(objGWADATA@fileInType == "METALSCRIPT") 		ls_objGWADATA <- GWADATA.getmetalscriptfiles(objGWADATA)
			else if(objGWADATA@fileInType == "FILELIST") 	ls_objGWADATA <- GWADATA.getfilelistfiles(objGWADATA)
			else if(grepl("*",objGWADATA@fileIn, fixed=T))	ls_objGWADATA <- GWADATA.getfiles(objGWADATA)
			else ls_objGWADATA <- list(objGWADATA)
			
			container.GWADATA = c(container.GWADATA , ls_objGWADATA)
			icount.GWADATA = icount.GWADATA + 1
		}
		# EASYMERGE
		if(strConfigFun == "EASYMERGE") {
			objEASYMERGE <- EASYMERGE.easymerge(objEASYMERGE.default, strConfigCommand, icount.GWADATA)	
			ls_objEASYMERGE <- list(objEASYMERGE)
			container.EASYMERGE = c(container.EASYMERGE , ls_objEASYMERGE)
			isEasyMerge = TRUE
		}
		if(strConfigFun == "EASYRBIND") {
			objEASYRBIND <- EASYRBIND.easyrbind(objEASYRBIND.default, strConfigCommand, icount.GWADATA)	
			ls_objEASYRBIND <- list(objEASYRBIND)
			container.EASYRBIND = c(container.EASYRBIND , ls_objEASYRBIND)
			isEasyRbind = TRUE
		}
	}
	
	if(length(aFileIn)>0) {
		for(iFileIn in 1:length(aFileIn)) {
			objGWADATA 	<- objGWADATA.default
			objGWADATA@fileIn <- aFileIn[iFileIn]
			container.GWADATA = c(container.GWADATA , list(objGWADATA))
			icount.GWADATA = icount.GWADATA + 1
		}
	}
	
	nGWA <- length(container.GWADATA)
	if(isEasyMerge) aiEasyMergeIDs = unlist(lapply(container.EASYMERGE, function(x) x@iMergeID))
	if(isEasyRbind) aiEasyRbindIDs = unlist(lapply(container.EASYRBIND, function(x) x@iRbindID))
	
	if(nGWA < 1) stop(paste("EASY ERROR:\n No EASYIN argument found in configuration part of ECF-file \n ",fileECFName,"\n PLease specify at least one input file !!!\n", sep=""))
	
	cat("Using:\n")
	for(i in 1:nGWA) {
		#stop()
		container.GWADATA[[i]] <- GWADATA.init(container.GWADATA[[i]])
		cat(paste("   + ",container.GWADATA[[i]]@fileIn,"\n"))
		if(isEasyMerge) {
			if(any(i == aiEasyMergeIDs)) {
				aiMatch = which(aiEasyMergeIDs == i)
				for(iMatch in aiMatch) {
					container.EASYMERGE[[iMatch]] <- EASYMERGE.init(container.EASYMERGE[[iMatch]])
					cat(paste("   + EASYMERGE ",container.EASYMERGE[[iMatch]]@fileIn,"\n"))
				}
			}
		}
		if(isEasyRbind) {
			if(any(i == aiEasyRbindIDs)) {
				aiMatch = which(aiEasyRbindIDs == i)
				for(iMatch in aiMatch) {
					if(container.EASYRBIND[[iMatch]]@blnHeader) container.EASYRBIND[[iMatch]] <- EASYRBIND.init(container.EASYRBIND[[iMatch]])
					else container.EASYRBIND[[iMatch]] <- EASYRBIND.init.nohead(container.EASYRBIND[[iMatch]], container.GWADATA[[i]])
					cat(paste("   + EASYRBIND ",container.EASYRBIND[[iMatch]]@fileIn,"\n"))
				}
			}
		}
	}
	
	### PathOut set to pathOut from first GWADATA
	pathOut <- ifelse(objGWADATA.default@pathOut != getwd(), objGWADATA.default@pathOut, container.GWADATA[[1]]@pathOut)
	
	#pathOut <- container.GWADATA[[1]]@pathOut
	# fileOutBody <- paste(pathOut,"/",fileECFName,sep="")
	# fileOutBody <- gsub(".ecf", "", fileOutBody)
	
	fileOutBodyIn <- paste(pathOut,"/",fileECFName,sep="")
	#fileOutBodyIn <- gsub(".ecf", "", fileOutBodyIn)
	if(substring(fileOutBodyIn, nchar(fileOutBodyIn)-3, nchar(fileOutBodyIn))==".ecf") 
		fileOutBodyIn <- substring(fileOutBodyIn,1,nchar(fileOutBodyIn)-4)
	
	cat("\n")
	cat("+++++\n")
	cat("Default output path is \n")
	cat(pathOut)
	cat("\n")
		
	#	return(list(objGWADATA.default, container.GWADATA, fileOutBody))

	#}
	
	
	########################
	################################################################################################################################################
	################################################################################################################################################
	
	cat("\n")
	cat("+++++\n")
	cat("Performing validity check on 10 rows from each file :\n")
	#stop()
	#aEqcScriptCommandsTmp=ls_aEqcScriptCommands[[1]]
	
	blnAddFileInTagToReport <- length(unique(unlist(lapply(container.GWADATA,function(x) x@fileInTag)))) > 1
	
	container.REPORT <- container.container.ReportPlots <- container.container.BReportPlots <- list()
	
	## each element refers to a command group
	
	################################################################################################################################################
	################################################################################################################################################
	########################

# stop()	
# iRun <- iMergeTag <- iCommandGroup <- iGWAmt <- iCommand <- 1


	for(iRun in 1:2) {
		
		if(iRun == 2 & blnValidityCheckOnly) break
		
		aMergeTags = unlist(lapply(container.GWADATA,function(x) x@fileInMergeTag))
		aMergeTagsUni = unique(aMergeTags)
		nMergeTags = length(aMergeTagsUni)
		
		container.Merge <- container.Rename <- container.CptidMap <- container.LoadInfo <- container.LiftOver <- list()
		
		iGWA <- 0
		nSubplots <- 0
		
		for(iMergeTag in 1:nMergeTags) {
			
			mergeTag = aMergeTagsUni[iMergeTag]
			
			## Reduce to GWADATA objects for current mergetag
			container.GWADATA.mergetag <- container.GWADATA[which(aMergeTags==mergeTag)]
			
			if(iRun == 1) {
				isValidScript = FALSE
			} else {
				cat("\n Passed validity check!\n")
				isValidScript = TRUE
			}
						
			nCommandGroups = length(ls_aEqcScriptCommands)
			
			idxMultiplotMerged <- idxMERGEMerged <- idxRENAMEMerged <- idxCPTIDMerged <- idxLOADINFOMerged <- idxLOMerged <- 0
			idxDeviceOffset <- idxDeviceOffsetIn
			
			for(iCommandGroup in 1:nCommandGroups) {
					
				aEqcScriptCommands = unlist(ls_aEqcScriptCommands[iCommandGroup])
				EqcScriptCodeType = aEqcScriptCodeType[iCommandGroup]
				
				if(isValidScript) {
					if(EqcScriptCodeType == "MERGED") {
						fileOutBody <- paste(fileOutBodyIn, ".merged",sep="")
					} else {
						fileOutBody <- fileOutBodyIn
					}
				} else {
					fileOutBody <- fileOutBodyIn
				}
					
				if(iMergeTag==1) {
					## initiate new report
					objREPORT <- REPORT(fileOutBody,objGWADATA.default@blnOverwriteResults)
				} else { 
					## add to existing report; there will be one per commandgroup! 
					objREPORT <- container.REPORT[[iCommandGroup]]
				}
				
				afileOutMultiPlot = c()
				container.BReportPlots = list()
				container.ReportPlots = list()
				
				container.Boxplots = list() ## [[1]] list of boxplotstat of boxplot 1
			
				if(EqcScriptCodeType=="GWA") {
					nGWAmt <- length(container.GWADATA.mergetag)
				} else {
					## MERGED
					nGWAmt <- 1
				}
				
				# isMergedSet = FALSE
				
				for(iGWAmt in 1:nGWAmt) {
					
					isGarbageCleaning <- ifelse(isValidScript, objGWA@blnGarbageCleaning, FALSE)
					isMemoryLog <- ifelse(isValidScript, objGWA@blnLogMemory, FALSE)
					isTimeLog <- ifelse(isValidScript, objGWA@blnLogTime, FALSE)
					
					idxBReportPlot <- idxReportPlot <- 0
					idxMultiplotGwaMt <- idxMERGEGwa <- idxRENAMEGwa <- idxCPTIDGwa <- idxLOADINFOGwa <- idxLOGwa <- 0
					
					if(EqcScriptCodeType=="GWA") {
						objGWA <- container.GWADATA.mergetag[[iGWAmt]]
						iGWA = iGWA + 1
					} else {
						## EqcScriptCodeType=="MERGED"
						objGWA <- objGWA.merged
					}
					if(isValidScript) {
						cat("\n")
						cat("+++++\n")
						cat(paste("Processing file:",objGWA@fileInShortName,"\n"))
					} else {
						cat("\n")
						cat(paste("   + ",objGWA@fileInShortName,"-> "))
					}
					.filetimeStart <- Sys.time()
					#if(!(objGWA@blnMergedStrat | objGWA@blnRbindTraits)) {
					# if(!(objGWA@blnMergedEasyin)) {
						# if(isValidScript) 	objGWA <- GWADATA.read(objGWA)
						# else				objGWA <- GWADATA.read.10rows(objGWA)
					# }
					if(EqcScriptCodeType=="GWA") {
						.funtimeStart <- Sys.time()
						if(isValidScript) {
							objGWA <- GWADATA.read(objGWA)
							if(objGWA@colFreq!="") {
								isKeep = !is.na(objGWA@tblGWA[[objGWA@colFreq]]) & (pmin(objGWA@tblGWA[[objGWA@colFreq]],1-objGWA@tblGWA[[objGWA@colFreq]]) >= objGWA@numMinMaf)
								objGWA = GWADATA.getrows(objGWA, which(isKeep))
							}
						} else {
							objGWA <- GWADATA.read.10rows(objGWA)
						}
						
						objGWA <- fnSetNumCol(objGWA)
						# objGWA <- fnSetFileInTagCol(objGWA)
						
						# objGWA <- fnSetTag(objGWA)

						# if(objGWA@astrSetNumCol[1]!="") {
							# for(strAddCol in objGWA@astrSetNumCol) {
								# ## strAddCol N=200
								# aAddColSplit = strsplit(strAddCol, "=", fixed=TRUE)[[1]]
								# newcolname=aAddColSplit[1]
								# newcolval=aAddColSplit[2]
								# objGWA <- GWADATA.cbind(objGWA, rep(as.numeric(newcolval),nrow(objGWA@tblGWA)), newcolname)
							# }
						# }
						
						.funtimeStop <- Sys.time()
						if(isTimeLog) 
							cat(paste("\t** Data reading time:", round(as.numeric(difftime(.funtimeStop,.funtimeStart,units="min")),2),"Minutes\n"))
							#cat(paste("\t** Function evaluation time:", round(as.numeric(difftime(.funtimeStop,.funtimeStart,units="min")),2),"Minutes\n"))
						if(isMemoryLog) {
							# cat(paste("\t** Current memory allocated:",system(paste("top -b -n 1 -p ",Sys.getpid()," | grep \'",paste(Sys.getpid(),Sys.getenv("USER")),"\' | awk '{print $5}'", sep=""),intern=TRUE),"\n",sep=" "))
							## By Lucho:
							cat(paste("\t** Current memory allocated:",system(paste("ps -o pid,user,pri,nice,vsz ",Sys.getpid()," | grep \'",paste(Sys.getpid(),Sys.getenv("USER")),"\' | awk '{print $5}'", sep=""),intern=TRUE),"\n",sep=" "))
						}
						### EASYMERGE
						if(isEasyMerge) {
							if(any(iGWAmt == aiEasyMergeIDs)) {
								.funtimeStart <- Sys.time()
								aiMatch = which(aiEasyMergeIDs == iGWAmt)
								for(iMatch in aiMatch) {
									objEASYMERGE <- container.EASYMERGE[[iMatch]]
									if(isValidScript) {
										objEASYMERGE <- EASYMERGE.read(objEASYMERGE)
										if(objEASYMERGE@colFreq!="") {
											isKeep = !is.na(objEASYMERGE@tblGWA[[objEASYMERGE@colFreq]]) & (pmin(objEASYMERGE@tblGWA[[objEASYMERGE@colFreq]],1-objEASYMERGE@tblGWA[[objEASYMERGE@colFreq]]) >= objEASYMERGE@numMinMaf)
											objEASYMERGE = GWADATA.getrows(objEASYMERGE, which(isKeep))
										}
									} else {
										objEASYMERGE <- EASYMERGE.read.10rows(objEASYMERGE)
									}
									objEASYMERGE <- fnSetNumCol(objEASYMERGE)
									
									EASYMERGE.GWADATA.valid(objEASYMERGE, objGWA)
									objGWA <- EASYMERGE.run(objEASYMERGE, objGWA)
								}
								.funtimeStop <- Sys.time()
								if(isTimeLog) 
									cat(paste("\t** Data merging time:", round(as.numeric(difftime(.funtimeStop,.funtimeStart,units="min")),2),"Minutes\n"))
								if(isMemoryLog) {
									# cat(paste("\t** Current memory allocated:",system(paste("top -b -n 1 -p ",Sys.getpid()," | grep \'",paste(Sys.getpid(),Sys.getenv("USER")),"\' | awk '{print $5}'", sep=""),intern=TRUE),"\n",sep=" "))	
									# By Lucho:
									cat(paste("\t** Current memory allocated:",system(paste("ps -o pid,user,pri,nice,vsz ",Sys.getpid()," | grep \'",paste(Sys.getpid(),Sys.getenv("USER")),"\' | awk '{print $5}'", sep=""),intern=TRUE),"\n",sep=" "))
								}
							}
						}
						### EASYRBIND
						if(isEasyRbind) {
							if(any(iGWAmt == aiEasyRbindIDs)) {
								.funtimeStart <- Sys.time()
								aiMatch = which(aiEasyRbindIDs == iGWAmt)
								for(iMatch in aiMatch) {
									objEASYRBIND <- container.EASYRBIND[[iMatch]]
									if(isValidScript) {
										objEASYRBIND <- EASYRBIND.read(objEASYRBIND)
										if(objEASYRBIND@colFreq!="") {
											isKeep = !is.na(objEASYRBIND@tblGWA[[objEASYRBIND@colFreq]]) & (pmin(objEASYRBIND@tblGWA[[objEASYRBIND@colFreq]],1-objEASYRBIND@tblGWA[[objEASYRBIND@colFreq]]) >= objEASYRBIND@numMinMaf)
											objEASYRBIND = GWADATA.getrows(objEASYRBIND, which(isKeep))
										}
									} else {
										objEASYRBIND <- EASYRBIND.read.10rows(objEASYRBIND)
									}
									objEASYRBIND <- fnSetNumCol(objEASYRBIND)
									objEASYRBIND <- fnSetFileInTagCol(objEASYRBIND)
									if(all(names(objGWA@tblGWA)!="fileInTag")) objGWA <- fnSetFileInTagCol(objGWA)
									
									objGWA <- EASYRBIND.run(objEASYRBIND, objGWA)
								}
								## clean up too much allocated space
								#if(isValidScript) stop()
								if(objGWA@numRowReadCount > 0 & objGWA@numRowReadCount < nrow(objGWA@tblGWA)) {
									objGWA@tblGWA <- objGWA@tblGWA[-((objGWA@numRowReadCount+1):nrow(objGWA@tblGWA)),]
								}
								.funtimeStop <- Sys.time()
								if(isTimeLog) 
									cat(paste("\t** Data rbind time:", round(as.numeric(difftime(.funtimeStop,.funtimeStart,units="min")),2),"Minutes\n"))
								if(isMemoryLog) {
									# cat(paste("\t** Current memory allocated:",system(paste("top -b -n 1 -p ",Sys.getpid()," | grep \'",paste(Sys.getpid(),Sys.getenv("USER")),"\' | awk '{print $5}'", sep=""),intern=TRUE),"\n",sep=" "))	
									# By Lucho:
									cat(paste("\t** Current memory allocated:",system(paste("ps -o pid,user,pri,nice,vsz ",Sys.getpid()," | grep \'",paste(Sys.getpid(),Sys.getenv("USER")),"\' | awk '{print $5}'", sep=""),intern=TRUE),"\n",sep=" "))
								}
							}
						}
					}
					
					###
					objREPORT	<-	REPORT.addval(objREPORT,"fileInShortName",objGWA@fileInShortName)
					if(blnAddFileInTagToReport) objREPORT	<-	REPORT.addval(objREPORT,"fileInTag",objGWA@fileInTag)
					#objREPORT	<-	REPORT.addval(objREPORT,"fileInTag",objGWA@fileInTag)
					objREPORT	<-	REPORT.addval(objREPORT,"numVarIn",dim(objGWA@tblGWA)[1])
					objREPORT	<-	REPORT.addval(objREPORT,"numVarOut","NA")
					
					#### Go through script commands
					
					nCommands = length(aEqcScriptCommands)
					
					for(iCommand in 1:nCommands) {
						
						# if(isValidScript & iGWAmt==2) stop()
						
						strCommand = aEqcScriptCommands[iCommand]
						strScriptFun = strsplit(strCommand," ")[[1]][1]
						
						.funtimeStart <- Sys.time()
						
						if(isValidScript) {
							cat(paste("   + ",gsub("--","\n      --",strCommand),sep=""))
							cat("\n")		
						}
						if(strScriptFun == "MERGEEASYIN") {
						
							astrSplitfileECFName = strsplit(fileECFName,".",fixed=T)[[1]]
							if(length(astrSplitfileECFName)>1) {
								fileOutShortName = paste(astrSplitfileECFName[-length(astrSplitfileECFName)],collapse=".")				
							} else {
								fileOutShortName = fileECFName
							}
							
							if(objGWA@fileInMergeTag != "") 
								fileOutShortName = paste(fileOutShortName, objGWA@fileInMergeTag, sep=".")
							
							objME <- MERGEEASYIN(strCommand, fileOutShortName)
							MERGEEASYIN.GWADATA.valid(objME, objGWA)
							if(objME@blnRbind) objGWA <- fnSetFileInTagCol(objGWA)
							
							if(iGWAmt == 1) {
							### start first
								objGWA.merged <- MERGEEASYIN.start(objME, objGWADATA.default, objGWA)
							} 
							
							if(iGWAmt > 1 & objGWA.merged@fileInMergeTag == objGWA@fileInMergeTag) {	### add 
								objGWA.merged <- MERGEEASYIN.run(objME, objGWA.merged, objGWA)
							}
						} 
						if(strScriptFun == "MIAMIPLOT") {
						
							objMIAMIPLOT <- MIAMIPLOT(strCommand)
							MIAMIPLOT.GWADATA.valid(objMIAMIPLOT, objGWA)
							
							if(isValidScript) {
								#fnOpenPng.Miami(objGWA, "miami")					
								fnOpenPlot(objGWA, objMIAMIPLOT)
								MIAMIPLOT.run(objMIAMIPLOT, objGWA)
								fnClosePlot()
							}
						}
						if(strScriptFun == "MHPLOT") {
							objMHPLOT <- MHPLOT(strCommand)
							MHPLOT.GWADATA.valid(objMHPLOT, objGWA)
							
							if(isValidScript) {
								#fnOpenPng.MH(objGWA, "mh")		
								fnOpenPlot(objGWA, objMHPLOT)
								MHPLOT.run(objMHPLOT, objGWA)						
								fnClosePlot()
							}
						}
						if(strScriptFun == "FANCYPLOT") {
							objFANCYPLOT <- FANCYPLOT(strCommand)
							FANCYPLOT.GWADATA.valid(objFANCYPLOT, objGWA)
							
							if(isValidScript) {
								#fnOpenPng.MH(objGWA, "mh")		
								fnOpenPlot(objGWA, objFANCYPLOT)
								FANCYPLOT.run(objFANCYPLOT, objGWA)						
								fnClosePlot()
							}
						}
						if(strScriptFun == "EFFECTPLOT") {
							objEP <- EFFECTPLOT(strCommand)
							EFFECTPLOT.GWADATA.valid(objEP, objGWA)
							
							if(isValidScript) {
								fnOpenPlot(objGWA, objEP)
								objREPORT <- EFFECTPLOT.run(objEP, objGWA,objREPORT)
								fnClosePlot()
								REPORT.write(objREPORT)
							}
						}
						if(strScriptFun == "FORESTPLOT") {
							
							objFP <- FORESTPLOT(strCommand)
							objFP <- FORESTPLOT.GWADATA.valid(objFP, objGWA)
							
							if(isValidScript) {
								fnOpenPlot(objGWA, objFP)
								# stop()
								FORESTPLOT.run(objFP, objGWA)
								fnClosePlot()
							}
						}
						if(strScriptFun == "ENRICHMENT") {
							objER <- ENRICHMENT(strCommand)
							ENRICHMENT.GWADATA.valid(objER, objGWA)
							
							if(isValidScript) {
								tblEnrichmentReport <- ENRICHMENT.run(objER,objGWA,objREPORT)

								pathOutER <- objGWA@pathOut
								if(strsplit(pathOutER,"")[[1]][nchar(pathOutER)] != "/") pathOutER <- paste(pathOutER,"/",sep="")	
								fileOutER <- paste(pathOutER,objGWA@fileInShortName,".enrichment.out",sep="")
								if(!file.exists(fileOutER)) {
									sink(fileOutER)
									cat(paste(names(tblEnrichmentReport),collapse="\t"))
									cat("\n")
									cat(paste(tblEnrichmentReport[1,],collapse="\t"))
									cat("\n")
								} else {
									sink(fileOutER,append=TRUE)
									cat(paste(tblEnrichmentReport[1,],collapse="\t"))
									cat("\n")
								}
								sink()
								
								# REPORT.write(objREPORT)
							}
						}
						if(strScriptFun == "ANNOTATE") {
							objANNOT 	<- ANNOTATE(strCommand)
							ANNOTATE.GWADATA.valid(objANNOT, objGWA)		
							objGWA 	<- ANNOTATE.run(objANNOT, objGWA)							
						}
						if(strScriptFun == "GWASCATALOG") {
							objGCAT 	<- GWASCATALOG(strCommand)
							GWASCATALOG.GWADATA.valid(objGCAT, objGWA)		
							objGWA 	<- GWASCATALOG.run(objGCAT, objGWA)
						}
						if(strScriptFun == "ADDGENE") {
							objADDGENE 	<- ADDGENE(strCommand)
							ADDGENE.GWADATA.valid(objADDGENE, objGWA)		
							objGWA 	<- ADDGENE.run(objADDGENE, objGWA)							
						}
						if(strScriptFun == "ADDGMAP") {
							
							objADDGMAP 	<- ADDGMAP(strCommand)
							ADDGMAP.GWADATA.valid(objADDGMAP, objGWA)	
							if(iGWAmt==1 & iMergeTag==1) {
								tblGmap <- ADDGMAP.read(objADDGMAP, blnReadAll = isValidScript)
							}
							objGWA 	<- ADDGMAP.run(objADDGMAP, objGWA, tblGmap)							
						}
						# if(strScriptFun == "MULTINDEP") {
							# # if(isValidScript) stop()
							# objMULTINDEP 	<- MULTINDEP(strCommand)
							# MULTINDEP.GWADATA.valid(objMULTINDEP, objGWA)				
							# objINDEP.tmp <- objMULTINDEP
							# for(i in 1:length(objMULTINDEP@acolMIndep)) {
								# objINDEP.tmp@rcdCriterion <- objMULTINDEP@arcdCriterion[i]
								# objINDEP.tmp@colIndep <- objMULTINDEP@acolMIndep[i]
								# objINDEP.tmp@strTag <- objMULTINDEP@astrMIndepTag[i]
								# objINDEP.tmp@blnAddIndepInfo <- TRUE
								
								# if(iGWAmt==1 & iMergeTag==1 & i==1) {
									# if(objINDEP.tmp@blnRecombRate) {
										# tblRRmi <- INDEP.read(objINDEP.tmp, blnReadAll = isValidScript)
									# } else {
										# tblRRmi <- data.frame()
									# }
								# }
															
								# lsOut 	<- INDEP.run(objINDEP.tmp, objGWA, objREPORT, tblRRmi, isValidScript)
								# objGWA <- lsOut[[1]]
								# # objGWA.tmp with columns about individuals indeps added!
								# ## strTag.aLociTag, strTag.aTopHit, strTag.aNumLocusSNPs
							# }
							# # if(isValidScript) stop()
 # # if(isValidScript) stop()
 # # # tTest = objGWA@tblGWA[!is.na(objGWA@tblGWA$crea.aTopHit)|!is.na(objGWA@tblGWA$cys.aTopHit),c(1,2,3,4,5,9,11,16,18)]
 # # write.table(objGWA@tblGWA, "tTestfull.txt",sep="\t",row.names=F,col.names=T,quote=F)

							# lsOut 	<- MULTINDEP.run(objMULTINDEP, objGWA, objREPORT, isValidScript)
							
							# objGWA <- lsOut[[1]]
							# objGWA.multindep <- lsOut[[2]]
							# objGWA.multindep.x.traittophits <- lsOut[[3]]
							# objGWA.multindep.x.crosstraittophits <- lsOut[[4]]
							# objREPORT <- lsOut[[5]]
							
							# rm(lsOut)
							# if(isValidScript) {
								# REPORT.write(objREPORT)
								# strTagIndep <- ifelse(objMULTINDEP@strTag != "", paste(".",objMULTINDEP@strTag,sep=""), objMULTINDEP@strTag)
								
								# if(nrow(objGWA.multindep@tblGWA)>0) GWADATA.write(objGWA.multindep, strSuffix = paste(strTagIndep,".multindep",sep=""))
								# if(nrow(objGWA.multindep.x.traittophits@tblGWA)>0) GWADATA.write(objGWA.multindep.x.traittophits, strSuffix = paste(strTagIndep,".multindepXtraittophits",sep=""))
								# if(nrow(objGWA.multindep.x.crosstraittophits@tblGWA)>0) GWADATA.write(objGWA.multindep.x.crosstraittophits, strSuffix = paste(strTagIndep,".multindepXcrosstraittophits",sep=""))
								# rm(objGWA.multindep)
								# rm(objGWA.multindep.x.traittophits)
								# rm(objGWA.multindep.x.crosstraittophits)
								# rm(strTagIndep)
							# } 
						# }
						if(strScriptFun == "INDEP") {
							
							# if(isValidScript) stop()
							
							objINDEP 	<- INDEP(strCommand)
							objINDEP = INDEP.GWADATA.valid(objINDEP, objGWA)	

							if(iGWAmt==1 & iMergeTag==1) {
								if(!objINDEP@fileRecombRate == "") {
									tblRRX <- INDEP.read(objINDEP, blnReadAll = isValidScript)
								} else {
									tblRRX <- data.frame()
								}
							}
							
							# if(isValidScript) stop()
							
							lsOut 	<- INDEP.run(objINDEP, objGWA, objREPORT, tblRRX, isValidScript)
							
							objGWA		 	<- lsOut[[1]]
							objGWA.indep 	<- lsOut[[2]]
							objGWA.indep.regionLeads 	<- lsOut[[3]]
							objGWA.indep.signalLeads 	<- lsOut[[4]]
							objGWA.indep.locusLeads 	<- lsOut[[5]]
							objGWA.indep.locusMiss 	<- lsOut[[6]]
							objREPORT 		<- lsOut[[7]]
							rm(lsOut)
							if(isValidScript) {
								REPORT.write(objREPORT)
								strTagIndep <- ifelse(objINDEP@strTag != "", paste(".",objINDEP@strTag,sep=""), objINDEP@strTag)
								GWADATA.write(objGWA.indep, strSuffix = paste(strTagIndep,".indep",sep=""))
								GWADATA.write(objGWA.indep.regionLeads, strSuffix = paste(strTagIndep,".indep.regionleads",sep=""))
								if(objINDEP@fileRecombRate != "" | objINDEP@colInGPos != "") GWADATA.write(objGWA.indep.signalLeads, strSuffix = paste(strTagIndep,".indep.signalleads",sep=""))
								if(any(!objINDEP@fileClumpBed=="")) {
									GWADATA.write(objGWA.indep.locusLeads, strSuffix = paste(strTagIndep,".indep.locusleads",sep=""))
									GWADATA.write(objGWA.indep.locusMiss, strSuffix = paste(strTagIndep,".indep.locusmiss",sep=""))
								}
								rm(objGWA.indep)
								rm(objGWA.indep.regionLeads)
								rm(objGWA.indep.signalLeads)
								rm(strTagIndep)
							} 
						}
						# if(strScriptFun == "INDEPX") {
							
							# # if(isValidScript) stop()
							
							# objINDEPX 	<- INDEPX(strCommand)
							# objINDEPX = INDEPX.GWADATA.valid(objINDEPX, objGWA)	

							# if(iGWAmt==1 & iMergeTag==1) {
								# if(!objINDEPX@fileRecombRate == "") {
									# tblRRX <- INDEPX.read(objINDEPX, blnReadAll = isValidScript)
								# } else {
									# tblRRX <- data.frame()
								# }
							# }
							
							# # if(isValidScript) stop()
							
							# lsOut 	<- INDEPX.run(objINDEPX, objGWA, objREPORT, tblRRX, isValidScript)
							
							# objGWA		 	<- lsOut[[1]]
							# objGWA.indep 	<- lsOut[[2]]
							# objGWA.indep.regionLeads 	<- lsOut[[3]]
							# objGWA.indep.signalLeads 	<- lsOut[[4]]
							# objREPORT 		<- lsOut[[5]]
							# rm(lsOut)
							# if(isValidScript) {
								# REPORT.write(objREPORT)
								# strTagIndep <- ifelse(objINDEPX@strTag != "", paste(".",objINDEPX@strTag,sep=""), objINDEPX@strTag)
								# GWADATA.write(objGWA.indep, strSuffix = paste(strTagIndep,".indep",sep=""))
								# GWADATA.write(objGWA.indep.regionLeads, strSuffix = paste(strTagIndep,".indep.regionleads",sep=""))
								# if(!objINDEPX@fileRecombRate == "") GWADATA.write(objGWA.indep.signalLeads, strSuffix = paste(strTagIndep,".indep.signalleads",sep=""))
								
								# if(objINDEPX@blnWriteExtendedRegion & nrow(objGWA.indep.regionLeads@tblGWA)>0) {
									
									# regionidcolname = names(objGWA.indep.regionLeads@tblGWA)[max(grep("regionId",names(objGWA.indep.regionLeads@tblGWA)))]
									# astrRegionId = GWADATA.getcol(objGWA.indep.regionLeads, regionidcolname)
									
									# extcolname = names(objGWA.indep.regionLeads@tblGWA)[max(grep("regionExtCoordinates",names(objGWA.indep.regionLeads@tblGWA)))]
									# astrRegionCoord = GWADATA.getcol(objGWA.indep.regionLeads, extcolname)
									# aChrR = unlist(lapply(strsplit(astrRegionCoord,":",fixed=T),function(x) x[1]))
									# aPosR = unlist(lapply(strsplit(astrRegionCoord,":",fixed=T),function(x) x[2]))
									# aPosR1 = as.integer(unlist(lapply(strsplit(aPosR,"_",fixed=T),function(x) x[1])))
									# aPosR2 = as.integer(unlist(lapply(strsplit(aPosR,"_",fixed=T),function(x) x[2])))
									
									# aChrGwa = GWADATA.getcol(objGWA, objINDEPX@colInChr)
									# aPosGwa = GWADATA.getcol(objGWA, objINDEPX@colInPos)
									
									# for(i in 1:length(aChrR)) {
										# rid = astrRegionId[i]
										# isRegion = aChrGwa==aChrR[i] & aPosGwa>=aPosR1[i] & aPosGwa<=aPosR2[i]
										# objGWA.regioni = GWADATA.getrows(objGWA, which(isRegion))
										
										# strReg = paste("chr",aChrR[i],"_",round(aPosR1[i]/1000000,1),"Mb_",round(aPosR2[i]/1000000,1),"Mb",sep="")
										
										# GWADATA.write(objGWA.regioni, strSuffix = paste(strTagIndep,".indep.region_",rid,".",strReg,sep=""))
									# }
									
								# }
								
								# rm(objGWA.indep)
								# rm(objGWA.indep.regionLeads)
								# rm(objGWA.indep.signalLeads)
								# rm(strTagIndep)
							# } 
						# }
						if(strScriptFun == "INDEP2") {
# stop()
							objINDEP2 	<- INDEP2(strCommand)
							objINDEP2 = INDEP2.GWADATA.valid(objINDEP2, objGWA)	

							if(iGWAmt==1 & iMergeTag==1) {
								if(!objINDEP2@fileRecombRate == "") {
									tblRR2 <- INDEP2.read(objINDEP2, blnReadAll = isValidScript)
								} else {
									tblRR2 <- data.frame()
								}
							}
							
							# if(isValidScript) stop()
							
							lsOut 	<- INDEP2.run(objINDEP2, objGWA, objREPORT, tblRR2, isValidScript)
							
							objGWA		 	<- lsOut[[1]]
							objGWA.indep 	<- lsOut[[2]]
							objGWA.indep.regionLeads 	<- lsOut[[3]]
							objGWA.indep.signalLeads 	<- lsOut[[4]]
							objREPORT 		<- lsOut[[5]]
							rm(lsOut)
							if(isValidScript) {
								REPORT.write(objREPORT)
								strTagIndep <- ifelse(objINDEP2@strTag != "", paste(".",objINDEP2@strTag,sep=""), objINDEP2@strTag)
								GWADATA.write(objGWA.indep, strSuffix = paste(strTagIndep,".indep",sep=""))
								GWADATA.write(objGWA.indep.regionLeads, strSuffix = paste(strTagIndep,".indep.regionleads",sep=""))
								if(!objINDEP2@fileRecombRate == "") GWADATA.write(objGWA.indep.signalLeads, strSuffix = paste(strTagIndep,".indep.signalleads",sep=""))
								
								rm(objGWA.indep)
								rm(objGWA.indep.regionLeads)
								rm(objGWA.indep.signalLeads)
								rm(strTagIndep)
							} 
						}
						# if(strScriptFun == "CLUMP") {
							# objCLUMP 	<- CLUMP(strCommand)
							# CLUMP.GWADATA.valid(objCLUMP, objGWA)							
							# #if(isValidScript) stop()
							# if(isValidScript) {
							
								# lsOut 	<- CLUMP.run(objCLUMP, objGWA, objREPORT)
								# objGWA		 			<- lsOut[[1]]
								# objGWA.clump 			<- lsOut[[2]]
								# objGWA.clump.x 			<- lsOut[[3]]
								# objGWA.clump.nomatch 	<- lsOut[[4]]
								# objREPORT 				<- lsOut[[5]]
								# rm(lsOut)
								# REPORT.write(objREPORT)
								# strTagClump <- ifelse(objCLUMP@strTag != "", paste(".",objCLUMP@strTag,sep=""), objCLUMP@strTag)
								# GWADATA.write(objGWA.clump, strSuffix = paste(strTagClump,".clump",sep=""))
								# GWADATA.write(objGWA.clump.x, strSuffix = paste(strTagClump,".clumpX",sep=""))
								# if(nrow(objGWA.clump.nomatch@tblGWA)>0) {
									# GWADATA.write(objGWA.clump.nomatch, strSuffix = paste(strTagClump,".clump_nomatch",sep=""))
								# }
								# rm(objGWA.clump)
								# rm(objGWA.clump.x)
								# rm(objGWA.clump.nomatch)
								# rm(strTagClump)
							# }
						# }
						# if(strScriptFun == "MULTCLUMP") {

							# objMULTCLUMP 	<- MULTCLUMP(strCommand)
							# MULTCLUMP.GWADATA.valid(objMULTCLUMP, objGWA)							

							# lsOut 	<- MULTCLUMP.run(objMULTCLUMP, objGWA, objREPORT,isValidScript)
							# objGWA		 		<- lsOut[[1]]
							# objGWA.multclump 	<- lsOut[[2]]
							# objGWA.multclump.x 	<- lsOut[[3]]
							# objGWA.multclump.nomatch 	<- lsOut[[4]]
							# objREPORT 			<- lsOut[[5]]
							
							# rm(lsOut)
							# if(isValidScript) {
								# REPORT.write(objREPORT)
								# strTagClump <- ifelse(objMULTCLUMP@strTag != "", paste(".",objMULTCLUMP@strTag,sep=""), objMULTCLUMP@strTag)
								# GWADATA.write(objGWA.multclump, strSuffix = paste(strTagClump,".multclump",sep=""))
								# GWADATA.write(objGWA.multclump.x, strSuffix = paste(strTagClump,".multclumpX",sep=""))
								# if(nrow(objGWA.multclump.nomatch@tblGWA)>0) {
									# GWADATA.write(objGWA.multclump.nomatch, strSuffix = paste(strTagClump,".multclump_nomatch",sep=""))
								# }
								# rm(objGWA.multclump)
								# rm(objGWA.multclump.x)
								# rm(objGWA.multclump.nomatch)
								# rm(strTagClump)
							# } 
						# }
						if(strScriptFun == "BONFERRONI") {
							objBF 	<- BONFERRONI(strCommand)
							lsOut 	<- BONFERRONI.run(objBF, objGWA, objREPORT,isValidScript)
							objGWA 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							rm(lsOut)
							if(isValidScript) {
								REPORT.write(objREPORT)
							}
						}
						if(strScriptFun == "FDR") {
							objFDR 	<- FDR(strCommand)				
							objGWA 	<- FDR.run(objFDR, objGWA)
						}
						if(strScriptFun == "PCA2STEP") {
							
							objPCA2STEP 	<- PCA2STEP(strCommand)
							objPCA2STEP 	<- PCA2STEP.GWADATA.valid(objPCA2STEP, objGWA)
							
							#if(isValidScript) stop()
							
							lsOut 		<- PCA2STEP.run(objPCA2STEP, objGWA, objREPORT, isValidScript)
							objGWA 		<- lsOut[[1]]
							objGWA.signif <- lsOut[[2]]
							objGWA.missbed <- lsOut[[3]]
							objREPORT 	<- lsOut[[4]]
							rm(lsOut)
							if(isValidScript) {
								if(nrow(objGWA.signif@tblGWA)>0) {
									GWADATA.write(objGWA.signif, strSuffix = paste0(".",objPCA2STEP@strTag,".signif"))
								}
								if(nrow(objGWA.missbed@tblGWA)>0) {
									GWADATA.write(objGWA.missbed, strSuffix = paste0(".",objPCA2STEP@strTag,".missing"))
								}
								REPORT.write(objREPORT)
							}
						}
						# if(strScriptFun == "CLUMPX") {
							# objCLUMPX 	<- CLUMPX(strCommand)
							# objCLUMPX 	<- CLUMPX.GWADATA.valid(objCLUMPX, objGWA)
							# # if(isValidScript) stop()
							# lsOut 		<- CLUMPX.run(objCLUMPX, objGWA, objREPORT,isValidScript)
							
							# objGWA		 	<- lsOut[[1]]
							# objGWA.clump 	<- lsOut[[2]]
							# objGWA.clump.locusLeads	<- lsOut[[3]]
							# objREPORT 		<- lsOut[[4]]
							# rm(lsOut)
							# if(isValidScript) {
								# REPORT.write(objREPORT)
								# strTagClump <- ifelse(objCLUMPX@strTag != "", paste(".",objCLUMPX@strTag,sep=""), objCLUMPX@strTag)
								# GWADATA.write(objGWA.clump, strSuffix = paste(strTagClump,".clump",sep=""))
								# GWADATA.write(objGWA.clump.locusLeads, strSuffix = paste(strTagClump,".clump.locusleads",sep=""))
							# }
						# }
						if(strScriptFun == "WEIGHTEDHYPOTHESIS") {
							objWH 	<- WEIGHTEDHYPOTHESIS(strCommand)
							WEIGHTEDHYPOTHESIS.GWADATA.valid(objWH,objGWA)
							
							#if(isValidScript) stop()
							
							lsOut 	<- WEIGHTEDHYPOTHESIS.run(objWH,objGWA,objREPORT)
						
							objGWA		 	<- lsOut[[1]]
							objGWA.wh.signif 	<- lsOut[[2]]
							objREPORT 		<- lsOut[[3]]
							rm(lsOut)
							if(isValidScript) {
								REPORT.write(objREPORT)
								GWADATA.write(objGWA.wh.signif, strSuffix = paste(objWH@strTag,".signif",sep=""))
								rm(objGWA.wh.signif)
							}
						}
						if(strScriptFun == "JOINTTEST") {
							objJT 	<- JOINTTEST(strCommand)				
							objGWA 	<- JOINTTEST.run(objJT, objGWA)
						}
						if(strScriptFun == "CALCPDIFF") {
							objPD 	<- CALCPDIFF(strCommand)				
							lsOut 	<- CALCPDIFF.run(objPD, objGWA, objREPORT)
							objGWA 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							rm(lsOut)
							gc()
							if(isValidScript) REPORT.write(objREPORT)
						}
						if(strScriptFun == "CALCPDIFFDIFF") {
							objPD 	<- CALCPDIFFDIFF(strCommand)				
							lsOut 	<- CALCPDIFFDIFF.run(objPD, objGWA, objREPORT)
							objGWA 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							rm(lsOut)
							gc()
							if(isValidScript) REPORT.write(objREPORT)
						}
						if(strScriptFun == "CALCPHET") {
							objCP 	<- CALCPHET(strCommand)				
							objGWA 	<- CALCPHET.run(objCP, objGWA)
						}
						if(strScriptFun == "METAANALYSIS") {
							objMA 	<- METAANALYSIS(strCommand)				
							objGWA 	<- METAANALYSIS.run(objMA, objGWA)
						}
						if(strScriptFun == "RENAMEMARKER") {
							
							if(EqcScriptCodeType=="GWA") {
								idxRENAMEGwa <- idxRENAMEGwa + 1
								idxRENAMEMerged <- idxRENAMEGwa
								idxRENAME = idxRENAMEGwa
							} else {
								## MERGED
								idxRENAMEMerged <- idxRENAMEMerged + 1
								idxRENAME = idxRENAMEMerged
							}
							if(iGWAmt==1 & iMergeTag==1) {
								## init new merge object; two times at each command group
								
								objRENAMEMARKER <- RENAMEMARKER(strCommand)
								objRENAMEMARKER <- RENAMEMARKER.read(objRENAMEMARKER, blnReadAll = isValidScript)
								
								container.Rename[[idxRENAME]] <- objRENAMEMARKER
							}
							objRENAMEMARKER <- container.Rename[[idxRENAME]]
							RENAMEMARKER.valid(objRENAMEMARKER)
							RENAMEMARKER.GWADATA.valid(objRENAMEMARKER, objGWA)
							## 
											
							lsOut <- RENAMEMARKER.run(objRENAMEMARKER, objGWA, objREPORT, isValidScript)

							objGWA 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							rm(lsOut)
							rm(objRENAMEMARKER)
						}
						if(strScriptFun == "CREATECPAID") {
							
							objCPA <- CREATECPAID(strCommand)
							CREATECPAID.GWADATA.valid(objCPA, objGWA)

							objGWA <- CREATECPAID.run(objCPA, objGWA, isValidScript)
							rm(objCPA)
						}
						# if(strScriptFun == "LIFTOVER") {
							
							# objLO <- LIFTOVER(strCommand)
							# objLO <- LIFTOVER.GWADATA.valid(objLO, objGWA)
# if(isValidScript) stop()
							# if() lsOut <- LIFTOVER.run(objLO, objGWA, isValidScript, objREPORT)
							# objGWA 		<- lsOut[[1]]
							# objREPORT 	<- lsOut[[2]]
							# tlift 		<- lsOut[[3]]
							# rm(objLO)
							# rm(lsOut)
						# }
						if(strScriptFun == "LIFTOVER") {
							
							objLO <- LIFTOVER(strCommand, objGWADATA.default)
							LIFTOVER.GWADATA.valid(objLO, objGWA)
							
							if(EqcScriptCodeType=="GWA") {
								idxLOGwa <- idxLOGwa + 1
								idxLOMerged <- idxLOGwa
								idxLO = idxLOGwa
							} else {
								## MERGED
								idxLOMerged <- idxLOMerged + 1
								idxLO = idxLOMerged
							}
							
							lsOut <- LIFTOVER.check(objLO, isValidScript, objGWA, objREPORT)
							objLO 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							blnLiftBackup = objLO@blnLift
							
							
							## init new merge object; two times at each command group
							if(iGWAmt==1 & iMergeTag==1) {
							
								lsOut <- LIFTOVER.read(objLO, isValidScript, objGWA, objREPORT)
								objLO 		<- lsOut[[1]]
								objREPORT 	<- lsOut[[2]]
								rm(lsOut)
								container.LiftOver[[idxLO]] <- objLO
							
							}
							objLO <- container.LiftOver[[idxLO]] 
							objLO@blnLift = blnLiftBackup
							
							lsOut <- LIFTOVER.run(objLO, objGWA, objREPORT, isValidScript)
							objGWA 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							rm(lsOut)
							rm(objLO)
							if(isValidScript) REPORT.write(objREPORT)
						}
						if(strScriptFun == "CREATECPTID") {
							
							objCPT <- CREATECPTID(strCommand, objGWADATA.default)
							CREATECPTID.GWADATA.valid(objCPT, objGWA)
							if(objCPT@blnUseMap) {
							
								if(EqcScriptCodeType=="GWA") {
									idxCPTIDGwa <- idxCPTIDGwa + 1
									idxCPTIDMerged <- idxCPTIDGwa
									idxCPTID = idxCPTIDGwa
								} else {
									## MERGED
									idxCPTIDMerged <- idxCPTIDMerged + 1
									idxCPTID = idxCPTIDMerged
								}
								if(iGWAmt==1 & iMergeTag==1) {
									## init new merge object; two times at each command group
									
									objCPT <- CREATECPTID.read(objCPT, blnReadAll = isValidScript, blnUseFastRead = objGWA@blnUseFastRead)
									container.CptidMap[[idxCPTID]] <- objCPT
									
								}
								objCPT <- container.CptidMap[[idxCPTID]]
							}
							
							lsOut <- CREATECPTID.run(objCPT, objGWA, objREPORT, isValidScript)
							objGWA 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							rm(lsOut)
							rm(objCPT)
							if(isValidScript) REPORT.write(objREPORT)
						}
						if(strScriptFun == "GETVARTYPE") {
							
							objGVT <- GETVARTYPE(strCommand)
							GETVARTYPE.GWADATA.valid(objGVT, objGWA)
							
							lsOut <- GETVARTYPE.run(objGVT, objGWA, objREPORT, isValidScript)
							objGWA 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							rm(lsOut)
							rm(objGVT)
							if(isValidScript) REPORT.write(objREPORT)
						}
						if(strScriptFun == "GETVAR") {
							# if(isValidScript) stop()
							objGV 	<- GETVAR(strCommand)
							
							objGV = GETVAR.GWADATA.valid(objGV, objGWA)
							
							# if(isValidScript & iGWAmt==2) stop()
							lsOut 	<- GETVAR.run(objGV, objGWA, objREPORT)
							
							objGWA.var <- lsOut[[1]]
							objREPORT <- lsOut[[2]]
							rm(lsOut)
							if(isValidScript) {
								
								
								REPORT.write(objREPORT)
								
								## write out for specific file: 
								strMarkerTag = gsub("\"","",objGV@strMarkerTag)
								strMarkerTag = gsub(":",".",strMarkerTag)
								if(nrow(objGWA.var@tblGWA)>0) GWADATA.write(objGWA.var, strSuffix = paste(".",strMarkerTag,sep=""),blnWrite10rows=FALSE)
								
								## write out combined file if requested 
								if(objGV@blnTry2Combine) {
									fileOutVar = paste(fileOutBody,".",strMarkerTag,".txt",sep="")
									## clean up if combined file already exists at first file
									if(iGWA == 1 & file.exists(fileOutVar)) file.remove(fileOutVar)	
									
									tblVarAdd = objGWA.var@tblGWA
									if(nrow(tblVarAdd) == 0) {
										anamestmp = names(tblVarAdd)
										tblVarAdd = as.data.frame(t(rep(NA,ncol(tblVarAdd))))
										names(tblVarAdd) = anamestmp
										rm(anamestmp)
									} 
									tblVarAdd = cbind(file=rep(objGWA.var@fileInShortName,nrow(tblVarAdd)),tblVarAdd,stringsAsFactors=FALSE)
								
									if(file.exists(fileOutVar)) {
										tblVarIn = read.table(fileOutVar, stringsAsFactors=FALSE, header=TRUE, sep="\t")
										if(identical(names(tblVarIn),names(tblVarAdd))) {
											write.table(tblVarAdd, fileOutVar, row.names=F, quote=F, sep="\t", col.names=F, append=TRUE)
										} else {
											cat(paste("EASY WARNING:",strScriptFun,"\n Appending of \n",objGWA.var@fileInShortName,"\n to the combined variant file not possible !!!\n", sep=""))
										}
										rm(tblVarIn)
									} else {
										write.table(tblVarAdd, fileOutVar, row.names=F, quote=F, sep="\t", col.names=T)
									}
									rm(tblVarAdd,fileOutVar)
								}
								rm(objGWA.var)
							}
						}
						if(strScriptFun == "LOADINFO") {
							
							if(EqcScriptCodeType=="GWA") {
								idxLOADINFOGwa <- idxLOADINFOGwa + 1
								idxLOADINFOMerged <- idxLOADINFOGwa
								idxLOADINFO = idxLOADINFOGwa
							} else {
								## MERGED
								idxLOADINFOMerged <- idxLOADINFOMerged + 1
								idxLOADINFO = idxLOADINFOMerged
							}
							if(iGWAmt==1 & iMergeTag==1) {
								## init new loadinfo object; two times at each command group
								objLI <- LOADINFO(strCommand, objGWADATA.default)
								LOADINFO.GWADATA.valid(objLI, objGWA)
								objLI <- LOADINFO.read(objLI, blnReadAll = isValidScript, blnUseFastRead = objGWA@blnUseFastRead)
								container.LoadInfo[[idxLOADINFO]] <- objLI
							}
							objLI <- container.LoadInfo[[idxLOADINFO]]
							lsOut <- LOADINFO.run(objLI, objGWA, objREPORT, isValidScript)
							objGWA 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							rm(lsOut)
							rm(objLI)
							if(isValidScript) REPORT.write(objREPORT)
						}
						if(strScriptFun == "HARMONIZEALLELES") {
							
							objHA <- HARMONIZEALLELES(strCommand)
							HARMONIZEALLELES.GWADATA.valid(objHA, objGWA)
							
							lsOut <- HARMONIZEALLELES.run(objHA, objGWA, objREPORT, isValidScript)

							objGWA 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							rm(lsOut)
							rm(objHA)
							if(isValidScript) REPORT.write(objREPORT)
						}
						if(strScriptFun == "MERGE") {
							if(EqcScriptCodeType=="GWA") {
								idxMERGEGwa <- idxMERGEGwa + 1
								idxMERGEMerged <- idxMERGEGwa
								idxMERGE = idxMERGEGwa
							} else {
								## MERGED
								idxMERGEMerged <- idxMERGEMerged + 1
								idxMERGE = idxMERGEMerged
							}
							if(iGWAmt==1 & iMergeTag==1) {
								## init new merge object; two times at each command group
								
								objMERGE <- MERGE(strCommand, objGWADATA.default, fileMerge)
								
								if(isValidScript) {
									objMERGE <- GWADATA.read(objMERGE)
								} else {
									objMERGE <- GWADATA.read.10rows(objMERGE)
								}
								
								container.Merge[[idxMERGE]] <- objMERGE
							}
							objMERGE <- container.Merge[[idxMERGE]]
							if(objMERGE@fileRefTag == "1" | objMERGE@fileRefTag == objGWA@fileInTag) {
								MERGE.GWADATA.valid(objMERGE, objGWA)
								lsOut <- MERGE.run(objMERGE, objGWA, objREPORT, isValidScript)
								objGWA 		<- lsOut[[1]]
								objREPORT 	<- lsOut[[2]]
								rm(lsOut)
								if(isValidScript) REPORT.write(objREPORT)
							}
							rm(objMERGE)
							
							
							# iMerge = iMerge + 1
							# if(iMerge>length(container.Merge)) {
							# # if(iGWA == 1) {
								# objMERGE <- MERGE(strCommand, objGWADATA.default)
								
								# if(isValidScript) objMERGE <- GWADATA.read(objMERGE)
								# else objMERGE <- GWADATA.read.10rows(objMERGE)
								
								# container.Merge[[iMerge]] <- objMERGE
							# }
							# objMERGE <- container.Merge[[iMerge]]
							# if(objMERGE@fileRefTag == "1" | objMERGE@fileRefTag == objGWA@fileInTag) {
								# MERGE.GWADATA.valid(objMERGE, objGWA)
								# lsOut <- MERGE.run(objMERGE, objGWA, objREPORT, isValidScript)
								# objGWA 		<- lsOut[[1]]
								# objREPORT 	<- lsOut[[2]]
								# rm(lsOut)
								# if(isValidScript) REPORT.write(objREPORT)
							# }
							# rm(objMERGE)
						}
						if(strScriptFun == "ADJUSTALLELESSIMPLE") {
							objAAS <- ADJUSTALLELESSIMPLE(strCommand)
							
							ADJUSTALLELESSIMPLE.GWADATA.valid(objAAS, objGWA)
							lsOut <- ADJUSTALLELESSIMPLE.run(objAAS, objGWA, objREPORT, isValidScript)
							
							objGWA 			<- lsOut[[1]]
							objREPORT 		<- lsOut[[2]]
							if(isValidScript) REPORT.write(objREPORT)
							
							rm(lsOut)
							rm(objAAS)
						}
						if(strScriptFun == "ADJUSTALLELES") {
							objAA <- ADJUSTALLELES(strCommand)
							
							ADJUSTALLELES.GWADATA.valid(objAA, objGWA)
							lsOut <- ADJUSTALLELES.run(objAA, objGWA, objREPORT, isValidScript)
							
							objGWA 			<- lsOut[[1]]
							objREPORT 		<- lsOut[[2]]
							if(isValidScript) REPORT.write(objREPORT)
							
							# objGWA.miss 	<- lsOut[[3]]
							# objGWA.invalid 	<- lsOut[[4]]
							rm(lsOut)
							rm(objAA)
						}
						if(strScriptFun == "AFCHECK") {
							objAC <- AFCHECK(strCommand)
							## objAC inherits from classes SPLOT and from ADJUSTALLELES (+ MERGE)
							SPLOT.GWADATA.valid(objAC, objGWA)
							
							if(isValidScript) {
								if(objAC@strMode == "singleplot") {
									fnOpenPlot(objGWA, objAC)
									SPLOT.run(objAC, objGWA)
									fnClosePlot()
								}
								if(objAC@strMode == "subplot") {
									
									if(EqcScriptCodeType=="GWA") {
										idxMultiplotGwaMt <- idxMultiplotGwaMt + 1
										idxMultiplotMerged <- idxMultiplotGwaMt
										idxAdd = idxMultiplotGwaMt
									} else {
										## MERGED
										idxMultiplotMerged <- idxMultiplotMerged + 1
										idxAdd = idxMultiplotMerged
									}

									# if(iGWAmt==1 & iMergeTag==1) {
									if(idxAdd > nSubplots) {
										## happens at each commandgroup ! 
										## init new multiplot
										nPlot = ifelse(EqcScriptCodeType=="GWA",nGWA,nMergeTags)
										fileOutMultiPlot = fnOpenMultiPlot(fileOutBody, objAC, nPlot, afileOutMultiPlot, objGWA@blnOverwriteResults)
										nSubplots <- idxAdd
										afileOutMultiPlot = c(afileOutMultiPlot,fileOutMultiPlot)
									} else {
										# add to existing plot; 
										fnAddPlot(idxAdd + idxDeviceOffset) # set device
									}
									## plot data
									SPLOT.run(objAC, objGWA)
									if((EqcScriptCodeType=="GWA" & iGWAmt == nGWAmt & iMergeTag == nMergeTags) | (EqcScriptCodeType=="MERGED" & iMergeTag == nMergeTags)) {
										## close multi plot
										fnClosePlot()
										idxDeviceOffset = idxDeviceOffset - 1
									}
								}
							}
							lsOut <- AFCHECK.run(objAC, objGWA, objREPORT)
							#objGWA.outlier <- lsOut[[1]]
							objGWA 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]] 
							rm(lsOut)
							# if(objAC@blnWriteOutlier & dim(objGWA.outlier@tblGWA)[1] > 0) {
								# GWADATA.write(objGWA.outlier, strSuffix = paste(".",objAC@strTag,".outlier",sep=""))
								# rm(objGWA.outlier)
							# }
							if(isValidScript) REPORT.write(objREPORT)
							#rm(objGWA.tmp)
							rm(objAC)
						
						}
						if(strScriptFun == "WRITE") {
							objWRITE 	<- WRITE(strCommand)			
							if(isValidScript) {
								objREPORT <- WRITE.run(objWRITE, objGWA, objREPORT)
								REPORT.write(objREPORT)
							}
						}
						if(strScriptFun == "ADDCOL") {
							objADDCOL 	<- ADDCOL(strCommand)				
							objGWA 		<- ADDCOL.run(objADDCOL, objGWA)
						}
						if(strScriptFun == "EDITCOL") {
							objEDITCOL 	<- EDITCOL(strCommand)				
							EDITCOL.GWADATA.valid(objEDITCOL, objGWA)
							objGWA 		<- EDITCOL.run(objEDITCOL, objGWA)
						}
						if(strScriptFun == "STRSPLITCOL") {
							objSTRSPLITCOL 	<- STRSPLITCOL(strCommand)				
							STRSPLITCOL.GWADATA.valid(objSTRSPLITCOL, objGWA)
							objGWA 		<- STRSPLITCOL.run(objSTRSPLITCOL, objGWA)
						}
						if(strScriptFun == "EXTRACTSNPS") {
							objES 	<- EXTRACTSNPS(strCommand)
							EXTRACTSNPS.GWADATA.valid(objES, objGWA)
							lsOut 	<- EXTRACTSNPS.run(objES, objGWA, objREPORT)
							tblOut 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							rm(lsOut)
							if(isValidScript) {				
								REPORT.write(objREPORT)
								strTagExtract <- ifelse(objES@strTag != "", paste(".",objES@strTag,sep=""), objES@strTag)
								if(strsplit(objGWA@pathOut,"")[[1]][nchar(objGWA@pathOut)] != "/") pathOutExtract <- paste(objGWA@pathOut,"/",sep="")
								else pathOutExtract<-objGWA@pathOut
								fileExtractOut = paste(pathOutExtract,objGWA@fileInShortName,strTagExtract,".extracted",sep="")
								write.table(tblOut, fileExtractOut, row.names=F, quote=F, sep="\t", na="NA")
								rm(tblOut)
								rm(strTagExtract)
								rm(pathOutExtract)
								rm(fileExtractOut)
							}
						}
						if(strScriptFun == "REMOVESNPS") {
							objRS 	<- REMOVESNPS(strCommand)
							REMOVESNPS.GWADATA.valid(objRS, objGWA)
							lsOut 	<- REMOVESNPS.run(objRS, objGWA, objREPORT)
							objGWA 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							objGWA.Removed 	<- lsOut[[3]]
							rm(lsOut)
							if(isValidScript) {
								REPORT.write(objREPORT)
								if(nrow(objGWA.Removed@tblGWA)>0) 
									GWADATA.write(objGWA.Removed, strSuffix = ".Removed")
								rm(objGWA.Removed)
							}
						}
						if(strScriptFun == "FLIPSTRAND") {
							objFS 	<- FLIPSTRAND(strCommand)
							FLIPSTRAND.GWADATA.valid(objFS, objGWA)
							lsOut 	<- FLIPSTRAND.run(objFS, objGWA, objREPORT)
							objGWA 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							rm(lsOut)
							if(isValidScript) REPORT.write(objREPORT)
						}
						if(strScriptFun == "EVALSTAT") {
							objEVALSTAT <- EVALSTAT(strCommand)
							objREPORT 	<- EVALSTAT.run(objEVALSTAT, objGWA, objREPORT)
							if(isValidScript) REPORT.write(objREPORT)
						}
						if(strScriptFun == "GETNUM") {
							objGN 		<- GETNUM(strCommand)
							objREPORT 	<- GETNUM.run(objGN, objGWA, objREPORT)
							if(isValidScript) REPORT.write(objREPORT)
						}
						if(strScriptFun == "CALCULATE") {
							objCALC 	<- CALCULATE(strCommand)				
							objREPORT 	<- CALCULATE.run(objCALC, objGWA, objREPORT)
							if(isValidScript) REPORT.write(objREPORT)
						}
						if(strScriptFun == "RADDCOL") {
							objRCALC 	<- RADDCOL(strCommand)				
							objREPORT 	<- RADDCOL.run(objRCALC, objREPORT)
							if(isValidScript) REPORT.write(objREPORT)
						}
						if(strScriptFun == "FILTER") {
							objFILTER 	<- FILTER(strCommand)
							lsOut 		<- FILTER.run(objFILTER, objGWA, objREPORT, !isValidScript)
							objGWA 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							rm(lsOut)
							if(isValidScript) REPORT.write(objREPORT)
						}
						if(strScriptFun == "CLEAN") {
							objCLEAN 	<- CLEAN(strCommand)
							lsOut 		<- CLEAN.run(objCLEAN, objGWA, objREPORT, !isValidScript)
							objGWA 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							objGWA.Cleaned 	<- lsOut[[3]]
							rm(lsOut)
							if(isValidScript) {
								REPORT.write(objREPORT)
								if(objCLEAN@blnWriteCleaned&nrow(objGWA.Cleaned@tblGWA)>0) 
									GWADATA.write(objGWA.Cleaned, 
													strSuffix = paste(".",objCLEAN@strCleanName,sep=""), 
													strMode = ifelse(nrow(objGWA.Cleaned@tblGWA)<10000, "txt","gz")
													)
								rm(objGWA.Cleaned)
							}
						}
						if(strScriptFun == "CRITERION") {
							objCRIT 	<- CRITERION(strCommand)
							lsOut 	<- CRITERION.run(objCRIT, objGWA, objREPORT)
							
							objGWA.Crit <- lsOut[[1]]
							objREPORT <- lsOut[[2]]
							rm(lsOut)
							if(isValidScript) {
								REPORT.write(objREPORT)
								#GWADATA.write(objGWA.Crit, strSuffix = ".criterion")
								if(objCRIT@blnWriteEmpty | nrow(objGWA.Crit@tblGWA)>0) 
									GWADATA.write(objGWA.Crit, strSuffix = paste(".",sub("numSNP_","",objCRIT@strCritName),sep=""))
								rm(objGWA.Crit)
							}
						}
						if(strScriptFun == "GETCOLS") {
							objGC 	<- GETCOLS(strCommand)
							objGWA 	<- GETCOLS.run(objGC, objGWA)
						}
						if(strScriptFun == "SORT") {
							objSORT 	<- SORT(strCommand)
							SORT.GWADATA.valid(objSORT, objGWA)
							objGWA 	<- SORT.run(objSORT, objGWA)
						}
						if(strScriptFun == "RENAMECOL") {
							objRC 		<- RENAMECOL(strCommand)
							lsOut 		<- RENAMECOL.run(objRC, objGWA, objREPORT)
							objGWA 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							rm(lsOut)
							if(isValidScript) REPORT.write(objREPORT)
						}
						if(strScriptFun == "REMOVECOL") {
							objRC 	<- REMOVECOL(strCommand)
							objGWA 	<- REMOVECOL.run(objRC, objGWA)
						}
						if(strScriptFun == "CLEANDUPLICATES") {
							objCD 		<- CLEANDUPLICATES(strCommand)
							CLEANDUPLICATES.GWADATA.valid(objCD, objGWA)
							if(isValidScript){
								lsOut 		<- CLEANDUPLICATES.run(objCD, objGWA, objREPORT)
								objGWA 		<- lsOut[[1]]
								objREPORT 	<- lsOut[[2]]
								rm(lsOut)
								REPORT.write(objREPORT)
							}
						}
						if(strScriptFun == "GC") {
							objGC 	<- GC(strCommand)
							GC.GWADATA.valid(objGC, objGWA)
							lsOut 	<- GC.run(objGC, objGWA, objREPORT)
							objGWA 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							rm(lsOut)
							if(isValidScript) REPORT.write(objREPORT)
						}
						if(strScriptFun == "QQPLOT") {
							
							objQQPLOT <- QQPLOT(strCommand)
							QQPLOT.GWADATA.valid(objQQPLOT, objGWA)
							
							if(isValidScript) {
								if(objQQPLOT@strMode == "singleplot") {
									fnOpenPlot(objGWA, objQQPLOT)
									objREPORT <- QQPLOT.run(objQQPLOT, objGWA, objREPORT)
									fnClosePlot()
								} 
								if(objQQPLOT@strMode == "subplot") {
									
									if(EqcScriptCodeType=="GWA") {
										idxMultiplotGwaMt <- idxMultiplotGwaMt + 1
										idxMultiplotMerged <- idxMultiplotGwaMt
										idxAdd = idxMultiplotGwaMt
									} else {
										## MERGED
										idxMultiplotMerged <- idxMultiplotMerged + 1
										idxAdd = idxMultiplotMerged
									}
									
									# if(iGWAmt==1 & iMergeTag==1) {
									if(idxAdd > nSubplots) {
										## happens at each commandgroup ! 
										## init new multiplot
										nPlot = ifelse(EqcScriptCodeType=="GWA",nGWA,nMergeTags)
										fileOutMultiPlot = fnOpenMultiPlot(fileOutBody, objQQPLOT, nPlot, afileOutMultiPlot, objGWA@blnOverwriteResults)
										nSubplots <- idxAdd
										afileOutMultiPlot = c(afileOutMultiPlot,fileOutMultiPlot)
									} else {
										# add to existing plot; 
										fnAddPlot(idxAdd + idxDeviceOffset) # set device
									}
									## plot data
									objREPORT <- QQPLOT.run(objQQPLOT, objGWA, objREPORT)
									if((EqcScriptCodeType=="GWA" & iGWAmt == nGWAmt & iMergeTag == nMergeTags) | (EqcScriptCodeType=="MERGED" & iMergeTag == nMergeTags)) {
										## close multi plot
										fnClosePlot()
										idxDeviceOffset = idxDeviceOffset - 1
									}
								}
							}
						}
						if(strScriptFun == "SPLOT") {
							
							objSPLOT <- SPLOT(strCommand)
							SPLOT.GWADATA.valid(objSPLOT, objGWA)
							
							if(isValidScript) {
								if(objSPLOT@strMode == "singleplot") {
									fnOpenPlot(objGWA, objSPLOT)
									SPLOT.run(objSPLOT, objGWA)
									fnClosePlot()
								}
								if(objSPLOT@strMode == "subplot") {
								
									if(EqcScriptCodeType=="GWA") {
										idxMultiplotGwaMt <- idxMultiplotGwaMt + 1
										idxMultiplotMerged <- idxMultiplotGwaMt
										idxAdd = idxMultiplotGwaMt
									} else {
										## MERGED
										idxMultiplotMerged <- idxMultiplotMerged + 1
										idxAdd = idxMultiplotMerged
									}

									# if(iGWAmt==1 & iMergeTag==1) {
									if(idxAdd > nSubplots) {
										## happens at each commandgroup ! 
										## init new multiplot
										nPlot = ifelse(EqcScriptCodeType=="GWA",nGWA,nMergeTags)
										fileOutMultiPlot = fnOpenMultiPlot(fileOutBody, objSPLOT, nPlot, afileOutMultiPlot, objGWA@blnOverwriteResults)
										nSubplots <- idxAdd
										afileOutMultiPlot = c(afileOutMultiPlot,fileOutMultiPlot)
									} else {
										# add to existing plot; 
										fnAddPlot(idxAdd + idxDeviceOffset) # set device
									}
									## plot data
									SPLOT.run(objSPLOT, objGWA)
									if((EqcScriptCodeType=="GWA" & iGWAmt == nGWAmt & iMergeTag == nMergeTags) | (EqcScriptCodeType=="MERGED" & iMergeTag == nMergeTags)) {
										## close multi plot
										fnClosePlot()
										idxDeviceOffset = idxDeviceOffset - 1
									}
								}
							}
						}
						if(strScriptFun == "PZPLOT") {
							objPZPLOT <- PZPLOT(strCommand)
							PZPLOT.GWADATA.valid(objPZPLOT, objGWA)
							
							if(isValidScript) {
								if(objPZPLOT@strMode == "singleplot") {
									#fnOpenPng(objGWA, "pz")
									#fnOpenPlot(objGWA, objPZPLOT, strSuffix = "pz")
									fnOpenPlot(objGWA, objPZPLOT)
									SPLOT.run(objPZPLOT, objGWA)
									fnClosePlot()
								} 
								if(objPZPLOT@strMode == "subplot") {
									
									if(EqcScriptCodeType=="GWA") {
										idxMultiplotGwaMt <- idxMultiplotGwaMt + 1
										idxMultiplotMerged <- idxMultiplotGwaMt
										idxAdd = idxMultiplotGwaMt
									} else {
										## MERGED
										idxMultiplotMerged <- idxMultiplotMerged + 1
										idxAdd = idxMultiplotMerged
									}

									# if(iGWAmt==1 & iMergeTag==1) {
									if(idxAdd > nSubplots) {
										## happens at each commandgroup ! 
										## init new multiplot
										nPlot = ifelse(EqcScriptCodeType=="GWA",nGWA,nMergeTags)
										fileOutMultiPlot = fnOpenMultiPlot(fileOutBody, objPZPLOT, nPlot, afileOutMultiPlot, objGWA@blnOverwriteResults)
										nSubplots <- idxAdd
										afileOutMultiPlot = c(afileOutMultiPlot,fileOutMultiPlot)
									} else {
										# add to existing plot; 
										fnAddPlot(idxAdd + idxDeviceOffset) # set device
									}
									## plot data
									SPLOT.run(objPZPLOT, objGWA)
									if((EqcScriptCodeType=="GWA" & iGWAmt == nGWAmt & iMergeTag == nMergeTags) | (EqcScriptCodeType=="MERGED" & iMergeTag == nMergeTags)) {
										## close multi plot
										fnClosePlot()
										idxDeviceOffset = idxDeviceOffset - 1
									}
								}
							}
						}
						if(strScriptFun == "BRPLOT") {
								
							#if(iGWA == nGWAmt) {
							objBRPLOT <- BRPLOT(strCommand)
							BRPLOT.REPORT.valid(objBRPLOT, objREPORT)
							
							if(isValidScript) {
								idxBReportPlot = idxBReportPlot + 1
								container.BReportPlots[[idxBReportPlot]]  <- objBRPLOT
								# fnOpenPlot(objREPORT, objRPLOT)
								# RPLOT.run(objRPLOT, objREPORT)
								# fnClosePlot()
							}
							#}
						}
						if(strScriptFun == "RPLOT") {
								
							#if(iGWA == nGWAmt) {
							objRPLOT <- RPLOT(strCommand)
							RPLOT.REPORT.valid(objRPLOT, objREPORT)
							
							if(isValidScript) {
								idxReportPlot = idxReportPlot + 1
								container.ReportPlots[[idxReportPlot]]  <- objRPLOT
								# fnOpenPlot(objREPORT, objRPLOT)
								# RPLOT.run(objRPLOT, objREPORT)
								# fnClosePlot()
							}
							#}
						}
					
						.funtimeStop <- Sys.time()
						if(isTimeLog) 
							cat(paste("\t** Function evaluation time:", round(as.numeric(difftime(.funtimeStop,.funtimeStart,units="min")),2),"Minutes\n"))
							#cat(paste("\t** Function evaluation time:", .funtimeStop-.funtimeStart,"\n"))
						if(isMemoryLog) {
							# cat(paste("\t** Current memory allocated:",system(paste("top -b -n 1 -p ",Sys.getpid()," | grep \'",paste(Sys.getpid(),Sys.getenv("USER")),"\' | awk '{print $5}'", sep=""),intern=TRUE),"\n",sep=" "))
							# By Lucho:
							cat(paste("\t** Current memory allocated:",system(paste("ps -o pid,user,pri,nice,vsz ",Sys.getpid()," | grep \'",paste(Sys.getpid(),Sys.getenv("USER")),"\' | awk '{print $5}'", sep=""),intern=TRUE),"\n",sep=" "))
						}
						if(isGarbageCleaning) gc()
										
						if(dim(objGWA@tblGWA)[1]==0) {
							
							if(isValidScript) {
								cat(paste("EASY WARNING:",strScriptFun,"\n After evaluation of \n",strScriptFun,"\n the GWA data set is empty!!!\n", sep=""))
								cat(paste("Skipping Easy2 evaluation for \n",objGWA@fileInShortName,"!!!\n", sep=""))
							}
							# HANDLE HERE OPEN PLOTS
							break
						}
					} ## iCommand

					if(REPORT.getval(objREPORT, "numVarOut")=="NA") objREPORT	<-	REPORT.setval(objREPORT,"numVarOut",dim(objGWA@tblGWA)[1])
					
					if(isValidScript) REPORT.write(objREPORT)
					
					objREPORT 	<- 	REPORT.newrow(objREPORT)
					
					if(!isValidScript) {
						cat("OK")
					} 
					.filetimeStop <- Sys.time()
					if(isTimeLog) 
						cat(paste("** File processing time:", round(as.numeric(difftime(.filetimeStop,.filetimeStart,units="min")),2),"Minutes\n"))
					
				} ## iGWAmt
				
				container.REPORT[[iCommandGroup]] <- objREPORT
				container.container.ReportPlots[[iCommandGroup]] <- container.ReportPlots
				container.container.BReportPlots[[iCommandGroup]] <- container.BReportPlots
				
			} ## iCommandGroup
		} ## iMergeTag
		
		for(iCommandGroup in 1:nCommandGroups) {
			
			objREPORT <- container.REPORT[[iCommandGroup]]
			container.ReportPlots <- container.container.ReportPlots[[iCommandGroup]]
			container.BReportPlots <- container.container.BReportPlots[[iCommandGroup]]
			
			## create collected RPLOTS:
			if(length(container.ReportPlots)>0) {
				for(iRP in 1:length(container.ReportPlots)) {
					objRP <- container.ReportPlots[[iRP]]
					fnOpenPlot(objREPORT, objRP)
					RPLOT.run(objRP, objREPORT)
					fnClosePlot()
				}
			}
			## create collected BRPLOTS:
			if(length(container.BReportPlots)>0) {
				for(iBRP in 1:length(container.BReportPlots)) {
					objBRP <- container.BReportPlots[[iBRP]]
					objBRP <- BRPLOT.REPORT.resize(objBRP, objREPORT)
					fnOpenPlot(objREPORT, objBRP)
					BRPLOT.run(objBRP, objREPORT)
					fnClosePlot()
				}
			}
		}
		
	} ## iRun
	
	########################
	################################################################################################################################################
	################################################################################################################################################
	
	cat("\n")
	cat("\n")
	cat("+++++\n")
	.timeStop <- Sys.time()
	cat(paste("Succesfully finished Easy2:", .timeStop,"\n"))
	#cat(paste("Elapsed time:", .timeStop-.timeStart,"\n"))
	cat(paste("Elapsed time:", round(as.numeric(difftime(.timeStop,.timeStart,units="min")),2),"Minutes \n"))
	aMemRows <- grep("Current memory allocated",scan(fileEcfOut,what="character",sep="\n",quiet=TRUE),value=TRUE)
	astrMem=unlist(lapply(strsplit(aMemRows," "),function(x) x[length(x)]))
	anumMem=rep(0,length(astrMem))
	strMaxMemAllocOut = ""
	astrMemUnit = c("m","g","t")
	for(strMemUnit in astrMemUnit) {
		isMemSize <- grepl(strMemUnit, astrMem)
		if(any(isMemSize)) {
			astrMemTmp = astrMem[isMemSize]
			astrMemTmp <- gsub(",",".",astrMemTmp,fixed=T)
			anumMemTmp = as.numeric(gsub(strMemUnit,"",astrMemTmp))
			strMaxMemAllocOut <- paste(max(anumMemTmp,na.rm=T), strMemUnit,sep="")
		}
	}
	if(strMaxMemAllocOut=="") strMaxMemAllocOut <- max(as.numeric(gsub(",",".",astrMem,fixed=T)),na.rm=T)
	cat(paste("Maximum memory allocated:", strMaxMemAllocOut,"\n"))
	cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
	
	#graphics.off()
	#closeAllConnections()
	
###############################################################################################################################################
##############################################################################################################################################	
# # # ######################## UNCOMMENT FROM HERE

	if(blnReturnGwadata & blnReturnReport) {
		return(list(TRUE, objGWA, objREPORT))
	} else if(!blnReturnGwadata & blnReturnReport) {
		return(list(TRUE, objREPORT))
	} else if(!blnReturnGwadata & blnReturnReport) {
		return(list(TRUE, objGWA))
	} else {
		return(TRUE)
	}
	
}
Easy2 <- function(fileECF,blnValidityCheckOnly=FALSE,blnReturnGwadata=FALSE,blnReturnReport=FALSE,aFileIn=c(),fileMerge=NA,pathOut=NA) { 
		
	### Wrapper for Easy2.run
	valout <-TRUE
	
	if(!file.exists(fileECF))
		stop(paste("EASY ERROR:\n ECF-file \n ",fileECF,"\n does not exist!!!\n", sep=""))
	
	#graphics.off()
	#closeAllConnections()
	options(show.error.messages = FALSE)
	options(warn=0)
	
	## get current device and connections
	idxDeviceIn <- dev.cur()
	tblConIn 	<- showConnections()
	
	fileLog <- paste(fileECF,".out",sep="")
	if(file.exists(fileLog)) file.remove(fileLog)
	
	#connection.Logfile <- file(fileLog,open='w')
	connection.Logfile <- file(fileLog,open='a')
	sink(connection.Logfile, split=TRUE)

	out <-  NULL
	out <- 	try(withCallingHandlers(
					Easy2.run(fileECF,blnValidityCheckOnly,blnReturnGwadata,blnReturnReport,aFileIn,fileMerge,pathOut), 
					warning=function(w) { cat(w$message); cat("\n") })
				)
	#out <- try(withCallingHandlers(Easy2.run.warning(fileECF,fileOutLog), warning=function(w) {cat(w$message); cat("\n")}))
	#withCallingHandlers({out <- try(Easy2.run.warning(fileECF,fileOutLog))}, warning=function(w) {cat(w$message); cat("\n")})
	
	if(class(out) == "try-error") {
		cat(out)
		valout <- FALSE
	} else {
		valout <- out
	}
	cat("\n")
	options(show.error.messages = TRUE)
	
	## close all open devices except those that were open at the start
	idxDeviceOut <- dev.cur()
	while(idxDeviceOut!=idxDeviceIn) {
		dev.off()
		idxDeviceOut <- dev.cur()
	}
	## close all open connections except those that were open at the start
	sink()
	close(connection.Logfile)
	
	#	Copy log file to output folder
	#if(paste(fileECF,".out",sep="") != paste(pathOut,"/",fileECFName,".out",sep="")) {
	aRowLog<-scan(fileLog,sep="\n",what="character",quiet=T)
	pathOut <- aRowLog[which(aRowLog=="Default output path is ")+1]
	
	fileEcfName <- strsplit(fileECF,"/")[[1]][length(strsplit(fileECF,"/")[[1]])]
	pathOut <- ifelse(pathOut==".",getwd(),pathOut)
	pathOut <- ifelse(substring(pathOut,nchar(pathOut),nchar(pathOut))!="/",paste(pathOut, "/" , sep=""),pathOut)
	fileLogOut <- paste(pathOut,fileEcfName,".out",sep="")
	fileLogName <- strsplit(fileLog,"/")[[1]][length(strsplit(fileLog,"/")[[1]])]
	if(fileLog != fileLogOut & paste(getwd(),"/",fileLogName,sep="") != fileLogOut) {
		file.copy(fileLog,fileLogOut,overwrite=TRUE)
		file.remove(fileLog)
	}
	
	#}
	# tblConOut 	<- showConnections()
	# while(nrow(tblConOut)>nrow(tblConIn)) {
		# sink()
		# tblConOut <- showConnections()
	# }
	#graphics.off()
	#closeAllConnections()
	return(valout)
}
EasyStrata2 <- function(fileECF,blnValidityCheckOnly=FALSE,blnReturnGwadata=FALSE,blnReturnReport=FALSE,aFileIn=c(),fileMerge=NA,pathOut=NA) { 


	return(Easy2(fileECF,blnValidityCheckOnly,blnReturnGwadata,blnReturnReport,aFileIn,fileMerge,pathOut))
}
EasyQC2 <- function(fileECF,blnValidityCheckOnly=FALSE,blnReturnGwadata=FALSE,blnReturnReport=FALSE,aFileIn=c(),fileMerge=NA,pathOut=NA) { 
	
	return(Easy2(fileECF,blnValidityCheckOnly,blnReturnGwadata,blnReturnReport,aFileIn,fileMerge,pathOut))
}
# #### For testing:
# EasyQC(fileECF,blnValidityCheckOnly,blnReturnGwadata,blnReturnReport,aFileIn)

