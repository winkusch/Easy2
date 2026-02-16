\name{EasyStrata2}
\alias{EasyStrata2}
\title{
Function to run EasyStrata2 (1.3.4.b42).
}
\description{
The Function takes the an EasyStrata2 config/script (ECF-) file and performs all steps defined in the ECF-file.
}
\usage{
EasyStrata2(fileECF,blnValidityCheckOnly,blnReturnGwadata,blnReturnReport,aFileIn,fileMerge,pathOut)
}
\arguments{
\item{fileECF}{
fileECF is an EasyStrata2 config/script (ECF) file.\cr
Please see 'EasyQC_9.0_Commands_140918.pdf' for more information on ECF functions and parameters.
}
\item{blnValidityCheckOnly}{
blnValidityCheckOnly is a logical value that indicates whether the full pipeline (blnValidityCheckOnly=FALSE, default) or only the validity check (blnValidityCheckOnly=TRUE) are supposed to be run.
}
\item{blnReturnGwadata}{
blnReturnGwadata is a logical value that indicates whether the most recent GWADATA object should be returned (default blnReturnGwadata=FALSE).
}
\item{blnReturnReport}{
blnReturnReport is a logical value that indicates whether the REPORT object should be returned (default blnReturnReport=FALSE).
}
\item{aFileIn}{
aFileIn is an array that pases paths to input files to the ecf (aFileIn=c(), default).
}
\item{fileMerge}{
fileMerge is the paths to a reference file that will be used by MERGE (fileMerge=NA, default; any value of --fileRef will be overwritten in case of fileMerge being used).
}
\item{pathOut}{
pathOut is the paths to store results (pathOut=NA, default; any value of --pathOut at DEFINE will be overwritten in case of pathOut being used).
}
}
\references{
http://www.genepi-regensburg.de/easystrata
}
\examples{
\donttest{installDir=system.file("extdata", package="Easy2")}
## Example Pipeline:
\donttest{ecfPipe=paste(extDataDir,"example_pipeline.ecf",sep="/")}
\donttest{EasyStrata2(ecfPipe)}
}
\author{
Thomas Winkler <thomas.winkler@klinik.uni-regensburg.de>
}