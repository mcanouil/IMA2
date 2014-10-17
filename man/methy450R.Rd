\name{IMA2.methy450R}
\alias{IMA2.methy450R}
\alias{exprmethy450}

\title{Load methylation 450k data}

\description{This function will load the methylation 450k data. The input information for the package consists basically of two files containing beta-value methylation data (including annotation) produced by BeadStudio or GenomeStudio, and sample phenotype data prepared by the user. \emph{IMA2.methy450R} function loads the input files with a single command described below and an \emph{exprmethy450} object will be created, which includes the following features: \emph{beta} value matrix, locus annotation, detection P-value and sample phenotype information. Besides, basic quality control information will be outputted in the QC.pdf, which include unsupervised sample clustering using all loci, boxplot for beta value of each sample, and barplot showing the percent of loci with detection P-value smaller than 1e-5 in each sample.}
\usage{
IMA2.methy450R(
    fileName,
    columnGrepPattern = list(beta = ".AVG_Beta", detectp = ".Detection.Pval"),
    groupfile,
    writePDF = FALSE
)
}

\arguments{
    \item{fileName}{This is the input file containing the 450k methylation data including beta value, detection p-value and loci annotation information. This file can be be produced from the illumina BeadStudio or GenomeStudio software.}
    \item{columnGrepPattern}{Specify the columns with corresponding characters within the 450k methylation data file produced by illumina BeadStudio or GenomeStudio software. These columns will be loaded into the function.}
    \item{groupfile}{This is the input file containing the phenotype data prepared by the users.}
    \item{writePDF}{Print control quality summary in a PDF file.}
    % \item{...}{Arguments in read.delim(), such as parameter "sep" could be set to sep="," if the GenoneStudio file is separated by comma instead of tab.}
}

\value{
This function will return an \emph{exprmethy450} class including these features: a beta value matrix, an annotation matrix, a detection p-value matrix and a list of sample phenotype information.
A QC.pdf file including the basic quality control information will be generated.
\item{bmatrix}{ the beta value matrix for each site in each sample}
\item{detectP}{ a matrix of detection p-value for each site in each sample}
\item{annot}{a matrix of annotation information for each targeted site}
\item{groupinfo}{a matrix of sample ID and phenotype information for each sample}
}

\author{
Dan Wang, Li Yan, Qiang Hu, Dominic J Smiraglia, Song Liu
Mickaël Canouil
}

\seealso{\code{\link{IMA2.methy450PP}}, \code{\link{IMA2.sitetest}}}

\examples{
\dontrun{
setwd(system.file("extdata", package = "IMA2"))
MethyFileName <- "SampleMethFinalReport.txt"
PhenoFileName <- "SamplePhenotype.txt"
data <- IMA2.methy450R(
    fileName = MethyFileName,
    columnGrepPattern = list(beta = ".AVG_Beta", detectp = ".Detection.Pval"),
    groupfile = PhenoFileName
)
}
}
