\name{IMA2-package}
\alias{IMA2-package}
\alias{IMA2}
\alias{IMA2-class}
\alias{snpsites-dataset}
\alias{transcript}
\docType{package}
\title{IMA2 (Illumina Methylation Analyzer 2)}
\description{
IMA2 is based on IMA package.

IMA2 is a package designed to automate the pipeline for exploratory analysis and summarization of site-level and region-level methylation changes in epigenetic studies utilizing the 450K DNA methylation microarray
IMA2 automates the tasks commonly required for the exploratory analysis and summarization of epigenetic data sets utilizing the 450K DNA methylation microarray. The package makes use of Illumina methylation annotation for region definition, as well as several Bioconductor packages for various preprocessing and differential testing steps. There are two major differences between IMA2 and existing packages for Infinium methylation microarray analysis. First, instead of analyzing CpG site only, IMA2 provide both site-level and region-level methylation analysis. Second, instead of manually calling individual R functions at the command line, IMA2 provides a pipeline which automate the tasks commonly required for the exploratory analysis and summarization of 450K microarray data. The user can either run the pipeline with default setting or specify optional routes in the parameter file of pipeline.

The main purpose of developing IMA2 package is to provide a range of commonly used analysis options for potential users to perform exploratory analysis and summarization of 450K microarray data in an automatic way. It is the best interest for the users to consult experienced bioinformatician/statistician about which specific analysis option should be chosen for their 450k microarray data. Written in open source R environment, it provides the flexibility for users to adopt, extend and customize the functionality for their specific needs. It can be used as an automatic pipeline to analyze specific regions as well as specific sites for downstream functional exploration and hypothesis generation.
}

\details{
    \tabular{ll}{
    Package: \tab IMA2\cr
    Title: \tab IMA2 (Illumina Methylation Analyzer 2)\cr
    Author: \tab Dan Wang, Li Yang, Qiang Hu, Dominic J Smiraglia, Song Liu\cr
    Contributor: \tab Mickaël Canouil\cr
    Maintainer: \tab Mickael Canouil <mickael.canouil at good.ibl.fr>\cr
    License: \tab GPL (>= 2)\cr
    Depends: \tab R (>= 3.0.0), methods, parallel, bioDist, graphics, grDevices, limma, preprocessCore, stats, utils\cr
    URL: \tab https://github.com/mcanouil/snpEnrichment\cr
    Encoding: \tab UTF-8\cr
    }
}

\author{
Dan Wang, Li Yang, Qiang Hu, Dominic J Smiraglia, Song Liu
Mickaël Canouil
}

\seealso{\code{\link{IMA2.methy450R}}, \code{\link{IMA2.methy450PP}}, \code{\link{IMA2.sitetest}}}
\examples{
\dontrun{
rm(list=ls())
options(stringAsfactors = FALSE)

### options in IMA2.methy450R ###
## Load data ##
# Specify the location of your R library
libPaths = "/home/R/myR" # Specify the location of your R library
# Specfiy the original methylation data produced by the GenomeStudio
MethyFileName = "../data/SampleMethFinalReport_Smiraglia.txt"
# Specify the phenotype for each sample
PhenoFileName = "../data/SamplePhenotype.txt"


### Preprocessing:IMA2.methy450PP ###
# The cutoff for sample-level detection Pvalue
samplefilterdetectP = 1e-5

# The percent of loci with detection Pvalue less than "samplefilterdetectP" in each sample
samplefilterperc = 0.75

# The cutoff for site-level detection Pvalue
sitefilterdetectP = 0.05

# The percent of samples with detection Pvalue less than "sitefilterdetectP" for each site
sitefilterperc = 0.5

# Remove the sites containing missing beta value
na.omit = TRUE

# Remove the sites on chromosome X
XYchrom = FALSE

# If TRUE, peak correction is performed
peakcorrection = FALSE

# If TRUE, quantile normalization performed
normalization = FALSE

# If FALSE, no transform is performed; if "arcsinsqr", arcsin square root
# transformation is performed; if "logit", logit transformation is performed
transfm = FALSE

# If FALSE, don't filter sites by the difference of group beta value. Otherwise,
# remove the sites with beta value difference smaller than the specified value
locidiff = FALSE

# Specify which two groups are considered to check the loci difference
#(if "locidiff" is not true)
locidiffgroup = c("g1","g2")

# If FALSE, keep the loci whose methylation level are measured by probes containing SNP(s) at/near
# the targeted CpG site; otherwise, filter out the list of SNP containing loci
# by specifying the snp file name and location
snpfilter = FALSE

# A list of SNP-containing probes (based on dbSNP v132) could be accessed by the command:
# snpfilter = system.file("extdata/snpsites.txt",package ="IMA")


### sitetest ###
# Other options of differential testing methods:
#"wilcox"/"pooled"/"satterthwaite" for the comparison between two group
testmethod = "limma"

# If "ON", covariates is continuous variable
concov = "OFF"

# Specify the case group index in the sample.txt file (if "concov" is "ON")
gcase = "g2"

# Specify the control group index in the sample.txt file (if "concov" is "ON")
gcontrol = "g1"

# Options for multiple testing correction.
#The user can choose the methods provided by p.adjust function of R stat package
Padj = "BH"

# Options for deriving an index of overall methylation value of each region.
#mean/median/tbrm: "tbrm" is Tukey's Biweight robust average
indexmethod = "mean"

# If ture, the differential test methods would change to the corresponding paired-test methods
paired = FALSE


### output the differential sites ###
# cut off for raw pvalue
rawpcut = NULL

# cut off for adjusted pvalue
adjustpcut = NULL

# cut off for beta value difference
betadiffcut = NULL


### Analysis Routes ###
# load the IMA2 package
library(IMA2)

# load the data
data <- IMA2.methy450R(
    fileName = MethyFileName,
    columnGrepPattern = list(beta = ".AVG_Beta", detectp = ".Detection.Pval"),
    groupfile = PhenoFileName
)

# QC filtering
dataf <- IMA2.methy450PP(
    data,
    na.omit = na.omit,
    normalization = normalization,
    peakcorrection = peakcorrection,
    transfm = transfm,
    samplefilterdetectP = samplefilterdetectP,
    samplefilterperc = samplefilterperc,
    sitefilterdetectP = sitefilterdetectP,
    locidiff = locidiff,
    locidiffgroup = locidiffgroup,
    XYchrom = XYchrom,
    snpfilter = snpfilter
)

# site-level testing with the "BH" adjustment
sitetest <- IMA2.sitetest(
    dataf,
    gcase = gcase,
    gcontrol = gcontrol,
    concov = concov,
    testmethod = testmethod,
    Padj = Padj,
    rawpcut = rawpcut,
    adjustpcut = adjustpcut,
    betadiffcut = betadiffcut,
    paired = paired
)
}
}

