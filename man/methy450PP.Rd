\name{IMA2.methy450PP}
\alias{IMA2.methy450PP}
\alias{methy450batch}

\title{Data preprocessing and quality control}

\description{It allows user to choose several filtering steps or modify filtering criteria for specific quality control purpose. These include whether or not to filter probes based on detection P-value; whether or not to remove the loci from the X or Y chromosome, or both; whether or not to perform peak transformation, whether or not to transfer the raw beta value using either arcsine square root or logit; whether or not to perform quantile normalization; whether or not to remove the loci containing missing beta values; whether or not to filter out loci whose methylation level are measured by probes containing SNP(s) at/near the targeted CpG site. The user can choose the preprocessing routes and corresponding cutoffs in the argument of this function.}

\usage{
IMA2.methy450PP(
    data,
    na.omit = TRUE,
    peakcorrection = FALSE,
    normalization = FALSE,
    transfm = c(FALSE, "arcsinsqr", "logit"),
    samplefilterdetectP = c(FALSE, 1e-05),
    samplefilterperc = 0.75,
    sitefilterdetectP = c(FALSE, 0.05),
    sitefilterperc = 0.75,
    locidiff = c(FALSE, 0.01),
    locidiffgroup = list("g1","g2"),
    XYchrom = c(FALSE,"X","Y", c("X","Y")),
    snpfilter = c(FALSE, "snpsites.txt")
)
}

\arguments{
    \item{data}{an exprmethy450 class returned by the IMA2.methy450R function}
    \item{na.omit}{if TRUE remove the sites containing missing value}
    \item{peakcorrection}{if TRUE, peak correction is performed based on the paper by sarah Dedeurwaerder et al.}
    \item{normalization}{if TRUE, quantile normalization performed}
    \item{transfm}{if FALSE, no transfm is performed, "arcsinsqr":arcsine square root transformation on beta value is performed, "logit":logit transformation on beta is performed}
    \item{samplefilterdetectP}{Default is false, i.e, no sample filtering by detection P-value. Otherwise, choose the cut off of detection P-value.}
    \item{samplefilterperc}{Keep the samples having at least specified percentage of sites with detection P-value less than the \bold{samplefilterdetectP}.}
    \item{sitefilterdetectP}{Default is false, i.e. no site filtering by detection p-value. Otherwise, choose the cut off of detection P-value.}
    \item{sitefilterperc}{Remove the sites having specified percentage of samples with detection P-value greater than \bold{sitefilterdetectP}.}
    \item{locidiff}{if FALSE, don't filter sites by the difference of group beta value. Otherwise, remove the sites with beta value difference greater than the specified value.}
    \item{locidiffgroup}{specify which two groups are considered to check the loci difference if locidiff is not true}
    \item{XYchrom}{if "X", remove the sites on chromosome X , if "Y", remove the sites on chromosome Y, if c("X","Y"), remove both on chromosome X and Y.}
    \item{snpfilter}{if FALSE, keep the loci whose methylation level are measured by probes containing SNP(s) at/near the targeted CpG site; otherwise filter out the list of snp-containing loci by specifying the snp file name and location}
}

\details{
    It allows user to choose several filtering steps or modify filtering criteria for specific quality control purpose.
    By default, IMA will filter out loci with missing beta value, from the X chromosome or with median detection P-value greater than 0.05.
    Users can choose to filter out loci whose methylation level are measured by probes containing SNP(s) at/near the targeted CpG site.
    The option for sample level quality control is also provided. Although the raw beta values will be analyzed as recommended by Illumina,
    users can choose Arcsine square root transformation when modeling the methylation level as the response in a linear model. Logit transformation
    is also available as an option. The default setting in IMA package for preprocessing is that no normalization will be performed.
    Although quantile normalization is available as an alternative preprocessing option, it should be pointed out that
    several literatures show that quantile normalization does not remove unwanted technical variation between samples in methylation analysis.
}

\value{
    This function will return a methy450batch class including:
    \item{bmatrix}{ a matrix of beta value for individual sites}
    \item{detectP}{ a matrix of detection p-value for individual sites}
    \item{annot}{a matrix of annotation information for individual sites}
    \item{groupinfo}{ a list of sample ID and phenotype of each sample}
    \item{TSS1500Ind}{ two lists of IDs - SID (site IDs) and PID (Position IDs) belonging to the TSS1500 region of each gene  }
    \item{TSS200Ind}{ two lists of IDs - SID (site IDs) and PID (Position IDs) belonging to the TSS200 region of each gene }
    \item{UTR5Ind}{ two lists of IDs - SID (site IDs) and PID (Position IDs) belonging to the 5' UTR region of each gene }
    \item{EXON1Ind}{two lists of IDs - SID (site IDs) and PID (Position IDs) belonging to the 1st EXON of each gene}
    \item{UTR3Ind}{two lists of IDs - SID (site IDs) and PID (Position IDs) belonging to the 3' UTR region of each gene}
    \item{GENEBODYInd}{two lists of IDs - SID (site IDs) and PID (Position IDs) belonging to the gene body region of each gene}
    \item{ISLANDInd}{two lists of IDS - SID (site IDs) and PID (Position IDs) belonging to the ISLAND region of each UCSC_CPG_ISLAND}
    \item{NSHOREInd}{two lists of IDs - SID (site IDs) and PID (Position IDs) belonging to the N Shore region of each UCSC_CPG_ISLAND}
    \item{SSHOREInd}{two lists of IDs - SID (site IDs) and PID (Position IDs) belonging to the S Shore region of each UCSC_CPG_ISLAND}
    \item{NSHELFInd}{two lists of IDs - SID (site IDs) and PID (Position IDs) belonging to the N Shelf region of each UCSC_CPG_ISLAND}
    \item{SSHELFInd}{two lists of IDs - SID (site IDs) and PID (Position IDs) belonging to the S Shelf region of each UCSC_CPG_ISLAND}
}

\author{
Dan Wang, Li Yan, Qiang Hu, Dominic J Smiraglia, Song Liu
Mickaël Canouil
}

\seealso{\code{\link{IMA2.methy450R}}, \code{\link{IMA2.sitetest}}}

\examples{
\dontrun{
setwd(system.file("extdata", package = "IMA2"))
MethyFileName <- "SampleMethFinalReport.txt"
PhenoFileName <- "SamplePhenotype.txt"
data <- IMA2.methy450R(
    file = MethyFileName,
    columnGrepPattern = list(beta = ".AVG_Beta", detectp = ".Detection.Pval"),
    groupfile = PhenoFileName
)
dataf <- IMA2.methy450PP(
    data,
    na.omit = TRUE,
    normalization = FALSE,
    transfm = FALSE,
    peakcorrection = TRUE,
    samplefilterdetectP = 1e-5,
    samplefilterperc =0.75,
    sitefilterdetectP = 0.05,
    sitefilterperc = 0.5,
    locidiff = FALSE,
    XYchrom = FALSE,
    snpfilter = FALSE
)
}
}
