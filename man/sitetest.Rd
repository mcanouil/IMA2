\name{IMA2.sitetest}
\alias{IMA2.sitetest}

\title{Site-level Differential Methylation Analysis}

\description{For each specific site, Wilcoxon rank-sum test (default), Student's t-test ("pooled" or "satterthwaite") and empirical Bayes statistics are available for inference in differential testing. Robust linear regression is available as an option to infer methylation change associated with continuous variable (e.g., age). A variety of multiple testing correction algorithms is available, including conservative Bonferroni correction and more liberal false discovery rate control. Users can specify the significance criteria in the parameter file.}

\usage{
IMA2.sitetest(
    dataf,
    gcase = "g2",
    gcontrol = "g1",
    testmethod = c("wilcox", "limma", "pooled", "satterthwaite"),
    Padj = "BH",
    concov = "OFF",
    rawpcut = NULL,
    adjustpcut = NULL,
    betadiffcut = NULL,
    paired = FALSE,
    nCPU = NULL,
    verbose = TRUE
)
}

\arguments{
    \item{dataf}{a \emph{methy450batch} class returned by the \emph{IMA2.methy450PP} function}
    \item{gcase}{the string names of case group/groups in the Phenotype.txt file}
    \item{gcontrol}{the string names of control group/groups in the Phenotype.txt file}
    \item{testmethod}{Wilcoxon rank-sum test (default), Student's t-test ("pooled" or "satterthwaite"), or empirical Bayes statistics.}
    \item{Padj}{The methods available for multiple testing correction. Users could choose any method provided in the \emph{p.adjust} function of R \emph{stat} package.}
    \item{concov}{if "ON", covariates is continuous variable, and the linear model would be used for testing the association between methylation level and continuous phenotype (e.g., age).}
    \item{rawpcut}{either null or the output cut-off for the raw p-value}
    \item{adjustpcut}{either null or the output cut-off for the adjusted p-value}
    \item{betadiffcut}{either null or the output cut-off for the beta value difference between two groups}
    \item{paired}{if TRUE, the test method would change to the corresponding paired-test method}
    \item{nCPU}{number of core to use.}
    \item{verbose}{if TRUE, print message.}
}

\value{return a matrix with "P-Value", "Adjust Pval", "beta-Difference" in separate columns for each site}

\author{
Dan Wang, Li Yan, Qiang Hu, Dominic J Smiraglia, Song Liu
Mickaël Canouil
}

\seealso{\code{\link{IMA2.methy450R}}, \code{\link{IMA2.methy450PP}}}
\examples{
\dontrun{
IMA2.sitetest(
    dataf,
    gcase = "g2",
    gcontrol = "g1",
    testmethod = c("wilcox", "limma", "pooled", "satterthwaite"),
    Padj = "BH",
    concov = "OFF",
    rawpcut = NULL,
    adjustpcut = NULL,
    betadiffcut = NULL,
    paired = FALSE
)
}
}
