IMA2.methy450R <- function (fileName, columnGrepPattern = list(beta = ".AVG_Beta", detectp = ".Detection.Pval"), groupfile, writePDF = FALSE, ...) {
    cat("................Reading data................\n")
    temp <- readLines(fileName, n = 20)
    nskip <- grep(columnGrepPattern$beta, temp, ignore.case = TRUE) - 1
    titleLine <- temp[nskip + 1]
    allcolname <- unlist(strsplit(titleLine, "\t"), use.names = FALSE)
    TargetID <- grep("TargetID", allcolname)
    if (length(TargetID)==0) {
        TargetID = 1
    } else {}

    betacol <- grep(columnGrepPattern$beta, allcolname)
    pvalcol <- grep(columnGrepPattern$detectp, allcolname)
    annotcol <- grep("ILMNID", allcolname):(length(allcolname))
    colClasses <- rep("NULL", length(allcolname))
    colClasses[TargetID] <- "character"
    colClasses[betacol] <- "numeric"
    colClasses[pvalcol] <- "numeric"
    colClasses[annotcol] <- "character"

    dataRead <- read.delim(file = fileName, header = TRUE, skip = nskip, colClasses = colClasses, sep = "\t",
                            na.string = c("NA", ""), check.names = FALSE, strip.white = TRUE, row.names = TargetID, ...)
    dataName <- colnames(dataRead)
    cat("......Extracting the beta value matrix......\n")
    betaCol <- grep(columnGrepPattern$beta, dataName)
    betamatrix <- dataRead[, betaCol]

    cat("........Extracting the pvalue matrix........\n")
    pvalCol <- grep(columnGrepPattern$detectp, dataName)
    detect_p <- dataRead[, pvalCol]

    cat("......Extracting the annotation matrix......\n")
    annotCol <- grep("ILMNID", dataName):(length(dataName))
    annotation <- dataRead[, annotCol]

    sname <- sub(columnGrepPattern$beta, "", colnames(betamatrix))
    colnames(betamatrix) <- sname
    colnames(detect_p) <- sname

    cat("...........Reading phenotype data...........\n")
    group <- data.frame(read.delim(groupfile, sep = "\t", header = TRUE))

    cat("Matching the orders of samples between phenotype data and beta value matrix.\n")
    index <- which(group[, 1]%in%colnames(betamatrix))
    if (sum(is.na(index)) > 0) {
        cat("ERROR:\nBelow samples couldn't be found in beta matrix:\n", group[is.na(index), 1], "\n.")
    } else {}

    betamatrix <- betamatrix[, as.character(group[index, 1])]
    detect_p <- detect_p[, as.character(group[index, 1])]
    rownames(group) <- group[index, 1]
    cat("Total CpG sites without any filtering are:", nrow(betamatrix), "\nTotal samples are:", ncol(betamatrix), "\n")

    x.methy450 <- new("exprmethy450", bmatrix = as.matrix(betamatrix), annot = as.matrix(annotation), detectP = as.matrix(detect_p), groupinfo = group)

    cat("..........Starting Quality Control..........\n")
    eset <- na.omit(betamatrix)
    samples <- paste(group[, 1], group[, 2], sep = "_")
    hc1 <- hclust(cor.dist(t(eset)), method = "average")
    if (writePDF) {
        pdf("./QC.pdf", width = 24)
            plot(hc1, samples, xlab = "Sample", main = "Clustering samples by all the CpG loci ", lwd = 2, font.axis = 2, font.lab = 2)
            boxplot(betamatrix, ylab = "beta Value")
            avgPval <- apply(detect_p, 2, function (x) { sum(x >= 1e-05) * 100/length(x) })
            barplot(avgPval, ylab = "% of detect pvalue >1e-5")
        dev.off()
    } else {}

    cat("An exprmethy450 class are created and the slotNames are:\n",
    slotNames(x.methy450), "\n")

    cat("Basic Quality Control information can be found in QC.pdf file\n")
    return(x.methy450)
}
