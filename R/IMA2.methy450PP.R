IMA2.methy450PP <- function (data, na.omit = TRUE, peakcorrection = FALSE, normalization = FALSE, transfm = c(FALSE, "arcsinsqr", "logit"),
                    samplefilterdetectP = c(FALSE, 1e-05), samplefilterperc = 0.75, sitefilterdetectP = c(FALSE, 0.05), sitefilterperc = 0.75, locidiff = c(FALSE, 0.01),
                    locidiffgroup = list("g1", "g2"), XYchrom = c(FALSE, "X", "Y", c("X", "Y")), snpfilter = c(FALSE, "snpsites.txt")) {
    Beta2M <- function (B) {
        return(log2(B / (1 - B)))
    }
    M2Beta <- function (M) {
        return((2 ^ M) / (2 ^ M + 1))
    }
    summits <- function (BetaValues) {
        d <- density(BetaValues)
        yneg <- d$y[1:which(d$x > M2Beta(0.0))[1]]
        ypos <- d$y[which(d$x > M2Beta(0.0))[1]:length(d$y)]
        sa <- d$x[which(d$y == max(yneg))]
        sb <- d$x[which(d$y == max(ypos))]
        return(c(Beta2M(sa), Beta2M(sb)))
    }
    correctI <- function (BetaValues, SI, SII) {
            return(BetaValues)
    }
    correctII <- function (BetaValues, SI, SII) {
            M <- Beta2M(BetaValues)
            sigma_u  <- SII[1] / SI[1]
            sigma_m  <- SII[2] / SI[2]
            M <- sapply(M, function (x) {
                if (x < 0) {
                    return(x / sigma_u)
                } else {
                    return(x / sigma_m)
                }
            })
            return(M2Beta(M))
    }
    peak.correction <- function (data, anno) {
        anno <- anno[row.names(data), ]
        TI  <- anno[,"INFINIUM_DESIGN_TYPE"] == 'I'
        TII <- anno[,"INFINIUM_DESIGN_TYPE"] == 'II'
        corrected.data <- apply(data, 2, function (B) {
            SI <- summits(B[TI])
            SII <- summits(B[TII])
            BI <- correctI(as.vector(B[TI]), SI, SII)
            BII <- correctII(as.vector(B[TII]), SI, SII)
            return(c(BI, BII))
        })
        row.names(corrected.data) <- c(row.names(data[TI, ]), row.names(data[TII, ]))
        return(corrected.data)
    }

    bmatrix <- data@bmatrix
    detect_p <- data@detectP
    annotation <- data@annot
    groupinfo <- data@groupinfo
    orignalrownm <- rownames(bmatrix)
    if (samplefilterdetectP) {
        goodsample <- colSums(detect_p <= samplefilterdetectP) >= samplefilterperc * nrow(detect_p)
        sample_names <- colnames(bmatrix)
        bmatrix <- bmatrix[, goodsample]
        detect_p <- detect_p[, goodsample]
        goodsample_names <- colnames(bmatrix)
        cat(abs(ncol(bmatrix) - length(goodsample)), "samples removed with at least",
            samplefilterperc * 100, "percentage sites having pvalue greater than",
            samplefilterdetectP, ".\n")
        tmp <- setdiff(sample_names, goodsample_names)
        if (length(tmp)!=0) {
            cat(paste0(tmp, collapse = ", "), "is/are removed.\n")
        } else {}
        groupinfo <- groupinfo[goodsample, ]
    }
    if (na.omit) {
        bmatrix <- na.omit(bmatrix)
        temp <- orignalrownm %in% rownames(bmatrix)
        detect_p <- detect_p[temp, ]
        annotation <- annotation[temp, ]
        temp <- nrow(data@bmatrix) - nrow(bmatrix)
        cat(temp, "sites contain missing value and are removed", ".\n")
    }
    if (snpfilter != FALSE) {
        snpsites <- read.delim(snpfilter, sep = "\t", stringsAsFactors = FALSE)[, "TargetID"]
        index <- rownames(bmatrix) %in% snpsites
        bmatrix <- bmatrix[!index, ]
        detect_p <- detect_p[!index, ]
        annotation <- annotation[!index, ]
        cat(sum(index), "sites contain snps and removed, \n")
    }
    if (XYchrom[1] != FALSE) {
        chr <- annotation[, "CHR"]
        index <- which(chr %in% XYchrom)
        good_chrom <- rownames(annotation)[-index]
        cat(length(index), "sites on chr", XYchrom, "are removed.\n")
    } else {
        good_chrom <- rownames(bmatrix)
    }
    if (sitefilterdetectP) {
        good_loci <- rownames(detect_p)[rowSums(detect_p <= sitefilterdetectP) >= sitefilterperc * ncol(detect_p)]
        cat(nrow(detect_p) - length(good_loci), "sites had at least",
            sitefilterperc * 100, "% samples with pvalue greater than",
            sitefilterdetectP, "and are removed.\n")
    } else {
        good_loci <- rownames(bmatrix)
    }
    if (locidiff) {
        c1 <- groupinfo[, 2] %in% locidiffgroup[[1]]
        c2 <- groupinfo[, 2] %in% locidiffgroup[[2]]
        con_mean <- rowMeans(bmatrix[, c1])
        trt_mean <- rowMeans(bmatrix[, c2])
        good_diff <- rownames(bmatrix)[abs(trt_mean - con_mean) >= locidiff]
        cat(length(good_diff), "sites had the beta difference between group great than",
            locidiff, "and are kept for the downstream analysis \n")
    } else {
        good_diff <- rownames(bmatrix)
    }
    all_good <- intersect(intersect(good_chrom, good_loci), good_diff)
    cat(length(all_good), "sites were retained from the original", length(orignalrownm), "sites.\n")
    bmatrix <- bmatrix[all_good, ]
    annotation <- annotation[all_good, ]
    detect_p <- detect_p[all_good, ]
    if (peakcorrection) {
        if (!na.omit) {
            cat("\tMissing value exist in the orignial data, \nPlease remove the missing value before peak correction, use na.omit = TRUE\n")
        }
        cat("Peak correction...\nThis part of code was provided by Matthieu Defrance <defrance@bigre.ulb.ac.be>\n")
        cat("Thanks for sharing the code with us.\n")
        cat("Dimension of beta matrix", dim(bmatrix), "\n")
        cat("Dimension of annotation", dim(annotation), "\n")
        bmatrix <- peak.correction(bmatrix, annotation)
        bmatrix <- bmatrix[rownames(annotation), ]
    } else {
        cat("\n")
        cat("\n")
        cat("Dimension of beta matrix", dim(bmatrix), "\n")
        cat("Dimension of annotation", dim(annotation), "\n")
    }
    if (normalization) {
        bmatrix <- normalize.quantiles(as.matrix(bmatrix))
        colnames(bmatrix) <- colnames(detect_p)
        rownames(bmatrix) <- rownames(detect_p)
        cat("Quantile normalization Performed\n")
    } else {}
    if (transfm == "arcsinsqr") {
        if (na.omit) {
            bmatrix <- asin(sqrt(bmatrix))
            cat("Transfer beta matrix by the arcsin square root\n")
        } else {
            cat("\tMissing value exist in the orignial data, \nPlease remove the missing value before transformation, use na.omit = TRUE\n")
            stop
        }
    } else {}
    if (transfm == "logit") {
        if (na.omit) {
            bmatrix[bmatrix == 0] <- min(bmatrix[bmatrix > 0], 0.001)/10
            bmatrix[bmatrix == 1] <- max(bmatrix[bmatrix < 1], 0.999) + (1 - max(bmatrix[bmatrix < 1], 0.999))/100
            bmatrix <- log(bmatrix/(1 - bmatrix))
            cat("Transfer beta matrix by the logit transformation \n")
        } else {
            cat("\tMissing value exist in the orignial data, \nPlease remove the missing value before transformation, use na.omit = TRUE\n")
            stop
        }
    } else {}
    cat(".......Split the annotation file to 11 annotated region categories.......\n\n")
    annot <- annotation
    name <- "UCSC_REFGENE_NAME"
    cpGsite <- as.character(annot[, 1])
    genelist <- strsplit(as.character(annot[, name]), ";")
    genelist[which(genelist == "character(0)")] = "NA"
    name <- "UCSC_REFGENE_GROUP"
    refgene <- strsplit(as.character(annot[, name]), ";")
    refgene[which(refgene == "character(0)")] <- "NA"
    listlength <- lapply(refgene, length)
    listlength[listlength == 0] <- 1
    col0 <- rep(1:nrow(annot), listlength)
    col1 <- rep(cpGsite, listlength)
    col2 <- unlist(genelist)
    col3 <- unlist(refgene)
    col4 <- rep(as.character(annotation[, "RELATION_TO_UCSC_CPG_ISLAND"]), listlength)
    col5 <- rep(as.character(annotation[, "UCSC_CPG_ISLANDS_NAME"]), listlength)
    splitToRegionlist <- function (grepname = c("TSS1500", "TSS200", "5'UTR", "1stExon", "Gene Body", "3'UTR")) {
        index <- col3 == grepname
        col1sub <- col1[index]
        col2sub <- col2[index]
        temp <- split(col1sub, col2sub)
        returnSID <- lapply(temp, unique)
        col0sub <- col0[index]
        temp <- split(col0sub, col2sub)
        returnPID <- lapply(temp, unique)
        return(Ind <- list(SID = returnSID, PID = returnPID))
    }
    TSS1500Ind <- splitToRegionlist(grepname = "TSS1500")
    TSS200Ind <- splitToRegionlist(grepname = "TSS200")
    UTR5Ind <- splitToRegionlist(grepname = "5'UTR")
    EXON1Ind <- splitToRegionlist(grepname = "1stExon")
    GENEBODYInd <- splitToRegionlist(grepname = "Body")
    UTR3Ind <- splitToRegionlist(grepname = "3'UTR")
    cat("TSS1500 region contains:   ", length(TSS1500Ind$SID), "UCSC REFGENE region \n")
    cat("TSS200 region contains:    ", length(TSS200Ind$SID), "UCSC REFGENE region\n")
    cat("5'UTR region contains:     ", length(UTR5Ind$SID), "UCSC REFGENE region\n")
    cat("1st Exon region contains:  ", length(EXON1Ind$SID), "UCSC REFGENE region\n")
    cat("Gene body region contains: ", length(GENEBODYInd$SID), "UCSC REFGENE region\n")
    cat("3'UTR region contains:     ", length(UTR3Ind$SID), "UCSC REFGENE region\n")
    splitToRegionlist2 <- function (grepname = c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf")) {
        index <- col4 == grepname
        col1sub <- col1[index]
        col5sub <- col5[index]
        temp <- split(col1sub, col5sub)
        returnSID <- lapply(temp, unique)
        col0sub <- col0[index]
        temp <- split(col0sub, col5sub)
        returnPID <- lapply(temp, unique)
        return(Ind <- list(SID = returnSID, PID = returnPID))
    }
    ISLANDInd <- splitToRegionlist2(grepname = "Island")
    NSHOREInd <- splitToRegionlist2(grepname = "N_Shore")
    SSHOREInd <- splitToRegionlist2(grepname = "S_Shore")
    NSHELFInd <- splitToRegionlist2(grepname = "N_Shelf")
    SSHELFInd <- splitToRegionlist2(grepname = "S_Shelf")
    cat("Island region contains:    ", length(ISLANDInd$SID), "UCSC CPG ISLAND region\n")
    cat("N_Shore region contains:   ", length(NSHOREInd$SID), "UCSC CPG ISLAND region\n")
    cat("S_Shore region contains:   ", length(SSHOREInd$SID), "UCSC CPG ISLAND region\n")
    cat("N_Shelf region contains:   ", length(NSHELFInd$SID), "UCSC CPG ISLAND region\n")
    cat("S_Shelf region contains:   ", length(SSHELFInd$SID), "UCSC CPG ISLAND region\n")
    # setClass("methy450batch", representation(bmatrix = "matrix",
        # annot = "matrix", detectP = "matrix", groupinfo = "data.frame",
        # TSS1500Ind = "list", TSS200Ind = "list", UTR5Ind = "list",
        # EXON1Ind = "list", GENEBODYInd = "list", UTR3Ind = "list",
        # ISLANDInd = "list", NSHOREInd = "list", SSHOREInd = "list",
        # NSHELFInd = "list", SSHELFInd = "list"), where = topenv(parent.frame()))
    x.methy450 <- new("methy450batch", bmatrix = as.matrix(bmatrix),
        annot = as.matrix(annotation), detectP = as.matrix(detect_p),
        groupinfo = groupinfo, TSS1500Ind = TSS1500Ind, TSS200Ind = TSS200Ind,
        UTR5Ind = UTR5Ind, EXON1Ind = EXON1Ind, GENEBODYInd = GENEBODYInd,
        UTR3Ind = UTR3Ind, ISLANDInd = ISLANDInd, NSHOREInd = NSHOREInd,
        SSHOREInd = SSHOREInd, NSHELFInd = NSHELFInd, SSHELFInd = SSHELFInd)
    cat("\nA methy450batch class is created and the slotNames are:\n", slotNames(x.methy450), "\n")

    # setMethod("show", "methy450batch",
        # function(object){
            # cat("*** Class methy450batch , method Show ***\n")
            # cat("* BetaMatrix (limited to matrix 5x5) = ", paste("(", paste(dim(object@bmatrix), collapse = "x"), ")", sep = ""), "\n")
            # nrowShow <- min(5 , nrow(object@bmatrix))
            # ncolShow <- min(5 , ncol(object@bmatrix))
            # if (nrow(object@bmatrix)!=0) {
                # print(formatC(object@bmatrix[1:nrowShow, 1:ncolShow]), quote = FALSE)
            # } else {}
            # cat("* ..... .....\n\n")

            # cat("* Dectect P-Value (limited to matrix 5x5) = ", paste("(", paste(dim(object@detectP), collapse = "x"), ")", sep = ""), "\n")
            # nrowShow <- min(5 , nrow(object@detectP))
            # ncolShow <- min(5 , ncol(object@detectP))
            # if (nrow(object@detectP)!=0) {
                # print(formatC(object@detectP[1:nrowShow, 1:ncolShow]), quote = FALSE)
            # } else {}
            # cat("* ..... .....\n\n")

            # cat("* Annotation (limited to matrix 5x5) = ", paste("(", paste(dim(object@annot), collapse = "x"), ")", sep = ""), "\n")
            # nrowShow <- min(5 , nrow(object@annot))
            # ncolShow <- min(5 , ncol(object@annot))
            # if (nrow(object@annot)!=0) {
                # print(formatC(object@annot[1:nrowShow, 1:ncolShow]), quote = FALSE)
            # } else {}
            # cat("* ..... .....\n\n")

            # cat("* Group Information (limited to matrix 5x5) = ", paste("(", paste(dim(object@groupinfo), collapse = "x"), ")", sep = ""), "\n")
            # nrowShow <- min(5 , nrow(object@groupinfo))
            # ncolShow <- min(5 , ncol(object@groupinfo))
            # if (nrow(object@groupinfo)!=0) {
                # print(formatC(as.matrix(object@groupinfo[1:nrowShow, 1:ncolShow])), quote = FALSE)
            # } else {}
            # cat("* ..... .....\n\n")

            # Regions <- c("TSS1500", "TSS200", "5\'UTR", "1stExon", "Gene Body", "3\'UTR", "Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf")
            # regions <- c("TSS1500Ind", "TSS200Ind", "UTR5Ind", "EXON1Ind", "GENEBODYInd", "UTR3Ind", "ISLANDInd", "NSHOREInd", "SSHOREInd", "NSHELFInd", "SSHELFInd")
            # for (iReg in 1:length(regions)) {
                # cat("* Regions", Regions[iReg],"(limited to 1st element) = ", paste("(", paste(length(eval(parse(text = paste("object@", regions[iReg], sep = "")))$SID), collapse = "x"), ")", sep = ""), "\n")
                # if (length(eval(parse(text = paste("object@", regions[iReg], sep = "")))[[1]])!=0) {
                    # cat("$SID"); print(eval(parse(text = paste("object@", regions[iReg], sep = "")))$SID[1], quote = FALSE)
                    # cat("$PID"); print(eval(parse(text = paste("object@", regions[iReg], sep = "")))$PID[1], quote = FALSE)
                # } else {}
                # cat("* ..... .....\n\n")
            # }

            # cat("******* End Show (methy450batch) *******\n")
        # }
    # )
    return(x.methy450)
}
