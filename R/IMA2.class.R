setClass("methy450batch",
    representation(
        bmatrix = "matrix",
        annot = "matrix",
        detectP = "matrix",
        groupinfo = "data.frame",
        TSS1500Ind = "list",
        TSS200Ind = "list",
        UTR5Ind = "list",
        EXON1Ind = "list",
        GENEBODYInd = "list",
        UTR3Ind = "list",
        ISLANDInd = "list",
        NSHOREInd = "list",
        SSHOREInd = "list",
        NSHELFInd = "list",
        SSHELFInd = "list"
    ),
    where = topenv(parent.frame())
)

setMethod("show", "methy450batch",
    function(object){
        cat("*** Class methy450batch , method Show ***\n")
        cat("* BetaMatrix (limited to matrix 5x5) = ", paste("(", paste(dim(object@bmatrix), collapse = "x"), ")", sep = ""), "\n")
        nrowShow <- min(5 , nrow(object@bmatrix))
        ncolShow <- min(5 , ncol(object@bmatrix))
        if (nrow(object@bmatrix)!=0) {
            print(formatC(object@bmatrix[1:nrowShow, 1:ncolShow]), quote = FALSE)
        } else {}
        cat("* ..... .....\n\n")

        cat("* Dectect P-Value (limited to matrix 5x5) = ", paste("(", paste(dim(object@detectP), collapse = "x"), ")", sep = ""), "\n")
        nrowShow <- min(5 , nrow(object@detectP))
        ncolShow <- min(5 , ncol(object@detectP))
        if (nrow(object@detectP)!=0) {
            print(formatC(object@detectP[1:nrowShow, 1:ncolShow]), quote = FALSE)
        } else {}
        cat("* ..... .....\n\n")

        cat("* Annotation (limited to matrix 5x5) = ", paste("(", paste(dim(object@annot), collapse = "x"), ")", sep = ""), "\n")
        nrowShow <- min(5 , nrow(object@annot))
        ncolShow <- min(5 , ncol(object@annot))
        if (nrow(object@annot)!=0) {
            print(formatC(object@annot[1:nrowShow, 1:ncolShow]), quote = FALSE)
        } else {}
        cat("* ..... .....\n\n")

        cat("* Group Information (limited to matrix 5x5) = ", paste("(", paste(dim(object@groupinfo), collapse = "x"), ")", sep = ""), "\n")
        nrowShow <- min(5 , nrow(object@groupinfo))
        ncolShow <- min(5 , ncol(object@groupinfo))
        if (nrow(object@groupinfo)!=0) {
            print(formatC(as.matrix(object@groupinfo[1:nrowShow, 1:ncolShow])), quote = FALSE)
        } else {}
        cat("* ..... .....\n\n")

        Regions <- c("TSS1500", "TSS200", "5\'UTR", "1stExon", "Gene Body", "3\'UTR", "Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf")
        regions <- c("TSS1500Ind", "TSS200Ind", "UTR5Ind", "EXON1Ind", "GENEBODYInd", "UTR3Ind", "ISLANDInd", "NSHOREInd", "SSHOREInd", "NSHELFInd", "SSHELFInd")
        for (iReg in 1:length(regions)) {
            cat("* Regions", Regions[iReg],"(limited to 1st element) = ", paste("(", paste(length(eval(parse(text = paste("object@", regions[iReg], sep = "")))$SID), collapse = "x"), ")", sep = ""), "\n")
            if (length(eval(parse(text = paste("object@", regions[iReg], sep = "")))[[1]])!=0) {
                cat("$SID"); print(eval(parse(text = paste("object@", regions[iReg], sep = "")))$SID[1], quote = FALSE)
                cat("$PID"); print(eval(parse(text = paste("object@", regions[iReg], sep = "")))$PID[1], quote = FALSE)
            } else {}
            cat("* ..... .....\n\n")
        }

        cat("******* End Show (methy450batch) *******\n")
    }
)

setClass("exprmethy450",
    representation(
        bmatrix = "matrix",
        annot = "matrix",
        detectP = "matrix",
        groupinfo = "data.frame"
    ),
    where = topenv(parent.frame())
)

setMethod("show", "exprmethy450",
    function(object){
        cat("*** Class exprmethy450 , method Show ***\n")
        cat("* BetaMatrix (limited to matrix 5x5) = ", paste("(", paste(dim(object@bmatrix), collapse = "x"), ")", sep = ""), "\n")
        nrowShow <- min(5 , nrow(object@bmatrix))
        ncolShow <- min(5 , ncol(object@bmatrix))
        if (nrow(object@bmatrix)!=0) {
            print(formatC(object@bmatrix[1:nrowShow, 1:ncolShow]), quote = FALSE)
        } else {}
        cat("* ..... .....\n\n")

        cat("* Detect P-Value (limited to matrix 5x5) = ", paste("(", paste(dim(object@detectP), collapse = "x"), ")", sep = ""), "\n")
        nrowShow <- min(5 , nrow(object@detectP))
        ncolShow <- min(5 , ncol(object@detectP))
        if (nrow(object@detectP)!=0) {
            print(formatC(object@detectP[1:nrowShow, 1:ncolShow]), quote = FALSE)
        } else {}
        cat("* ..... .....\n\n")

        cat("* Annotation (limited to matrix 5x5) = ", paste("(", paste(dim(object@annot), collapse = "x"), ")", sep = ""), "\n")
        nrowShow <- min(5 , nrow(object@annot))
        ncolShow <- min(5 , ncol(object@annot))
        if (nrow(object@annot)!=0) {
            print(formatC(object@annot[1:nrowShow, 1:ncolShow]), quote = FALSE)
        } else {}
        cat("* ..... .....\n\n")

        cat("* Group Information (limited to matrix 5x5) = ", paste("(", paste(dim(object@groupinfo), collapse = "x"), ")", sep = ""), "\n")
        nrowShow <- min(5 , nrow(object@groupinfo))
        ncolShow <- min(5 , ncol(object@groupinfo))
        if (nrow(object@groupinfo)!=0) {
            print(formatC(as.matrix(object@groupinfo[1:nrowShow, 1:ncolShow])), quote = FALSE)
        } else {}
        cat("* ..... .....\n\n")

        cat("******* End Show (exprmethy450) *******\n")
    }
)