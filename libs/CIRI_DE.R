#!/usr/bin/env Rscript
rm(list=ls())
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(optparse))

# Parse parameters
parser <- OptionParser()
parser <- add_option(parser, c("--lib"), action="store", default=NA, type='character',
                     help="library csv matrix")
parser <- add_option(parser, c("--bsj"), action="store", default=NA, type='character',
                     help="bsj csv matrix")
# parser <- add_option(parser, c("--circ"), action="store", default=NA, type='character',
#                      help="circ info matrix")
parser <- add_option(parser, c("--gene"), action="store", default=NA, type='character',
                     help="gene count matrix")
parser <- add_option(parser, c("--out"), action="store", default=NA, type='character',
                     help="output file")

# opt <- parse_args(parser, args = c("--lib=AD_lib.csv",
#                                    "--bsj=AD_bsj.csv",
#                                    "--gene=AD_gene_count_matrix.csv",
#                                    "--out=AD_de.csv"))
opt <- parse_args(parser)

# main point of program is here, do this whether or not "verbose" is set
if (is.na(opt$lib) || is.na(opt$bsj) || is.na(opt$gene) || is.na(opt$out)) {
    cat("Please specify --lib/--bsj/--gene/--out, refer to the manual for detailed instruction!\n",
        file=stderr())
    quit()
}

# Load data
lib_mtx <- read.csv(opt$lib, row.names = 1)
gene_mtx <- read.csv(opt$gene, row.names = 1)
gene_mtx <- gene_mtx[,rownames(lib_mtx)]
bsj_mtx <- read.csv(opt$bsj, row.names = 1)

gene_DGE <- DGEList(counts = gene_mtx, group = lib_mtx$Group)
gene_idx <- filterByExpr(gene_DGE)
gene_DGE <- gene_DGE[gene_idx, keep.lib.sizes=FALSE]
gene_DGE <- calcNormFactors(gene_DGE)


if ("Subject" %in% colnames(lib_mtx)) {
    subject <- factor(lib_mtx$Subject)
    treat <- factor(lib_mtx$Group, levels=c("C", "T"))

    design <- model.matrix(~subject + treat)
} else {
    treat <- factor(lib_mtx$Group, levels=c("C", "T"))
    design <- model.matrix(~treat)
}

# design <- model.matrix(~factor(lib_mtx$Group))
# gene_DGE <- estimateDisp(gene_DGE, design, robust = TRUE)
# gene_fit <- glmFit(gene_DGE, design)
# gene_lrt <- glmLRT(gene_fit)
#
# gene_df <- gene_lrt$table
# gene_order <- order(gene_lrt$table$PValue)
# gene_df$DE <- decideTestsDGE(gene_lrt)
# gene_df <- gene_df[gene_order, ]
# gene_df$FDR <- p.adjust(gene_df$PValue, method="fdr")

circ_DGE <- DGEList(counts = bsj_mtx,
                    group = lib_mtx$Group,
                    lib.size = gene_DGE$samples[, "lib.size"],
                    norm.factors = gene_DGE$samples[, "norm.factors"])

# circ_idx <- filterByExpr(circ_DGE, min.count=)
# circ_DGE <- circ_DGE[circ_idx, , keep.lib.sizes=TRUE]
# head(circ_df[rowSums(bsj_mtx >= 2) >= nrow(lib_mtx) / 2,])

circ_DGE <- estimateDisp(circ_DGE, design, robust = TRUE)
circ_fit <- glmFit(circ_DGE, design)
circ_lrt <- glmLRT(circ_fit)

circ_df <- circ_lrt$table
circ_order <- order(circ_lrt$table$PValue)
circ_df$DE <- decideTestsDGE(circ_lrt)
circ_df <- circ_df[circ_order, ]
circ_df$FDR <- p.adjust(circ_df$PValue, method="fdr")
write.csv(circ_df, file=opt$out, quote = FALSE)
