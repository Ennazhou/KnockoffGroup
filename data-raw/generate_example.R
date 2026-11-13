# Generate the packaged example entirely from random draws.
# Run from the package root with Rscript data-raw/generate_example.R.
set.seed(20261112)
ngenes <- 3L
snps_per_gene <- 4L
ntrios <- 60L
nindividuals <- 180L
nvariants <- ngenes * snps_per_gene

simulate_haplotypes <- function(n) {
    haplotypes <- matrix(0L, n, nvariants)
    for (gene in seq_len(ngenes)) {
        columns <- (gene - 1L) * snps_per_gene + seq_len(snps_per_gene)
        haplotypes[, columns[1]] <- rbinom(n, 1, 0.3)
        for (column in columns[-1]) {
            copy_previous <- rbinom(n, 1, 0.75)
            new_allele <- rbinom(n, 1, 0.3)
            haplotypes[, column] <- ifelse(copy_previous == 1L,
                haplotypes[, column - 1L], new_allele)
        }
    }
    haplotypes
}

founders <- simulate_haplotypes(4L * ntrios)
dat.hap <- matrix(0L, 6L * ntrios, nvariants)
for (family in seq_len(ntrios)) {
    parental_rows <- (family - 1L) * 4L + seq_len(4L)
    family_rows <- (family - 1L) * 6L + seq_len(6L)
    dat.hap[family_rows[1:4], ] <- founders[parental_rows, ]
    dat.hap[family_rows[5], ] <- founders[parental_rows[sample.int(2L, 1L)], ]
    dat.hap[family_rows[6], ] <- founders[parental_rows[2L + sample.int(2L, 1L)], ]
}
dat <- dat.hap[seq(1L, nrow(dat.hap), 2L), ] + dat.hap[seq(2L, nrow(dat.hap), 2L), ]
population_haplotypes <- simulate_haplotypes(2L * nindividuals)
dat.pop <- population_haplotypes[seq(1L, 2L * nindividuals, 2L), ] +
    population_haplotypes[seq(2L, 2L * nindividuals, 2L), ]
y.b <- rbinom(nindividuals, 1L, plogis(-0.8 + 0.7 * dat.pop[, 1]))
y.c <- 0.7 * dat.pop[, 1] + rnorm(nindividuals)
start <- seq(1000L, by = 1000L, length.out = ngenes)
end <- start + 301L
pos <- as.vector(t(outer(start, seq(0L, 300L, 100L), '+')))
gene.clusters <- rep(seq_len(ngenes), each = snps_per_gene)
geneinfo <- data.frame(start, end, nsnp = snps_per_gene, gene_id = seq_len(ngenes))
window_matrix <- geneinfo[, c('start', 'end', 'gene_id')]
window.bed <- data.frame(chr = 1L, start, end)
KnockoffGroup.example <- list(dat = dat, dat.hap = dat.hap, dat.pop = dat.pop,
    pos = pos, geneinfo = geneinfo, window_matrix = window_matrix,
    window.bed = window.bed, gene.clusters = gene.clusters, y.b = y.b, y.c = y.c)
stopifnot(all(dat %in% 0:2), all(dat.hap %in% 0:1), length(unique(y.b)) == 2L)
save(KnockoffGroup.example, file = 'data/KnockoffGroup.example.rda', compress = 'xz', version = 2)
