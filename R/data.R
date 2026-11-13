#' Synthetic genotype and haplotype example
#'
#' Simulated inputs for population and trio analyses, generated with seed
#' 20261112 by `data-raw/generate_example.R`. All coordinates and genes are
#' artificial. The population individuals are sampled independently of the trios.
#'
#' @format A named list containing
#' \describe{
#'   \item{dat}{A 180 by 12 genotype matrix for 60 trios. Rows within each trio are father, mother, offspring.}
#'   \item{dat.hap}{A 360 by 12 binary haplotype matrix. Each trio has six rows, two per individual. The offspring's first and second haplotypes are paternal and maternal.}
#'   \item{dat.pop}{A 180 by 12 genotype matrix for unrelated population individuals.}
#'   \item{pos}{A length-12 vector of artificial variant positions.}
#'   \item{geneinfo}{Three gene records with start, end, nsnp, and gene_id.}
#'   \item{window_matrix}{Gene boundaries with start, end, and gene_id.}
#'   \item{window.bed}{The same boundaries with chr, start, and end. Ends are exclusive for the population association function.}
#'   \item{gene.clusters}{A length-12 integer vector assigning four variants to each gene.}
#'   \item{y.b}{A length-180 binary outcome for dat.pop.}
#'   \item{y.c}{A length-180 continuous outcome for dat.pop.}
#' }
#' @source Generated entirely from random draws. The generation script is included in the repository.
"KnockoffGroup.example"
