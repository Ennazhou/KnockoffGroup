test_that("dispersion probabilities match the manuscript Davies calculation", {
  probability <- KnockoffGroup:::Get.p.SKAT(
    c(3, 2), diag(c(1, 2)), matrix(c(1, 1), ncol = 1), c(1, 1), NULL)
  expect_equal(as.numeric(probability), 0.01673880723180532, tolerance = 1e-12)
})

test_that("multiple trio genes use their own variant columns", {
  data(KnockoffGroup.example)
  x <- KnockoffGroup.example
  set.seed(15)
  ko <- create_knockoff(x$window_matrix, x$gene.clusters, x$dat.hap, x$pos, M = 2)
  together <- KnockoffGroup_trio(x$geneinfo, x$dat, ko, x$pos)
  separately <- do.call(rbind, lapply(seq_len(nrow(x$geneinfo)), function(gene) {
    columns <- which(x$gene.clusters == gene)
    KnockoffGroup_trio(x$geneinfo[gene, , drop = FALSE], x$dat[, columns, drop = FALSE],
      ko[, columns, , drop = FALSE], x$pos[columns])
  }))
  expect_equal(together, separately)
})

test_that("q-values handle tied scores and preserve input order", {
  expect_equal(knockoff_qvalues(c(0, 0, 1), c(3, 2, 1), 10), c(.05, .05, 1))
  expect_equal(knockoff_qvalues(c(0, 1), c(3, 3), 1), c(1, 1))
  expect_equal(knockoff_qvalues(c(0, NA, 1), c(0, 2, 1), 10), c(1, 1, 1))
})

test_that("population analysis retains dimensions for one SNP or one gene", {
  data(KnockoffGroup.example)
  x <- KnockoffGroup.example
  for (nvariant in c(1L, 4L)) {
    columns <- seq_len(nvariant)
    dat <- x$dat.pop[, columns, drop = FALSE]
    for (M in c(1L, 2L)) {
      set.seed(100)
      ko <- create_knockoff(x$window_matrix[1, , drop = FALSE],
        x$gene.clusters[columns], dat, x$pos[columns], design = "population", M = M)
      for (trait in c("C", "D")) {
        null <- KS.prelim(if (trait == "C") x$y.c else x$y.b, out_type = trait)
        result <- KS.chr_ko(null, 1, ko, dat, x$pos[columns],
          as.matrix(x$window.bed[1, , drop = FALSE]), M = M)
        expect_equal(nrow(result$result.window), 1L)
        expect_true(all(is.finite(result$result.window[, c("kappa", "tau", "P_KS")])))
      }
    }
  }
})
