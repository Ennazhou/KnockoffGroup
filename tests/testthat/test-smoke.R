test_that("example data creates trio and population knockoffs", {
  data(KnockoffGroup.example)
  idx_trio <- seq_len(30)
  idx_hap <- seq_len(60)
  idx_pop <- seq_len(30)

  set.seed(1)
  dat_ko_trio <- create_knockoff(
    KnockoffGroup.example$window_matrix,
    KnockoffGroup.example$gene.clusters,
    KnockoffGroup.example$dat.hap[idx_hap, , drop = FALSE],
    KnockoffGroup.example$pos,
    design = "trio",
    M = 2
  )

  set.seed(1)
  dat_ko_pop <- create_knockoff(
    KnockoffGroup.example$window_matrix,
    KnockoffGroup.example$gene.clusters,
    KnockoffGroup.example$dat.pop[idx_pop, , drop = FALSE],
    KnockoffGroup.example$pos,
    design = "population",
    M = 2
  )

  expect_equal(dim(dat_ko_trio), c(length(idx_trio), ncol(KnockoffGroup.example$dat), 2))
  expect_equal(dim(dat_ko_pop), c(length(idx_pop), ncol(KnockoffGroup.example$dat), 2))
})

test_that("example data runs exported analysis entry points", {
  data(KnockoffGroup.example)
  idx_dat <- seq_len(30)
  idx_hap <- seq_len(60)

  set.seed(2)
  dat_ko_trio <- create_knockoff(
    KnockoffGroup.example$window_matrix,
    KnockoffGroup.example$gene.clusters,
    KnockoffGroup.example$dat.hap[idx_hap, , drop = FALSE],
    KnockoffGroup.example$pos,
    design = "trio",
    M = 2
  )

  trio_result <- KnockoffGroup_trio(
    KnockoffGroup.example$geneinfo,
    KnockoffGroup.example$dat[idx_dat, , drop = FALSE],
    dat_ko_trio,
    KnockoffGroup.example$pos,
    chr = "1"
  )

  set.seed(3)
  y <- rnorm(length(idx_dat))
  prelim <- KS.prelim(y, out_type = "C")
  dat_ko_pop <- create_knockoff(
    KnockoffGroup.example$window_matrix,
    KnockoffGroup.example$gene.clusters,
    KnockoffGroup.example$dat.pop[idx_dat, , drop = FALSE],
    KnockoffGroup.example$pos,
    design = "population",
    M = 2
  )

  pop_result <- KS.chr_ko(
    prelim,
    chr = 1,
    dat_ko = dat_ko_pop,
    dat = KnockoffGroup.example$dat.pop[idx_dat, , drop = FALSE],
    pos = KnockoffGroup.example$pos,
    window.bed = KnockoffGroup.example$window.bed,
    M = 2,
    thres.single = 0,
    thres.ultrarare = 0
  )

  expect_s3_class(trio_result, "data.frame")
  expect_true(all(c("result.single", "result.window") %in% names(pop_result)))
})
