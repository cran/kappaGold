# mkuhn, 2024-09-11
# test utils-functions from this package


test_that("Jackknife SE and bias", {
  
  # Expected results from boostrap::jackknife
  
  # mean
  mean_f <- function(x) mean(datasets::mtcars$mpg[x])
  expect_equal(mean_f(1:32), expected = 20.090625, tolerance = 1e-7)
  
  # jackknife of mean
  mean_jk <- victorinox(mean_f, 1:32)
  expect_type(mean_jk, type = "list")
  expect_named(mean_jk, expected = c("bias_j", "se_j"))
  expect_identical(mean_jk[["bias_j"]], 0)
  expect_equal(mean_jk[["se_j"]], 1.06542395937282, tolerance = 1e-7)
  
  
  # Pearson correlation
  cor_f <- function(x) cor(datasets::mtcars$mpg[x], datasets::mtcars$disp[x])
  expect_equal(cor_f(1:32), expected = -0.847551379262479, tolerance = 1e-9)
  
  # jackknife of Pearson correlation
  cor_jk <- victorinox(cor_f, 1:32)
  expect_type(cor_jk, type = "list")
  expect_named(cor_jk, expected = c("bias_j", "se_j"))
  expect_equal(cor_jk[["bias_j"]], -0.000945434421600888, tolerance = 1e-9)
  expect_equal(cor_jk[["se_j"]], 0.0395197382974496, tolerance = 1e-7)
  
})

