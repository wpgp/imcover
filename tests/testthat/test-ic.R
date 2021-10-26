context("ic core")

# test data
df <- data.frame(region = c(1,1,1),
                 iso = c("A","B","C"),
                 time = c(2020, 2020, 2020),
                 vax = c('DTP1', 'DTP1', 'DTP1'),
                 coverage = c(.8, .7, .99),
                 source = c('admin', 'admin', 'admin'))


test_that("we can create ic from data.frame", {
  icdf <- as.ic_data(df)
  expect_identical(class(icdf), c("ic.df","data.frame"))
})

test_that("we can subset ic object", {
  icdf <- as.ic_data(df)

  expect_equal(icdf[[2]], df$iso)
  expect_equal(icdf[,2], icdf)
  expect_equal(icdf[,2, drop = T], df$iso)
  expect_equal(nrow(icdf[1,]), 1)
})

test_that("rbind works on ic object", {
  ic1 <- as.ic_data(df)
  ic2 <- rbind(ic1, ic1)

  expect_equal(nrow(ic2), 2*nrow(ic1))
  expect_equal(ncol(ic2), ncol(ic1))
  expect_equal(names(ic2), names(ic1))
})

test_that("merge works on ic object", {
  ic1 <- as.ic_data(df)
  ic2 <- ic1
  ic2$extra <- c(7, 8, 9)
  x <- merge(ic1, ic2)

  expect_equal(x, ic2)

  names(ic2)[2] <- 'country'
  x <- merge(ic1, ic2, sort = F)
  expect_equal(nrow(x), nrow(ic1))
  expect_equal(ncol(x), ncol(ic2))

  x <- merge(ic1, ic2, attr.x = FALSE, sort = F)
  expect_equal(x, ic2)
})
