test.helpers <- function() {
  expect_equal(last(c(1:3)), 3)
  expect_equal(first(c(1:3)), 2)
  expect_equal(dropLast(c(1:3)), c(1,2))
}
test.helpers()
