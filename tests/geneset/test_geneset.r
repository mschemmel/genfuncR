test.helpers <- function() {
  expect_equal(last(c(1:3)), 3)
  expect_equal(first(c(1:3)), 1)
  expect_equal(dropLast(c(1:3)), c(1, 2))
  expect_true(inRange(3, 0, 10))
  expect_true(inRange(3, 3, 10))
  expect_false(inRange(3, 4, 10))
  expect_true(inRange(10, 3, 10))
  expect_false(inRange(11, 3, 10))

}
test.helpers()

test.layout <- function() {
  expect_equal(getAnnoYScale(numeric(0)), c(0, 0.5, 1))
  expect_equal(getAnnoYScale(2), pretty(c(0, 2)))
  expect_equal(getAnnoYScale(c(1:10), c(2:5)), pretty(c(2:5)))

  expect_equal(getAnnoYBreaks(numeric(0)), c(0, 0.5, 1))
  expect_equal(getAnnoYBreaks(c(1:3)), c(0, 0.5, 1))
  expect_equal(getAnnoYBreaks(c(1:2)), c(0, 0.5, 1))

  expect_equal(relativePosition(100, 50, 150), 0.5)
}
test.layout()

test.data <- function() {
  expect_warning(checkChromosomes(c("Chr1", "Chr2", "Chr3")))
  expect_equal(checkChromosomes(c("Chr1", "Chr2", "Chr3")), "Chr1")
  expect_equal(checkChromosomes(c("Chr1")), "Chr1")
}
test.data()
