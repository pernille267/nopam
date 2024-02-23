library(data.table)
library(nopam.smoothing)
library(np, quietly = TRUE)

# Load data
alkaline <- fread(file = "~/nopam_warnings/alkaline-phosphatase-from-2021-01-01-to-2023-01-01-raw.csv") |> expand_data()
test_data_1 <- alkaline[`Instrument Code` == "A2-B" & `Laboratory Code` == "MID_2" & `Measured At` >= "2021-11-30" & `Measured At` <= "2022-01-08"]
dates <- as.numeric(test_data_1$`Measured At` - min(test_data_1$`Measured At`, na.rm = TRUE)) + 1
medians <- test_data_1$Median
normalized_medians <- (test_data_1$Median - 81.3302) / 1.8043

actual_lc_1 <- laks(date = dates, median = medians, bandwidth = 11, average = 81.3302, standard_deviation = 1.8043, approximate = FALSE, tol = 0) |> as.data.table()
actual_lc_2 <- laks(date = dates, median = medians, bandwidth = 11, average = 81.3302, standard_deviation = 1.8043, approximate = FALSE, tol = 0.01) |> as.data.table()
actual_lc_3 <- laks(date = dates, median = medians, bandwidth = 11, average = 81.3302, standard_deviation = 1.8043, approximate = TRUE, tol = 0) |> as.data.table()
actual_lc_4 <- laks(date = dates, median = medians, bandwidth = 11, average = 81.3302, standard_deviation = 1.8043, approximate = TRUE, tol = 0.01) |> as.data.table()
actual_ll_1 <- llks(date = dates, median = medians, bandwidth = 11, average = 81.3302, standard_deviation = 1.8043, approximate = FALSE, tol = 0) |> as.data.table()
actual_ll_2 <- llks(date = dates, median = medians, bandwidth = 11, average = 81.3302, standard_deviation = 1.8043, approximate = FALSE, tol = 0.01) |> as.data.table()
actual_ll_3 <- llks(date = dates, median = medians, bandwidth = 11, average = 81.3302, standard_deviation = 1.8043, approximate = TRUE, tol = 0) |> as.data.table()
actual_ll_4 <- llks(date = dates, median = medians, bandwidth = 11, average = 81.3302, standard_deviation = 1.8043, approximate = TRUE, tol = 0.01) |> as.data.table()

expected_lc_11 <- npreg(bws = 11, txdat = dates[!is.na(medians)], tydat = medians[!is.na(medians)], regtype = "lc", gradients = TRUE)
expected_lc_12 <- npreg(bws = 11, txdat = dates[!is.na(medians)], tydat = normalized_medians[!is.na(medians)], regtype = "lc", gradients = TRUE)
expected_ll_11 <- npreg(bws = 11, txdat = dates[!is.na(medians)], tydat = medians[!is.na(medians)], regtype = "ll", gradients = TRUE)
expected_ll_12 <- npreg(bws = 11, txdat = dates[!is.na(medians)], tydat = normalized_medians[!is.na(medians)], regtype = "ll", gradients = TRUE)

expected_lc_1_and_2 <- data.table(date = dates[!is.na(medians)],
                                  smoothed_median = expected_lc_11$mean,
                                  smoothed_normalized_median = expected_lc_12$mean,
                                  gradients = atan(expected_lc_12$grad[,]) * 180 / pi)

expected_ll_1_and_2 <- data.table(date = dates[!is.na(medians)],
                                  smoothed_median = expected_ll_11$mean,
                                  smoothed_normalized_median = expected_ll_12$mean,
                                  gradients = atan(expected_ll_12$grad[,]) * 180 / pi)


test_that(desc = "Test direct corresponse with npreg", code = {
  expect_equal(object = actual_lc_1$smoothed_median, expected = expected_lc_1_and_2$smoothed_median, tolerance = 1e-6)
  expect_equal(object = actual_ll_1$smoothed_median, expected = expected_ll_1_and_2$smoothed_median, tolerance = 1e-6)
  expect_equal(object = actual_lc_1$gradients, expected = expected_lc_1_and_2$gradients, tolerance = 1e-4)
  expect_equal(object = actual_ll_1$gradients, expected = expected_ll_1_and_2$gradients, tolerance = 1e-6)
})

test_that(desc = "Test approximate corresponse with npreg", code = {
  expect_equal(object = actual_lc_2$smoothed_median, expected = expected_lc_1_and_2$smoothed_median, tolerance = 1e-3)
  expect_equal(object = actual_ll_2$smoothed_median, expected = expected_ll_1_and_2$smoothed_median, tolerance = 1e-3)
  expect_equal(object = actual_lc_2$gradients, expected = expected_lc_1_and_2$gradients, tolerance = 1e-1)
  expect_equal(object = actual_ll_2$gradients, expected = expected_ll_1_and_2$gradients, tolerance = 1e-1)
})

sd_gradients_lc <- atan(c(0, diff(expected_lc_1_and_2$smoothed_normalized_median, lag = 2) / 2, 0)) * 180 / pi
sd_gradients_ll <- atan(c(0, diff(expected_ll_1_and_2$smoothed_normalized_median, lag = 2) / 2, 0)) * 180 / pi

test_that(desc = "Test second differences slopes with npreg", code = {
  expect_equal(object = actual_lc_3$gradients, expected = sd_gradients_lc, tolerance = 1e-6)
  expect_equal(object = actual_ll_3$gradients, expected = sd_gradients_ll, tolerance = 1e-6)
  expect_equal(object = actual_lc_4$gradients, expected = sd_gradients_lc, tolerance = 1e-1)
  expect_equal(object = actual_ll_4$gradients, expected = sd_gradients_ll, tolerance = 1e-1)
})
