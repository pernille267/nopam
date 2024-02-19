# Load packages
library(data.table)
library(nopam.smoothing)


# Load test data
alkaline <- fread(file = "~/nopam_warnings/alkaline-phosphatase-from-2021-01-01-to-2023-01-01-raw.csv")
creatinine <- fread(file = "~/nopam_warnings/creatinine-from-2021-01-01-to-2023-01-01-raw.csv")
glucose <- fread(file = "~/nopam_warnings/glucose-from-2021-01-01-to-2023-01-01-raw.csv")
potassium <- fread(file = "~/nopam_warnings/potassium-from-2021-01-01-to-2023-01-01-raw.csv")
tsh <- fread(file = "~/nopam_warnings/tsh-from-2021-01-01-to-2023-01-01-raw.csv")

# Fill missing dates with NA
#-------------------------------------------------------------------------------
id_cols <- c("Analyte Name", "Country", "Manufacturer Name", "Instrument Model Name", "Laboratory Code", "Instrument Id", "Instrument Code")
date_col <- "Measured At"
measurement_cols <- c("Median", "Hyper Percentage", "Hypo Percentage")

alkaline_list <- split(x = alkaline, by = id_cols, keep.by = TRUE, sorted = FALSE)
creatinine_list <- split(x = creatinine, by = id_cols, keep.by = TRUE, sorted = FALSE)
glucose_list <- split(x = glucose, by = id_cols, keep.by = TRUE, sorted = FALSE)
potassium_list <- split(x = potassium, by = id_cols, keep.by = TRUE, sorted = FALSE)
tsh_list <- split(x = tsh, by = id_cols, keep.by = TRUE, sorted = FALSE)

expand_group_data <- function(data, id_cols) {
  # Generate date sequence for this group
  all_dates <- seq(min(data[[date_col]]), max(data[[date_col]]), by = "1 day")

  # Number of rows for each id combination
  n_rows <- length(all_dates)

  # Create expanded data for this group
  expanded_group_data <- CJ("Measured At" = all_dates)

  for (id in id_cols) {
    expanded_group_data[[id]] <- rep(data[[id]][1], n_rows)
  }

  # Merge with original data
  out <- merge(expanded_group_data, data, by = c("Measured At", id_cols), all.x = TRUE)
  out$`Is Weekend` <- ifelse(weekdays(out$`Measured At`) == "Saturday" | weekdays(out$`Measured At`) == "Sunday", TRUE, FALSE)
  out$`Weekday` <- weekdays(out$`Measured At`)
  return(out)
}

alkaline <- lapply(X = alkaline_list, FUN = function(x) expand_group_data(x, id_cols = id_cols)) |> rbindlist(idcol = NULL)
creatinine <- lapply(X = creatinine_list, FUN = function(x) expand_group_data(x, id_cols = id_cols)) |> rbindlist(idcol = NULL)
glucose <- lapply(X = glucose_list, FUN = function(x) expand_group_data(x, id_cols = id_cols)) |> rbindlist(idcol = NULL)
potassium <- lapply(X = potassium_list, FUN = function(x) expand_group_data(x, id_cols = id_cols)) |> rbindlist(idcol = NULL)
tsh <- lapply(X = tsh_list, FUN = function(x) expand_group_data(x, id_cols = id_cols)) |> rbindlist(idcol = NULL)
#-------------------------------------------------------------------------------

test_data_1 <- creatinine[`Laboratory Code` == "LOW_2" & `Instrument Code` == "VITROS_2"]
bias_warning_1 <- xwarning_process0(data = test_data_1, pg_data = creatinine, from = "2022-01-01", to = "2022-12-31", measure = "median", warning = "bias", tol = 5, bw2 = 31, bw3 = 366)

test_that("Correctness", {
  expect_equal(bias_warning_1[`Measured At` == "2022-02-16"]$`Monthly Number of Results`, 27)
  expect_equal(bias_warning_1[`Measured At` == "2022-02-16"]$`Yearly Number of Results`, 352)
  expect_equal(bias_warning_1[`Measured At` == "2022-02-16"]$`Monthly Median`, 73.373411, tolerance = 1e-3)
  expect_equal(bias_warning_1[`Measured At` == "2022-02-16"]$`Monthly Start Date`, as.IDate("2022-01-17"))
})
