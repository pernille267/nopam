library(data.table)
library(nopam.smoothing)

test_data_1 <- fread(file = "~/nopam_warnings/lc_ll_demo/test_data_1.csv")[`Measured At` <= "2022-01-05"]
dates <- as.numeric(test_data_1$`Measured At` - min(test_data_1$`Measured At`, na.rm = TRUE)) + 1
medians <- test_data_1$Median

laks(date = dates, median = medians, bandwidth = 11, average = 83.85682, standard_deviation = 19.49721, approximate = TRUE, tol = 0.01)
