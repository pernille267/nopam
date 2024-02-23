#' Capitalize a character value
#'
#' @param x A \code{character} string to be capitalized.
#'
#' @return The input string, but where the first letter is uppercase.
#' @export
#'
#' @examples capitalize("norway")
capitalize <- function(x){
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
}

#' Get Date Breaks Using The Range of Dates
#'
#' @param date_range An \code{integer}.
#'
#' @return The date breaks
#' @export
#'
#' @examples get_date_tick_distance(39)
get_date_tick_distance <- function(date_range){
  if(date_range <= 20){
    return("4 days")
  }
  else if(date_range <= 40){
    return("1 week")
  }
  else if(date_range <= 80){
    return("2 weeks")
  }
  else if(date_range <= 160){
    return("1 month")
  }
  else if(date_range <= 320){
    return("2 months")
  }
  else if(date_range <= 480){
    return("3 months")
  }
  else if(date_range <= 640){
    return("4 months")
  }
  else{
    return("6 months")
  }
  return("6 months")
}

#' Convert From Default to Their Corresponding Values
#'
#' @param plotting_arguments A named \code{list} of plotting arguments
#' @param date_range An \code{integer}.
#' @param measure A \code{character} signifying which measure to smooth.
#' @param warning A \code{character} signifying which warning process that is of interest to simulate.
#'
#' @return A named \code{list} with the actual plotting arguments.
#' @export
#'
#' @examples get_default_plotting_arguments(list("warning_point_shape" = 3))
get_default_plotting_arguments <- function(plotting_arguments, date_range, measure = c("median", "hyper", "hypo"), warning = c("slope", "bias", "peer_group")){
  warning <- warning[1]
  measure <- measure[1]

  provided_plotting_arguments <- names(plotting_arguments)

  if(any("title" == provided_plotting_arguments)){
    if(plotting_arguments$title == "default" | is.na(plotting_arguments$title)){
      if(warning == "peer_group"){
        plotting_arguments$title <- paste0(capitalize(measure), " peer group warning process")
      }
      else{
        plotting_arguments$title <- paste0(capitalize(measure), " ", warning, " warning process")
      }

    }
  }
  else{
    plotting_arguments$title <- "default"
  }

  if(any("subtitle" == provided_plotting_arguments)){
    if(plotting_arguments$subtitle == "default" | is.na(plotting_arguments$subtitle)){
      if(warning == "slope"){
        plotting_arguments$subtitle <- paste0("Line segments below the smoothed ", measure, " value curve indicate warnings")
      }
      else if(warning == "bias"){
        plotting_arguments$subtitle <- paste0("Line segments inbetween monthly and yearly median curves of the ", measure, " values indicate warnings")
      }
      else if(warning == "peer_group"){
        plotting_arguments$subtitle <- paste0("Line segments inbetween monthly instrument and peer group median curves of the ", measure, " values indicate warnings")
      }
    }
  }
  else{
    plotting_arguments$subtitle <- "default"
  }

  if(any("xlab" == provided_plotting_arguments)){
    if(plotting_arguments$xlab == "default" | is.na(plotting_arguments$xlab)){
      plotting_arguments$xlab <- "Date (YYYY-MM-DD)"
    }
  }
  else{
    plotting_arguments$xlab <- "default"
  }

  if(any("ylab" == provided_plotting_arguments)){
    if(plotting_arguments$ylab == "default" | is.na(plotting_arguments$ylab)){
      if(warning == "slope"){
        plotting_arguments$ylab <- capitalize(measure)
      }
      else if(warning == "bias"){
        plotting_arguments$ylab <- paste0("Monthly and yearly median of ", measure, " values")
      }
      else if(warning == "peer_group"){
        plotting_arguments$ylab <- paste0("Monthly median of ", measure, " values")
      }
    }
  }
  else{
    plotting_arguments$ylab <- "default"
  }

  if(any("legend_show" == provided_plotting_arguments)){
    if(plotting_arguments$legend_show == "default" | is.na(plotting_arguments$legend_show)){
      plotting_arguments$legend_show <- FALSE
    }
  }
  else{
    plotting_arguments$legend_show <- "default"
  }

  if(any("legend_position" == provided_plotting_arguments)){
    if(plotting_arguments$legend_position == "default" | is.na(plotting_arguments$legend_position)){
      plotting_arguments$legend_position <- "bottom"
    }
  }
  else{
    plotting_arguments$legend_position <- "default"
  }

  if(any("date_tick_distance" == provided_plotting_arguments)){
    if(plotting_arguments$date_tick_distance == "default" | is.na(plotting_arguments$date_tick_distance)){
      plotting_arguments$date_tick_distance <- get_date_tick_distance(date_range)
    }
  }
  else{
    plotting_arguments$date_tick_distance <- "default"
  }

  if(any("unit" == provided_plotting_arguments)){
    if(plotting_arguments$unit == "default" | is.na(plotting_arguments$unit)){
      if(measure == "median"){
        plotting_arguments$unit <- ""
      }
      else if(any(measure %in% c("hyper", "hypo"))){
        plotting_arguments$unit <- "%"
      }
    }
  }
  else{
    plotting_arguments$unit <- "default"
  }

  if(any("raw_data_shape" == provided_plotting_arguments)){
    if(plotting_arguments$raw_data_shape == "default" | is.na(plotting_arguments$raw_data_shape)){
      plotting_arguments$raw_data_shape <- 3
    }
  }
  else{
    plotting_arguments$raw_data_shape <- "default"
  }

  if(any("raw_data_alpha" == provided_plotting_arguments)){
    if(plotting_arguments$raw_data_alpha == "default" | is.na(plotting_arguments$raw_data_alpha)){
      plotting_arguments$raw_data_alpha <- 0.3
    }
  }
  else{
    plotting_arguments$raw_data_alpha <- "default"
  }

  if(any("raw_data_size" == provided_plotting_arguments)){
    if(plotting_arguments$raw_data_size == "default" | is.na(plotting_arguments$raw_data_size)){
      plotting_arguments$raw_data_size <- 0.3
    }
  }
  else{
    plotting_arguments$raw_data_size <- "default"
  }

  if(any("smooth_data_linewidth" == provided_plotting_arguments)){
    if(plotting_arguments$smooth_data_linewidth == "default" | is.na(plotting_arguments$smooth_data_linewidth)){
      plotting_arguments$smooth_data_linewidth <- 0.5
    }
  }
  else{
    plotting_arguments$smooth_data_linewidth <- "default"
  }

  if(any("warning_segment_color" == provided_plotting_arguments)){
    if(plotting_arguments$warning_segment_color == "default" | is.na(plotting_arguments$warning_segment_color)){
      plotting_arguments$warning_segment_color <- "gray"
    }
  }
  else{
    plotting_arguments$warning_segment_color <- "default"
  }

  if(any("warning_segment_linewidth" == provided_plotting_arguments)){
    if(plotting_arguments$warning_segment_linewidth == "default" | is.na(plotting_arguments$warning_segment_linewidth)){
      plotting_arguments$warning_segment_linewidth <- 0.25
    }
  }
  else{
    plotting_arguments$warning_segment_linewidth <- "default"
  }

  if(any("warning_point_shape" == provided_plotting_arguments)){
    if(plotting_arguments$warning_point_shape == "default" | is.na(plotting_arguments$warning_point_shape)){
      plotting_arguments$warning_point_shape <- 4
    }
  }
  else{
    plotting_arguments$warning_point_shape <- "default"
  }

  if(any("warning_point_color" == provided_plotting_arguments)){
    if(plotting_arguments$warning_point_color == "default" | is.na(plotting_arguments$warning_point_color)){
      plotting_arguments$warning_point_color <- "black"
    }
  }
  else{
    plotting_arguments$warning_point_color <- "default"
  }

  if(any("warning_point_size" == provided_plotting_arguments)){
    if(plotting_arguments$warning_point_size == "default" | is.na(plotting_arguments$warning_point_size)){
      plotting_arguments$warning_point_size <- 0.5
    }
  }
  else{
    plotting_arguments$warning_point_size <- "default"
  }

  return(plotting_arguments)
}

#' Plot a Warning Process
#'
#' @param data A \code{data.table} containing necessary warning data.
#' @param from A \code{IDate} object that is the first date where a warning is calculated.
#' @param to A \code{IDate} object that is the last date where a warning is calculated
#' @param measure A \code{character} signifying which measure to smooth.
#' @param warning A \code{character} signifying which warning process that is of interest to simulate.
#' @param strip_percentage A \code{double} between 0 and 100. Data outside are stripped from the plot.
#' @param plot A \code{logical} value. Set to \code{TRUE} to plot the plot right away. Set to \code{FALSE} to store the plot.
#' @param plotting_arguments A \code{list} of optional plotting arguments.
#'
#' @return A \code{ggplot2} object.
#' @export
#'
#' @examples print(1)
warning_process_plot <- function(data, from, to, measure = c("median", "hyper", "hypo"), warning = c("slope", "bias", "peer_group"), strip_percentage = 1, plot = FALSE, plotting_arguments = list("title" = "default", "subtitle" = "default", "xlab" = "default", "ylab" = "default", "legend_show" = "default", "legend_position" = "default", "date_tick_distance" = "default", "unit" = "default", "raw_data_shape" = "default", "raw_data_alpha" = "default", "raw_data_size" = "default", "smooth_data_linewidth" = "default", "warning_segment_color" = "default", "warning_segment_linewidth" = "default", "warning_point_shape" = "default", "warning_point_color" = "default", "warning_point_size" = "default")){

  m <- `Is Warning` <- `Measured At` <- `Instrument Code` <- `sm` <- `Monthly Median` <- `Yearly Median` <- `Peer Group Monthly Median` <- NULL
  warning <- warning[1]
  measure <- measure[1]
  date_range <- diff(range(data$`Measured At`, na.rm = TRUE))
  plotting_arguments <- get_default_plotting_arguments(plotting_arguments, date_range = date_range, measure = measure, warning = warning)
  plotting_arguments <- get_default_plotting_arguments(plotting_arguments, date_range = date_range, measure = measure, warning = warning)
  n_instruments <- length(unique(data$`Instrument Code`))
  plot_grid_nrow <- c(1, 2, 3, 2, 3, 3, 3, 3, 3, 4, 4, 4)
  plot_grid_ncol <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3)
  if(n_instruments < 0){
    stop("error")
  }
  else if(n_instruments <= 12){
    plot_grid_nrow <- plot_grid_nrow[n_instruments]
    plot_grid_ncol <- plot_grid_ncol[n_instruments]
  }
  else if(n_instruments > 12){
    plot_grid_nrow <- 4
  }
  if(warning == "slope"){

    if(measure == "median"){
      data$`m` <- data$`Median`
      data$`sm` <- data$`Smoothed Median`
    }
    else if(measure == "hyper"){
      data$`m` <- data$`Hyper Percentage`
      data$`sm` <- data$`Smoothed Hyper Percentage`
    }
    else if(measure == "hypo"){
      data$`m` <- data$`Hypo Percentage`
      data$`sm` <- data$`Smoothed Hypo Percentage`
    }

    upper_and_lower <- quantile(x = data$m, probs = c(strip_percentage / 200, 1 - strip_percentage / 200), na.rm = TRUE, names = FALSE)

    plot_object <- ggplot() +
      geom_vline(xintercept = data[is.na(m)][`Is Warning` == TRUE]$`Measured At`, color = plotting_arguments$warning_segment_color, linewidth = plotting_arguments$warning_segment_linewidth, alpha = 0.5) +
      geom_vline(xintercept = as.IDate(from)) +
      geom_point(data = data, mapping = aes(x = `Measured At`, y = m, color = `Instrument Code`), shape = plotting_arguments$raw_data_shape, alpha = plotting_arguments$raw_data_alpha, size = plotting_arguments$raw_data_size, show.legend = plotting_arguments$legend_show, na.rm = TRUE) +
      geom_line(data = data[!is.na(m)], mapping = aes(x = `Measured At`, y = sm, color = `Instrument Code`), alpha = 0.5, linewidth = plotting_arguments$smooth_data_linewidth / 2, show.legend = plotting_arguments$legend_show, na.rm = TRUE) +
      geom_line(data = data, mapping = aes(x = `Measured At`, y = sm, color = `Instrument Code`), linewidth = plotting_arguments$smooth_data_linewidth, show.legend = plotting_arguments$legend_show, na.rm = TRUE) +
      geom_segment(data = data, mapping = aes(x = `Measured At`, xend = `Measured At`, y = upper_and_lower[1], yend = sm, alpha = `Is Warning`), color = plotting_arguments$warning_segment_color, linewidth = plotting_arguments$warning_segment_linewidth, show.legend = plotting_arguments$legend_show, na.rm = TRUE) +
      geom_point(data = data, mapping = aes(x = `Measured At`, y = sm, alpha = `Is Warning`), shape = plotting_arguments$warning_point_shape, color = plotting_arguments$warning_point_color, size = plotting_arguments$warning_point_size, show.legend = plotting_arguments$legend_show, na.rm = TRUE) +
      facet_wrap(facets = . ~ `Instrument Code`, nrow = plot_grid_nrow, ncol = plot_grid_ncol, scales = "free_x") +
      labs(title = plotting_arguments$title, subtitle = plotting_arguments$subtitle) +
      scale_x_date(name = plotting_arguments$xlab, breaks = plotting_arguments$date_tick_distance, limits = c(as.IDate(from) - 38, as.IDate(to))) +
      scale_y_continuous(name = plotting_arguments$ylab, n.breaks = 10, limits = upper_and_lower) +
      scale_alpha_manual(values = c("FALSE" = 0, "TRUE" = 1)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 20, vjust = 0.9),
            axis.ticks.length.x = unit(x = 3, units = "mm"),
            legend.position = plotting_arguments$legend_position,
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_rect(fill = "black"),
            strip.text = element_text(face = "bold", color = "white"))

    if(plot){
      plot(plot_object)
    }
    else{
      return(plot_object)
    }
  }
  if(warning == "bias"){
    plot_object <- ggplot(data = data) +
      geom_vline(xintercept = as.IDate(from)) +
      geom_segment(mapping = aes(x = `Measured At`, xend = `Measured At`, y = `Monthly Median`, yend = `Yearly Median`, alpha = `Is Warning`), color = plotting_arguments$warning_segment_color, linewidth = plotting_arguments$warning_segment_linewidth, show.legend = plotting_arguments$legend_show) +
      geom_step(mapping = aes(x = `Measured At`, y = `Monthly Median`, color = `Instrument Code`), linewidth = plotting_arguments$smooth_data_linewidth, show.legend = plotting_arguments$legend_show) +
      geom_step(mapping = aes(x = `Measured At`, y = `Yearly Median`), color = "black", linewidth = plotting_arguments$smooth_data_linewidth, show.legend = plotting_arguments$legend_show) +
      geom_point(mapping = aes(x = `Measured At`, y = `Monthly Median`, alpha = `Is Warning`), shape = plotting_arguments$warning_point_shape, color = plotting_arguments$warning_point_color, size = plotting_arguments$warning_point_size, show.legend = plotting_arguments$legend_show) +
      facet_wrap(facets = . ~ `Instrument Code`, nrow = plot_grid_nrow, ncol = plot_grid_ncol, scales = "free_x") +
      labs(title = plotting_arguments$title, subtitle = plotting_arguments$subtitle) +
      scale_x_date(name = plotting_arguments$xlab, breaks = plotting_arguments$date_tick_distance, limits = c(as.IDate(from), as.IDate(to))) +
      scale_y_continuous(name = plotting_arguments$ylab, n.breaks = 10) +
      scale_alpha_manual(values = c("FALSE" = 0, "TRUE" = 1)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 20, vjust = 0.9),
            axis.ticks.length.x = unit(x = 3, units = "mm"),
            legend.position = plotting_arguments$legend_position,
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_rect(fill = "black"),
            strip.text = element_text(face = "bold", color = "white"))

    if(plot){
      plot(plot_object)
    }
    else{
      return(plot_object)
    }
  }
  if(warning == "peer_group"){
    plot_object <- ggplot(data = data) +
      geom_vline(xintercept = as.IDate(from)) +
      geom_segment(mapping = aes(x = `Measured At`, xend = `Measured At`, y = `Monthly Median`, yend = `Peer Group Monthly Median`, alpha = `Is Warning`), color = plotting_arguments$warning_segment_color, linewidth = plotting_arguments$warning_segment_linewidth, show.legend = plotting_arguments$legend_show) +
      geom_step(mapping = aes(x = `Measured At`, y = `Monthly Median`, color = `Instrument Code`), linewidth = plotting_arguments$smooth_data_linewidth, show.legend = plotting_arguments$legend_show) +
      geom_step(mapping = aes(x = `Measured At`, y = `Peer Group Monthly Median`), color = "black", linewidth = plotting_arguments$smooth_data_linewidth, show.legend = plotting_arguments$legend_show) +
      geom_point(mapping = aes(x = `Measured At`, y = `Monthly Median`, alpha = `Is Warning`), shape = plotting_arguments$warning_point_shape, color = plotting_arguments$warning_point_color, size = plotting_arguments$warning_point_size, show.legend = plotting_arguments$legend_show) +
      facet_wrap(facets = . ~ `Instrument Code`, nrow = plot_grid_nrow, ncol = plot_grid_ncol, scales = "free_x") +
      labs(title = plotting_arguments$title, subtitle = plotting_arguments$subtitle) +
      scale_x_date(name = plotting_arguments$xlab, breaks = plotting_arguments$date_tick_distance, limits = c(as.IDate(from), as.IDate(to))) +
      scale_y_continuous(name = plotting_arguments$ylab, n.breaks = 10) +
      scale_alpha_manual(values = c("FALSE" = 0, "TRUE" = 1)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 20, vjust = 0.9),
            axis.ticks.length.x = unit(x = 3, units = "mm"),
            legend.position = plotting_arguments$legend_position,
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_rect(fill = "black"),
            strip.text = element_text(face = "bold", color = "white"))

    if(plot){
      plot(plot_object)
    }
    else{
      return(plot_object)
    }
  }
}
