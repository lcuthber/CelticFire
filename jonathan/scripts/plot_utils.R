cropped_ggsave <- function(save_path, plot=ggplot2::last_plot(), ...) {
  ggplot2::ggsave(save_path, plot=plot, ...)
  knitr::plot_crop(save_path)
}