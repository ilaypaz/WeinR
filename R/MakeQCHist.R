#' GC content histogram (average per read)
#'
#' Produces a histogram of per-read GC content values from a BAM file.
#'
#' @param prGCContent BAM file path or numeric vector/list of per-read GC content values.
#' @param file_name Character title for the plot.
#'
#' @return A ggplot object, or NULL if no finite GC values are available.
#' @keywords internal
make_gc_plot <- function(prGCContent, file_name) {
  
  if (base::length(prGCContent) == 1L &&
      base::is.character(prGCContent) &&
      base::grepl("\\.bam$", prGCContent, ignore.case = TRUE)) {
    
    bam_data <- Rsamtools::scanBam(
      prGCContent,
      param = Rsamtools::ScanBamParam(
        what = "seq"
      )
    )
    
    sequences <- toupper(as.character(bam_data[[1]]$seq))
    
    gc_count <- lengths(regmatches(sequences, gregexpr("[GC]", sequences)))
    valid_count <- lengths(regmatches(sequences, gregexpr("[ATGC]", sequences)))
    
    prGCContent <- gc_count / valid_count
    prGCContent[valid_count == 0] <- NA_real_
  }
  
  prGCContent <- as.numeric(base::unlist(prGCContent, use.names = FALSE))
  prGCContent <- prGCContent[base::is.finite(prGCContent)]
  
  if (base::length(prGCContent) == 0L) {
    return(NULL)
  }
  
  ggplot2::ggplot(
    data = data.frame(gc_content = prGCContent),
    ggplot2::aes(x = gc_content)
  ) +
    ggplot2::geom_histogram(
      bins = 40,
      fill = "#CC79A7",
      color = "white",
      alpha = 0.9
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      expand = c(0, 0),
      breaks = scales::pretty_breaks(n = 10)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = c(0, 0),
      breaks = scales::pretty_breaks(n = 10)
    ) +
    ggplot2::labs(
      x = "GC Content (proportion)",
      y = "Read Count",
      title = file_name
    ) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.6),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.6),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    )
}

make_gc_plot(
  "/Users/ilaypaz/Downloads/capstone/PBM60157_1pct.bam",
  "PBM60157_1pct GC Content"
)
