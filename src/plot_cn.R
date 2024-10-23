library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(gridExtra)

# Function to read and format coverage/clipping data
read_window_data <- function(coverage_file, left_clip_file, right_clip_file) {
  # Read files
  coverage <- read.table(coverage_file, 
                        col.names = c("chrom", "start", "end", "coverage"))
  left_clips <- read.table(left_clip_file, 
                          col.names = c("chrom", "start", "end", "left_clips"))
  right_clips <- read.table(right_clip_file, 
                           col.names = c("chrom", "start", "end", "right_clips"))
  
  # Merge data
  window_data <- coverage %>%
    left_join(left_clips, by = c("chrom", "start", "end")) %>%
    left_join(right_clips, by = c("chrom", "start", "end"))
  
  return(window_data)
}

# Function to read CNV calls
read_cnv_calls <- function(cnv_file) {
  cnv_calls <- read.table(cnv_file, 
                         col.names = c("chrom", "start", "end", 
                                     "copy_number", "confidence"))
  return(cnv_calls)
}

# Function to plot genomic region
plot_genomic_region <- function(window_data, cnv_calls, 
                              chrom, start_pos, end_pos,
                              highlight_regions = NULL) {
  # Filter data for region
  region_data <- window_data %>%
    filter(chrom == !!chrom,
           start >= start_pos,
           end <= end_pos)
  
  region_calls <- cnv_calls %>%
    filter(chrom == !!chrom,
           start >= start_pos,
           end <= end_pos)
  
  # Calculate plotting ranges
  max_coverage <- max(region_data$coverage) * 1.1
  max_clips <- max(c(region_data$left_clips, region_data$right_clips)) * 1.1
  
  # Create coverage plot
  p_coverage <- ggplot() +
    # Coverage data
    geom_line(data = region_data,
              aes(x = start, y = coverage),
              color = "grey50", size = 0.5) +
    # Left clips
    geom_line(data = region_data,
              aes(x = start, y = left_clips, color = "Left Clips"),
              size = 0.5) +
    # Right clips
    geom_line(data = region_data,
              aes(x = start, y = right_clips, color = "Right Clips"),
              size = 0.5) +
    # CNV calls as step function
    geom_step(data = region_calls,
              aes(x = start, y = copy_number * max_coverage/4,
                  color = "Copy Number"),
              size = 1) +
    # Optional highlighted regions
    {if (!is.null(highlight_regions))
      geom_rect(data = highlight_regions,
                aes(xmin = start, xmax = end,
                    ymin = -Inf, ymax = Inf),
                fill = "yellow", alpha = 0.2)
    } +
    # Aesthetics
    scale_color_manual(values = c("Left Clips" = "blue",
                                 "Right Clips" = "red",
                                 "Copy Number" = "darkgreen")) +
    scale_y_continuous(name = "Coverage",
                      sec.axis = sec_axis(~.* 4/max_coverage,
                                        name = "Copy Number")) +
    scale_x_continuous(labels = comma) +
    theme_bw() +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank()) +
    labs(x = paste("Position on", chrom),
         color = "Data Type")
  
  return(p_coverage)
}

# Function to create confidence track
plot_confidence_track <- function(cnv_calls, chrom, start_pos, end_pos) {
  region_calls <- cnv_calls %>%
    filter(chrom == !!chrom,
           start >= start_pos,
           end <= end_pos)
  
  p_conf <- ggplot(region_calls) +
    geom_segment(aes(x = start, xend = end,
                     y = confidence, yend = confidence),
                 color = "purple", size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(x = NULL,
         y = "Confidence")
  
  return(p_conf)
}

# Example usage function
plot_cnv_region <- function(window_data, cnv_calls,
                          chrom, start_pos, end_pos,
                          highlight_regions = NULL,
                          output_file = NULL) {
  # Create main plot
  p_main <- plot_genomic_region(window_data, cnv_calls,
                               chrom, start_pos, end_pos,
                               highlight_regions)
  
  # Create confidence track
  p_conf <- plot_confidence_track(cnv_calls, chrom, start_pos, end_pos)
  
  # Combine plots
  combined_plot <- grid.arrange(
    p_main, p_conf,
    heights = c(4, 1),
    ncol = 1
  )
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, combined_plot,
           width = 12, height = 8)
  }
  
  return(combined_plot)
}

# Example usage:
# window_data <- read_window_data("coverage.bed", "left_clips.bed", "right_clips.bed")
# cnv_calls <- read_cnv_calls("cnv_calls.bed")
# 
# # Plot specific region
# plot_cnv_region(
#   window_data = window_data,
#   cnv_calls = cnv_calls,
#   chrom = "chr1",
#   start_pos = 1000000,
#   end_pos = 2000000,
#   highlight_regions = data.frame(
#     start = c(1200000),
#     end = c(1300000)
#   ),
#   output_file = "cnv_plot.pdf"
# )