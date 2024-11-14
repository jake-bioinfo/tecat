datasetInput <- function() {
  result.df <- read.csv("./data/example_data/output_results.csv")
  return(result.df)
}
result.df <- read.csv("./data/example_data/output_results.csv")


require(tidyverse)
require(grid)
source("./scripts/def_plotting.R")

plot_num_reads <- function(input.df = result.df) {
  bam_read_count <- ddply(input.df, .(sample, ref_name, chromEnd), summarize,
                          reads.per.chr = length(unique(telo_name)),
                          tel.length = mean(telomere_length)
  )
  
  
  plot_bam_read_count <- bam_read_count[order(bam_read_count$chromEnd, decreasing = TRUE), ]
  pbrc <- plot_bam_read_count %>%
    slice(gtools::mixedorder(ref_name))
  pbrc$ref_name <- factor(pbrc$ref_name, levels = unique(pbrc$ref_name))
  
  title <- "Number of Reads by End"
  y.lab <- "Telomere Read Count"
  x.lab <- "Chromosome"
  color.lab <- "Sample Type"
  
  bar_pl <- def.bar(
    pbrc,
    title, x.lab, y.lab,
    color.lab
  ) +
    def_th + theme(legend.key.size = unit(1, "lines")) +
    theme(legend.title = element_text(size = 16)) +
    theme(
      legend.position = "bottom",
      legend.justification = "right"
    )
  return(bar_pl)
}

#---- Density Plot
# Ampls for result.dd
get_dens_plot <- function(input.df = result.df) {
  modes <- ddply(input.df, .(sample),
                 summarize,
                 mode.1 = amps(telomere_length)$Peaks[, 1]
                 [which(amps(telomere_length)$Peaks[, 2] >=
                          sort(
                            amps(telomere_length)$Peaks[, 2],
                            TRUE
                          )[1])][1]
  )
  result.tel.stats <- get_mean_sd(input.df)
  
  h.title <- "B. Telomere Length Distribution"
  h.x <- "Telomere Length (kb)"
  h.y <- "Density (f/kb)"
  h.color <- "Sample Type"
  h.fill <- "Sample Type"
  grob <- grobTree(textGrob("Mode: dashed; Mean: solid",
                            x = 0.0, y = 0.90, hjust = 0,
                            gp = gpar(col = "red", fontsize = 13, fontface = "italic")
  ))
  
  # Generate histogram
  hist <- ggplot(
    input.df,
    aes(x = telomere_length, color = sample, fill = sample)
  ) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 60,
                   position = "identity",
                   alpha = 0.10
    ) +
    geom_vline(
      data = modes,
      aes(
        xintercept = mode.1,
        color = sample
      ),
      linetype = "dashed",
      alpha = 1,
      linewidth = 1
    ) +
    geom_density(alpha = 0.40) +
    geom_vline(
      data = result.tel.stats,
      aes(
        xintercept = mean.l,
        color = sample
      ),
      linetype = "solid",
      alpha = 1,
      linewidth = 1
    ) +
    geom_density(alpha = 0.40) +
    annotation_custom(grob)
  
  # Placing labels and theme
  densplot <- hist + hist_lb(
    h.title, h.x,
    h.y, h.color,
    h.fill
  ) +
    def_th + theme(
      legend.key.size = unit(4, "lines"), legend.position = "top",
      legend.text = element_text(size = 38),
      legend.spacing.x = unit(1.0, "cm")
    ) +
    cl_pal
  
  return(densplot)
}

viola_plot <- function(input.df = result.df, prime = "3'") {
  # Select
  if (prime == "3'") {
    select <- "3'"
  } else {
    select <- "5'"
  }
  
  # Reformatting data
  samples <- unique(input.df$sample)
  df <- filter(input.df, chromEnd == select)
  df <- df %>%
    slice(gtools::mixedorder(ref_name))
  df$ref_name <- factor(df$ref_name, levels = unique(df$ref_name))
  # Set defined colors
  color_ls <- c('#ed6b65', '#93bbff', '#5cd27f')
  for (i in 1:length(unique(df$sample))){
    df$colorls[df$sample == unique(df$sample)[i]] <- color_ls[i]
  }
  
  # Generating plot
  viola <- plot_ly(
    type = "violin", side = "negative", legendgroup = samples,
    scalegroup = samples, name = samples, width = 1600,
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    ), 
    points = "all",
    colors = "Set2"
  ) %>%
    layout(
      yaxis = list(title = list(text = "Telomere Length (kb)", font = list(size = 17))), xaxis = list(title = list(text = "Chromosome", font = list(size = 17))),
    title = list(text = paste0("Length by Chromosome End: ", select), font = list(size = 20), y = 1,x = .5)
    )
  
  sidecheck <- "positive"
  i = 1
  for (item in samples) {
    
    if (sidecheck == "positive") {
      viola <- viola %>%
        add_trace(
          x = df$ref_name[df$sample == item], y = df$telomere_length[df$sample == item],
          color = I(color_ls[i]),
          hoveron = "points+kde",
          legendgroup = item,
          scalegroup = item,
          name = item,
          side = "positive"
        )
      sidecheck <- "negative"
      # }
    } else {
      viola <- viola %>%
        add_trace(
          x = df$ref_name[df$sample == item], y = df$telomere_length[df$sample == item],
          color = I(color_ls[i]),
          hoveron = "points+kde",
          legendgroup = item,
          scalegroup = item,
          name = item,
          side = "negative"
        )
      sidecheck <- "positive"
    }
    i = i+1
  }
  return(viola)
}

barpldffun <- function(resultdf) {
  barpldf <- result.df %>%
    group_by(seqnames, s.name, chr.end) %>%
    summarise(`Telomere Read Count` = n())
  return(barpldf)
}