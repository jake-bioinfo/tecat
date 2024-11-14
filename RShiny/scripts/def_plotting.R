# Setting default theme for plotting
tx_fam <- "HersheySans"
def_th <- theme(
  
  # Hide panel borders and remove grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  panel.background = element_blank(),
  plot.margin = margin(t = 0, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
      
  # Adding title dimensions
  text = element_text(size = 28,
                      face = "bold", 
                      family = tx_fam),
  plot.title = element_text(size = 42, 
                            face = "plain", 
                            hjust = 0,
                            family = tx_fam,
                            margin = margin(t = 0, r = 0, b = 1, 
                                            l = 0, unit = "cm")),
  title = element_text(size = 35, 
                       face = "plain", 
                       hjust = 0.5,
                       family = tx_fam), 
  
  # Adding axis title margins
  axis.text.x = element_text(angle = 90, margin = margin(t = 0.3, unit = "cm")),
  axis.title.x = element_text(margin = margin(0.5, 0, 0, 0, unit = "cm"), 
                              face = "plain", size = 36),
  axis.title.y = element_text(margin = margin(0, 0.5, 0, 0, unit = "cm"), 
                              face = "plain", size = 36))
  
# Default scatter plot labels
plot_labels <- function(scat_title, y.lab, x.lab, color.lab, shape.lab) {
  plot_labels <- labs(title = scat_title, 
                      y = y.lab, 
                      x = x.lab, 
                      color = color.lab,
                      shape = shape.lab)
  return(plot_labels)
}

# Histogram labels
hist_lb <- function(hist_title, hist.x, hist.y, 
                    hist.color, hist.fill) { 
  hist_lb <- labs(title = hist_title , 
                  x = hist.x, 
                  y = hist.y, 
                  color = hist.color,
                  fill = hist.fill)
  return(hist_lb)
}

# Figure labels
fig_lb <- function(fig_title) {
  fig_lb <- labs(title = fig_title)
  return(fig_lb)
}

# default color palette
cl_pal <- scale_colour_brewer(palette = "Set1")

# def ggsave function
def.ggsave <- function(filename, plot = plot, dpi = 320, 
                       scale = 2, width = 16, 
                       height = 8, units = "in") {
  ggsave(filename=filename, plot = plot,
         height=height, width=width, 
         dpi=dpi, scale = scale, 
         units = "in", limitsize = FALSE)
}

# Scatter plot function
def.scatter <- function(df, scat_title, x.lab, y.lab,
                        color.lab, shape.lab,
                        read.length = read.length, 
                        tel.length = tel.length) {
  scatter <- ggplot(df, aes(x = read.length, y = tel.length)) +
    stat_density_2d(aes(fill = stat(nlevel)), 
                    geom = "polygon", position = "identity") +
    scale_fill_viridis_c()
  
  scatter <- scatter + plot_labels(scat_title, y.lab, x.lab, color.lab, shape.lab) +
    def_th +
    guides(colour = guide_legend(override.aes = list(size=10)),
           shape = guide_legend(override.aes = list(size=10))) +
    cl_pal

  return(scatter)
}

# Bar graph plot function
def.bar <- function(df, title, x.lab, y.lab,
                    color.lab) {
  bar <- ggplot(transform(df,
                          chr.end=factor(chromEnd, levels = c("5'", "3'"))), 
                   aes(x = ref_name, y = reads.per.chr, fill = sample)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_x_discrete(labels = as.character(df$ref_name), breaks = df$ref_name) +
  facet_wrap( ~ chr.end) 
  
  bar <- bar + hist_lb(title, x.lab, y.lab, color.lab, color.lab) +
    def_th +
    guides(colour = guide_legend(override.aes = list(size=10)),
           shape = guide_legend(override.aes = list(size=10))) +
    cl_pal

  return(bar)
}

# Retrieve mean and standard deviation from data sets
get_mean_sd <- function(input.df = result.df) {

  rts <- data.frame(sample = unique(input.df$sample), mean.l = 1, sd.l = 1)

  for (sample_types in unique(input.df$sample)) {
    samp_mean <- mean(filter(input.df, sample == sample_types)$telomere_length)
    rts <- rts %>% mutate(mean.l = case_when(sample== sample_types ~ samp_mean, T ~ mean.l))
  }
  for (sample_types in unique(input.df$sample)) {
    samp_sd <- sd(filter(input.df, sample == sample_types)$telomere_length)
    rts <- rts %>% mutate(sd.l = case_when(sample == sample_types ~ samp_sd, T ~ sd.l))
  }
  return(rts)
}
