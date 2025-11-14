#### load packages ####
library(tidyverse)
library(cowplot)
library(truncnorm)
#### snakemake files ####

#### testing files ####
root_dir <- "~/Desktop/github/unimod_human"
figure_dir <- file.path(root_dir, "outputs/simulation/figures/steric_hindrance")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
#### end of parsing arguments ####

#### set up parameters ####
cell_num <- 20000
total_sites <- 2001

k <- 50
ksd <- c(0, 5, 10, 15, 25)
kmin <- 17
kmax <- 200

zeta <- 2000
zeta_sd <- 1000
zeta_min <- 500
zeta_max <- 4000

scale_factor <- 1e-3 # for quick visualization

theme_set(cowplot::theme_cowplot())

set.seed(20220816)

# plot varied pause sites
truncated_normal <- function(cell_num, k, ksd, kmin, kmax) {
  y <- vector()
  while (length(y) < cell_num) {
    x <- rnorm(cell_num * 2, mean = k, sd = ksd)
    x <- x[(x >= kmin) & (x <= kmax)]
    y <- c(y, x)
  }
  y <- round(y[1:cell_num]) 
  return(y)
}

plot_varied_k <- function(cell_num, k, ksd, kmin, kmax, add_sd_label = FALSE) {
 
  y <- truncated_normal(cell_num, k, ksd, kmin, kmax)
  
  p <- tibble(pause_site = y) %>%
    ggplot(aes(x = pause_site)) +
    geom_histogram(aes(y = ..density..), binwidth = 2) +
    # https://stackoverflow.com/questions/27644550/plotting-the-poisson-distribution-using-ggplot2s-stat-function
    stat_function(geom="line", n=kmax, fun=dtruncnorm, args=list(kmin, kmax, k, ksd),
                  color = "blue", linetype = "dashed", lwd = 1) +
    coord_cartesian(xlim = c(1, kmax)) +
    labs(x = "Position of Pause Sites", y = "Density")
  
  if (add_sd_label) {
  p <- p + annotate(geom = 'text', label = paste("sd =", ksd),
                    x = 150, y = Inf, vjust = 3)
  }
  
  return(p)
}

p_ls <- map(ksd[2:5], plot_varied_k, cell_num = cell_num, k = 50, kmin = kmin, kmax = kmax, add_sd_label = TRUE)
p <- plot_grid(plotlist = p_ls, nrow = 2, align = "vh")
ggsave(file.path(figure_dir, "histogram_varied_pause_sites_k50.png"), plot = p,
       width = 10, height = 6)

p <- plot_varied_k(cell_num = cell_num, ksd = 25, k = 50, kmin = kmin, kmax = kmax, add_sd_label = FALSE)
ggsave(file.path(figure_dir, "histogram_varied_pause_sites_k50sd25.pdf"), plot = p,
       width = 6, height = 3)

p_ls <- map(ksd[2:5], plot_varied_k, cell_num = cell_num, k = 70, kmin = kmin, kmax = kmax, add_sd_label = TRUE)
p <- plot_grid(plotlist = p_ls, nrow = 2, align = "vh")
ggsave(file.path(figure_dir, "histogram_varied_pause_sites_k70.png"), plot = p,
       width = 10, height = 6)

## get some fractions for correcting beta ##
# the fraction of cells when new initiation of a second RNAP is possible
pause_sites <- truncated_normal(cell_num = cell_num, k = 50, ksd = 25, kmin = kmin, kmax = kmax)
f33 <- mean(pause_sites > 33)
f50 <- mean(pause_sites > 50)
f70 <- mean(pause_sites > 70)
# the fraction of cells when a second RNAP is initiated and pause at position 21 or above
pause_sites <- truncated_normal(cell_num = cell_num, k = 50, ksd = 25, kmin = kmin, kmax = kmax)
f33 <- mean(pause_sites > 33 + 20)
f50 <- mean(pause_sites > 50 + 20)
f70 <- mean(pause_sites > 70 + 20)

# plot varied elongation rates
zv <- vector()
while (length(zv) < total_sites * cell_num * scale_factor) {
  x <- rnorm(total_sites * cell_num * scale_factor * 2, mean = zeta, sd = zeta_sd)
  x <- x[(x >= zeta_min) & (x <= zeta_max)]
  zv <- c(zv, x)
}
zv <- zv[1:(total_sites * cell_num * scale_factor)] 

p <- tibble(zeta = zv) %>%
  ggplot(aes(x = zeta)) +
  geom_histogram(aes(y = ..density..), bins = 50) +
  # https://stackoverflow.com/questions/27644550/plotting-the-poisson-distribution-using-ggplot2s-stat-function
  stat_function(geom="line", n=zeta_max, fun=dtruncnorm, args=list(zeta_min, zeta_max, zeta, zeta_sd),
                color = "blue", linetype = "dashed", lwd = 1) +
  coord_cartesian(xlim = c(1, zeta_max)) +
  labs(x = expression("Elongation Rate "*zeta))

ggsave(file.path(figure_dir, "histogram_varied_elongation_rates.png"), plot = p,
       width = 6, height = 3)