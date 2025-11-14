#### description ####
# compare the original version and the c version of SimPol

#### log file ####
log <- file(snakemake@log[[1]], open="wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### snakemake files ####

#### load packages ####
library(tidyverse)
library(rtracklayer)

#### testing files ####
root_dir <- "~/Desktop/github/unimod_human"

figure_dir <- file.path(root_dir, "outputs/simulation/figures/compare_simpol")

meta_in <- file.path(root_dir, "metadata/simulation_params_lrt.csv")

#### end of parsing arguments ####
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
theme_set(cowplot::theme_cowplot())

kmax <- 200

meta <- read_csv(meta_in, show_col_types = FALSE)

meta$simPolv1 <- 
  map(file.path(root_dir, "outputs/simulation/data_lrt", paste0(meta$param_id, ".plus.bw")),
      ~ import.bw(.x)$score)

meta$simPolv2 <- 
  map(file.path(root_dir, "outputs/simulation/data_lrt", meta$param_id, "combined_cell_data.csv"),
      ~ as.integer(readLines(.x)[-1]))

simPolv1_tbl <- bind_rows(map(meta$simPolv1,
                              ~ c(
                                "tss_v1" = sum(.x[1:kmax]), "gb_v1" = sum(.x[(kmax + 1):length(.x)])
                              )))

simPolv2_tbl <- bind_rows(map(meta$simPolv2,
                              ~ c(
                                "tss_v2" = sum(.x[1:kmax]), "gb_v2" = sum(.x[(kmax + 1):length(.x)])
                              )))

meta <- meta %>% bind_cols(simPolv1_tbl, simPolv2_tbl)

rate_levels <- c(0.1, 0.2, 0.5, 0.8, 1, 1.2, 2, 5, 10)
meta <- meta %>%
  mutate(a_range = factor(a_range, levels = rate_levels),
         b_range = factor(b_range, levels = rate_levels))

p <- meta %>%
   ggplot(aes(x = tss_v1, y = tss_v2, color = a_range, shape = b_range)) +
   geom_point(size = 4) +
   scale_shape_manual(values=seq(0,15)) +
   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
   labs(x = "SimPol v1", y = "SimPol v2",
        color = expression("True "*alpha*zeta),
        shape = expression("True "*beta*zeta),
        title = "RNAP # within TSS region") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(figure_dir, "RNAP_number_within_TSS_region.png"), plot = p,
       width = 8, height = 6)

p <- meta %>%
  ggplot(aes(x = gb_v1, y = gb_v2, color = a_range, shape = b_range)) +
  geom_point(size = 4) +
  scale_shape_manual(values= seq(0, 15)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "SimPol v1", y = "SimPol v2",
       color = expression("True "*alpha*zeta),
       shape = expression("True "*beta*zeta),
       title = "RNAP # within gene body") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(figure_dir, "RNAP_number_within_gene_body.png"), plot = p,
       width = 8, height = 6)  

# pick one example to visualize read counts within TSS region
plot_tbl <- meta %>%
  filter(a_range == 1, b_range == 1)

p <- plot_tbl %>%
  select(simPolv1, simPolv2) %>%
  unnest(cols = c(simPolv1, simPolv2)) %>%
  mutate(position = row_number()) %>%
  filter(position <= 300) %>%
  pivot_longer(cols = c(simPolv1, simPolv2)) %>% 
  ggplot(aes(x = position, y = value)) +
  geom_col() +
  geom_vline(xintercept = c(17, 200), linetype = "dashed", color = "gray") +
  facet_grid(name ~ .) +
  labs(x = "Distance from TSS", y = "Number of RNAP")

ggsave(file.path(figure_dir, "RNAP_number_within_TSS_region_per_nuc.png"), plot = p,
       width = 12, height = 6)  
