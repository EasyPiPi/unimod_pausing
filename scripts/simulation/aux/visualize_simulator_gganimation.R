#### Simulation of PRO-seq reads ####
library(S4Vectors)
library(tidyverse)
library(gganimate)
library(profvis)
library(microbenchmark)

root_dir <- "~/Desktop/github/unimod_human"
result_dir <- file.path(root_dir, "outputs/simulation")

# key parameters
k <- 50
m <- 250
l <- 20000 - k -m 

# k <- 5
# l <- 20
# m <- 5

alpha <- 1 # initiation rate, maybe 2.7 polymerases per cell per min?
beta <- 0.2 # pause release rate
gamma <- 1 # termination rate
zeta <- 50 # elongation rate 

cell_num <- 2e3
s <- 2 # Pol II size

delta_t <- 0.01
t <- 500
steps <- t / delta_t

total_sites <- k + l + m + 1
  
#### case for one gene in a cell ####
# probability vector to control transition from state to stat
# may need to accommodate varied zeta later
prob_init <- Rle(values = c(alpha * delta_t, delta_t, beta * delta_t, delta_t, gamma * delta_t) * zeta,
            lengths = c(1, k - 1, 1, l + m - 1, 1))

# initialize a vector to holding Pol II presence and absence in 
pos_init <- Rle(values = c(1, 0), lengths = c(1, k + l + m))

# assign the positions to a cell, get probability for comparison with draws to 
# determine which polymerases could be moved 
pos <- pos_init

# initialize a list to record polymerase positions
m_pos <- List()

# record start time
t1 <- Sys.time()

for (i in 1:steps) {
  # calculate the transition probability
  prob <- prob_init * pos
  # find which positions has probability larger than 0
  prob_pos <- prob > 0

  # determine whether polymerase can move or not
  # criteria 1, probability larger than random draw
  pos_pending <- which(prob_pos)
  
  # randomly draw numbers for comparison
  draws <- runif(n = length(pos_pending), min = 0, max = 1)
  
  c1 <- prob[pos_pending] > draws
  
  # criteria 2, enough space ahead to let polymerase go
  # test_rle <- Rle(values = c(1, 0, 1, 0, 1), lengths = c(1, 1, 1, 3, 1))
  # dist <- test_rle == 0
  
  if (length(c1) > 1) {
    c2 <- diff(pos_pending) > s 
    c2 <- c(c2, TRUE) } else {
      c2 <- c(TRUE)  
    }
  
  # decide which polymerase can eventually move
  pos_moving <- pos_pending[as.vector(c1) & c2]
  
  # advancing the polymerases to the next position
  if (length(pos_moving) > 0) {
    # set original spot as 0
    pos[pos_moving] <- 0
    pos[1] <- 1
    if (pos_moving[length(pos_moving)] == total_sites) {
      if (length(pos_moving) != 1) {
        pos[pos_moving[1:(length(pos_moving) - 1)] + 1] <- 1
      }} else {
        pos[pos_moving + 1] <- 1
      }
  }
  # record polymerase position
  m_pos[[i]] <- pos
}

# record end time
t2 <- Sys.time()

# make string to record params for figure and animation names
name_params <- paste0("k", k, "l", l, "m", m, "s", s,
                      "a", alpha, "b", beta, "g", gamma, "z", zeta,
                      "t", delta_t)
name_params <- str_remove_all(name_params, "\\.")

save.image(file = file.path(result_dir, paste0(name_params, ".RData")))

#### visualize polymerase position ####
# convert S4 vector List to a normal matrix
# m_pos_mx <- Reduce(rbind, lapply(m_pos, as.vector))
m_pos_mx <- Reduce(rbind, lapply(m_pos[(length(m_pos) - 999) : length(m_pos)], as.vector))

# exclude the first position which represents the free state
m_pos_mx <- m_pos_mx[, 2:NCOL(m_pos_mx)]
# rename columns to represent polymerase position on DNA
colnames(m_pos_mx) <- seq(1, NCOL(m_pos_mx))
# convert matrix to tibble for ggplot2
m_pos_plt <- m_pos_mx %>%
  as_tibble() %>%
  mutate(step = row_number()) %>%
  pivot_longer(cols = !step, names_to = "position", values_to = "polymerase") %>%
  mutate(position = factor(position, levels = seq(1, length(base::unique(position)))))

# visualize positions by static figure
p <- m_pos_plt %>%
  filter(polymerase > 0, step <= 100, position %in% c(1:70, 19971:20000)) %>%
  ggplot(aes(x = position, y = step)) +
  geom_point() +
  labs(y = 'Time point',
       x = 'Polymerase position') +
  ylim(c(100, 1))

ggsave(file.path(result_dir, 
                 paste0("polymerase_position_", name_params ,".png")), plot = p,
       width = 20, height = 12)

# visualize positions by animation
m_pos_ani <- m_pos_plt %>%
  # filter(polymerase > 0, step <= 5000) %>%
  filter(polymerase > 0, step <= 100, position %in% c(1:70)) %>%
  ggplot(aes(x = position, y = polymerase)) +
  # geom_point(size = 10) +
  geom_point(size = 2) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = 20)) + 
  labs(title = 'Time point: {frame_time}',
       x = 'Polymerase position') + 
  transition_time(step)

# configure animation for output
# m_pos_ani_out <- animate(m_pos_ani, height = 200, width = 1200, fps = 10, nframes = steps)

m_pos_ani_out <- animate(m_pos_ani, height = 200, width = 2000, fps = 10, nframes = 100)

anim_save(file.path(result_dir,
                    paste0("polymerase_position_", name_params ,".gif")),
          animation = m_pos_ani_out)
