###################################################
###################################################
## Disentangling the Effects of Intercropping on ##
## Vector-Borne Plant-Virus Dynamics             ##
###################################################
###################################################

# The code below was used to generate all numerical 
# simulations within the aforementioned manuscript

############################################
## Fig. 5: Maximum Amplification Envelope ##
############################################

library(expm)
library(tidyverse)
library(patchwork)
library(ggnewscale)

# parameters
P1   <- 7397
P2   <- 9443
VT   <- 342
mu1  <- 0.56
mu2  <- 0.6
muV  <- 1.26
Ptot <- P1 + P2
largec1 <- 6.51
larged1 <- 7.62
smallc1 <- 0.18
smalld1 <- 0.64

# interval chosen for r_T: [0,5]
Tmax   <- 5
t_vals <- seq(0, Tmax, length.out = 50)

# transmission rate variables
c_vals <- seq(0.1, 9, length.out = 50)
d_vals <- seq(0.1, 9, length.out = 50)

# jacobian
compute_J <- function(c2, d2, c1, d1) {
  matrix(c(
    -mu1,          0,         d1 * P1 / Ptot,
    0,        -mu2,         d2 * P2 / Ptot,
    VT * c1/Ptot, VT * c2/Ptot,  -muV
  ), 3, 3, byrow=TRUE)
}

#basic reproduction number, R0
compute_R0 <- function(c2, d2, c1, d1) {
  sum_term <- (c1 * d1 * P1 / mu1) + (c2 * d2 * P2 / mu2)
  sqrt( VT / (muV * Ptot^2) * sum_term )
}

# hamiltonian of jacobian
compute_RJ <- function(c2, d2, c1, d1) {
  J <- compute_J(c2, d2, c1, d1)
  H <- (J + t(J)) / 2
  max(Re(eigen(H, only.values=TRUE)$values))
}

# max. amp. envelope, r_T
compute_rT <- function(c2, d2, c1, d1) {
  R0  <- compute_R0(c2, d2, c1, d1)
  RJ  <- compute_RJ(c2, d2, c1, d1)
  # only compute r_T if R0<1 (stable DFE) AND RJ>0 (reactive)
  if (!(R0 < 1 && RJ > 0)) return(NA_real_)
  J <- compute_J(c2, d2, c1, d1)
  max(sapply(t_vals, function(t) norm(expm(J * t), type="2")), na.rm=TRUE)
}

# low-high scenarios
scenarios <- tibble(
  label = c("Low", "High"),
  c1    = c(smallc1, largec1),
  d1    = c(smalld1, larged1)
)

# RJ: reactivity
# when RJ<0 = initial resilience

# plot r_T gradient, where RJ<0 and where R0>1
plots <- map(scenarios %>% split(.$label), function(scn) {
  df <- expand.grid(c2 = c_vals, d2 = d_vals) %>%
    as_tibble() %>%
    rowwise() %>%
    mutate(
      R0  = compute_R0(c2, d2, scn$c1, scn$d1),
      RJ  = compute_RJ(c2, d2, scn$c1, scn$d1),
      r_T = compute_rT(c2, d2, scn$c1, scn$d1)
    ) %>%
    ungroup()
  
  # midpoint at mid of r_T range
  midpt <- mean(range(df$r_T, na.rm=TRUE))
  
  df <- df %>%
    mutate(
      overlay = case_when(
        R0 > 1     ~ "Outbreak",
        RJ < 0     ~ "Initial Resilience",
        TRUE       ~ NA_character_
      )
    )
  
  ggplot(df, aes(x = c2, y = d2)) +
    # base r_T heatmap (only where R0<1 and RJ>0)
    geom_raster(data = filter(df, is.na(overlay)), aes(fill = r_T)) +
    
    # r_T color scale
    scale_fill_gradient2(
      low      = "#FAF18A", 
      mid      = "#F0EAD6", 
      high     = "#4E95D9",
      midpoint = midpt,
      na.value = "grey90",
      name     = expression(r[T])
    ) +
    
    ggnewscale::new_scale_fill() +
    geom_raster(
      data = filter(df, !is.na(overlay)),
      aes(fill = overlay)
    ) +
    scale_fill_manual(
      values = c(
        "Outbreak" = "darkblue",
        "Initial Resilience" = "#FFC20A"
      ),
      guide = "none"
    ) +
    
    labs(
      title = scn$label,
      x     = expression(c[2]),
      y     = expression(d[2])
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    coord_fixed()
  
})

# combine scenarios
plots$Low + plots$High +
  plot_layout(nrow=1, guides="keep")

############################
## Fig.6 and 7: Delta R_0 ##
############################

## The code below is for generating Fig. 7

## Rerun this code with relevant parameter values
## to generate Fig. 6

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpattern)

# Fig. 7 parameters
params <- list(
  c1 = 3.65, 
  c2 = 0.93, 
  d1 = 9.02, 
  d2 = 2.37, mu1 = 1.25,
  mu2 = 2.31, 
  V_T = 500, 
  muV = 1.78
)

# R0 for monoculture system
R0_crop <- function(x, y, p) {
  with(p, sqrt((V_T / (muV * x^2)) * (c1 * d1 * x / mu1 + y*0)))
}

# R0 for intercropped system
R0_intercrop <- function(x, y, p) {
  with(p, sqrt(
    (V_T / (muV * (x + y)^2)) * 
      ((c1 * d1 * x / mu1) + (c2 * d2 * y / mu2))
  ))
}

# Define a grid of total population size variables
x_vals <- seq(1, 20000, length.out = 500)
y_vals <- seq(1, 20000, length.out = 500)

# calculate R0 for different systems
# as well as Delta R0
# and split into various scenarios/regimes
grid <- expand.grid(x = x_vals, y = y_vals) %>%
  mutate(
    R0_crop = R0_crop(x, y, params),
    R0_intercrop = R0_intercrop(x, y, params),
    rel_change = (R0_intercrop - R0_crop) / R0_crop,
    regime = case_when(
      R0_crop < 1 & R0_intercrop < 1 ~ "Both Stable",
      R0_crop > 1 & R0_intercrop > 1 ~ "Both Outbreaks",
      R0_crop < 1 & R0_intercrop > 1 ~ "Stable to Outbreak",
      R0_crop > 1 & R0_intercrop < 1 ~ "Outbreak to Stable"
    )
  )

# add in signs of Delta R0 for plotting
grid_rel <- grid %>%
  mutate(
    R0_crop = R0_crop(x, y, params),
    R0_intercrop = R0_intercrop(x, y, params),
    rel_change = (R0_intercrop - R0_crop) / R0_crop,
    regime = case_when(
      R0_crop < 1 & R0_intercrop < 1 ~ "Both Stable",
      R0_crop > 1 & R0_intercrop > 1 ~ "Both Outbreaks",
      R0_crop < 1 & R0_intercrop > 1 ~ "Stable to Outbreak",
      R0_crop > 1 & R0_intercrop < 1 ~ "Outbreak to Stable"
    ),
    rel_sign = ifelse(rel_change > 0, "+", "-"),
    regime_sign = paste0(regime, " (", rel_sign, ")")
  )

# filter out NA or Inf values
grid_clean <- grid_rel %>%
  filter(!is.na(rel_change), !is.infinite(rel_change))

fill_colors <- c(
  "Both Stable (+)"       = "lightgreen",
  "Both Stable (-)"       = "forestgreen",
  
  "Both Outbreaks (+)"    = "yellow", 
  "Both Outbreaks (-)"    = "gold", 
  
  "Stable to Outbreak (+)"= "#FF7F00",
  "Stable to Outbreak (-)"= "#FDB96E",
  
  "Outbreak to Stable (+)"= "lightblue",
  "Outbreak to Stable (-)"= "darkblue"
)

# plot
ggplot(grid_clean, aes(x = x, y = y, fill = regime_sign)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_manual(values = fill_colors, name = expression("Regime" ~ (Delta*R[0]))) +
  labs(
    x = expression(P[1]),
    y = expression(P[2]),
    title = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

##############################
## Fig. 8 (right): Simplex ###
##############################

# The following code was used to generate
# the simplex triangle plot in Fig. 8 (right)

# The same code, with parameters changed, generates
# plots similar to that of Fig. 8 (left), where 
# only one plant dominates in terms of total population size

library(ggtern)
library(dplyr)

# total effort allowed
C <- round(runif(1,1,10))
# total population sizes, where two are equal
P <- c(90, 90, 20)
names(P) <- c("alpha1", "alpha2", "alpha3")

V_T <- round(runif(1,97,300))

# grid of alpha values over simplex
step <- 0.01
alpha1 <- seq(0, C, by = step)
alpha2 <- seq(0, C, by = step)
grid <- expand.grid(alpha1 = alpha1, alpha2 = alpha2) %>%
  mutate(
    # alpha3 computed so that alphas sum to 1
    alpha3 = C - alpha1 - alpha2
  ) %>%
  filter(alpha3 >= 0) %>%
  mutate(
    fval = alpha1 * P[1] + alpha2 * P[2] + alpha3 * P[3]
  )

# find optimal point, i.e. all effort on max P
opt_index <- which.max(P)
opt_point <- c(0, 0, 0)
opt_point[opt_index] <- C
opt_df <- as.data.frame(as.list(opt_point))
names(opt_df) <- c("alpha1", "alpha2", "alpha3")

# compute max R0
P_sum <- sum(P)
P_max <- max(P)
R0_max <- sqrt((V_T / (P_sum^2)) * C * P_max)

triangle_border <- data.frame(
  alpha1 = c(C, 0, 0),
  alpha2 = c(0, C, 0),
  alpha3 = c(0, 0, C)
)

# simplex plot of f for varying alpha values,
# showing where optima lie 
ggtern(data = grid, aes(x = alpha1, y = alpha2, z = alpha3)) +
  geom_point(aes(color = fval), linewidth = 1, alpha = 1) +
  geom_polygon(data = triangle_border, aes(x = alpha1, y = alpha2, z = alpha3),
               fill = NA, color = "black", size = 1) +
  scale_color_gradientn(colors = c("gold", "white", "darkblue")) +
  labs(
    color = expression(f) 
  ) +
  theme(
    legend.position = "right"
  ) +
  labs(x = expression(alpha[1]), y = expression(alpha[2]), z = expression(alpha[3]))

########################################
## Supp. Material: Infected Bias Plot ##
########################################

library(deSolve)
library(ggplot2)
library(dplyr)
library(patchwork)

# parameters
P1 <- 10000; P2 <- 5000; VT <- 500
d1 <- 2.08; d2 <- 3.17
c1 <- 7.21; c2 <- 0.59
mu1 <- 0.16; mu2 <- 1.85
muV <- 1.5
Tmax <- 100
initcond <- c(A = 24, B = 14, C = 21)

# Note: The model below (model_reduced) can also be used to 
# simulate individual system trajectories by plotting output 
# from desolve() for various initial conditions,
# i.e. for Supp. Material section "Simulation Scenarios" 

model_reduced <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    dA <- (d1 * C * (P1 - A)) / (P1 + (v1 - 1) * A + P2 + (v2 - 1) * B) - mu1 * A
    dB <- (d2 * C * (P2 - B)) / (P1 + (v1 - 1) * A + P2 + (v2 - 1) * B) - mu2 * B
    dC <- (VT - C) * ((c1 * A + c2 * B) / (P1 + (v1 - 1) * A + P2 + (v2 - 1) * B)) - muV * C
    list(c(dA, dB, dC))
  })
}

# bias values
v_vals <- seq(0.1, 2, length.out = 30)
param_grid <- expand.grid(v1 = v_vals, v2 = v_vals)

library(purrr)

results <- param_grid %>%
  pmap_dfr(function(v1, v2) {
    
    parms <- list(
      P1 = P1, P2 = P2, VT = VT,
      d1 = d1, d2 = d2,
      c1 = c1, c2 = c2,
      mu1 = mu1, mu2 = mu2, muV = muV,
      v1 = v1, v2 = v2
    )
    
    out <- ode(
      y = initcond,
      times = seq(0, Tmax, by = 1),
      func = model_reduced,
      parms = parms,
      method = "ode45"
    )
    
    #take last 10 iterations
    out_df <- as.data.frame(out)
    tail_window <- out_df %>% 
      filter(time >= Tmax - 10)
    
    tibble(
      v1 = v1,
      v2 = v2,
      # get mean densities/abundances over time interval
      P1I_T = mean(tail_window$A),
      P2I_T = mean(tail_window$B),
      VI_T  = mean(tail_window$C)
    )
  })

plot_contour <- function(df, zvar, title) {
  ggplot(df, aes(x = v1, y = v2)) +
    geom_raster(aes_string(fill = zvar)) +
    #geom_contour(aes_string(z = zvar), color = "white", alpha = 0.6) +
    scale_fill_gradient(low = "darkblue", high = "gold", name = NULL) +
    labs(
      title = title,
      x = expression(v[1]),
      y = expression(v[2])
    ) +
    coord_fixed(ratio = 3.5) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      axis.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5)
    )
}

p1 <- plot_contour(results, "P1I_T", "Crop")
p2 <- plot_contour(results, "P2I_T", "Intercrop")
p3 <- plot_contour(results, "VI_T",  "Vector")

# combine all three variables
(p1 | p2 | p3)

###########################################
## Supp. Material: Transient Bounds Plot ##
###########################################

library(expm)
library(ggplot2)
library(tidyr)
library(dplyr)

# parameters
d1 <- 8.4
d2 <- 3.35
P1 <- 1000
P2 <- 200
v1 <- 1.88
v2 <- 1.96
VT <- 50
c1 <- 6.54
c2 <- 12.14
mu1 <- 2.63
mu2 <- 5.93
muV <- 4.44
Ptot <- P1 + P2

# jacobian
J <- matrix(c(
  -mu1,       0,        d1 * P1 / Ptot,
  0,     -mu2,       d2 * P2 / Ptot,
  VT*c1/Ptot, VT*c2/Ptot,    -muV
), 3, 3, byrow = TRUE)

# next-generation matrix
Fmat <- matrix(c(
  0,   0,     d1 * P1 / Ptot,
  0,   0,     d2 * P2 / Ptot,
  VT*c1/Ptot, VT*c2/Ptot, 0
), 3, 3, byrow = TRUE)

Vmat <- diag(c(mu1, mu2, muV))
NGM  <- Fmat %*% solve(Vmat)
R0   <- max(abs(eigen(NGM, only.values = TRUE)$values))

# initial condition vector y0 (so that sums to 1)
x <- c(0.0129, 0.165)
w_raw <- c(x[1], x[2] - x[1], 1 - x[2])
w_rounded <- round(w_raw, 3)
diff <- 1 - sum(w_rounded)
i_max <- which.max(w_rounded)
w_rounded[i_max] <- w_rounded[i_max] + diff
y0 <- w_rounded

# times
times <- seq(0, 3, length.out = 1095)

# compute linearised trajectory at each time 
# i.e. matrix exponential times initial condition
Y_linear <- sapply(times, function(t) expm(J * t) %*% y0)
out <- as.data.frame(t(Y_linear))
names(out) <- c("P1I", "P2I", "VI")
out$time <- times

# calculate the norm of the trajectory vector at each time
out$N <- apply(out[, c("P1I", "P2I", "VI")], 1, 
               function(row) norm(as.matrix(row), type = "2"))

# compute reactivity and bounds
H <- (J + t(J)) / 2
Rreact <- max(eigen(H)$values)

# get norm of matrix exponential at each time
r_t <- sapply(times, function(t) norm(expm(J * t), type = "2"))
r_max <- max(r_t)

# compute the Kreiss bound
f_r1 <- function(r) {
  if (r <= 0) return(-Inf)
  M <- try(solve(r * diag(3) - J), silent = TRUE)
  if (inherits(M, "try-error")) return(-Inf)
  r * norm(M, type = "2")
}
opt <- optimize(f_r1, interval = c(1e-6, 20), maximum = TRUE)
Kstar <- opt$objective

out_long <- out %>%
  select(time, P1I, P2I, VI, N) %>%
  pivot_longer(cols = -time, names_to = "variable", values_to = "value")

bounds <- data.frame(
  name = c("K", "R", "r"),
  value = c(Kstar, Rreact, r_max),
  linetype = c("solid", "dotted", "dashed"),
  color = c("black", "black", "black")
)

bounds_long <- bounds %>%
  mutate(variable = name) %>%
  select(-name) %>%
  group_by(variable, value, linetype, color) %>%
  do(data.frame(time = out$time, value = .$value)) %>%
  ungroup()

# combine data
combined_df <- bind_rows(
  out_long %>% mutate(linetype = "solid", color = variable),
  bounds_long
) %>%
  mutate(variable = factor(variable, 
                           levels = c("P1I", "P2I", "VI", "N", "K", "R", "r")))

# plot all bounds and linearised dynamics
ggplot(combined_df, aes(x = time, 
                        y = value, 
                        color = variable, 
                        linetype = variable)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c(
    "P1I" = "darkblue",
    "P2I" = "blue",
    "VI" = "gold",
    "N" = "yellow",
    "K" = "black",
    "R" = "black",
    "r" = "grey"
  ),
  labels = c(
    "P1I" = "Crop",
    "P2I" = "Intercrop",
    "VI" = "Vector",
    "N" = "Norm",
    "K" = expression(K[J]),
    "R" = expression(R[J]),
    "r" = expression(r[T])
  )) +
  scale_linetype_manual(values = c(
    "P1I" = "solid",
    "P2I" = "solid",
    "VI" = "solid",
    "N" = "solid",
    "K" = "solid",
    "R" = "dotted",
    "r" = "solid"
  ), labels = c(
    "P1I" = "Crop",
    "P2I" = "Intercrop",
    "VI" = "Vector",
    "N" = "Norm",
    "K" = expression(K[J]),
    "R" = expression(R[J]),
    "r" = expression(r[T])
  )) +
  labs(x = "Time", y = "Value", color = "", linetype = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right", 
        legend.box = "vertical",
        aspect.ratio = 1) +
  coord_cartesian(ylim = c(0, max(combined_df$value) * 1.1))
