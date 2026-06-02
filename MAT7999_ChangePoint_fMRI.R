library(oro.nifti)
library(changepoint)
library(igraph)
library(gSeg)
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)

################
# PART 0 — SETUP
################

#### 0.1 Choose Subject, Task, and Run

sub_id <- "sub-01" # 1-20 subjects
task_name <- "nofeedback" # nofeedback, auditoryfeedback, visualfeedback
run_id <- "01"


#### 0.2 Inspect Available Experimental Runs

####Available functional runs in this dataset:####
# - faceexplocalizer: functional localizer used to identify pSTS ROI
# - auditoryfeedback: imagery task with auditory neurofeedback
# - visualfeedback: imagery task with visual neurofeedback
# - nofeedback: imagery task without feedback


#### 0.3 Select Representative Task Condition

# The nofeedback run is selected for initial analysis because it preserves
# the full facial-expression imagery task structure while excluding added
# auditory and visual feedback mechanisms. This provides the clearest
# baseline setting for evaluating changepoint behavior before extending
# analysis to more the complex feedback-based conditions.


#### 0.4 Define file paths

# Define path to selected subject's BOLD scan for analysis
bold_path <- file.path(
  "data_raw", sub_id, "ses-01", "func",
  paste0(sub_id, "_ses-01_task-", task_name, "_run-", run_id, "_bold.nii.gz")
)


# Define path to selected subject's event timing file for experimental comparison
events_path <- file.path(
  "data_raw", sub_id, "ses-01", "func",
  paste0(sub_id, "_ses-01_task-", task_name, "_run-", run_id, "_events.tsv")
)

#### 0.5 Verify that the required raw data files exist before analysis begins

if (!file.exists(bold_path)) {
  stop("BOLD file not found. Check bold_path.")
}

if (!file.exists(events_path)) {
  stop("Event file not found. Check events_path.")
}

#### 0.6 Build function to extract the mst

# Build mst between all time points
build_mst_from_dist <- function(D) { #D is the pairwise distance matrix
  g_full <- graph_from_adjacency_matrix( 
    D, #makes a full weighted graph from matrix D
    mode = "undirected", # i -> j <=> j -> i
    weighted = TRUE, # defines matrix entries as edge weights
    diag = FALSE #invalid function when i = j
  )
  mst(g_full, weights = E(g_full)$weight) #extract mst from g_full
}

##########################################
# PART 1 — SINGLE-SUBJECT DATA PREPARATION
##########################################

#### 1.1 Read Raw fMRI Data

# Read the subject's raw BOLD scan in NIfTI format(scan object; not raw numbers)
nii <- readNIfTI(bold_path, reorient = FALSE)

# Extract the 4-dim image array from the NIfTI object
img <- nii@.Data

# Check for 64x64x33x300
dim(img)

#### 1.2 Construct Brain Mask

mask <- img[,,,1] != 0 #voxels with non-zero intensity on first brain volume measured
num_masked_voxels <- sum(mask) #count retained voxels
num_masked_voxels #should be less than or equal to 64*64*33 = 135,168

#### 1.3 Reshape to Time-by-Voxel Matrix
T_len <- dim(img)[4] #store number of time points in the fMRI scan
V_len <- num_masked_voxels

# Extract all masked voxel values by repeating the mask across all time 
# points and reshape into a voxel-by-time matrix
X_vox <- matrix(img[rep(mask, T_len)], nrow = V_len, ncol = T_len)
dim(X_vox) #voxels x 300

X <- t(X_vox) # time-by-voxel matrix for analysis
dim(X) #300 x voxels

#### 1.4 Inspect Example Voxel Signal

summary(X[,1]) #first voxel over all time

plot(
  X[,1],
  type = "l",
  main = "First kept voxel over time",
  xlab = "Time point",
  ylab = "Intensity"
)#should have moderate fluctuations, continuous, finite trends

#### 1.5 Remove Unusable Voxels(most likely not needed)

vox_sd <- apply(X, 2, sd) #sd of each voxel over time
vox_mean <- colMeans(X) #average voxel intensity over time

# Keep positive, finite, and non-constant voxels
keep <- is.finite(vox_sd) & vox_sd > 0 &
  is.finite(vox_mean) & vox_mean > 0 

X2 <- X[, keep, drop = FALSE] #keep final good voxels
dim(X2) #300 x less than or equal to num_masked_voxels

##########################################
# PART 1B — 1-D MEAN INTENSITY TIME SERIES
##########################################

mean_intensity <- rowMeans(X2)

intensity_df <- data.frame(
  time = 1:length(mean_intensity),
  mean_intensity = mean_intensity
)

# Change points from graph-based methods for sub-01 nofeedback
cp_df <- data.frame(
  cp = c(267, 150),
  method = c("Original Graph CP", "Weighted Graph CP")
)

fig_mean_intensity <- ggplot(intensity_df, aes(x = time, y = mean_intensity)) +
  geom_line(linewidth = 0.8, color = "black") +
  
  geom_vline(
    data = cp_df,
    aes(xintercept = cp, linetype = method),
    color = "red",
    linewidth = 1
  ) +
  
  scale_linetype_manual(
    values = c(
      "Original Graph CP" = "solid",
      "Weighted Graph CP" = "dashed"
    )
  ) +
  
  labs(
    title = "Mean fMRI Intensity Over Time",
    x = "Time Point",
    y = "Mean Intensity",
    linetype = "Method"
  ) +
  
  theme(
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

fig_mean_intensity

ggsave(
  filename = "figures/mean_intensity_time_series.png",
  plot = fig_mean_intensity,
  width = 7,
  height = 5,
  dpi = 300
)

##########################################
# PART 2A — DIMENSIONALITY REDUCTION (PCA)
##########################################

# Apply PCA to the cleaned voxel-by-time matrix.
# Rows = time points, Columns = retained voxels.

pca_fit <- prcomp(X2, center = TRUE, scale. = TRUE)


#### 2A.1 Variance Explained


# Standard deviations of PCs
pca_sd <- pca_fit$sdev

# Eigenvalues (variance of each PC)
pca_var <- pca_sd^2

# Proportion of variance explained
pca_pve <- pca_var / sum(pca_var)

# Cumulative variance explained
cum_var <- cumsum(pca_pve)


#### 2A.2 Choose Number of PCs

# Smallest number of PCs explaining at least 80%
k <- which(cum_var >= 0.80)[1]
k


############################
# 2A.3 PCA Diagnostic Figures
############################

pca_df <- data.frame(
  PC = 1:length(pca_pve),
  PVE = pca_pve,
  CumVar = cum_var
)

# Scree Plot: First 30 PCs

# Visual elbow chosen for display
elbow_pc <- 19

fig_scree <- ggplot(pca_df[1:30, ], aes(x = PC, y = PVE)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = elbow_pc, linetype = 2) +
  annotate(
    "text",
    x = elbow_pc + 1,
    y = max(pca_df$PVE[1:20]) * 0.75,
    label = paste0("levels out\naround PC ", elbow_pc),
    hjust = 0,
    size = 3.5
  ) +
  labs(
    title = "Scree Plot of First 30 Principal Components",
    x = "Principal Component",
    y = "Proportion of Variance Explained"
  )

fig_scree

ggsave(
  filename = "figures/pca_screeplot.png",
  plot = fig_scree,
  width = 7,
  height = 5,
  dpi = 300
)


# Cumulative Variance Plot
fig_pca_cumvar <- ggplot(pca_df[1:300, ], aes(x = PC, y = CumVar)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.80, linetype = 2) +
  geom_vline(xintercept = k, linetype = 2) +
  annotate(
    "text",
    x = k - 70,
    y = 0.87,
    label = paste0("80% threshold\nk = ", k),
    size = 3.8
  ) +
  labs(
    title = "Cumulative Variance Explained by PCA",
    x = "Number of Principal Components",
    y = "Cumulative Proportion of Variance Explained"
  )

fig_pca_cumvar

ggsave(
  filename = "figures/pca_cumulative_variance.png",
  plot = fig_pca_cumvar,
  width = 7,
  height = 5,
  dpi = 300
)


#### 2A.4 Build Reduced Score Matrix

Z <- pca_fit$x[, 1:k, drop = FALSE]
dim(Z)   # expected: 300 x k


#### 2A.5 Selected PCs Over Time

# Compare selected PCs representing early, elbow-level, and retained components
pc_select <- c(1, 4, 19, k)
pc_labels <- paste0("PC", pc_select)

pc_selected_df <- data.frame(
  Time = 1:nrow(pca_fit$x)
)

for (pc in pc_select) {
  pc_selected_df[[paste0("PC", pc)]] <- pca_fit$x[, pc]
}

pc_selected_long <- pc_selected_df %>%
  tidyr::pivot_longer(
    cols = -Time,
    names_to = "Component",
    values_to = "Score"
  )

# Force legend order
pc_selected_long$Component <- factor(
  pc_selected_long$Component,
  levels = pc_labels
)

fig_selected_pcs <- ggplot(
  pc_selected_long,
  aes(x = Time, y = Score, color = Component)
) +
  geom_line(linewidth = 0.5) +
  labs(
    title = "Selected Principal Components Over Time",
    x = "Time Point",
    y = "PC Score",
    color = "Component"
  ) +
  theme(
    legend.position = "top"
  )

fig_selected_pcs

ggsave(
  filename = "figures/selected_pcs_time.png",
  plot = fig_selected_pcs,
  width = 7,
  height = 5,
  dpi = 300
)

##########################################
# PART 2B — BENCHMARK CHANGE POINT METHODS
##########################################
# Benchmark methods are used as comparisons to the graph-based method.
# PC1 is used for classical univariate methods.
# The first 5 PCs are used for multivariate benchmark methods.

pc1 <- Z[, 1]
Z_bench <- Z[, 1:min(5, ncol(Z)), drop = FALSE]
n_bench <- nrow(Z_bench)

#### 1. PELT

# PELT on PC1 tended to over-segment this dataset.
# We report the number of detected changes rather than listing all of them.

pelt_result <- cpt.mean(
  pc1,
  method = "PELT",
  penalty = "Manual",
  pen.value = "10*log(n)"
)

pelt_cpts <- cpts(pelt_result)
length(pelt_cpts)
cpts(pelt_result)


#### 2. Binary Segmentation

binseg_result <- cpt.mean(
  pc1,
  method = "BinSeg",
  penalty = "BIC",
  Q = 5
)

binseg_cpts <- cpts(binseg_result)
binseg_cpts


#### 3. Energy Method

library(ecp)

set.seed(123)

energy_result <- e.divisive(
  Z_bench,
  R = 499,
  sig.lvl = 0.05
)

energy_cpts <- energy_result$estimates

# Remove endpoints if included
energy_cpts_interior <- energy_cpts[
  energy_cpts != 1 & energy_cpts != (n_bench + 1)
]

energy_cpts
energy_cpts_interior


#### 4. Kernel / MMD Scan Method

# Compute Gaussian-kernel MMD scan statistic on first 5 PCs.

Z_mmd <- Z_bench
n_mmd <- nrow(Z_mmd)

D2 <- as.matrix(dist(Z_mmd))^2

# Median heuristic for Gaussian kernel bandwidth
sigma2 <- median(D2[D2 > 0])

K <- exp(-D2 / (2 * sigma2))

mmd_stat <- function(t, K) {
  idx1 <- 1:t
  idx2 <- (t + 1):nrow(K)
  
  K11 <- K[idx1, idx1, drop = FALSE]
  K22 <- K[idx2, idx2, drop = FALSE]
  K12 <- K[idx1, idx2, drop = FALSE]
  
  mean(K11) + mean(K22) - 2 * mean(K12)
}

# Avoid very small segments
candidate_t <- 30:(n_mmd - 30)
mmd_vals <- sapply(candidate_t, mmd_stat, K = K)
mmd_tau <- candidate_t[which.max(mmd_vals)]
mmd_tau


#### 5. Benchmark Summary Table

benchmark_results <- data.frame(
  method = c(
    "PELT",
    "Binary Segmentation",
    "Energy",
    "Kernel / MMD"
  ),
  input = c(
    "PC1",
    "PC1",
    "First 5 PCs",
    "First 5 PCs"
  ),
  result = c(
    paste0(
      "Over-segmented; ",
      length(pelt_cpts),
      " changepoints with 10*log(n) penalty"
    ),
    paste(binseg_cpts, collapse = ", "),
    paste(energy_cpts_interior, collapse = ", "),
    as.character(mmd_tau)
  )
)

benchmark_results


################################################
# PART 3 — GRAPH-BASED CHANGE POINT CONSTRUCTION
################################################

# 3.1 Define Original Feature-Based Distance

# Compute correlation between time points in PCA space
C <- cor(t(Z))

# Convert correlation to a distance measure
D_orig <- 1 - C

# Set diagonal to zero
diag(D_orig) <- 0

# Check dimensions of original distance matrix
dim(D_orig)

# 3.2 Define Temporal Distance

# Number of time points
n <- nrow(Z)

# Time index
time_index <- 1:n

# Compute normalized temporal separation |i - j| / n
D_time <- as.matrix(dist(time_index, method = "manhattan")) / n

# Check dimensions of temporal distance matrix
dim(D_time)

# 3.3 Define Temporally Weighted Distance

# Choose temporal weighting parameter
lambda <- 0.25

# Combine feature distance with temporal penalty
D_temp <- D_orig + lambda * D_time

# Set diagonal to zero
diag(D_temp) <- 0

# Check dimensions of temporally weighted distance matrix
dim(D_temp)

summary(as.vector(D_orig)) #check [0,2]
summary(as.vector(D_time)) #check [0,1]
summary(as.vector(D_temp))

#### 3.4 Build mst from D_orig and D_temp

#construct MST from original feature-based distance
g_orig <- build_mst_from_dist(D_orig)

#construct MST from temporally weighted distance
g_temp <- build_mst_from_dist(D_temp)

#check graph size is 300 vertices and 299 edges
vcount(g_orig)
ecount(g_orig)
vcount(g_temp)
ecount(g_temp)

#### 3.5 Convert MSTs to Edge Matrices

# Convert original MST to a two-column edge matrix
E_orig <- as_edgelist(g_orig, names = FALSE)
E_orig <- matrix(as.integer(E_orig), ncol = 2)

# Convert temporally weighted MST to a two-column edge matrix
E_temp <- as_edgelist(g_temp, names = FALSE)
E_temp <- matrix(as.integer(E_temp), ncol = 2)

# Inspect edge matrices are 299 x 2
head(E_orig)
dim(E_orig)
head(E_temp)
dim(E_temp)

#### 3.6 Run Graph-Based Change Point Detection

#apply graph-based scan statistic to original MST
res_orig <- gseg1(n, E_orig)

#apply graph-based scan statistic to temporally weighted MST
res_temp <- gseg1(n, E_temp)

#extract estimated change-point locations
res_orig$scanZ$generalized$tauhat
res_temp$scanZ$generalized$tauhat

#extract approximate p-values for generalized scan statistic
res_orig$pval.appr$generalized
res_temp$pval.appr$generalized

E_orig_sorted <- t(apply(E_orig, 1, sort))
E_temp_sorted <- t(apply(E_temp, 1, sort))

orig_labels <- apply(E_orig_sorted, 1, paste, collapse = "-")
temp_labels <- apply(E_temp_sorted, 1, paste, collapse = "-")

sum(orig_labels %in% temp_labels)
sum(!(temp_labels %in% orig_labels))

#### 3.7 Plot Original Graph-Based Scan Statistic

# Estimated change point from original graph method
tau_orig <- res_orig$scanZ$generalized$tauhat

# Extract generalized scan statistic values
orig_stat <- res_orig$scanZ$generalized$S

# Check that scan statistic was extracted correctly
length(orig_stat)
head(orig_stat)

# Build data frame for plotting
orig_df <- data.frame(
  time = seq_along(orig_stat),
  scan_statistic = orig_stat
)

# Create original graph-based scan statistic plot
fig_orig_scan <- ggplot(orig_df, aes(x = time, y = scan_statistic)) +
  geom_line(linewidth = 1) +
  geom_vline(
    xintercept = tau_orig,
    color = "red",
    linewidth = 1
  ) +
  labs(
    title = "Original Graph-Based Scan Statistic",
    x = "Time Point",
    y = "Generalized Scan Statistic"
  )

fig_orig_scan

# Save figure
ggsave(
  filename = "figures/original_graph_scan_statistic.png",
  plot = fig_orig_scan,
  width = 7,
  height = 5,
  dpi = 300
)

# Show figure of mst
# Use a fixed layout so the saved figure is consistent
lay <- layout_with_fr(g_orig)

png(
  filename = "figures/mst_change_point_edge_trimmed.png",
  width = 900,
  height = 700,
  res = 150
)

par(mar = c(0, 0, 0, 0))

plot(
  g_orig,
  layout = lay,
  edge.color = edge_colors,
  vertex.label = NA,
  vertex.size = 4,
  vertex.color = "orange",
  vertex.frame.color = "black",
  edge.width = 1.2,
  main = "",
  margin = -0.15
)

dev.off()

change_point <- res_orig$scanZ$generalized$tauhat
print(change_point)
######################################
# PART 4 — LAMBDA SENSITIVITY ANALYSIS
######################################

#### 4.1 Define Lambda Grid

lambda_vals <- c(0, 0.1, 0.25, 0.5, 1, 2, 5, 10)

#create results table for lambda sensitivity study
results_lambda <- data.frame(
  lambda = lambda_vals,
  tauhat = NA,
  pval = NA
)

# 4.2 Compute Change Point Estimates Across Lambda

for (i in seq_along(lambda_vals)) {
  
  lambda <- lambda_vals[i]
  
  # Construct temporally weighted distance matrix
  D_lambda <- D_orig + lambda * D_time
  diag(D_lambda) <- 0
  
  # Build MST from weighted distance
  g_lambda <- build_mst_from_dist(D_lambda)
  
  # Convert MST to edge matrix
  E_lambda <- as_edgelist(g_lambda, names = FALSE)
  E_lambda <- matrix(as.integer(E_lambda), ncol = 2)
  
  # Run graph-based changepoint detection
  res_lambda <- suppressMessages(
    suppressWarnings(
      gseg1(n, E_lambda)
    )
  )
  
  # Store outputs
  results_lambda$tauhat[i] <- res_lambda$scanZ$generalized$tauhat
  results_lambda$pval[i] <- res_lambda$pval.appr$generalized
}

#### 4.3 Check Results

results_lambda

#plot lambda vs estimated chage point
ggplot(results_lambda, aes(x = lambda, y = tauhat)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_y_continuous(breaks = c(153,154)) +
  labs(
    title = "Estimated Change-Point vs Temporal Weighting",
    x = expression(lambda),
    y = expression(hat(tau))
  ) +
  theme_minimal()


#######################################
# Part 5: Connect Results to Experiment
#######################################

#### 5.1 Read Repetition Time from JSON

library(jsonlite)

#path to metadata JSON
json_path <- file.path(
  "data_raw", sub_id, "ses-01", "func",
  paste0(sub_id, "_ses-01_task-", task_name, "_run-01_bold.json")
)

#read metadata
scan_info <- fromJSON(json_path)

#extract TR
TR <- scan_info$RepetitionTime
TR

#convert tau into seconds
tau_seconds <- (154 - 1) * TR
tau_seconds

events <- read.delim(events_path)
events[
  tau_seconds >= events$onset &
    tau_seconds < (events$onset + events$duration),
]

nofeedback_events


##############################################
# PART — CONNECT CHANGE POINTS TO EVENTS FILE
##############################################



TR <- 2
bold_delay <- 6

events <- read.delim(events_path)

cp_context <- data.frame(
  estimate = c("Original graph", "Weighted graph", "Stabilized weighted"),
  tauhat = c(267, 150, 154)
) %>%
  mutate(
    scan_time_sec = (tauhat - 1) * TR,
    delay_adjusted_sec = scan_time_sec - bold_delay
  )

cp_event_context <- cp_context %>%
  rowwise() %>%
  mutate(
    event_row = list(
      events %>%
        filter(
          delay_adjusted_sec >= onset,
          delay_adjusted_sec < onset + duration
        )
    )
  ) %>%
  tidyr::unnest(event_row)

cp_event_context

##############################################
# Part 6 - Compare Three Tasks For One Subject
##############################################

task_curves <- data.frame(
  lambda = rep(c(0,0.1,0.25,0.5,1,2,5,10), 3),
  tauhat = c(
    c(267,267,210,153,153,151,150,150),   # nofeedback
    c(70,70,70,195,151,150,150,150),      # visual
    c(162,162,139,139,145,145,150,150)    # auditory
  ),
  task = rep(c("No Feedback","Visual","Auditory"), each = 8)
)

ggplot(task_curves, aes(x = lambda, y = tauhat, color = task)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "Sub-01: Change-Point vs Temporal Weighting Across Tasks",
    x = expression(lambda),
    y = expression(hat(tau)),
    color = "Task"
  ) +
  theme_minimal()


#########################################
# PART 7 — MULTI-SUBJECT ANALYSIS BY TASK
#########################################

# Tasks to compare
task_list <- c("nofeedback", "visualfeedback", "auditoryfeedback")

# Analysis settings
k_multi <- 220
lambda_multi <- 5

# Empty table to store all results
results_all_tasks <- data.frame()

for (task_cur in task_list) {
  
  for (sid in subject_ids) {
    
    cat("Processing", sid, "for task", task_cur, "\n")
    
    subj_result <- tryCatch(
      analyze_subject_graph_cp(
        sub_id = sid,
        task_name = task_cur,
        k = k_multi,
        lambda = lambda_multi
      ),
      error = function(e) {
        message("Skipping ", sid, " for task ", task_cur, ": ", e$message)
        return(NULL)
      }
    )
    
    if (!is.null(subj_result)) {
      results_all_tasks <- rbind(results_all_tasks, subj_result)
    }
  }
}

results_all_tasks
saveRDS(results_all_tasks, "data_work/results_all_tasks.rds")
View(results_all_tasks)


results_multi <- results_all_tasks %>%
  filter(task == "nofeedback")

##############################################
# PART 8B — CROSS-TASK CHANGE POINT FIGURE
##############################################

# Clean task labels for plotting
results_all_tasks$task_label <- factor(
  results_all_tasks$task,
  levels = c("nofeedback", "visualfeedback", "auditoryfeedback"),
  labels = c("No Feedback", "Visual Feedback", "Auditory Feedback")
)

# Boxplot of all tasks for each subject
plot_all_tasks <- results_all_tasks %>%
  select(subject, task, tau_orig, tau_temp) %>%
  pivot_longer(
    cols = c(tau_orig, tau_temp),
    names_to = "method",
    values_to = "tau"
  ) %>%
  mutate(
    method = recode(
      method,
      tau_orig = "Original",
      tau_temp = "Weighted"
    ),
    task = factor(
      task,
      levels = c("nofeedback", "visualfeedback", "auditoryfeedback"),
      labels = c("No Feedback", "Visual Feedback", "Auditory Feedback")
    )
  )

fig_all_tasks <- ggplot(plot_all_tasks, aes(x = method, y = tau)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.7, size = 2) +
  geom_hline(yintercept = 150, linetype = "dashed") +
  facet_wrap(~ task, nrow = 1) +
  labs(
    title = "Graph-Based Change-Point Estimates Across Subjects and Tasks",
    x = "Method",
    y = expression(hat(tau))
  ) +
  theme_bw(base_size = 14)

fig_all_tasks

ggsave(
  filename = "figures/all_tasks_method_comparison.png",
  plot = fig_all_tasks,
  width = 11,
  height = 4.5,
  dpi = 300
)

##############################################
# PART 9 — NO-FEEDBACK MULTI-SUBJECT COMPARISON
##############################################

# Paired differences
results_multi$diff_tau <- results_multi$tau_temp - results_multi$tau_orig

# Basic summaries
summary(results_multi$tau_orig)
summary(results_multi$tau_temp)
summary(results_multi$diff_tau)

mean(results_multi$tau_orig)
mean(results_multi$tau_temp)
mean(results_multi$diff_tau)

sd(results_multi$tau_orig)
sd(results_multi$tau_temp)
sd(results_multi$diff_tau)

median(results_multi$tau_orig)
median(results_multi$tau_temp)
median(results_multi$diff_tau)

IQR(results_multi$tau_orig)
IQR(results_multi$tau_temp)
IQR(results_multi$diff_tau)


# Variance comparison
var.test(results_multi$tau_orig, results_multi$tau_temp)

##################
# Part 10 - Figures
##################

# Figure data preperation
results_long <- results_multi %>%
  select(subject, tau_orig, tau_temp) %>%
  pivot_longer(
    cols = c(tau_orig, tau_temp),
    names_to = "method",
    values_to = "tauhat"
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c("tau_orig", "tau_temp"),
      labels = c("Original", "Temporally weighted")
    )
  )
results_multi$tau_temp_factor <- factor(
  results_multi$tau_temp,
  levels = 147:153
)

# Boxplot
fig_boxplot <- ggplot(results_long, aes(x = method, y = tauhat)) +
  geom_boxplot() +
  labs(
    title = "Estimated Change Points Across Subjects",
    x = "Method",
    y = expression(hat(tau))
  )

fig_boxplot

ggsave(
  filename = "figures/multi_subject_boxplot.png",
  plot = fig_boxplot,
  width = 7,
  height = 5,
  dpi = 300
)

# Scatterplot
fig_scatter <- ggplot(results_multi, aes(x = tau_orig, y = tau_temp)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(
    title = "Original vs Temporally Weighted Change Point Estimates",
    x = expression(hat(tau)[orig]),
    y = expression(hat(tau)[temp])
  )

fig_scatter

ggsave(
  filename = "figures/orig_vs_weighted_scatter.png",
  plot = fig_scatter,
  width = 7,
  height = 5,
  dpi = 300
)

# Bar chart of weighted estimates
fig_weighted_bar <- ggplot(results_multi, aes(x = tau_temp_factor)) +
  geom_bar() +
  labs(
    title = "Frequency of Weighted Change Point Estimates",
    x = expression(hat(tau)[temp]),
    y = "Count"
  )

fig_weighted_bar

ggsave(
  filename = "figures/weighted_cp_frequency.png",
  plot = fig_weighted_bar,
  width = 7,
  height = 5,
  dpi = 300
)

# Brain image
image(nii, z = 15, plot.type = "single")
png("brain_slice.png", width = 600, height = 600)
image(nii, z = 15, plot.type = "single")
dev.off()

#########################################
# OPTIONAL: REPRESENTATIVE SUBJECT SEARCH
#########################################
# This section was used to identify subjects
# near cohort central tendency. It is not
# required for the main report results.

median_temp <- median(results_multi$tau_temp)

results_multi$dist_to_med_temp <-
  abs(results_multi$tau_temp - median_temp)

results_multi[order(results_multi$dist_to_med_temp), ]

