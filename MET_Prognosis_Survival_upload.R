#MET Prognosis
#Analysis of the data using a negative binomial survival model

#Set working directory
#Loading in data
df <- read.csv("prognosis_data_cleaned.csv", stringsAsFactors = FALSE)

#Libraries
library(rstan)
library(loo)
library(rethinking)
library(tidyr)
library(cowplot)
library(ggplot2)
library(ggdist)
library(gtsummary)
library(gt)
library(dplyr)

# Stan optimization
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

HDILow<- function(x, HDI=0.90) {
  sortedPts = sort( x)
  ciIdxInc = ceiling( HDI * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc 
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ] }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ] 
  return( HDImin)
}

HDIHigh<- function(x, HDI=0.90) {
  sortedPts = sort( x)
  ciIdxInc = ceiling( HDI * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc 
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ] }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ] 
  return( HDImax)
}

# Reproducibility
set.seed(123)

  ####Setting up the data####
colnames(df)
df2 <- df[c(3:5, 8, 60, 63, 64)]

#Removing participants who failed due to msk injury
df2 <- df[c(3:5, 8, 10,60, 63, 64)]
df3 <- df2[!(df2$msk_injury_fail=="fail_msk_injury"),]

#Creating a censored column
df3$censored <- ifelse(is.na(df3$length_of_recovery), 1, 0)

#Adding in the last known day that censored values were still injured
df3$clinic_mc_date_clean <- as.Date(df3$clinic_mc_date_clean)
df3$date_conc <- as.Date(df3$date_conc)

df3$length_of_recovery[df3$censored == 1] <- as.numeric(df3$clinic_mc_date_clean[df3$censored == 1] - 
               df3$date_conc[df3$censored == 1])

df_analysis <- df3[c(4, 5, 7:9)]

#####Survival NBM Model#####
  ####Loading the model####
survival2_nbm <- stan_model("survival2_nbm.stan")

  ####Setting up the data####
met_data <- list(
  N = nrow(df_analysis),
  G = 2,
  y = df_analysis$length_of_recovery,
  censored =  df_analysis$censored,
  group = df_analysis$met_pass_fail
)

str(met_data)

  ####Running model####
m_met2 <- sampling(survival2_nbm, data = met_data,
                  iter = 2000, chains = 4, seed = 987)
precis(m_met2,depth=2)
saveRDS(m_met2, "m_met2.rds")
m_met2 <- readRDS("m_met2.rds")

#Extracting samples
post_met <- rstan::extract(m_met2) 
days <- 0:185
  ####Creating Posterior plots####
#Summarize posterior survival per group
summarize_survival_base <- function(R_mat) {
  day <- 0:(ncol(R_mat) - 1)
  mean <- apply(R_mat, 2, mean)
  lower <- apply(R_mat, 2, quantile, probs = 0.05)
  upper <- apply(R_mat, 2, quantile, probs = 0.95)
  data.frame(day = day, mean = mean, lower = lower, upper = upper)
}

#Putting into an array
R_array <- list(100-(post_met$R_group1 * 100),(100-(post_met$R_group2*100)))
R_group1 <- (100-(post_met$R_group1 * 100))
R_group2 <- (100-(post_met$R_group2 * 100))
R_group <- array(NA,dim=c(4000,2,201))
R_group[,1,] <- R_group1
R_group[,2,] <- R_group2
G <- dim(R_group)[2]

#Creating colours and labels
group_colors <- rainbow(G)
group_labels <-c("Fail", "Pass")

#Summarizing ribbons for each group
summary_list <- lapply(1:G, function(g) {
  R_g <- R_array[[g]]
  summarize_survival_base(R_g)
})

#Plot setup
plot(NULL, xlim = c(0, 200), ylim = c(0, 100),
     xlab = "Days", ylab = "Probability Recovered by Time t (%)",
     main = "Posterior Predictive Recovery Curves by MET Result")

# Plot ribbons and medians
for (g in 1:G) {
  df <- summary_list[[g]]
  polygon(c(df$day, rev(df$day)),
          c(df$lower, rev(df$upper)),
          col = adjustcolor(group_colors[g], alpha.f = 0.3), border = NA)
  lines(df$day, df$mean, col = group_colors[g], lwd = 2)
  abline(v = c(28), col = "black", lty = 2, lwd = 2)
}

legend("bottomright", legend = group_labels, col = group_colors, lwd = 2, bty = "n")

  ####Posterior Predictive Checks####
#Empirical survival for uncensored observations
compute_empirical_survival <- function(y) {
  sapply(0:185, function(d) sum(d > y) / length(y))
}

#New plot
plot(NULL, xlim = c(0, 185), ylim = c(0, 100),
     xlab = "Days", ylab = "Probability Recovered by Time t (%)",
     main = "Posterior Survival Ribbons + Observed Data")

#Looping over groups
for (g in 1:G) {
  df <- summary_list[[g]]
  
  # Posterior
  polygon(c(df$day, rev(df$day)),
          c(df$lower, rev(df$upper)),
          col = adjustcolor(group_colors[g], alpha.f = 0.3), border = NA)
  lines(df$day, df$mean, col = group_colors[g], lwd = 2)
  
  # Observed (dots for uncensored recovery)
  y_obs <- df_analysis$length_of_recovery[df_analysis$met_pass_fail == g & 
                                            df_analysis$censored == 0]
  if (length(y_obs) > 0) {
    emp <- compute_empirical_survival(y_obs)
    points(0:185, emp * 100, col = group_colors[g], pch = 16)
  }
}

legend("bottomright",
       legend = c(c("Fail (Model)", "Pass (Model)"),
                  c("Fail (Observed)", "Pass (Observed)")),
       col = rep(group_colors, 2),
       pch = c(rep(NA, G), rep(16, G)),
       lty = c(rep(1, G), rep(NA, G)),
       lwd = c(rep(2, G), rep(NA, G)),
       bty = "n")

  ####Checking values####
df_fail_rec <- summary_list[[1]] #Fail
df_pass_rec <- summary_list[[2]] #Pass

group_df <- as.data.frame(rstan::extract(m_met2))
group_df$fail <- exp(group_df$alpha.1)
group_df$pass <- exp(group_df$alpha.2)
group_df$contrast_fp <- (group_df$fail-group_df$pass)
apply(group_df[c(531:533)],2,mean)
apply(group_df[c(531:533)],2,HDILow)
apply(group_df[c(531:533)],2,HDIHigh)
length(which(group_df$contrast_fp > 0))/length(group_df$contrast_fp)

  ####Extracting survival predictions#####
post_group  <- rstan::extract(m_met2)
y_sim_group <- post_group$y_sim

#Empirical cumulative recovery (only uncensored data)
y_obs <- met_data$y[met_data$censored == 0]
empirical_cum <- sapply(0:185, function(d) mean(y_obs <= d))
df_empirical <- data.frame(day = 0:185, recovery = empirical_cum * 100)

#Posterior predictive cumulative recovery from y_sim
get_cum_recovery_draws <- function(y_sim) {
  sapply(0:185, function(d) rowMeans(y_sim <= d))
}
group_cum_mat <- get_cum_recovery_draws(y_sim_group) * 100

#Summarizing posterior predictions
summarize_curve <- function(cum_mat) {
  data.frame(
    day = 0:185,
    median = apply(cum_mat, 2, median),
    lower  = apply(cum_mat, 2, quantile, 0.05),
    upper  = apply(cum_mat, 2, quantile, 0.95)
  )
}
df_group_cum  <- summarize_curve(group_cum_mat)

#Cumulative Recovery Curves
plot(NULL, xlim = c(0, 185), ylim = c(0, 100),
     xlab = "Days Since Injury", ylab = "Probability Recovered by Time t (%)",
     main = "Cumulative Recovery: Posterior vs. Observed")

#Posterior ribbons
polygon(c(df_group_cum$day, rev(df_group_cum$day)),
        c(df_group_cum$lower, rev(df_group_cum$upper)),
        col = adjustcolor("blue", alpha.f = 0.2), border = NA)
lines(df_group_cum$day, df_group_cum$median, col = "blue", lwd = 2)

#Empirical cumulative recovery
lines(df_empirical$day, df_empirical$recovery, col = "black", lwd = 2, lty = 2)

legend("bottomright", legend = c("MET Model", "Empirical"),
       col = c("blue", "black"), lwd = 2, lty = c(1, 1, 2), bty = "n")

####Graph for Manuscript####
colnames(group_df)
plot_met_post_px <- rbind.data.frame(gather(group_df[c(531:533)]))
plot_met_post_px$stage <- c(rep("MET Fail",4000),rep("MET Pass",4000),rep("Contrast (Difference)",4000))
plot_met_post_px$group <- c(rep("MET Fail",4000),rep("MET Pass",4000),rep("Contrast (Difference)",4000))

fig_met_post_px <- plot_met_post_px %>%
  ggplot(aes(y = stage,x = value,fill = group)) + 
  theme_classic()+
  theme(strip.text.x = element_text(face = 'bold',size = 10),
        axis.title.x = element_text(size = 9,face = 'bold'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = 'bold',size = 9),
        legend.position = 'none',
        axis.text.x = element_text(face = 'bold',size = 8),
        title = element_text(face = 'bold'))+
  scale_fill_manual("Group",values = c("midnightblue","#FF0000","#00FFFF")) +xlim(-40,120)+
  stat_halfeye(.width = c(0.70,0.90),alpha = 0.9) +
  xlab("Days to Recovery")
fig_met_post_px

fig_met_post_px_2 <- plot_met_post_px %>%
  ggplot(aes(x = value,fill = group)) + 
  theme_classic()+
  theme(strip.text.x = element_text(face = 'bold',size = 10),
        axis.title.x = element_text(size = 9,face = 'bold'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = 'bold',size = 9),
        legend.position = 'none',
        axis.text.x = element_text(face = 'bold',size = 8),
        title = element_text(face = 'bold'))+
  geom_histogram(position = "identity",
                 alpha = 0.8,
                 bins = 500) +
  facet_wrap(~ stage, ncol = 1, scales = "free_y") +
  scale_fill_manual("Group",values = c("midnightblue","#FF0000","#00FFFF")) +xlim(-20,100)+
  xlab("Days to Recovery")
fig_met_post_px_2

#Using the below figure for manuscript
fig_met_post_px_3 <- plot_met_post_px %>%
  ggplot(aes(y = stage,x = value,fill = group)) +
  stat_histinterval(orientation = "horizontal",breaks = 200,                    
    scale = 1,alpha = 0.9,slab_color = NA) +
  theme_classic() +
  theme(strip.text.x = element_text(face = 'bold', size = 10),
    axis.title.x = element_text(size = 9, face = 'bold'),
    axis.title.y = element_blank(),
    axis.text.y = element_text(face = 'bold', size = 9),
    legend.position = 'none',
    axis.text.x = element_text(face = 'bold', size = 8),
    title = element_text(face = 'bold')) +
  scale_fill_manual("Group",values = c("midnightblue", "#FF0000", "#6699FF")) +
  xlim(-20, 100) + geom_vline(xintercept = 0, colour = "black", linetype = 2) +
  xlab("Days to Recovery")
fig_met_post_px_3

#Second figure for manuscript
post <- rstan::extract(m_met2)
R_array <- list(100-(post$R_group1 * 100),(100-(post$R_group2*100)))
days <- 0:200
G <- length(R_array)

#Choosing colours and labels
group_colors <- c("#FF0000", "#6699FF")
group_labels <-c("Fail", "Pass")

#Summarizing ribbons for each group
summary_list <- lapply(1:G, function(g) {
  R_g <- R_array[[g]]
  summarize_survival_base(R_g)
})
plot(NULL, xlim = c(0, 185), ylim = c(0, 100),
       xlab = "Days", ylab = "Probability Recovered by Time t (%)")
  
#Plotting ribbons and medians
for (g in 1:G) {
    df <- summary_list[[g]]
    polygon(c(df$day, rev(df$day)),
            c(df$lower, rev(df$upper)),
            col = adjustcolor(group_colors[g], alpha.f = 0.3), border = NA)
    lines(df$day, df$mean, col = group_colors[g], lwd = 2)
}
  
legend("bottomright", legend = group_labels, col = group_colors, lwd = 2, bty = "n")

plot_data <- do.call(rbind, lapply(1:G, function(g) {
  summary_list[[g]] %>% mutate(Group = paste0("Group ", g))
}))

plot_data$Group <- ifelse(plot_data$Group == "Group 1","Fail","Pass")

#Creating a ggplot2 plot instead
rec_curve_plot <- ggplot(plot_data, aes(x = day)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Group), alpha = 0.3) +
  geom_line(aes(y = mean, color = Group), linewidth = 1) +
  theme_classic() +
  theme(strip.text.x = element_text(face = 'bold',size = 10),
        axis.title.x = element_text(size = 9,face = 'bold'),
        axis.text.y = element_text(face = 'bold',size = 9),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.position = c(0.95,.1),
        legend.spacing.y = unit(0, "pt"),
        legend.spacing.x = unit(0, "pt"),
        legend.box.margin=margin(-0.2,-0.2,-0.2,-0.2),
        axis.text.x = element_text(face = 'bold',size = 8),
        title = element_text(face = 'bold'))+
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  xlab("Days to Recovery") + ylab("Probability Recovered (%)")
rec_curve_plot

#For publication
progn_fig  <- plot_grid(fig_met_post_px_3,rec_curve_plot,labels = "AUTO",ncol=2)
progn_fig
ggsave("prognosis_fig.jpg",progn_fig,dpi=600)

