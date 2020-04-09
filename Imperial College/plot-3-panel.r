library(ggplot2)
library(tidyr)
library(dplyr)
library(rstan)
library(data.table)
library(lubridate)
library(gdata)
library(EnvStats)
library(matrixStats)
library(scales)
library(gridExtra)
library(ggpubr)
library(bayesplot)
library(cowplot)

source("utils/geom-stepribbon.r")
#---------------------------------------------------------------------------
make_three_pannel_plot <- function(){
  
  args <- commandArgs(trailingOnly = TRUE)
  
  filename2 <- args[1]
  load(paste0("results/", filename2))
  print(sprintf("loading: %s",paste0("results/",filename2)))
  data_interventions <- read.csv("data/interventions.csv", 
                                 stringsAsFactors = FALSE)
  covariates <- data_interventions[1:11, c(1,2,3,4,5,6, 7, 8)]
  for(i in 1:11){
    print(i)
    N <- length(dates[[i]])
    country <- countries[[i]]
    
    predicted_casos <- colMeans(prediction[,1:N,i])
    predicted_casos_li <- colQuantiles(prediction[,1:N,i], probs=.025)
    predicted_casos_ui <- colQuantiles(prediction[,1:N,i], probs=.975)
    predicted_casos_li2 <- colQuantiles(prediction[,1:N,i], probs=.25)
    predicted_casos_ui2 <- colQuantiles(prediction[,1:N,i], probs=.75)
    
    
    estimated_mortes <- colMeans(estimated.mortes[,1:N,i])
    estimated_mortes_li <- colQuantiles(estimated.mortes[,1:N,i], probs=.025)
    estimated_mortes_ui <- colQuantiles(estimated.mortes[,1:N,i], probs=.975)
    estimated_mortes_li2 <- colQuantiles(estimated.mortes[,1:N,i], probs=.25)
    estimated_mortes_ui2 <- colQuantiles(estimated.mortes[,1:N,i], probs=.75)
    
    rt <- colMeans(out$Rt[,1:N,i])
    rt_li <- colQuantiles(out$Rt[,1:N,i],probs=.025)
    rt_ui <- colQuantiles(out$Rt[,1:N,i],probs=.975)
    rt_li2 <- colQuantiles(out$Rt[,1:N,i],probs=.25)
    rt_ui2 <- colQuantiles(out$Rt[,1:N,i],probs=.75)
    
    
    # delete these 2 lines
    covariates_country <- covariates[which(covariates$Country == country), 2:8]   
    
    # Remove sport
    covariates_country$sport = NULL 
    covariates_country$travel_restrictions = NULL 
    covariates_country_long <- gather(covariates_country[], key = "key", 
                                      value = "value")
    covariates_country_long$x <- rep(NULL, length(covariates_country_long$key))
    un_dates <- unique(covariates_country_long$value)
    
    for (k in 1:length(un_dates)){
      idxs <- which(covariates_country_long$value == un_dates[k])
      max_val <- round(max(rt_ui)) + 0.3
      for (j in idxs){
        covariates_country_long$x[j] <- max_val
        max_val <- max_val - 0.3
      }
    }
    
    
    covariates_country_long$value <- as_date(covariates_country_long$value) 
    covariates_country_long$country <- rep(country, 
                                           length(covariates_country_long$value))
    
    data_country <- data.frame("time" = as_date(as.character(dates[[i]])),
                               "country" = rep(country, length(dates[[i]])),
                               "reported_casos" = reported_casos[[i]], 
                               "reported_casos_c" = cumsum(reported_casos[[i]]), 
                               "predicted_casos_c" = cumsum(predicted_casos),
                               "predicted_min_c" = cumsum(predicted_casos_li),
                               "predicted_max_c" = cumsum(predicted_casos_ui),
                               "predicted_casos" = predicted_casos,
                               "predicted_min" = predicted_casos_li,
                               "predicted_max" = predicted_casos_ui,
                               "predicted_min2" = predicted_casos_li2,
                               "predicted_max2" = predicted_casos_ui2,
                               "mortes" = mortes_por_pais[[i]],
                               "mortes_c" = cumsum(mortes_por_pais[[i]]),
                               "estimated_mortes_c" =  cumsum(estimated_mortes),
                               "death_min_c" = cumsum(estimated_mortes_li),
                               "death_max_c"= cumsum(estimated_mortes_ui),
                               "estimated_mortes" = estimated_mortes,
                               "death_min" = estimated_mortes_li,
                               "death_max"= estimated_mortes_ui,
                               "death_min2" = estimated_mortes_li2,
                               "death_max2"= estimated_mortes_ui2,
                               "rt" = rt,
                               "rt_min" = rt_li,
                               "rt_max" = rt_ui,
                               "rt_min2" = rt_li2,
                               "rt_max2" = rt_ui2)
    save(data_country, country, covariates_country_long, file = "debug.RData")
    make_plots(data_country = data_country, 
               covariates_country_long = covariates_country_long,
               filename2 = filename2,
               country = country)
    
  }
}

#---------------------------------------------------------------------------
make_plots <- function(data_country, covariates_country_long, 
                       filename2, country){
  
  data_casos_95 <- data.frame(data_country$time, data_country$predicted_min, 
                              data_country$predicted_max)
  names(data_casos_95) <- c("time", "casos_min", "casos_max")
  data_casos_95$key <- rep("nintyfive", length(data_casos_95$time))
  data_casos_50 <- data.frame(data_country$time, data_country$predicted_min2, 
                              data_country$predicted_max2)
  names(data_casos_50) <- c("time", "casos_min", "casos_max")
  data_casos_50$key <- rep("fifty", length(data_casos_50$time))
  data_casos <- rbind(data_casos_95, data_casos_50)
  levels(data_casos$key) <- c("ninetyfive", "fifty")
  
  p1 <- ggplot(data_country) +
    geom_bar(data = data_country, aes(x = time, y = reported_casos), 
             fill = "coral4", stat='identity', alpha=0.5) + 
    geom_ribbon(data = data_casos, 
                aes(x = time, ymin = casos_min, ymax = casos_max, fill = key)) +
    xlab("") +
    ylab("Numero diario de infeccoes") +
    scale_x_date(date_breaks = "semanas", labels = date_format("%e %b")) + 
    scale_fill_manual(name = "", labels = c("50%", "95%"),
                      values = c(alpha("deepskyblue4", 0.55), 
                                 alpha("deepskyblue4", 0.45))) + 
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "None") + 
    guides(fill=guide_legend(ncol=1))
  
  data_mortes_95 <- data.frame(data_country$time, data_country$death_min, 
                               data_country$death_max)
  names(data_mortes_95) <- c("time", "death_min", "death_max")
  data_mortes_95$key <- rep("nintyfive", length(data_mortes_95$time))
  data_mortes_50 <- data.frame(data_country$time, data_country$death_min2, 
                               data_country$death_max2)
  names(data_mortes_50) <- c("time", "death_min", "death_max")
  data_mortes_50$key <- rep("fifty", length(data_mortes_50$time))
  data_mortes <- rbind(data_mortes_95, data_mortes_50)
  levels(data_mortes$key) <- c("ninetyfive", "fifty")
  
  
  p2 <-   ggplot(data_country, aes(x = time)) +
    geom_bar(data = data_country, aes(y = mortes, fill = "reported"),
             fill = "coral4", stat='identity', alpha=0.5) +
    geom_ribbon(
      data = data_mortes,
      aes(ymin = death_min, ymax = death_max, fill = key)) +
    scale_x_date(date_breaks = "semanas", labels = date_format("%e %b")) +
    scale_fill_manual(name = "", labels = c("50%", "95%"),
                      values = c(alpha("deepskyblue4", 0.55), 
                                 alpha("deepskyblue4", 0.45))) + 
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.position = "None") + 
    guides(fill=guide_legend(ncol=1))
  
  
  plot_labels <- c("Bloqueio completo", 
                   "Eventos públicos proibidos",
                   "Fechamento de escolas",
                   "Autoisolamento",
                   "Distanciamento social")
  
  # Plotting interventions
  data_rt_95 <- data.frame(data_country$time, 
                           data_country$rt_min, data_country$rt_max)
  names(data_rt_95) <- c("time", "rt_min", "rt_max")
  data_rt_95$key <- rep("nintyfive", length(data_rt_95$time))
  data_rt_50 <- data.frame(data_country$time, data_country$rt_min2, 
                           data_country$rt_max2)
  names(data_rt_50) <- c("time", "rt_min", "rt_max")
  data_rt_50$key <- rep("fifty", length(data_rt_50$time))
  data_rt <- rbind(data_rt_95, data_rt_50)
  levels(data_rt$key) <- c("ninetyfive", "fifth")
  
  p3 <- ggplot(data_country) +
    geom_stepribbon(data = data_rt, aes(x = time, ymin = rt_min, ymax = rt_max, 
                                        group = key,
                                        fill = key)) +
    geom_hline(yintercept = 1, color = 'black', size = 0.1) + 
    geom_segment(data = covariates_country_long,
                 aes(x = value, y = 0, xend = value, yend = max(x)), 
                 linetype = "dashed", colour = "grey", alpha = 0.75) +
    geom_point(data = covariates_country_long, aes(x = value, 
                                                   y = x, 
                                                   group = key, 
                                                   shape = key, 
                                                   col = key), size = 2) +
    xlab("") +
    ylab(expression(R[t])) +
    scale_fill_manual(name = "", labels = c("50%", "95%"),
                      values = c(alpha("seagreen", 0.75), alpha("seagreen", 0.5))) + 
    scale_shape_manual(name = "Intervencoes", labels = plot_labels,
                       values = c(21, 22, 23, 24, 25, 12)) + 
    scale_colour_discrete(name = "Intervencoes", labels = plot_labels) + 
    scale_x_date(date_breaks = "semanas", labels = date_format("%e %b"), 
                 limits = c(data_country$time[1], 
                            data_country$time[length(data_country$time)])) + 
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position="right")
  save(p1,p2,p3, file="debug2.RData")
  p <- plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(1, 1, 2))
  save_plot(filename = paste0("figures/", country, "_three_pannel_", filename2, ".pdf"), 
            p, base_width = 14)
}

#-----------------------------------------------------------------------------------------------
#filename <- "base-joint-1236305.pbs.Rdata"
# make_three_pannel_plot(filename)

make_three_pannel_plot()
