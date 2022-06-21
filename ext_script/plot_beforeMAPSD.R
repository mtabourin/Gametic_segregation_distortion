#!/bin/R

library(ggplot2)
library(reshape)
library(plyr)

args = commandArgs(trailingOnly=TRUE)

input_name = paste0('/shared/projects/gametic_segregation_distortion/results/table/', args[1], '.table.out.toPlot.txt')
output_name = paste0('/shared/projects/gametic_segregation_distortion/results/plots/', args[1], '.table.out.toPlot.pdf')

df = read.table(input_name, header = T)

df = melt(df, id = c("scaffold", "mean_window"))

df_1 <- df[ df$variable == c("somatic_mean", "germline_mean"), ]
df_2 <- df[ df$variable == c("baseline", "difference_mean"), ]

df_1$variable <- revalue(df_1$variable, c("somatic_mean" = "Somatic", "germline_mean" = "Germline"))
df_2$variable <- revalue(df_2$variable, c("baseline" = "Somatic", "difference_mean" = "Germline"))

df_plot_1 <- ggplot(df_1, aes(x = mean_window, y = value, color = variable, linetype = variable) ) +
  geom_line() + expand_limits(y=c(0.35,0.65)) +
  facet_wrap(scaffold ~ ., ncol = 4, nrow = 2, scales = 'free_x' ) +
  scale_color_manual(values = c("red", "red")) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  ylab("Raw ancestry ratio") +
  xlab("Mean window position (bp)")


df_plot_2 <- ggplot(df_2, aes(x = mean_window, y = value, color = variable, linetype = variable) ) +
  geom_line() + expand_limits(y=c(-0.15,0.10)) +
  facet_wrap(scaffold ~ ., ncol = 4, nrow = 2, scales = 'free_x' ) +
  scale_color_manual(values = c("blue", "blue")) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  ylab("Ancestry ratio difference") +
  xlab("Mean window position (bp)")

pdf(output_name)

df_plot_1
df_plot_2

dev.off()


