#!/bin/R

library(ggplot2)
library(reshape)
library(plyr)

args = commandArgs(trailingOnly=TRUE)

main = paste0('/shared/projects/gametic_segregation_distortion/results/table/', args[1], '_MAPSD.main.txt')
bootstrap = paste0('/shared/projects/gametic_segregation_distortion/results/table/', args[1], '_MAPSD.bootstrap.txt')
windows = paste0('/shared/projects/gametic_segregation_distortion/results/table/', args[1], '_MAPSD.windows.txt')

output = paste0('/shared/projects/gametic_segregation_distortion/results/plots/', args[1], '_likelihood.pdf')


df_windows = read.table(windows, header = F)
colnames(df_windows) <- c('scaffold', 'mean_window', 'somatic', 'germline', 'lnL_somatic', 'lnL_germline', 'lnL_diff')


df_main = read.table(main, header = F)
colnames(df_main) <- c('scaffold', 'position', 'somatic', 'germline', 'lnL_somatic', 'lnL_germline')
df_main$lnL_diff = - 2 * ( df_main$lnL_somatic - df_main$lnL_germline )


df_bootstrap = read.table(bootstrap, header = F)
colnames(df_bootstrap) <- c('scaffold', 'num_bootstrap', 'position', 'somatic', 'germline', 'lnL_somatic', 'lnL_germline')
inf_position_bootstrap <- ddply(df_bootstrap, .(scaffold), summarize, inf = as.numeric(quantile(position, c(0.05))))
sup_position_bootstrap <- ddply(df_bootstrap, .(scaffold), summarize, sup = as.numeric(quantile(position, c(0.95))))


df_windows  = melt(df_windows , id = c('scaffold', 'mean_window'))

df_windows_1 <- df_windows[ df_windows$variable == c("somatic", "germline"), ]

df_windows_plot_1 <- ggplot( df_windows_1, aes(x = mean_window, y = value, color = variable, linetype = variable) ) +
	geom_line() +
	facet_grid(scaffold ~ ., scales = 'free_x' ) +
	scale_color_manual(values = c("red", "red")) +
	scale_linetype_manual(values = c("dotted", "solid")) +
	ylab("Raw ancestry ratio") +
	xlab("Mean window position (bp)")


df_windows_2 <- df_windows[ df_windows$variable == c("lnL_diff"), ]
df_windows_2$value <- - 2 * df_windows_2$value

df_windows_plot_2 <- ggplot( df_windows_2 ) +
	geom_line( aes(x = mean_window, y = value) ) +
	geom_point( data = df_main, aes(x = position, y = lnL_diff), color = 'red' ) +
	geom_vline( data = inf_position_bootstrap, aes(xintercept = inf), color = 'green') +
	geom_vline( data = sup_position_bootstrap, aes(xintercept = sup), color = 'green') +
	facet_grid(scaffold ~ ., scales = 'free_y' ) 


pdf(output)

	df_windows_plot_2
	df_windows_plot_1

dev.off()

