## NEW RECALL
## load libraries
library(dplyr)
library(ggplot2)
library(colorspace)
library(data.table)
library(tidyr)

## INPUTS
methods = snakemake@params$methods %>% strsplit(" ") %>% unlist()
ER_curve_tables = snakemake@params$ERCurve_tables %>% strsplit(" ") %>% unlist()
# color palette
cpFilePlotting = snakemake@input$colorPalette
cpList = readRDS(cpFilePlotting)
# methods config
methods_config = read.table(file=snakemake@params$methods_config, header=TRUE, sep="\t")
# output
out_file = snakemake@output$er_combined

### FORMAT DATA
## aggregate and process ERCurve tables
# process first table
temp = read.table(file = ER_curve_tables[1],header = TRUE, fill = TRUE) %>% drop_na()
methods_temp = dplyr::filter(methods_config, method==methods[1])
temp$method = methods[1]
temp$pred_name_long = methods_temp$pred_name_long[1]
temp$nPoints = nrow(temp)
ER_curve_all = temp

if (length(methods) > 1) {
  for (i in 2:length(methods)) {
    temp = read.table(file = ER_curve_tables[i],header = TRUE, fill = TRUE) %>% drop_na()
    methods_temp = dplyr::filter(methods_config, method==methods[i])
    temp$method = methods[i]
    temp$pred_name_long = methods_temp$pred_name_long[1]
    temp$nPoints = nrow(temp)
    ER_curve_all = rbind(ER_curve_all, temp)
    }
}

df = ER_curve_all
df[is.na(df)] = 0
df = df[order(df$threshold),]

df_binary = dplyr::filter(df, nPoints==2, threshold==1)

df = dplyr::filter(df, nPoints!=2)


# data for binary predictors to be plotted as points


g=ggplot(data=df, aes(x=recall.LCL, y=enrichment.LCL, color=pred_name_long)) +
  geom_line(linewidth=1) +
  geom_point(data=df_binary, aes(x=recall.LCL, y=enrichment.LCL, color=pred_name_long), size=4) +
  scale_color_manual(values=unlist(cpList)) +
  ylab("Enrichment (LCL eQTLs/GM12878 predictions)") + xlab("Recall (LCL eQTLs in GM12878 predictions)") +
  labs(col="Predictor") +
  ylim(c(0, 50)) +
  theme_classic() + theme(axis.text = element_text(size = 7), axis.title = element_text(size = 8))

ggsave(out_file, g, width=8, height=6)

