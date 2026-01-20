# ======================================================================
# == load libraries ==
# ======================================================================

library(tidyverse)

# ======================================================================
# == Snakemake integation ==
# ======================================================================

# define the meta tables
# meta_tables = expand(rules.runSamplePreprocessing.output.meta,sample_name=SAMPLES.keys())
# meta_tables <- 
meta_tables <- snakemake@input$meta_tables

message("input metadata: ", meta_tables)

# define the output
# out_id_object <-  "results/Seurat/object/01_connect_5k_pbmc_NGSC3_ch1_gex_1_obj_preQC.rds"
# out_id_meta <- "results/Seurat/table/01_connect_5k_pbmc_NGSC3_ch1_gex_1_meta_preQC.rds"
out_plot_mito <- snakemake@output$plot_mito
out_plot_feature <- snakemake@output$plot_feature
out_LUT_QC <- snakemake@output$LUT_QC
out_meta_total <- snakemake@output$meta_total

message("output plot mito: ", out_plot_mito)
message("output plot feature: ", out_plot_feature)
message("output table meta total before QC: ", out_LUT_QC)
message("output table LUT QC: ", out_LUT_QC)

# ======================================================================
# == load input ==
# ======================================================================

# define the basename
file <- basename(meta_tables)
message("input sample names: ", file)

sample_name <- file %>%
  str_remove(pattern = "_meta_preQC.tsv")

meta_total <- pmap(list(meta_tables,sample_name),function(x,x_name){
  read_tsv(x) %>%
    mutate(dataset = x_name)
}) %>%
  bind_rows()

# define some plotting parameters
# define number of rows in the panel
test <- length(sample_name)

message("number of samples processed: ", test)

# test <- 10
panel_row <- test %>% sqrt() %>% round(digits = 0)
panel_col <- (test / panel_row) %>% ceiling()

# ======================================================================
# == run the standard processing ==
# ======================================================================

# print head of the full metadata and its dimensions
head(meta_total)
dim(meta_total)

# plot mito prop
p_mito <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "percent.mt") %>%
  ggplot(aes(x=value))+geom_histogram(binwidth = 0.05) +
  facet_wrap(orig.ident~var,scales = "free",ncol = panel_col) +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = c(0,1,2,5,10,20,40,60,80,100)) +
  geom_vline(xintercept = c(10),col="red",linetype="dashed") +
  annotate("rect", xmin=0, xmax=10, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank())

# plot features 
p_feature <- meta_total %>%
  gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.globin)) %>%
  filter(var == "nFeature_RNA") %>%
  ggplot(aes(x=value)) + 
  geom_histogram(binwidth = 0.05) +
  facet_wrap(orig.ident~var,scales = "free",ncol = panel_col) +
  theme_bw() +
  scale_x_continuous(trans = "log1p",breaks = c(200,500,1000,2000,4000,6000,8000,10000,20000)) +
  geom_vline(xintercept = c(1000,6000),col="red",linetype="dashed") +
  annotate("rect", xmin=1000, xmax=6000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))

# define the table to fill the QC parameters
# after the lates FG we decided to define some all inclusive aprameters
# LUT_df <- data.frame(sample_name = sample_name,
#                      featureLow_thr = NA,
#                      featureHigh_thr = NA,
#                      mito_thr = NA)

LUT_df <- data.frame(sample_name = sample_name,
                     featureLow_thr = 0,
                     featureHigh_thr = Inf,
                     mito_thr = Inf)

# ======================================================================
# == save output ==
# ======================================================================

# save plot histo mito prop reads
ggsave(plot = p_mito,
       out_plot_mito,
       width = 4 * panel_col,
       height = 4 * panel_row)

# save plot histo feature reads
ggsave(plot = p_feature,
       out_plot_feature,
       width = 4 * panel_col,
       height = 4 * panel_row)

# save tables
write_tsv(meta_total,out_meta_total)
write_csv(LUT_df,out_LUT_QC)

