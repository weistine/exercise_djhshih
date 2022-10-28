library(tidyverse)
library(data.table)
library(ggplot2)
library(patchwork)
library(ggsignif)
library(ggsci)
library(data.table)
library(plyr)

setwd("D:/Desktop/exercise-r-main/exercise-r-main")

# Q1.Mutation data analysis --------
mut_dat <- readRDS("data/mutations_sclc_ucologne_2015.rds")

## Q1.1 -------
q1.1_top10_gene <- mut_dat %>% 
  group_by(gene) %>% 
  dplyr::summarise(freq = n()) %>%
  arrange(.,desc(freq)) %>%
  .$gene %>% .[1:10]

iden_samp <- mut_dat %>%
  group_by(sample_id) %>%
  dplyr::summarise(num = n())
quan <- quantile(iden_samp$num,c(0.8,0.9))
q1.1_samp_filt <- iden_samp %>% 
  filter(num >= quan[1] & num <= quan[2])

## Q1.2 ---------
mut_ef <- fread("data/mutation_effects.tsv") %>%
  mutate(variant_class = tolower(variant_class))
mut_dat$effect <- mapvalues(tolower(mut_dat$variant_class),mut_ef$variant_class,mut_ef$effect)
q1.2_eff_sum <- mut_dat %>%
  filter(effect %in% c("loss_of_function","neutral")) %>%
  group_by(gene,effect) %>%
  dplyr::summarise(freq = n()) %>%
  pivot_wider(names_from = effect,values_from = freq,values_fill = 0)

## Q1.3 -------
# using Fisher's exact test
los <- sum(q1.2_eff_sum$loss_of_function)
non_los <- sum(q1.2_eff_sum$neutral)
test_data <- q1.2_eff_sum %>%
  transmute(glos = loss_of_function,
            nglos = los-loss_of_function,
            gneu = neutral,
            ng = non_los-neutral)
test_los <- lapply(1:nrow(test_data), function(i){
  t <- fisher.test(matrix(t(test_data[i,2:5]),ncol = 2))
  d <- test_data[i,] %>%
    dplyr::mutate(pvalue = t$p.value,
                  OR = t$estimate[[1]],
                  OR.lower95 = t$conf.int[1],
                  OR.upper95 = t$conf.int[2])
  return(d)
  }) %>% Reduce(rbind,.)
q1.3_high_los <- test_los %>% filter(pvalue < 0.05 & OR >1)

## Q1.4 -------
#assume the genes with the significantly higher proportion of loss-of-function mutations are the candidate tumour suppressor genes
pvalue <- q1.3_high_los$pvalue
names(pvalue) <- q1.3_high_los$gene
q1.3_high_los$q.value <- p.adjust(sort(pvalue),"BH")
q1.4_adj_los <- q1.3_high_los %>% transmute(.,gene = gene,
                                  p.value=pvalue,
                                  q.value=q.value,
                                  OR=OR,
                                  OR.lower95 = OR.lower95,
                                  OR.upper95 = OR.upper95) %>% 
  filter(q.value < 0.05)
write.csv(q1.4_adj_los,file = "result/Q1.4_candidate_tumour_suppressor_genes.csv")

rm(list = ls()[!str_detect(ls(),"q1.")])
save.image(file="result/Q1_Mutation_data_analysis_result.Rda")


# Q2.Transcriptomic data normalization -------
## Q2.1 -------
exp_dat <- readRDS("data/expr_sclc_ucologne_2015.rds")

#judge whether need to transformation 
qx <- as.numeric(quantile(exp_dat, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC){
  q2.1_exprset <- log2(exp_dat+1)
  print("log2 transform finished")
}

## Q2.2 -------
#with the data in Q2.2_test.txt to test the median polish algorithm from scartch
#the result in Q2.2_test_scratch_output.txt

## Q2.3 -------
q2.3_med_scr_test <- read.table("result/Q2.2_test_scratch_output.txt") %>%.$V1 %>% matrix(.,ncol = 5) %>% t()

q2.3_med_stat_test <- read.table("result/Q2.2_test.txt") %>%.$V1 %>% matrix(.,ncol = 5) %>% t() %>% 
  medpolish(.) %>% .$residuals

#draw heatmap to compare the residuals of median polish algorithm in scratch and stat
mat <- q2.3_med_scr_test
before_exp <- Heatmap(mat,
                      name = "scratch_exp",
                      col = colorRampPalette(c("#4575B4", "white", "#D73027"))(50),
                      #top_annotation = ann.top,
                      show_column_names = F,
                      show_row_names = F,
                      cluster_columns = T,
                      cluster_rows = T,
                      column_title = "scratch_expr",
                      #row_labels = gene,
                      #row_names_gp = gpar(fontsize = 3.2),
                      #row_names_side = "left",
                      width = ncol(mat) * unit(8, "mm"), # 格子的宽度
                      height = nrow(mat) * unit(8, "mm"),
                      #column_split = heat_dat[2, ], # 设置分类向量来分割热图，故会把热图根据组名进行分组
                      #column_title_rot = 90, # 分块之后每块的注释用column_title进行修改名字
                      #column_title_gp = gpar(fontsize = 6, fontface = "bold"),
                      #column_gap = unit(0.4, "mm"),
                      heatmap_legend_param = list(
                        labels_gp = gpar(fontsize = 7),
                        title_gp = gpar(fontsize = 8, fontface = "bold"),
                        legend_height = unit(1.5, "cm"),
                        grid_width = unit(0.3, "cm")
                      ))
mat <- q2.3_med_stat_test
after_exp <- Heatmap(mat,
                     name = "stat_exp",
                     col = colorRampPalette(c("#4575B4", "white", "#D73027"))(50),
                     #top_annotation = ann.top,
                     show_column_names = F,
                     show_row_names = F,
                     cluster_columns = T,
                     cluster_rows = T,
                     column_title = "stat_expr",
                     #row_labels = gene,
                     #row_names_gp = gpar(fontsize = 3.2),
                     #row_names_side = "left",
                     width = ncol(mat) * unit(8, "mm"), # 格子的宽度
                     height = nrow(mat) * unit(8, "mm"),
                     #column_split = heat_dat[2, ], # 设置分类向量来分割热图，故会把热图根据组名进行分组
                     #column_title_rot = 90, # 分块之后每块的注释用column_title进行修改名字
                     #column_title_gp = gpar(fontsize = 6, fontface = "bold"),
                     #column_gap = unit(0.4, "mm"),
                     heatmap_legend_param = list(
                       labels_gp = gpar(fontsize = 7),
                       title_gp = gpar(fontsize = 8, fontface = "bold"),
                       legend_height = unit(1.5, "cm"),
                       grid_width = unit(0.3, "cm")
                     ))
p1 <- before_exp+after_exp
pdf(file = "result/Q2.2_heatmap_of_residuals_of_scratch_and_stat_median_polish.pdf", width = 12, height = 8, family = "GB1")
p1
dev.off()

# the result is same.Taking into account time cost,because scratch takes a lot of running time to get residuals,we use median polish algorithm in stat

q2.3_med_exp <- medpolish(q2.1_exprset) %>% .$residuals
save(q2.3_med_exp,file = "result/Q2.3_gene_expression_matrix_with_medpolish_result.Rda")

## Q2.4 --------
library(ComplexHeatmap)

quantile(q2.1_exprset,seq(0,1,0.1))
mat <- q2.1_exprset
mat[mat > 5] <- 5
before_exp <- Heatmap(mat,
               name = "bef_exp",
               col = colorRampPalette(c("#4575B4", "white", "#D73027"))(50),
               #top_annotation = ann.top,
               show_column_names = F,
               show_row_names = F,
               cluster_columns = T,
               cluster_rows = T,
               column_title = "before_expression",
               #row_labels = gene,
               #row_names_gp = gpar(fontsize = 3.2),
               #row_names_side = "left",
               width = ncol(mat) * unit(1, "mm"), # 格子的宽度
               height = nrow(mat) * unit(0.005, "mm"),
               #column_split = heat_dat[2, ], # 设置分类向量来分割热图，故会把热图根据组名进行分组
               #column_title_rot = 90, # 分块之后每块的注释用column_title进行修改名字
               #column_title_gp = gpar(fontsize = 6, fontface = "bold"),
               #column_gap = unit(0.4, "mm"),
               heatmap_legend_param = list(
                 labels_gp = gpar(fontsize = 7),
                 title_gp = gpar(fontsize = 8, fontface = "bold"),
                 legend_height = unit(1.5, "cm"),
                 grid_width = unit(0.3, "cm")
               ))

quantile(q2.3_med_exp,seq(0,1,0.1))
mat <- q2.3_med_exp
mat[mat > 1] <- 1
mat[mat < -1] <- (-1)
after_exp <- Heatmap(mat,
               name = "aft_exp",
               col = colorRampPalette(c("#4575B4", "white", "#D73027"))(50),
               #top_annotation = ann.top,
               show_column_names = F,
               show_row_names = F,
               cluster_columns = T,
               cluster_rows = T,
               column_title = "after_expression",
               #row_labels = gene,
               #row_names_gp = gpar(fontsize = 3.2),
               #row_names_side = "left",
               width = ncol(mat) * unit(1, "mm"), # 格子的宽度
               height = nrow(mat) * unit(0.005, "mm"),
               #column_split = heat_dat[2, ], # 设置分类向量来分割热图，故会把热图根据组名进行分组
               #column_title_rot = 90, # 分块之后每块的注释用column_title进行修改名字
               #column_title_gp = gpar(fontsize = 6, fontface = "bold"),
               #column_gap = unit(0.4, "mm"),
               heatmap_legend_param = list(
                 labels_gp = gpar(fontsize = 7),
                 title_gp = gpar(fontsize = 8, fontface = "bold"),
                 legend_height = unit(1.5, "cm"),
                 grid_width = unit(0.3, "cm")
               ))
p1 <- before_exp+after_exp
pdf(file = "result/Q2.4_heatmap_of_results_before_and_after_median_polish.pdf", width = 12, height = 8, family = "GB1")
p1
dev.off()

## Q2.5 -------
#the result is q2.3_med_exp

rm(list = ls()[!str_detect(ls(),"q2.")])
save.image(file="result/Q2_Transcriptomic_data_normalization_result.Rda")


# Q3.Differential gene expression analysis --------
## Q3.1 ---------
ph_dat <- fread("data/pheno_sclc_ucologne_2015.tsv") 
ph_dat[ph_dat == ""] <- NA
q3.1_ph_dat <- ph_dat %>% drop_na(uicc_tumor_stage) %>%
  mutate(stage = str_replace(uicc_tumor_stage,"a|b|B",""),
         group = ifelse(stage == "I"|stage == "II","early","advanced"))

## Q3.2-------

#Use the R package "limma" to identify differentially expressed genes
#When fitting the linear model, age and sex were included as covariates to remove the influence 
#differentially expressed genes satisfy the standard:|logFC|>1 & P.Value < 0.05 

library(limma)

load("result/Q2.3_gene_expression_matrix_with_medpolish_result.Rda")

samp_inter <- intersect(colnames(q2.3_med_exp),q3.1_ph_dat$patient_id)
design_dat <- q3.1_ph_dat %>%
  filter(patient_id %in% samp_inter) %>%
  transmute(sample = patient_id,
            age = age,
            sex = ifelse(sex == "Male",1,2),
            group = group)

expr_dat <- q2.3_med_exp %>% as.data.frame() %>% dplyr::select(samp_inter) %>% as.matrix()

design <- model.matrix(~0+factor(design_dat$group)+design_dat$age+design_dat$sex)  
rownames(design) <- design_dat$sample
colnames(design) <- c("case","control","age","sex")
contrast.matrix <- makeContrasts(case-control,levels = design)

fit <- lmFit(expr_dat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)
q3.2_DEG <- topTable(fit2, coef=1, number = nrow(expr_dat),adjust="BH") %>%
  filter(abs(logFC)>1 & P.Value < 0.05 )

rm(list = ls()[!str_detect(ls(),"q3.")])
save.image(file="result/Q3_Differential_gene_expression_analysis_result.Rda")


# Q4.Integrative analysis -------
## Q4.1 -------
load("result/Q2.3_gene_expression_matrix_with_medpolish_result.Rda")

var_dat <- fread("data/sv_sclc_ucologne_2015.tsv")
pv_exp <- q2.3_med_exp %>% as.data.frame() %>%
  rownames_to_column(var = "symbol") %>%
  pivot_longer(-symbol,names_to = "sample_id")

q4.1_exp_var <- var_dat %>%
  dplyr::select(sample_id,site1_hugo_symbol,site2_hugo_symbol,site2_effect_on_frame,event_info) %>%
  pivot_longer(-c(sample_id,site2_effect_on_frame,event_info),values_to = "symbol") %>%
  inner_join(.,pv_exp,by=c("sample_id","symbol"))

## Q4.2 --------
#the SV matrix only contain two samples involved in structural variants,whether the SV matrix is not complete?
#when identifying SVs,the invovled genes in samples with the SV only has one sample,so i could not use t test,
#so i compare the expression level in the sample with the median of other samples' expression level

highly_gene <- c()
for (i in 1:nrow(q4.1_exp_var)) {
  oth_med <- pv_exp %>% 
    filter(symbol == q4.1_exp_var[i,]$symbol & sample_id != q4.1_exp_var[i,]$sample_id) %>% 
    .$value %>% median()
  if(q4.1_exp_var[i,]$value > oth_med) 
    highly_gene  <- c(q4.1_exp_var[i,]$symbol,highly_gene)
}

q4.2_sel_SV <- var_dat %>%
  filter(site1_hugo_symbol %in% highly_gene & site2_hugo_symbol %in% highly_gene & site2_effect_on_frame == "in-frame")

rm(list = ls()[!str_detect(ls(),"q4.")])
save.image(file="result/Q4_Integrative_analysis_result.Rda")
