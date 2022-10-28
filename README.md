# Programming Exercises for R

** all analysis script in `script/script.R` **
* also with some questions when I analyzed the data *

## Mutation data analysis

The mutation data is in `data/mutations_sclc_ucologne_2015.rds`.

1. Identify the top 10 most frequently mutated genes.
   Identify samples whose mutation count is in the 80 to 90 percentile.
2. Categorize variants based on their expected effects using the
   `data/mutation_effects.tsv` table.
   Generate a count matrix containing the numbers of loss-of-function
   and neutral mutations for each gene.
3. Implement a statistical test that determines whether a gene has a
   significantly higher proportion of loss-of-function mutations
   (excluding mutations with uncertain effects),
   compared to other genes.
4. Identify candidate tumour suppressor genes using this statistical test,
   adjusting for multiple hypothesis testing.
   The output table should contain:
     - gene symbol
     - an estimate of effect size
     - p value
     - q value
5. Perform a literature search and explain the function of each candidate
   gene in the context of cancer, as well as specifically in small cell 
   lung cancer.
   
** Q1-Q4 answer in `result/Q1_Mutation_data_analysis_result.Rda` **    
** Q5 answer in `result/Q1.5_candidate_tumour_suppressor_genes.xlsx` **

## Transcriptomic data normalization

The transcriptomic data is in `data/expr_sclc_ucologne_2015.rds`.

1. Perform an appropriate log transformation on the data.
2. Implement a median polish algorithm from scartch.  
    ** `result\Q2.2_test.txt` is the test data for the algorithm **   
    ** scartch program in `result\Q2.2_scratch_median_polish.sb3` **   
    ** scartch output in `result\Q2.2_test_scratch_output.txt` **   
3. Compare the residuals of your algorithm and `stats::medpolish`.
    ** answer in `result\Q2.4_heatmap_of_results_before_and_after_median_polish.pdf` **    
    ** `result\Q2.3_gene_expression_matrix_with_medpolish_result.Rda' is the residuals of expression matrix with `stats::medpolish` **   
4. Plot heatmaps of the results before and after median polish.
    ** answer in `result\Q2.4_heatmap_of_results_before_and_after_median_polish.pdf` **   
5. Output the median polished residual matrix as the normalized transcriptomic data.

** other answer in `result\Q2_Transcriptomic_data_normalization_result.Rda` **


## Differential gene expression analysis

The clinical data is in `data/pheno_sclc_ucologne_2015.tsv`.

1. Define two groups of tumours as early stage (stages I-II) vs. advanced stage
   tumours (stages III-IV), while excluding samples missing stage information.
2. Identify genes that differentially expressed in early vs. advanced stage
   tumours using an appropriate R package.
   
 ** answer in `result\Q3_Differential_gene_expression_analysis_result.Rda` **

## Integrative analysis

The structural variant data is in `data/sv_sclc_ucologne_2015.tsv`.

1. For each gene involved in a structural variant, determine the expression 
   level of the gene in the sample that harbours the structural variant.
2. Identify structural variants that satisfy the following critiera:
      - The involved pair of genes both have elevated expression levels
        in the affected samples compared to unaffected samples.
      - The second gene in the pair is in frame.

** answer in `result\Q4_Integrative_analysis_result.Rda` **


