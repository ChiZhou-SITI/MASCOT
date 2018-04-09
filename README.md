# MASCOT: Model-based Analysis of Single-cell CRISPR knockOuT screening

* **MASCOT** is the first one-step applicable pipeline based on topic model to analyze single-cell CRISPR screening data, which could help to prioritize the knockout gene impact in a cellular heterogeneity level.
* **MASCOT** is an integrated pipeline for model-based analysis of single cell CRISPR knockout screening data. **MASCOT** consists of three steps: **data preprocessing**, **model building** and **knockout effect prioritizing**:<br> 
    * **Data preprocessing**: Besides the routine quality control and data normalization applied in single-cell RNA-seq analysis, **MASCOT** addresses several considerations that should be taken into account for such a novel data type: (1) Filtering out cells with an extremely large proportion of zero knockout expression among the control cells, as these cells are inherently noisy or meaningless. (2) Reducing the false positive rate of gene knockout and (3) Filtering out knockout cells without a sufficient number of cells to capture the corresponding perturbation phenotype;<br>
    * **Model building**: **MASCOT** builds an analytical model based on Topic Models to handle single-cell CRISPR screening data. The concept of topic models was initially presented in machine learning community and has been successfully applied to gene expression data analysis. A key feature of topic model is that it allows each knockout sample to process a proportion of membership in each functional topic rather than to categorize the sample into a discrete cluster. Such a topic profile, which is derived from large-scale cell-to-cell different knockout samples, allows for a quantitative description of the biologic function of cells under specific gene knockout conditions. **MASCOT** addresses several specific issues when applying the topic model to this specific data type: (1) Differences between sgRNA knockout efficiencies is are considered and rectified in the evaluation of the knockout effect on the corresponding cells. (2) The distribution of topics between cases and controls is affected by the ratio of their sample numbers, and such a sample imbalance issue is addressed when evaluating the knockout effect. (3) The optimal topic number is automatically selected by MASCOT in a data-driven manner.<br>
    * **Knockout effect prioritizing**: Based on the model-based perturbation analysis, **MASCOT** can quantitatively estimate and prioritize the individual gene knockout effect on cell phenotypes from two different perspectives, i.e., prioritizing the gene knockout effect as an overall perturbation effect, or in a functional topic-specific way.<br>    
* **Input data format description**:  
* For illustration purposepurposes, we only took the unstimulated cell data applied in CROP-seq as a an example plot.<br>
    * Install: you can install the **MASCOT** package from Github using **devtools** packages.<br>
    ```r
    require(devtools)
    install_github("BinDuan/MASCOT")
    library(MASCOT)
    ```
    * The data **crop_unstimulated.rda** is a dataset containing "crop_unstimulated_expression_profile", "crop_unstimulated_sample_info_gene" and "crop_unstimulated_sample_info_sgRNA".
    ```r
    data("crop_unstimulated",MASCOT)
    # expression_profile
    dim(crop_unstimulated_expression_profile)
    ```
    ```
    ## [1] 36722  2646
    ```
    ```r
    crop_unstimulated_expression_profile[1:3,1:3]
    ```
    ```
    ##          GCAGTCCTTCTN ACGTAGGGGTAN AAACAACCGAAN
    ## A1BG                0            0            0
    ## A1BG-AS1            0            0            0
    ## A1CF                0            0            0
    ```
    ```r
    # sample_info_gene
    length(crop_unstimulated_sample_info_gene)
    ```
    ```
    ## [1] 2646
    ```
    ```r
    class(crop_unstimulated_sample_info_gene)
    head(crop_unstimalated_sample_info_gene)
    ```
    ```
    ## [1] "character"
    
    ## GCAGTCCTTCTN ACGTAGGGGTAN AAACAACCGAAN TCAGTGGCTTCT AGTATTCTCACN TTATAGCATGCA 
    ##      "NR4A1"     "NFATC2"       "CTRL"       "CTRL"       "CTRL"       "CTRL"
    ```
    ```r
    # sample_info_sgRNA
    length(crop_unstimulated_sample_info_sgRNA)
    ```
    ```
    ## [1] 2646
    ```
    ```r
    class(crop_unstimulated_sample_info_sgRNA)
    head(crop_unstimulated_sample_info_sgRNA)
    ```
    ```
    ##[1] "character"
    
    ##         GCAGTCCTTCTN          ACGTAGGGGTAN    AAACAACCGAAN   TCAGTGGCTTCT    AGTATTCTCACN    TTATAGCATGCA  
    ## "Tcrlibrary_NR4A1_1"  "Tcrlibrary_NFATC2_1"    "CTRL00698"    "CTRL00320"     "CTRL00087"     "CTRL00640" 
    ```
    
    * The first step: data preprocessing.
    ```r
    # integrate the input data and filter filter mitochondrial ribosomal protein(^MRP) and ribosomal protein(^RP)
    crop_seq_list<-InputAndPreprocess(crop_unstimulated_expression_profile,crop_unstimulated_sample_info_gene,crop_unstimulated_sample_info_sgRNA,sample_info_batch=NULL)
    
    # quality control
    crop_seq_qc<-singleCellCRISPRscreen_qc(crop_seq_list$expression_profile,crop_seq_list$sample_info_gene,gene_low=500,species="Hs",plot=T)
    
    # other filterings, including "zero_ratio", "sgRNA efficiency" and "phenotype capture".
    crop_seq_filtered<-cellFilteringAndKOefficiencyCalculating(crop_seq_qc$expression_profile_qc,crop_seq_qc$sample_info_gene_qc,crop_seq_list$sample_info_sgRNA,nonzero=0.01,grna_cell_num=10,fold_change=0.5,plot=T)
    
    # plot cells with different component.
    component<-plot_filtering_information_component(crop_seq_list$sample_info_gene,crop_seq_qc$sample_info_gene_qc,crop_seq_filtered$nonzeroRatio,crop_seq_filtered$sample_info_gene_qc_zr_se,crop_seq_filtered$sample_info_gene_qc_zr_se_pc)
    ```
    * The second step: model building
    ```r
    # obtain high dispersion different genes.
    crop_seq_vargene<-getHighDispersionDifferenceGenes(crop_seq_filtered$expression_profile_qc_zr_se_pc,crop_seq_filtered$sample_info_gene_qc_zr_se_pc,plot=T)
    
    # get topics and select topic number automatically.
   optimalTopics<-getTopics(crop_seq_vargene,crop_seq_filtered$sample_info_gene_qc_zr_se_pc,plot=T)
    
    # plot heatmap between cells and topics.
   plot_cellAndTopic(optimalTopics)
   
    # annotate each topic's functions. Hs(homo sapiens) or Mm(mus musculus) are available.
   topic_enrichment<-topic_functionAnnotation(optimalTopics,species="Hs")
   
    # get offtarget information. This step won't affect the final ranking result, but give you offtarget information. If you don't want to consider this factor, you can skip this step. 
    library(CRISPRseek)
    library("BSgenome.Hsapiens.UCSC.hg38")
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
    gRNAFilePath<-"~/grna.fa"
    crop_results <- offTargetAnalysis(inputFilePath = gRNAFilePath, findgRNAs = FALSE,findgRNAsWithREcutOnly = FALSE,findPairedgRNAOnly = FALSE, BSgenomeName = Hsapiens,txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,min.score=1,scoring.method = "CFDscore",orgAnn = org.Hs.egSYMBOL, max.mismatch = 3,outputDir=getwd(), overwrite = TRUE)
    # then, check the offtargets further.
    offTarget_Info<-getOffTargetInfo(crop_results,crop_seq_filtered$expression_profile_qc_zr_se_pc,crop_seq_filtered$sample_info_gene_qc_zr_se_pc,crop_seq_list$sample_info_sgRNA)
    
    ```
    * The third step: knockout effect prioritizing
    ```r
    # knockout effect prioritizing
    
    # calculate topic distribution for each cell
    ratioDiff<-getRatioAndDiff(optimalTopics,crop_seq_filtered$sample_info_gene_qc_zr_se_pc,KO_efficiency=crop_seq_filtered$KO_efficiency)
    
    # calculate the overall perturbation effect ranking list with offtarget_info calculated
    rankTotal_result<-rankOverall(ratioDiff,offTarget_hash=offTarget_info)
    # calculate the overall perturbation effect ranking list without offtarget_info calculated
    rankOverall_result<-rankOverall(ratioDiff)
    
    # calculate the topic-specific ranking list 
    rankTopicSpecific_result<-rankTopicSpecific(rankOverall_result$rankOverall_result_detail)
    
    ```
    * output
    ```r
    # output
    
    # output the overall perturbation effect ranking list with concision and detailed styles.
    write.table(rankOverall_result$rankOverall_result_concision,"~/rankOverall_result_concision.txt",col.names=T,row.names=F,quote=F,sep="\t")
    write.table(rankOverall_result$rankOverall_result_detail,"~/rankOverall_result_detail.txt",col.names=T,row.names=T,quote=F,sep="\t")
    
    # output the topic-specific ranking list with concision and detailed styles.
    write.table(rankTopicSpecific_result$rankTopicSpecific_result_concision,"~/rankTopicSpecific_result_concision.txt",col.names=T,row.names=T,quote=F,sep="\t")
    write.table(rankTopicSpecific_result$rankTopicSpecific_result_detail,"/rankTopicSpecific_result_detail.txt",col.names=T,row.names=T,quote=F,sep="\t")
    
 
