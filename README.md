#MASCOT: **M**odel-based **A**nalysis of **S**ingle-cell **C**RISPR knock**O**u**T** screening

- **MASCOT** is the first one-step applicable pipeline based on topic model to analyze single-cell CRISPR screening data, which could help to prioritize the knockout gene impact in a cellular heterogeneity level.

- **MASCOT** is an integrated pipeline for model-based analysis of single cell CRISPR knockout screening data. **MASCOT** consists of three steps: data preprocessing, model building and knockout effect prioritizing: ##data preprocessing **(1)** In the first step, besides the routine quality control and data normalization applied in single-cell RNA-seq analysis, **MASCOT** addresses several considerations that should be taken into account for such a novel data type: (a) Filtering out cells with an extremely large proportion of zero knockout expression among the control cells, as these cells are inherently noisy or meaningless. (b) Reducing the false positive rate of gene knockout and (c) Filtering out knockout cells without a sufficient number of cells to capture the corresponding perturbation phenotype; **(2)** Second, **MASCOT** builds an analytical model based on Topic Models to handle single-cell CRISPR screening data. The concept of topic models was initially presented in machine learning community and has been successfully applied to gene expression data analysis. A key feature of topic model is that it allows each knockout sample to process a proportion of membership in each functional topic rather than to categorize the sample into a discrete cluster. Such a topic profile, which is derived from large-scale cell-to-cell different knockout samples, allows for a quantitative description of the biologic function of cells under specific gene knockout conditions. **MASCOT** addresses several specific issues are faced and addressed in MASCOT when applying the topic model to this specific data type: (1) Difference Differences between sgRNA knockout efficiencies is are considered and rectified in the evaluation of the knockout impact effect on the corresponding cells. (2) The distribution of topics between case cases and control controls is affected by the ratio of their sample numbers, and such a sample imbalance issue is addressed when evaluating the knockout effect (Supplementary Fig. 2b and see Online Methods). (3) The optimal topic number is automatically selected by MASCOT in a data-driven manner (see Online Methods).

Finally, based on the topic-model based perturbation analysis, MASCOT can quantitatively estimate and prioritize the individual gene knockout impact on cell phenotypes in two different perspectives (Fig. 1, Fig. 3, Supplementary Fig. 1, and see Online Methods), i.e., prioritizing the gene knockout effect either as an overall perturbation impact, or in a functional topic-specific way (Fig. 3b and see Online Methods).  

- As an example, we presented our in-house oral microbial data from 39 human samples (collected by the Stomatological Hospital Affiliated to Tongji University, Shanghai,China). These sample data are 16s-rRNA-sequencing-based, which are targeted on a clinically classical oral diease, i.e. the oral lichen planus (OLP), a chronic oral disease without clear pathological mechanism and effective treatment. These samples can be grouped as the controls (**NOT_OLP**) and two OLP disease sub-types - **OLP_non-erosive** and **OLP_erosive**. Using **MetaTopics**, the relationship between topics and a specific disease is deciphered, which provides new clues to understand the pathological mechanism of OLP. The following examples present an comprehensive illustration of using **MetaTopics** in the analysis of these samples.
- Install: you can install the **MetaTopics** package from Github using **devtools** packages.


```r
require(devtools)
install_github("bm2-lab/MetaTopics")
library(MetaTopics)
```

- The data **meta_counts** is a matrix containing the read counts of the microbe in each individual sample mapped to a certain microbial reference. The data **genus_2_phylum** is a data frame containing the annotation for the microbe in data meta_counts


```r
dim(meta_counts)
```

```
## [1]  39 129
```

```r
colnames(meta_counts)[1:10]
```

```
##  [1] "Abiotrophia"     "Acholeplasma"    "Achromobacter"  
##  [4] "Acidaminococcus" "Acidovorax"      "Acinetobacter"  
##  [7] "Actinobacillus"  "Actinomyces"     "Aggregatibacter"
## [10] "Alloscardovia"
```

```r
head(genus_2_phylum)
```

```
##                           genus         phylum presence  colour
## Abiotrophia         Abiotrophia     Firmicutes        1 #DECBE4
## Acholeplasma       Acholeplasma        Unknown        0 #E5D8BD
## Achromobacter     Achromobacter Proteobacteria        1 #CCEBC5
## Acidaminococcus Acidaminococcus        Unknown        0 #E5D8BD
## Acidovorax           Acidovorax        Unknown        0 #E5D8BD
## Acinetobacter     Acinetobacter Proteobacteria        1 #CCEBC5
```

The column phylum has the corresponding phylum level annotation for each microbiome. The column colour is used for abundance.plot as an identification of the phylum level. In our example, the data has 39 samples and the data **lls** contains the disease type annotation information for each sample.


```r
rownames(meta_counts)
```

```
##  [1] "NO2"  "NO6"  "NO8"  "NO9"  "NO13" "NO14" "NO15" "NO16" "NO18" "NO19"
## [11] "NO22" "NO23" "NO24" "NO25" "NO28" "NO30" "NO31" "NO32" "NO34" "NO35"
## [21] "NO36" "NO37" "NO38" "NO39" "NO40" "NO41" "NO42" "NO43" "NO44" "NO45"
## [31] "NO46" "NO48" "NO49" "NO50" "NO52" "NO53" "NO54" "NO55" "NO56"
```

```r
table(lls)
```

```
## lls
##        NOT_OLP    OLP_erosive OLP_no-erosive 
##             16             14              9
```

Users can explore the abundance distribution profile of the data:


```r
meta_abundance <- micro.abundance(meta_counts,1)
genus_2_phylum=genus_2_phylum[colnames(meta_abundance),]
abundance.plot(meta_abundance,classification = genus_2_phylum$phylum,col=genus_2_phylum$colour)
```

![](Readme_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

- In order to use topic model, some odd samples and microbe need to be filtered. Users can use the function **noise.removal** in package **BiotypeR** to perform this procedure.


```r
library(BiotypeR)
data.denoized=noise.removal(t(meta_counts), percent=0.01)
samples=colnames(data.denoized)
bacteria=rownames(data.denoized)
data.final=meta_counts[samples,bacteria]
dim(data.final)
```

```
## [1] 39 88
```

After the pre-processing, the final data contains 39 samples with 88 microbe taxons.

- With the processed data in-hand, users can use cross-validation to find the appropriate topic number for topic model. The function **selectK** could be used to select the appropriate topic number and the function **plot_perplexity** helps to visualize the returned perplexity and likelihood in the topic number selection.


```r
library(slam)
library(topicmodels)
dtm=as.simple_triplet_matrix(data.final)
seed_num=2014
fold_num=5
kv_num = c(2:30)
sp=smp(cross=fold_num,n=nrow(dtm),seed=seed_num)
control = list(seed = seed_num, burnin = 1000,thin = 100, iter = 1000)
#not run: system.time((ctmK=selectK(dtm=dtm,kv=kv_num,SEED=seed_num,cross=fold_num,sp=sp,method='Gibbs',control=control)))
#not run: plot_perplexity(ctmK,kv_num)
```

- If users specify the topic number, function **LDA** in package **topicmodels** can be used to build the model. Here is an example with the topic number 10 specified.


```r
Gibbs_model_example = LDA(dtm, k = 10, method = "Gibbs",
            control = list(seed = seed_num, burnin = 1000,thin = 100, iter = 1000))
```

- The returned model is a S4 Object. The element **beta** in the model contains the estimation of the probability of each microbe in each topic.


```r
dim(Gibbs_model_example@beta)
```

```
## [1] 10 88
```

```r
apply(exp(Gibbs_model_example@beta),1,sum)
```

```
##  [1] 1 1 1 1 1 1 1 1 1 1
```

The function **plot_beta** is a visualization way to view the matrix result. This plot is based on ggplot2. The parameter **prob** is a cutoff used to restrict the number of points on the plots. If a microbe has a probability smaller than the cutoff in a topic, it will not be shown in the visualizd result.


```r
library(ggplot2)
plot_beta(Gibbs_model_example,prob=0.01)
```

![](Readme_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

The element **gamma** in the model containes the estimation of the probability of each topic in each sample. This function will act as a visualization representation of the result matrix. The parameter聽prob聽is a cutoff used to restrict the number of points on the plots. If a topic has a probability smaller than the cutoff in an individual, it will not be shown in the visualizd result.


```r
plot_gamma(Gibbs_model_example,lls,prob=0.05)
```

![](Readme_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

In order to interpret the relationships between the sub-communities and disease, Quetelet Index is introduced to estimate the relative change of the observation/occurence frequency of a latent sub-community among all the samples compared to that among the samples with a certain disease statue. Function **qindex** is used to compute the Quelete Index from the topic model. Quelete Index quantitatively describes the degree of the influence of the specific topic on certain disease. Parameter **prob** is a probability cutoff used to identify a meaningful sub-community observation. For a certain individual, the topics with probability no smaller than prob will be thought as a meaningful observation in this individual.


```r
Q_values <- qindex(Gibbs_model_example,lls,0.05)
head(Q_values)
```

```
##     label variable      quelete
## 1 NOT_OLP   Topic1  0.083333333
## 2 NOT_OLP   Topic2  0.772727273
## 3 NOT_OLP   Topic3  0.027667984
## 4 NOT_OLP   Topic4 -0.129186603
## 5 NOT_OLP   Topic5  0.083333333
## 6 NOT_OLP   Topic6 -0.004784689
```

```r
Q_values$Sign = factor(ifelse(Q_values$quelete<0,'Negative','Positive'),levels=c('Positive','Negative'))
Q_values$Abs = abs(Q_values$quelete)
p<-ggplot(Q_values,aes(label,variable))
p+geom_point(aes(cex=Abs,colour=Sign))+
  theme_bw(base_size = 12,base_family = "")+
  labs(x='Disease')
```

![](Readme_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

- Package LDAvis, build by ***Carson Sievert*** is a great tool to visualize and interpre the topics. It provides an web-based interaction application to show the relationship among topics and rank the terms, i.e. bacteria here, in each topic based on your topic model.


```r
library(LDAvis)
json <- with(TwentyNewsgroups,
             createJSON(phi, theta, doc.length, vocab, term.frequency))

json <- createJSON(phi = exp(Gibbs_model_example@beta),
                   theta = Gibbs_model_example@gamma,
                   doc.length = apply(data.final,1,sum),
                   vocab = Gibbs_model_example@terms,
                   term.frequency = apply(data.final,2,sum))
#not run
#serVis(json)
serVis(json, out.dir='olp_html',open.browser = FALSE)
```

```
## Warning in dir.create(out.dir): 'olp_html'已存在
```

![](LDAvis_webpage.png)

