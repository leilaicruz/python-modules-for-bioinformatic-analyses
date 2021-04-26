---
title: "Bioinformatic-analyses"
output: pdf_document
author: Leila Inigo de la Cruz
classoption: onecolumn
pdf_document:
latex_engine: pdflatex
toc: true
lof: true
numberSections: true
highlight: tango
sectionsDepth: 3
chapters: True
figPrefix:
  - "Fig."
  - "Figs."
secPrefix:
  - "Section"
  - "Sections"
fontsize: 12pt
geometry: margin=0.5in
autoEqnLabels: true
cref: true
crossref: true
---



All of the functions and scripts here developed are aimed to answer some simple questions that require the data that is already publicly published. 
The main source of data are the databases:

- [YeastMine](https://yeastmine.yeastgenome.org/yeastmine/begin.do)
- [BioGRid](https://thebiogrid.org/)
- [SGD](http://sgd-archive.yeastgenome.org/curation/literature/)


## Questions guiding bioinformatic studies

1. ***Are the number of genetic interactors of essential genes more than for non essential genes?*** 

    - [GO TO THE SCRIPT](../src(source-code)/script_interactors-of-essential-genes.py)

![Essential genes are slightly more connected](../output_images/essential-and-not-essential-genes-number-of-interactors.png)

- Conclusion:
    - Essential genes are slightly more connected than non essential genes. 
    - There are genes low connected that are also essential genes. 


2. ***Is there any correlation with the common go terms and common interactors with the type of genetic interaction a pair of genes has?*** 

    - [GO TO SCRIPT](../src(source-code)/script_big-loop-to-know-interactors-and-go-terms-common.py)

![There is no apparent correlation between the common go terms and common interactors with the type of genetic interaction all gene pairs share. ](../output_images/no-correlation-based-on-the-common-interact-and-go-terms-type.png)

- Conclusion:
    - There is no apparent correlation between common go terms and common interactors with the type of genetic  interaction. 
    - It is interesting that the SL curve is in the tail for both measurements. So, generally SL pairs share more interactors among them and also more common go terms , which could mean that they in general belong to the same functions/modules. 


3. ***Are paralogs functionally divergent?***

- [GO TO SCRIPT](../src(source-code)/script_common_go_for_paralogs.py)

![Evidence that paralogs have functionally diverged](../output_images/functional-diversification-of-paralogs.png)

- Conclusion:
    - According to the figure, most of the paralogs have functionally diverged into different functions because most of them **do not share common go terms**. 

4. ***How many paralogs are also synthetic lethal pairs? ***

- [GO TO SCRIPT](../src(source-code)/script_paralogs_and_SL_relationship.py)

![A quarter of paralogs are also synthetic lethals](../output_images/one-quarter-of-paralogs-are-SL.png)

- Conclusion:
    - A quarter of the paralogs so far founded in budding yeast are also synthetic lethal.

5. ***Are the interaction scores from Constanzo et al 2016 fitness data from SGA experiments following the multiplicative model?***

- [GO TO SCRIPT](../src(source-code)/script_multiplicative-model-go-terms-ocurrence.py)

![Scores from multiplicative model checked from SGA experiments](../output_images/BEM1_data_from_constanzo-check-of-the-scores.png)

- Conclusion:
    - Yes, the scores follows a multiplicative model , because you can see that the points above have positive scores and the points below have negative scores.

6. ***How are the fitness map from SATAY vs SGA of dpl1 gene? ***

- [GO TO SCRIPT](../src(source-code)/script_fitness-map-Constanzo-vs-SATAY.py)

![SATAY gives much more insight from the fitness map upon a gene deletion than SGA.](../output_images/constanzo-vs-satay-dpl1-fitness-map.png)

- Conclusion:
    - With SATAY we will have much more data to fill the whole fitness map compared to the existing available SGA data. 





 