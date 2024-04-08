# Tissue Informative Markers for cfDNA
Code for the analysis of tissue informative markers in cfDNA data

## TIMs 

TIMs were called using publicly available reference data from [ENCODE](https://www.encodeproject.org/) and [Blueprint](http://blueprint-data.bsc.es/#!/). 

#### dependencies 
- Python3
- Python packages: bottleneck, NumPy

## Methylation Pipeline 

To map reads and call methylation, we used [BsBolt](https://github.com/NuttyLogic/BSBolt). Duplicate reads were removed using [UMI-tools](https://github.com/CGATOxford/UMI-tools). 

#### dependencies 
- BsBolt
- Samtools 
- UMI-tools 

## Disease prediction models

Models to predict ALS disease status and phenotypes were developed using [BigStatsR](https://cran.r-project.org/web/packages/bigstatsr/index.html), a statistical package optimized for genomic data. 

#### dependencies 
- R 
- BigStatsR
- data.table
- dplyr
- pROC



