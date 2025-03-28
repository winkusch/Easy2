# Easy2 (= EasyQC + EasyStrata)

R package for GWAS quality control (= EasyQC) and for evaluation of stratified GWAS (= EasyStrata) . 

## Description

The EasyX R package combines the full functionality of two other R packages: 
- EasyQC (Quality control for GWAS, [Winkler et al. Nat Protoc 2014](https://pubmed.ncbi.nlm.nih.gov/24762786/))
- EasyStrata (Stratified GWAS evaluation, [Winkler et al. Bioinformatics 2015](https://pubmed.ncbi.nlm.nih.gov/25260699/)) 
</ul>
Any functions exported from EasyQC and EasyStrata are also available through Easy2.    
More information as well as template scripts and reference files can be found at [www.genepi-regensburg.de/software](https://www.genepi-regensburg.de/software). 

### Dependencies

Easy2 depends on the following other R packages:  
- Cairo, plotrix, data.table, forestplot 
</ul>
These need to be installed prior to the installation of EasyX. 
Using parallelization for functions INDEP and PCA2STEP, you will further require R packages bigsnpr, foreach and doParallel. 

### Installing

We recommend to use 'devtools' in R to install the package directly from github:  
```
library(devtools)
install_github("winkusch/Easy2")
```
Alternatively, you can download the tarball and install it in R: 
```
install.packages("EasyX_VERSION.tar.gz")
```

### Documentation

Please see the [wiki](https://github.com/winkusch/Easy2/wiki) for Easy2 functions and parameters. 

### Executing program

To start the program, you will have to load the package in R and call an ecf-script by the EasyX function: 
```
library(Easy2)
Easy2("/path2script/template.ecf")
```  
Please check our website [www.genepi-regensburg.de/software](https://www.genepi-regensburg.de/software) for ecf-script templates and reference files. 
## Help

Please contact thomas.winkler@klinik.uni-regensburg.de if you need assistance. 

## Authors

Thomas Winkler, University of Regensburg, Germany

## Version History

* 1.2.4 first git release

## Citation

If you are using the package for GWAS quality control, please cite [Winkler et al. Nat Protoc 2014](https://pubmed.ncbi.nlm.nih.gov/24762786/).   
If you are using the package for (stratified) GWAS results evaluation, please cite [Winkler et al. Bioinformatics 2015](https://pubmed.ncbi.nlm.nih.gov/25260699/). 
