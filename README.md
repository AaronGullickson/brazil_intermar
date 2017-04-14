# Patterns of Racial and Educational Assortative Mating in Brazil

The files in this repository will reproduce the analysis used in the paper "Patterns of Racial and Educational Assortative Mating in Brazil" by Aaron Gullickson and Florencia Torche, published in [Demography](http://link.springer.com/article/10.1007/s13524-014-0300-2). The analysis was performed in Stata, and R, on Brazilian Census data from [IPUMS](http://www.ipums.org).

## The Data

The data used here were compiled from 6% random sample of the 2000 Brazilian Census. The data were extracted from the [IPUMS](http://www.ipums.org) system. Under the terms of the standard IPUMS contract, We are not currently allowed to re-distribute the data. We do, however, provide in our base directory the couple-level data that resulted from merging partners together from the original data. This dataset includes information about the union type, the race of each spouse, and the education of each spouse. It is saved as a g-zipped STATA data file.

Users should note that there is a mislabeling in Figure 1 of the published paper. The points that correspond to Brazil White/Brown and Brazil Brown/Black should be reversed. That error has been corrected in this repository.

## Running the Programs

The repository includes a bash script entitled `runeverything.sh` that will run the entire analysis. However, in order to use the script properly, users will need to edit it to provide the proper path to a stata executable. Alternatively, users can just run each script separately. Here is the order in which scripts should be run as well as a description of each script:

1.	 `getaggdata.do` - This STATA script will read in the couple-level data file `brazil.dta` in the `input` diretory, perform some data restrictions and then will aggregate this data into tables that will be exported as csv files to the `output` directory. Note that users will need to unzip the raw data file before running this script.

2. `analysis.R` - This R script will read in the tables in the `output` directory and perform the entire analysis, producing multiple log-linear models. A description of these models is provided in the comments at the beginning of this file.

3. `figures.Rmd` - This R markdown script will produce the figures from the paper as an html document entitled `figures.html`. It requires the models from `analysis.R` and will source this script in its start up.

There is also a `runfull.R` script but this is only used by the `runeverything.sh` bash script.
