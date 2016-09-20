################################
# runfull.R 
# this will run the full anaysis in R
# currently it just knits the R markdown file 
# since that file internally calls the analysis
################################

library(rmarkdown)
render("figures.Rmd")
