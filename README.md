# R Package {nphDesign}
Steps to install this package on your R console. The required R version is 4.0.3 or later. 

(1);   
install.packages("devtools") #install devtools package if not installed yet.

(2); 
install.packages("roxygen2") #install roxygen2 package if not installed yet.

(3); 
rm(list = ls()) #clean enviornment

(4); 
remove.packages("nphDesign") #remove previous unsuccessfully installed package or earlier version of the package.

(5); 
devtools::install_github("philiphe2018/nphDesign") #install the package from github repository.

(6) Test installation; 
library(nphDesign) #call nphDesign package

(7); 
?wlr  #Call the R documentation of function wlr()
help(package="nphDesign") #Overview of the package functions


If it doesn't work somehow, close R, then try again. 
