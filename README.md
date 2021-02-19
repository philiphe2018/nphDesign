# R Package {nphDesign}
#Steps to install this package on your R console. The required R version is 4.0.3 or later. 

#(1) install devtools package if not installed yet;

install.packages("devtools")

#(2) install roxygen2 package if not installed yet;

install.packages("roxygen2") 

#(3) #clean enviornment;

rm(list = ls()) 

#(4) remove previous unsuccessfully installed package or earlier version of the package.; 

remove.packages("nphDesign") 

#(5)install the package from github repository.; 

devtools::install_github("philiphe2018/nphDesign") 

#(6) Test installation; 

library(nphDesign) 

#(7) Overview of the package functions;

help(package="nphDesign") 

#If it doesn't work somehow, close R, then try again. 
