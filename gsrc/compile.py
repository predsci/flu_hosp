#!/usr/bin/python
import os
import string
os.system("make")
os.system("R CMD SHLIB detsirh.o") # ignbin.o toms343.o utils.o")
os.system("R CMD SHLIB stosirh.o ignbin.o")
os.system("R CMD SHLIB mcmc.o detsirh.o") 
#os.system("R CMD SHLIB gridcovid.o ignbin.o toms343.o utils.o")
os.system("/bin/rm -rf *.o")
