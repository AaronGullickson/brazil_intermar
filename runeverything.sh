#!/bin/bash

#Set your stata path
stata=/Applications/Stata/StataMP.app/Contents/MacOS//stata-mp

#get rid of old dat
rm output/*
rm logs/*

#unzip the input data
gunzip input/brazil.dta.gz

#run the analysis
$stata -b do getaggdata.do
rm getaggdata.log
R CMD BATCH runfull.R
rm runfull.Rout

#zip the input data again
gzip input/brazil.dta
