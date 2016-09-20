/*****************************************
* getaggdata.do
*
* This program will read in individual-level data on 
* Brazilian couples, apply sample restrictions, and 
* aggregate the data for the log-linear analyisis. 
*
* The input data are based on Brazilian census data 
* downloaded from IPUMS that have been merged together
* to get male and female partners on the same line of data. 
******************************************/

log using logs/log_getaggdata.txt, replace text

use input/brazil, clear

*identify cohabitors
generate cohab = marsth==220

*age restriction: husband 25-35 and wife 14-60
keep if agem>=25 & agem<=35 & agef>=14 & agef<=60

*drop if missing values
drop if racem==. | racef==. | edfm==. | edff==.

*collapse educational categories into a smaller set
gen edf = 1
replace edf = 2 if edff==2
replace edf = 3 if edff==3 | edff==4
replace edf = 4 if edff==5
replace edf = 5 if edff==6
table edff edf

gen edm = 1
replace edm = 2 if edfm==2
replace edm = 3 if edfm==3 | edfm==4
replace edm = 4 if edfm==5
replace edm = 5 if edfm==6
table edfm edm
*we do not put labels on these so that we can treat them as scalars
*in R easier, when necessary

*rename the key variables to be consistent with model presentation
rename racem HR
rename racef WR
rename edm HE
rename edf WE

*shorter race labels
label drop racedef
label define racedef 1 "W" 2 "Br" 3 "Bl"
label values HR racedef
label values WR racedef

*adjust weights so that they correctly sum to sample size
egen wtmean = mean(wtperm)
replace wtperm = wtperm/wtmean

*save this data temporarily, so I can aggregate it different ways
save temp, replace

*aggregate the data
collapse (sum) wtperm, by(HR WR HE WE cohab)
rename wtperm Freq

*export it to csv
export delimited HR WR HE WE Freq using output/table_marriages.csv if !cohab, replace
export delimited HR WR HE WE Freq using output/table_cohabs.csv if cohab, replace

*now do the same but collapse to white and black
use temp, replace

replace HR=2 if HR==3
replace WR=2 if WR==3

label drop racedef
label define racedef 1 "W" 2 "Bl"
label values HR racedef
label values WR racedef

collapse (sum) wtperm, by(HR WR HE WE cohab)
rename wtperm Freq

*export it to csv
export delimited HR WR HE WE Freq using output/tableWB_marriages.csv if !cohab, replace
export delimited HR WR HE WE Freq using output/tableWB_cohabs.csv if cohab, replace

erase temp.dta

log close
