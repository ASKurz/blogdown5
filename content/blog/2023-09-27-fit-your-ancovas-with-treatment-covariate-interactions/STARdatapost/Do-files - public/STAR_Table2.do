capture log close
local d Z
local directory DIRECTORY HERE
log using "`directory'\Table2.log", replace

****************************************************************
* PROGRAM: STAR_Table2
* PROGRAMMER: Bruno Ferman (based on program by written Simone Schaner)
* PURPOSE: Makes Table 2 of Angrist, Lang and Oreopoulos (2008) 
*****************************************************************


clear
cd "`directory'"
set mem 50m
set more off
set linesize 200

use STAR_public_use, clear	

keep if noshow==0
	
g byte group1= 1
g byte group2= sex=="M"
g byte group3= sex=="F"

* CONTROLS
local basic i.sex i.mtongue i.hsgroup i.numcourses_nov1
local all   i.sex i.mtongue i.hsgroup  i.numcourses_nov1 i.lastmin i.mom_edn i.dad_edn

mat t2=J(30,6,.)
local i=0
local j=1

* DUMMY FOR HAS FALL GRADES
g byte hasfall= grade_20059_fall!=.

foreach group in 1 2 3 {
	foreach var in totcredits_year1 mathsci hasfall {
	g var`i'=`var'
	sum var`i' if cont & group`group'
		mat t2[`j'+1,`i'*2+1]=r(mean)
		mat t2[`j'+2,`i'*2+1]=r(sd)
	xi: reg var`i' ssp sfp sfsp `basic' if group`group', r
		mat t2[`j'+3,`i'*2+1]=_b[ssp]
		mat t2[`j'+4,`i'*2+1]=_se[ssp]
		mat t2[`j'+5,`i'*2+1]=_b[sfp]
		mat t2[`j'+6,`i'*2+1]=_se[sfp]
		mat t2[`j'+7,`i'*2+1]=_b[sfsp]
		mat t2[`j'+8,`i'*2+1]=_se[sfsp]
		mat t2[`j'+9,`i'*2+1]=e(N)
	xi: reg var`i' ssp sfp sfsp `all' if group`group', r
		mat t2[`j'+3,`i'*2+2]=_b[ssp]
		mat t2[`j'+4,`i'*2+2]=_se[ssp]
		mat t2[`j'+5,`i'*2+2]=_b[sfp]
		mat t2[`j'+6,`i'*2+2]=_se[sfp]
		mat t2[`j'+7,`i'*2+2]=_b[sfsp]
		mat t2[`j'+8,`i'*2+2]=_se[sfsp]
		mat t2[`j'+9,`i'*2+2]=e(N)
	drop var`i'
	local ++i
	}
	local i=0
	local j=`j'+10
	}
	
	mat colnames t2= numc_bas numc_all msci_bas msci_all csrv_bas csvr_all 
	mat rownames t2= all Cont_Mean "x" ssp "x" sfp "x" sfsp "x" obs ///
		males Cont_Mean "x" ssp "x" sfp "x" sfsp "x" obs ///
		females Cont_Mean "x" ssp "x" sfp "x" sfsp "x" obs 
	mat li t2
	
log close

