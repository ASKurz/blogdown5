capture log close
local d Z
local directory DIRECTORY HERE
log using "`directory'\Table6.log", replace

****************************************************************
* PROGRAM: STAR_Table6
* PROGRAMMER: Bruno Ferman (based on program written by Simone Schaner)
* PURPOSE: Makes Table 6 of Angrist, Lang and Oreopoulos (2008) 
*****************************************************************

clear
cd "`directory'"
set mem 50m
set more off
set linesize 200

use STAR_public_use, clear	

keep if noshow==0

* ADDITIONAL SAMPLE SELECTION
foreach var in prob_year goodstanding_year credits_earned {
	foreach num in 1 2 {
	replace `var'`num'=. if prob_year1==.
	}
	}

* SAMPLE IN ALL CONTROLS SPEC
keep if sex~="" & mtongue~="" & hsgroup~=. & numcourses_nov1~=. & lastmin~=. & mom_edn~=. & dad_edn~=.

* CONTROLS
local all    i.sex i.mtongue i.hsgroup i.numcourses_nov1 i.lastmin i.mom_edn i.dad_edn

mat t6=J(40,6,.)

local i=0
local j=1

g byte group1= 1
g byte group2= sex=="M"
g byte group3= sex=="F"

foreach var in GPA_year1 prob_year1 goodstanding_year1 credits_earned1  {
	foreach group in 1 2 3 {
	g var`i'=`var'
	xi: reg var`i' ssp sfp sfsp `all' if group`group', r
		mat t6[`j'+3,`group']=_b[ssp]
		mat t6[`j'+4,`group']=_se[ssp]
		mat t6[`j'+5,`group']=_b[sfp]
		mat t6[`j'+6,`group']=_se[sfp]
		mat t6[`j'+7,`group']=_b[sfsp]
		mat t6[`j'+8,`group']=_se[sfsp]
		mat t6[`j'+9,`group']=e(N)
	qui sum var`i' if e(sample) & cont
		mat t6[`j'+1,`group']=r(mean)
		mat t6[`j'+2,`group']=r(sd)
	drop var`i'
	local ++i
	}	
	local j=`j'+10
	}
	
local i=0
local j=1

foreach var in GPA_year2 prob_year2 goodstanding_year2 credits_earned2  {
	foreach group in 1 2 3 {
		g var`i'=`var'
	xi: reg var`i' ssp sfp sfsp `all' if group`group', r
		mat t6[`j'+3,`group'+3]=_b[ssp]
		mat t6[`j'+4,`group'+3]=_se[ssp]
		mat t6[`j'+5,`group'+3]=_b[sfp]
		mat t6[`j'+6,`group'+3]=_se[sfp]
		mat t6[`j'+7,`group'+3]=_b[sfsp]
		mat t6[`j'+8,`group'+3]=_se[sfsp]
		mat t6[`j'+9,`group'+3]=e(N)
	qui sum var`i' if e(sample) & cont
		mat t6[`j'+1,`group'+3]=r(mean)
		mat t6[`j'+2,`group'+3]=r(sd)
	drop var`i'
	local ++i
	}	
	local j=`j'+10
	}
	

mat rownames t6= A_GPA cmean "x" ssp "x"  sfp "x" sfsp "x" obs ///
				  B_Prob cmean "x" ssp "x"  sfp "x" sfsp "x" obs	///
				  C_GoodSt cmean "x" ssp "x"  sfp "x" sfsp "x" obs ///
				  D_CredErn cmean "x" ssp "x"  sfp "x" sfsp "x" obs
				  				  
mat li t6

log close
				  
