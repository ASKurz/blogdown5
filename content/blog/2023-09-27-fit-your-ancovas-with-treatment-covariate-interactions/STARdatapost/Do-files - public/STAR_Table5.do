capture log close
local d Z
local directory DIRECTORY HERE
log using "`directory'\Table5.log", replace

****************************************************************
* PROGRAM: STAR_Table5
* PROGRAMMER: Bruno Ferman (based on program written by Simone Schaner)
* PURPOSE: Makes Table 5 of Angrist, Lang and Oreopoulos (2008) 
*****************************************************************

clear
cd "`directory'"
set mem 50m
set more off
set linesize 200

use STAR_public_use, clear	

keep if noshow==0

local basic  i.sex i.mtongue i.hsgroup i.numcourses_nov1
local all   i.sex i.mtongue i.hsgroup i.numcourses_nov1 i.lastmin i.mom_edn i.dad_edn

mat t5=J(24,6,.)

local i=0
local j=1

g byte group1= 1
g byte group2= sex=="M"
g byte group3= sex=="F"
g byte sspany= ssp | sfsp

replace GPA_year1=. if grade_20059_fall==. /*ONLY DO GPA FOR FALL GRADE SAMPLE*/
replace grade_20059_fall=. if GPA_year1==. /*ONLY DO FALL GRADE FOR GPA SAMPLE*/

foreach var in grade_20059_fall GPA_year1 {
	foreach group in 1 2 3 {
	g var`i'=`var'
	xi: reg var`i' ssp sfp sfsp `all' if group`group', r
		mat t5[`j'+3,`group']=_b[ssp]
		mat t5[`j'+4,`group']=_se[ssp]
		mat t5[`j'+7,`group']=_b[sfp]
		mat t5[`j'+8,`group']=_se[sfp]
		mat t5[`j'+9,`group']=_b[sfsp]
		mat t5[`j'+10,`group']=_se[sfsp]
		mat t5[`j'+11,`group']=e(N)
	qui sum var`i' if e(sample) & cont
		mat t5[`j'+1,`group']=r(mean)
		mat t5[`j'+2,`group']=r(sd)
		
	xi: reg var`i' ssp sfpany `all' if group`group', r
		mat t5[`j'+3,`group'+3]=_b[ssp]
		mat t5[`j'+4,`group'+3]=_se[ssp]
		mat t5[`j'+5,`group'+3]=_b[sfp]
		mat t5[`j'+6,`group'+3]=_se[sfp]
		mat t5[`j'+11,`group'+3]=e(N)
	qui sum var`i' if e(sample) & cont
		mat t5[`j'+1,`group'+3]=r(mean)
		mat t5[`j'+2,`group'+3]=r(sd)
		
		drop var`i'
	local ++i
	}	
	local j=`j'+12
	}
	
mat colnames t5= all males females all males females 
mat rownames t5= A_fallgr cmean "x" ssp "x" sfpany "x" sfp "x" sfsp "x" obs ///
				  B_GPAyr1 cmean "x" ssp "x" sfpany "x" sfp "x" sfsp "x" obs	

mat li t5
log close
