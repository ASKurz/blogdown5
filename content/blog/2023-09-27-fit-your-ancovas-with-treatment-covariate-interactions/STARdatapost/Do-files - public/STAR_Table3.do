capture log close
local d Z
local directory DIRECTORY HERE
log using "`directory'\Table3.log", replace

****************************************************************
* PROGRAM: STAR_Table3
* PROGRAMMER: Bruno Ferman (based on program written by Simone Schaner)
* PURPOSE: Makes Table 3 of Angrist, Lang and Oreopoulos (2008) 
*****************************************************************

clear
cd "`directory'"
set mem 50m
set more off
set linesize 200

use STAR_public_use, clear	

keep if noshow==0

* CONTROLS
local basic  i.sex i.mtongue i.hsgroup i.numcourses_nov1
local all   i.sex i.sex i.mtongue i.hsgroup i.numcourses_nov1 i.lastmin i.mom_edn i.dad_edn

mat t3=J(26,8,.)

local i=0
foreach var in signup {

g var`i'=`var'

	xi: reg var`i' ssp sfp sfsp `basic', r
		mat t3[2,`i'*2+1]=_b[ssp]
		mat t3[3,`i'*2+1]=_se[ssp]
		mat t3[4,`i'*2+1]=_b[sfp]
		mat t3[5,`i'*2+1]=_se[sfp]
		mat t3[6,`i'*2+1]=_b[sfsp]
		mat t3[7,`i'*2+1]=_se[sfsp]
		mat t3[8,`i'*2+1]=e(N)
	qui sum var`i' if e(sample) & cont 
		mat t3[1,`i'*2+1]=r(mean)
		
	xi: reg var`i' ssp sfp sfsp `all', r
		mat t3[2,`i'*2+2]=_b[ssp]
		mat t3[3,`i'*2+2]=_se[ssp]
		mat t3[4,`i'*2+2]=_b[sfp]
		mat t3[5,`i'*2+2]=_se[sfp]
		mat t3[6,`i'*2+2]=_b[sfsp]
		mat t3[7,`i'*2+2]=_se[sfsp]
		mat t3[8,`i'*2+2]=e(N)
	qui sum var`i' if e(sample) & cont
		mat t3[1,`i'*2+2]=r(mean)
		
	xi: reg var`i' ssp sfp sfsp `basic' if sex=="M", r
		mat t3[11,`i'*2+1]=_b[ssp]
		mat t3[12,`i'*2+1]=_se[ssp]
		mat t3[13,`i'*2+1]=_b[sfp]
		mat t3[14,`i'*2+1]=_se[sfp]
		mat t3[15,`i'*2+1]=_b[sfsp]
		mat t3[16,`i'*2+1]=_se[sfsp]
		mat t3[17,`i'*2+1]=e(N)
	qui sum var`i' if e(sample) & cont
		mat t3[10,`i'*2+1]=r(mean)
	
	xi: reg var`i' ssp sfp sfsp `all' if sex=="M", r
		mat t3[11,`i'*2+2]=_b[ssp]
		mat t3[12,`i'*2+2]=_se[ssp]
		mat t3[13,`i'*2+2]=_b[sfp]
		mat t3[14,`i'*2+2]=_se[sfp]
		mat t3[15,`i'*2+2]=_b[sfsp]
		mat t3[16,`i'*2+2]=_se[sfsp]
		mat t3[17,`i'*2+2]=e(N)
		qui sum var`i' if e(sample) & cont
		mat t3[10,`i'*2+2]=r(mean)

	xi: reg var`i' ssp sfp sfsp `basic' if sex=="F", r
		mat t3[20,`i'*2+1]=_b[ssp]
		mat t3[21,`i'*2+1]=_se[ssp]
		mat t3[22,`i'*2+1]=_b[sfp]
		mat t3[23,`i'*2+1]=_se[sfp]
		mat t3[24,`i'*2+1]=_b[sfsp]
		mat t3[25,`i'*2+1]=_se[sfsp]
		mat t3[26,`i'*2+1]=e(N)
	qui sum var`i' if e(sample) & cont
		mat t3[19,`i'*2+1]=r(mean)

	xi: reg var`i' ssp sfp sfsp `all' if sex=="F", r
		mat t3[20,`i'*2+2]=_b[ssp]
		mat t3[21,`i'*2+2]=_se[ssp]
		mat t3[22,`i'*2+2]=_b[sfp]
		mat t3[23,`i'*2+2]=_se[sfp]
		mat t3[24,`i'*2+2]=_b[sfsp]
		mat t3[25,`i'*2+2]=_se[sfsp]
		mat t3[26,`i'*2+2]=e(N)
	qui sum var`i' if e(sample) & cont
		mat t3[19,`i'*2+2]=r(mean)
		drop var`i'
		local ++i
	}
	

foreach var in used_ssp used_adv used_fsg {

g var`i'=`var'

	xi: reg var`i' ssp sfp sfsp `basic', r
		mat t3[2,`i'*2+1]=_b[ssp]
		mat t3[3,`i'*2+1]=_se[ssp]
		mat t3[6,`i'*2+1]=_b[sfsp]
		mat t3[7,`i'*2+1]=_se[sfsp]
		mat t3[8,`i'*2+1]=e(N)
	qui sum var`i' if e(sample) & cont 
		mat t3[1,`i'*2+1]=r(mean)
		
	xi: reg var`i' ssp sfp sfsp `all', r
		mat t3[2,`i'*2+2]=_b[ssp]
		mat t3[3,`i'*2+2]=_se[ssp]
		mat t3[6,`i'*2+2]=_b[sfsp]
		mat t3[7,`i'*2+2]=_se[sfsp]
		mat t3[8,`i'*2+2]=e(N)
	qui sum var`i' if e(sample) & cont
		mat t3[1,`i'*2+2]=r(mean)
		
	xi: reg var`i' ssp sfp sfsp `basic' if sex=="M", r
		mat t3[11,`i'*2+1]=_b[ssp]
		mat t3[12,`i'*2+1]=_se[ssp]
		mat t3[15,`i'*2+1]=_b[sfsp]
		mat t3[16,`i'*2+1]=_se[sfsp]
		mat t3[17,`i'*2+1]=e(N)
	qui sum var`i' if e(sample) & cont
		mat t3[10,`i'*2+1]=r(mean)
	
	xi: reg var`i' ssp sfp sfsp `all' if sex=="M", r
		mat t3[11,`i'*2+2]=_b[ssp]
		mat t3[12,`i'*2+2]=_se[ssp]
		mat t3[15,`i'*2+2]=_b[sfsp]
		mat t3[16,`i'*2+2]=_se[sfsp]
		mat t3[17,`i'*2+2]=e(N)
		qui sum var`i' if e(sample) & cont
		mat t3[10,`i'*2+2]=r(mean)

	xi: reg var`i' ssp sfp sfsp `basic' if sex=="F", r
		mat t3[20,`i'*2+1]=_b[ssp]
		mat t3[21,`i'*2+1]=_se[ssp]
		mat t3[24,`i'*2+1]=_b[sfsp]
		mat t3[25,`i'*2+1]=_se[sfsp]
		mat t3[26,`i'*2+1]=e(N)
	qui sum var`i' if e(sample) & cont
		mat t3[19,`i'*2+1]=r(mean)

	xi: reg var`i' ssp sfp sfsp `all' if sex=="F", r
		mat t3[20,`i'*2+2]=_b[ssp]
		mat t3[21,`i'*2+2]=_se[ssp]
		mat t3[24,`i'*2+2]=_b[sfsp]
		mat t3[25,`i'*2+2]=_se[sfsp]
		mat t3[26,`i'*2+2]=e(N)
	qui sum var`i' if e(sample) & cont
		mat t3[19,`i'*2+2]=r(mean)
		drop var`i'
		local ++i
	}


	
mat colnames t3= sign_bas sign_all ssp_bas ssp_all adv_bas adv_all fsg_bas fsg_all
mat rownames t3= contmean ssp "x" sfp "x" sfsp "x" obs males contmean ssp "x" sfp "x" sfsp ///
				"x" obs females contmean ssp "x" sfp "x" sfsp "x" obs


mat li t3
log close
