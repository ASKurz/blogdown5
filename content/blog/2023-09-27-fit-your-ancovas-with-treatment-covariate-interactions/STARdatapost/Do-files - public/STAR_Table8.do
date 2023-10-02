capture log close
local d Z
local directory DIRECTORY HERE
log using "`directory'\Table8.log", replace

****************************************************************
* PROGRAM: STAR_Table8
* PROGRAMMER: Bruno Ferman (based on program written by Simone Schaner)
* PURPOSE: Makes Table 8 of Angrist, Lang and Oreopoulos (2008) 
*****************************************************************

clear
cd "`directory'"
set mem 50m
set more off
set linesize 200

use STAR_public_use, clear	

keep if noshow==0
keep if sex=="F" /*This table uses only women*/

* SAMPLE IN ALL CONTROLS SPEC
keep if sex~="" & mtongue~="" & hsgroup~=. & numcourses_nov1~=. & lastmin~=. & mom_edn~=. & dad_edn~=.

* CONTROL SETS
local all    i.sex i.mtongue i.hsgroup i.numcourses_nov1 i.lastmin i.mom_edn i.dad_edn i.year


foreach var in prob_year goodstanding_year credits_earned {
	foreach num in 1 2 {
	replace `var'`num'=. if prob_year1==.
	}
	}

keep   sex mtongue hsgroup numcourses_nov1 lastmin mom_edn dad_edn GPA_year1 GPA_year2 ///
       ssp sfp sfsp control prob_year1 prob_year2 credits_earned* ssp_p sfp_p sfsp_p sfpany_p


gen id=_n
	
* STACK THE DATA
reshape long GPA_year goodstanding_year prob_year credits_earned, i(id) j(year)

local i=0
mat t8=J(13,6,.)

foreach var in GPA_year prob_year credits_earned {
		g var`i'=`var'
	xi: ivreg var`i' (ssp_p sfp_p sfsp_p = ssp sfsp sfp) `all' , r clu(id)
		mat t8[3,`i'*2+1]=_b[ssp]
		mat t8[4,`i'*2+1]=_se[ssp]
		mat t8[5,`i'*2+1]=_b[sfp]
		mat t8[6,`i'*2+1]=_se[sfp]
		mat t8[7,`i'*2+1]=_b[sfsp]
		mat t8[8,`i'*2+1]=_se[sfsp]
		mat t8[13,`i'*2+1]=e(N)
	qui sum var`i' if e(sample) & cont
		mat t8[1,`i'*2+1]=r(mean)
		mat t8[2,`i'*2+1]=r(sd)
	xi: ivreg2 var`i' (ssp_p sfpany_p= ssp sfsp sfp) `all'  , r clu(id) partial(`all')
		mat t8[3,`i'*2+2]=_b[ssp]
		mat t8[4,`i'*2+2]=_se[ssp]
		mat t8[9,`i'*2+2]=_b[sfp]
		mat t8[10,`i'*2+2]=_se[sfp]
		mat t8[13,`i'*2+2]=e(N)
	qui sum var`i' if e(sample) & cont
		mat t8[1,`i'*2+2]=r(mean)
		mat t8[2,`i'*2+2]=r(sd)		
* OVERID TEST
		mat t8[11,`i'*2+2]=e(j) /*e(N)*e(r2)*/ 
		mat t8[12,`i'*2+2]=e(jp) /*pval of test*/
	drop var`i'
	local ++i
	}	
	
mat colnames t8 = GPA GPA Prob Prob CredErn CredErn
mat rownames t8 = ContMean "x" SSP "x" SFP "x" SFSP "x" SFSP_any "x" overid pvalue obs

mat li t8

log close

