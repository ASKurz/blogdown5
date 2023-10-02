capture log close
local d Z
local directory DIRECTORY HERE
log using "`directory'\Table1.log", replace

****************************************************************
* PROGRAM: STAR_Table1
* PROGRAMMER: Bruno Ferman (based on program written by Simone Schaner)
* PURPOSE: Makes Table 1 of Angrist, Lang and Oreopoulos (2008) 
*****************************************************************

clear
cd "`directory'"
set mem 50m
set more off
set linesize 200

use STAR_public_use, clear	
	
mat t1_all=J(36,6,.)
mat t1_female=J(36,6,.)
mat t1_male=J(36,6,.)

local i=0
foreach var in numcourses_nov1 noshow compsurv {
	g var`i'=`var'
	sum var`i' if control
		mat t1_all[`i'*2+1,1]=r(mean)
		mat t1_all[`i'*2+2,1]=r(sd)

	reg var`i' ssp sfp sfsp
		mat t1_all[`i'*2+1,2]=_b[ssp]
		mat t1_all[`i'*2+2,2]=_se[ssp]

		mat t1_all[`i'*2+1,3]=_b[sfp]
		mat t1_all[`i'*2+2,3]=_se[sfp]
		
		mat t1_all[`i'*2+1,4]=_b[sfsp]
		mat t1_all[`i'*2+2,4]=_se[sfsp]
		mat t1_all[`i'*2+1,6]=e(N)
		
		test ssp sfp sfsp
		mat t1_all[`i'*2+1,5]=r(F)
		mat t1_all[`i'*2+2,5]=r(p)

	local ++i		
}


foreach var in female gpa0 age english hcom chooseUTM work1 mom1 mom2 /// 
	  dad1 dad2 lm_rarely lm_never graddeg finish4 {
	g var`i'=`var'
	sum var`i' if control&noshow==0
		mat t1_all[`i'*2+1,1]=r(mean)
		mat t1_all[`i'*2+2,1]=r(sd)

	reg var`i' ssp sfp sfsp if noshow==0
		mat t1_all[`i'*2+1,2]=_b[ssp]
		mat t1_all[`i'*2+2,2]=_se[ssp]

		mat t1_all[`i'*2+1,3]=_b[sfp]
		mat t1_all[`i'*2+2,3]=_se[sfp]
		
		mat t1_all[`i'*2+1,4]=_b[sfsp]
		mat t1_all[`i'*2+2,4]=_se[sfsp]
		mat t1_all[`i'*2+1,6]=e(N)
		
		test ssp sfp sfsp
		mat t1_all[`i'*2+1,5]=r(F)
		mat t1_all[`i'*2+2,5]=r(p)

	local ++i		
}


mat colnames t1_all=Cont_Mean SSPvCont SFPvCont SFSPvCont Fstat Obs
mat rownames t1_all=numcourses_nov1 "x" noshow "x" compsurv "x" female "x" gpa0 "x" age "x" english "x" ///
	LiveHome x chooseUTM "x" PWork "x" mom1 "x" mom2 "x" dad1 "x" dad2 "x" lm_rarely "x" lm_never "x" ///
	graddeg x finish4 x


mat li t1_all

log close

