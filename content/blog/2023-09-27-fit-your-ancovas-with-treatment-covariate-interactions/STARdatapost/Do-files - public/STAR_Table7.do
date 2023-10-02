capture log close
local d Z
local directory DIRECTORY HERE
log using "`directory'\Table7.log", replace

****************************************************************
* PROGRAM: STAR_Table7
* PROGRAMMER: Bruno Ferman (based on program written by Simone Schaner)
* PURPOSE: Makes Table 7 of Angrist, Lang and Oreopoulos (2008) 
*****************************************************************

clear
cd "`directory'"
set mem 50m
set more off
set linesize 200

use STAR_public_use, clear	

* SAMPLE IN ALL CONTROLS SPEC
keep if sex~="" & mtongue~="" & hsgroup~=. & numcourses_nov1~=. & lastmin~=. & mom_edn~=. & dad_edn~=.


* CONTROL SETS
local all i.sex i.mtongue i.hsgroup i.numcourses_nov1 i.lastmin i.mom_edn i.dad_edn i.year


keep if noshow==0

keep   sex mtongue hsgroup numcourses_nov1 lastmin mom_edn dad_edn GPA_year1 GPA_year2 ///
       ssp sfp sfsp control

gen id=_n

* STACK THE DATA
reshape long GPA_year, i(id) j(year)


mat t7=J(30,6,.)

local i=1
local j=1

g byte group1= 1
g byte group2= sex=="M"
g byte group3= sex=="F"

foreach group in 2 3 {
	
	xi: reg GPA_year ssp sfp sfsp `all'  if group`group', clu(id)
		mat t7[`j'+3,`i']=_b[ssp]
		mat t7[`j'+4,`i']=_se[ssp]
		mat t7[`j'+7,`i']=_b[sfp]
		mat t7[`j'+8,`i']=_se[sfp]
		mat t7[`j'+9,`i']=_b[sfsp]
		mat t7[`j'+10,`i']=_se[sfsp]
		mat t7[`j'+11,`i']=e(N)
		
		qui sum GPA_year if e(sample) & cont, detail
				mat t7[`j'+1,`i']=r(mean)
				mat t7[`j'+2,`i']=r(sd)
		local ++i
		
		foreach quantile in .1 .25 .5 .75 .9 {

	xi: bootstrap _b[ssp] _b[sfp] _b[sfsp],  reps(500) cluster(id): qreg GPA_year ssp sfp sfsp `all' if group`group', q(`quantile')
		mat t7[`j'+3,`i']=_b[_bs_1]
		mat t7[`j'+4,`i']=_se[_bs_1]
		mat t7[`j'+7,`i']=_b[_bs_2]
		mat t7[`j'+8,`i']=_se[_bs_2]
		mat t7[`j'+9,`i']=_b[_bs_3]
		mat t7[`j'+10,`i']=_se[_bs_3]
		mat t7[`j'+11,`i']=e(N)
		
	qui sum GPA_year if e(sample) & cont, detail
		local val=`quantile'*100
		
		mat t7[`j'+1,`i']=r(p`val')
		*mat t7[`j'+2,`i']=r(sd)
	local ++i
	}	
	local j=`j'+12
	local i=1
	}

local i=1
local j=26

xi: reg GPA_year ssp sfp sfsp i.year i.hsgroup if group3, clu(id)
	mat t7[`j',`i']=_b[sfsp]
	mat t7[`j'+1,`i']=_se[sfsp]

xi: reg GPA_year ssp sfp sfsp i.year   if group3, clu(id)
	mat t7[`j'+2,`i']=_b[sfsp]
	mat t7[`j'+3,`i']=_se[sfsp]
	mat t7[`j'+4,`i']=e(N)

local ++i

foreach quantile in .1 .25 .5 .75 .9 {

xi: bootstrap _b[ssp] _b[sfp] _b[sfsp],  reps(500) cluster(id): qreg GPA_year ssp sfp sfsp i.year i.hsgroup if group3, q(`quantile')
	mat t7[`j',`i']=_b[_bs_3]
	mat t7[`j'+1,`i']=_se[_bs_3]

xi: bootstrap _b[ssp] _b[sfp] _b[sfsp],  reps(500) cluster(id): qreg GPA_year ssp sfp sfsp i.year if group3, q(`quantile')
	mat t7[`j'+2,`i']=_b[_bs_3]
	mat t7[`j'+3,`i']=_se[_bs_3]
	mat t7[`j'+4,`i']=e(N)
	local ++i
	}	

mat colnames t7= yr1_ols yr1_q10 yr1_q25 yr1_q50 yr1_q75 yr1_q90 

mat rownames t7= A_GPA_male cmean "x" ssp "x" sfpany "x" sfp "x" sfsp "x" obs ///
		     B_GPA_fem cmean "x" ssp "x" sfpany "x" sfp "x" sfsp "x" obs ///
                 C_Limited_Cov SFSP_year "x" SFSP_year_HS "x" obs
				  
mat li t7

log close
