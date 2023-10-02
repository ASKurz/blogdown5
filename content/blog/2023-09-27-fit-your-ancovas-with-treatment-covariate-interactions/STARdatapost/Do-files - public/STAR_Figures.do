capture log close
local d Z
local directory DIRECTORY HERE
log using "`directory'\Figures.log", replace

****************************************************************
* PROGRAM: STAR_Figures
* PROGRAMMER: Bruno Ferman (based on program written by Simone Schaner)
* PURPOSE: Makes figures of Angrist, Lang and Oreopoulos (2008) 
*****************************************************************

clear
cd "`directory'"
set mem 50m
set more off
set linesize 200

use STAR_public_use, clear	

keep if noshow==0

gen group=0 if control==1
replace group=1 if ssp==1
replace group=2 if sfp==1
replace group=3 if sfsp==1

gen star=ssp+sfp+sfsp

* WHETHER OR NOT TO INCLUDE XLINES
local xline xline(2.7 3.0 3.3, lstyle(grid))

* GET K-S P VALUES
preserve
keep if group<=1 & sex=="M"
ksmirnov GPA_year1, by(star)
	local M_ssp=round(r(p_cor),.001)
restore
preserve
keep if group==0 | group==2 & sex=="M"
ksmirnov GPA_year1, by(star)
	local M_sfp=round(r(p_cor),.001)
restore
preserve
keep if group==0 | group==3 & sex=="M"
ksmirnov GPA_year1, by(star)
	local M_sfsp=round(r(p_cor),.001)
restore

preserve
keep if group<=1 & sex=="F"
ksmirnov GPA_year1, by(star)
	local F_ssp=round(r(p_cor),.001)
restore
preserve
keep if group==0 | group==2 & sex=="F"
ksmirnov GPA_year1, by(star)
	local F_sfp=round(r(p_cor),.001)
restore
preserve
keep if group==0 | group==3 & sex=="F"
ksmirnov GPA_year1, by(star)
	local F_sfsp=round(r(p_cor),.001)
restore

* "NORMALIZE" GPA VALUES
 replace GPA_year1= GPA_year1+.3 if hsgroup==1
 replace GPA_year1= GPA_year1-.3 if hsgroup==3

 
foreach sex in M F {
kdensity GPA_year1 if group==0 & sex=="`sex'" & GPA_year1>=0, g(xaxis01`sex' dens01`sex') gauss clwidth(thick)  nodraw
kdensity GPA_year1 if (group==1 & sex=="`sex'") & GPA_year1>=0,  g(xaxis02`sex' dens02`sex') gauss clwidth(thick) nodraw

label var dens01`sex' "Control"
label var dens02`sex' "SSP"

twoway (connected dens01`sex' xaxis01`sex' if xaxis01`sex'>=0 & xaxis01`sex'<=4, sort msymbol(none) clcolor(lavender) clwidth(vthick)) /// 
	   (connected dens02`sex' xaxis02`sex' if xaxis02`sex'>=0 & xaxis02`sex'<=4, sort msymbol(none) clcolor(dknavy) clpat(solid) clwidth(thick)), ///
		xlab(0(4)4 .3 .7 1 1.3 1.7 2 2.3 2.7 3 3.3 3.7 4) xti("First Year GPA") name(Fig1PanelA_`sex', replace) ///
		title("Panel A: Control vs. SSP", size(medium)) `xline' /// 
		text(.02 2.7 "1000", place(nw) si(small)) text(.02 2.7 "cutoff", place(sw) si(small)) ///
		text(.02 3 "2500", place(n) si(small)) text(.02 3 "cutoff", place(s) si(small)) ///
		text(.02 3.3 "5000", place(ne) si(small)) text(.02 3.3 "cutoff", place(se) si(small)) nodraw
*		note("K-S p-value: 0``sex'_ssp'", pos(6) ring(5) just(center) size(medsmall)) nodraw		
drop dens01`sex' dens02`sex' xaxis01`sex' xaxis02`sex'
	
kdensity GPA_year1 if group==0 & sex=="`sex'" & GPA_year1>=0, g(xaxis01`sex' dens01`sex') gauss clwidth(thick) nodraw
kdensity GPA_year1 if (group==2 & sex=="`sex'") & GPA_year1>=0, g(xaxis02`sex' dens02`sex') gauss clwidth(thick) nodraw

label var dens01`sex' "Control"
label var dens02`sex' "SFP"

twoway (connected dens01`sex' xaxis01`sex', sort msymbol(none) clcolor(lavender) clwidth(vthick))  ///
		(connected dens02`sex' xaxis02`sex', sort msymbol(none) clcolor(dknavy) clpat(solid) clwidth(thick)), /// 
		xlab(0(4)4 .3 .7 1 1.3 1.7 2 2.3 2.7 3 3.3 3.7 4) xti("First Year GPA") ///
		name(Fig1PanelB_`sex', replace) title("Panel B: Control vs. SFP", size(medium)) ///
		text(.02 2.7 "1000", place(nw) si(small)) text(.02 2.7 "cutoff", place(sw) si(small)) ///
		text(.02 3 "2500", place(n) si(small)) text(.02 3 "cutoff", place(s) si(small)) `xline' ///
		text(.02 3.3 "5000", place(ne) si(small)) text(.02 3.3 "cutoff", place(se) si(small)) nodraw
*		note("K-S p-value: 0``sex'_sfp'", pos(6) ring(5) just(center) size(medsmall))
		
drop dens01`sex' dens02`sex' xaxis01`sex' xaxis02`sex'

kdensity GPA_year1 if group==0 & sex=="`sex'" & GPA_year1>=0, g(xaxis01`sex' dens01`sex') gauss clwidth(thick) nodraw
kdensity GPA_year1 if (group==3 & sex=="`sex'") & GPA_year1>=0, g(xaxis02`sex' dens02`sex') gauss clwidth(thick) nodraw

label var dens01`sex' "Control"
label var dens02`sex' "SFSP"

twoway (connected dens01`sex' xaxis01`sex', sort msymbol(none) clcolor(lavender) clwidth(vthick))  ///
		(connected dens02`sex' xaxis02`sex', sort msymbol(none) clcolor(dknavy) clpat(solid) clwidth(thick)), /// 
		xlab(0(4)4 .3 .7 1 1.3 1.7 2 2.3 2.7 3 3.3 3.7 4) xti("First Year GPA") ///
		name(Fig1PanelC_`sex', replace) title("Panel B: Control vs. SFSP", size(medium)) ///
		text(.02 2.7 "1000", place(nw) si(small)) text(.02 2.7 "cutoff", place(sw) si(small)) ///
		text(.02 3 "2500", place(n) si(small)) text(.02 3 "cutoff", place(s) si(small)) ///
		text(.02 3.3 "5000", place(ne) si(small)) text(.02 3.3 "cutoff", place(se) si(small)) `xline' nodraw
*		note("K-S p-value: 0``sex'_sfsp'", pos(6) ring(5) just(center) size(medsmall)) 
		
drop dens01`sex' dens02`sex' xaxis01`sex' xaxis02`sex'

}

* COMBINE AND SAVE GRAPHS
gr combine Fig1PanelA_M Fig1PanelB_M Fig1PanelC_M, rows(1) col(3) ysize(6) xsize(13) saving(fig1a_males3wy_v4, replace) ///
	title("Figure 1a. Males' Normalized First-year GPA", size(medium)) ///
	note("Notes: These figures plot the smoothed kernel densities of first year GPA." ///
		"The K-S p-value is a test for equality of distributions.")
gr combine Fig1PanelA_F Fig1PanelB_F Fig1PanelC_F, rows(1) col(3) ysize(6) xsize(13) saving(fig1b_females3wy_v4, replace) ///
	title("Figure 1b. Females' Normalized First-year GPA", size(medium)) ///
		note("Notes: These figures plot the smoothed kernel densities of first year GPA.")
		
log close

gr use fig1a_males3wy_v4
gr use fig1b_females3wy_v4
