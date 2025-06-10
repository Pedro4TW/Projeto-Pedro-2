************* SET PATH *************	
global user = 1 // 1 = Mac; 2 = PC	

* set ROOT to the corresponding path where the replicationi folder is saved

if $user==1{
	cd "ROOT/codes/"
	global raw_calibration = "ROOT/codes/country grids/summary_tables/"
	global load_gravity = "ROOT/codes/country grids/gravity"
	global save_draft = "ROOT/draft/final tables/"
	global save_draft2 = "ROOT/draft/"
}
if $user==2{
	cd "ROOT\codes\"
	global raw_calibration = "ROOT\codes\country grids\summary_tables\"
	global load_gravity = "ROOT\codes\country grids\gravity"
	global save_draft = "ROOT\draft\final tables\"
	global save_draft2 = "ROOT\draft\"
}

import excel using "$raw_calibration\all_cfactuals.xls",sheet("all_results") firstrow clear

********** keep cong = "on" **********
keep if cong=="on"

********** keep cfac **********
keep if cfac_type=="exp_eng" || cfac_type=="exp_til"  || cfac_type=="mis_eng" 

* define mobility as 0/1
g temp=Mobility
drop Mobility
g Mobility=1 if temp=="on"
replace Mobility=0 if temp=="off"
drop temp

* define some outcome variables
g y_calib = Y_calib/L_calib
g y_exp = Y_exp/L_exp

g c_calib = C_calib/L_calib
g c_exp = C_exp/L_exp

g XMY_exp =  XY_exp+ MY_exp
g XMY_calib =  XY_calib+ MY_calib

* define parameter space
g low_beta = (beta==0.13)
g convex = (beta>gamma)

* unique identifier of type of counterfactual
egen cfactual_type_id = group( convex cfac_type Ngoods city_allocation Mobility )
isid cfactual_type_id Country location
bys cfactual_type_id Country location: g rep=_N
drop cfactual_type_id

* check number of counterfactuals for each case
egen tag_cfac = tag( cfactual_id ICC )

table Country cfac_type Ngoods if convex==0 & Mobility==0 & city_allocation=="largest cities", c(sum tag_cfac)
table Country cfac_type Ngoods if convex==1 & Mobility==0 & city_allocation=="largest cities", c(sum tag_cfac)
table Country cfac_type Ngoods if convex==0 & Mobility==1 & city_allocation=="largest cities", c(sum tag_cfac)
table Country cfac_type Ngoods if convex==1 & Mobility==1 & city_allocation=="largest cities", c(sum tag_cfac)

* these are complete
table Country cfac_type Ngoods if convex==1 & Mobility==0 & city_allocation=="nuts", c(sum tag_cfac)
table Country cfac_type Ngoods if convex==1 & Mobility==1 & city_allocation=="nuts", c(sum tag_cfac)

* order and sort
order cfactual_id ICC location
sort cfactual_id ICC location
isid cfactual_id ICC location

* benchmark
g benchmark = convex & Mobility & Ngoods==10 & city_allocation=="largest cities"

* verify that there are 24 counterfactuals for each benchmark case:
table benchmark cfac_type if benchmark, c(sum tag_cfac)

* normalize calibrated Z, H and c
g h_calib = H_calib/L_calib
	
bysort cfactual_id ICC: egen avZ = mean(Z_calib)  // average Z within country
by cfactual_id ICC: egen avh = mean(h_calib)    // average h within country
by cfactual_id ICC: egen avH = mean(H_calib)    // average h within country
by cfactual_id ICC: egen avc_calib = mean(c_calib) // average c within country
by cfactual_id ICC: egen avc_exp = mean(c_exp)  // average c_exp within country

replace Z_calib=Z_calib/avZ
replace h_calib=h_calib/avh
replace H_calib=H_calib/avH
replace c_calib=c_calib/avc_calib
replace c_exp=c_exp/avc_exp
	
* express welfare gain in % points
replace welfare_gain = ( welfare_gain-1 )*100	

local vars = " XY_calib MY_calib XY_exp MY_exp XMY_exp XMY_calib Z_calib L_calib L_exp L_data Y_calib Y_exp Y_data C_calib C_exp y_calib y_exp I_obs I_exp c_exp c_calib h_calib H_calib"
foreach var in `vars'{
	g l`var'=log(`var')
}

local vars = "L c y XY MY XMY"
foreach var in `vars'{
	g g`var' = ( `var'_exp-`var'_calib )/( 0.5*( `var'_exp+`var'_calib ) )
}

local vars = "L c y XY MY XMY"
foreach var in `vars'{
	g g`var'2 = l`var'_exp-l`var'_calib 
}

g dI = I_exp-I_obs
g gI = ( I_exp-I_obs )/( 0.5*( I_exp+I_obs ) )

************* calibrated d0

sum d0 if convex==1 & Mobility==1 & Ngoods==10 & city_allocation=="largest cities"
g temp1 = string(r(mean),"%12.5fc")
listtex temp1 in 1 using "$save_draft2/d0_mobil_convex.tex", replace
drop temp1

sum d0 if convex==1 & Mobility==0 & Ngoods==10 & city_allocation=="largest cities"
g temp2 = string(r(mean),"%12.5fc")
listtex temp2 in 1 using "$save_draft2/d0_fixed_convex.tex", replace
drop temp2

sum d0 if convex==0 & Mobility==1 & Ngoods==10 & city_allocation=="largest cities"
g temp3 = string(r(mean),"%12.5fc")
listtex temp3 in 1 using "$save_draft2/d0_mobil_nonconvex.tex", replace
drop temp3

sum d0 if convex==0 & Mobility==0 & Ngoods==10 & city_allocation=="largest cities"
g temp4 = string(r(mean),"%12.5fc")
listtex temp4 in 1 using "$save_draft2/d0_fixed_nonconvex.tex", replace
drop temp4

************* intra-regional trade share

destring intra_reg_calib, replace
sum intra_reg_calib if benchmark
g temp1 = string(r(mean)*100,"%12.0fc")
g temp2 = string(r(sd)*100,"%12.0fc")
listtex temp1 in 1 using "$save_draft2/intra_reg.tex", replace
listtex temp2 in 1 using "$save_draft2/intra_reg_sd.tex", replace
drop temp1 temp2

************* Paper table: placement of infrastructure and population

preserve

	keep if benchmark	

	****TABLE 1-a

	*column1
	xi: reg gI lL_calib ly_calib lI_obs i.ICC if cfac_type=="mis_eng", vce(cluster ICC)
	eststo table_1a_1

	*column2
	xi: reg gI lL_calib ly_calib lI_obs i.ICC if cfac_type=="exp_eng", vce(cluster ICC)
	eststo table_1a_2

	*column3
	xi: reg gI lL_calib ly_calib lI_obs i.ICC if cfac_type=="exp_til", vce(cluster ICC)
	eststo table_1a_3

	esttab table_1a* using "$save_draft2/table_1a.tex", b(%9.3f) keep(lL_calib ly_calib lI_obs) coeflabels(lL_calib "Population"  ly_calib "Tradeable Income per Capita" lI_obs "Infrastructure") ///
		mtitle( "Reallocation" "Expansion (GEO)" "Expansion (FOC)") noleg nonotes /// 
		replace not ///
		s(N r2_a, label("Observations" "Adjusted R-squared") fmt(%9.0f %9.2f) ) 


	****TABLE 1-b
	
	*column1
	xi: reg gL lL_calib ly_calib lc_calib lI_obs gI diff_producer i.ICC if cfac_type=="mis_eng", vce(cluster ICC)
	eststo table_1b_1

	*column2
	xi: reg gL lL_calib ly_calib lc_calib lI_obs gI diff_producer i.ICC if cfac_type=="exp_eng", vce(cluster ICC)
	eststo table_1b_2

	*column3
	xi: reg gL lL_calib ly_calib lc_calib lI_obs gI diff_producer i.ICC if cfac_type=="exp_til", vce(cluster ICC)
	eststo table_1b_3

	esttab table_1b* using "$save_draft2/table_1b.tex", b(%9.3f) keep(lL_calib ly_calib lc_calib lI_obs gI diff_producer) ///
	                                                               coeflabels(lL_calib "Population"  ly_calib "Tradeable Income per Capita" lc_calib "Consumption per Capita" ///
																              lI_obs "Infrastructure" gI "Infrastructure Growth" diff_producer "Differentiated Producer") ///
		mtitle( "Reallocation" "Expansion (GEO)" "Expansion (FOC)") noleg nonotes /// 
		replace not ///
		s(N r2_a, label("Observations" "Adjusted R-squared") fmt(%9.0f %9.2f)) 
		
restore

************* Appendix Table A.4

preserve

	keep if convex & Mobility & cfac_type=="exp_eng" & city_allocation=="largest cities"
	
	
	*column1
	xi: reg gI lL_calib ly_calib lI_obs i.ICC if Ngoods==5, vce(cluster ICC)
	eststo table_app_1
	qui	estadd local dep "Investment"
	qui	estadd local ngoods "5"
	
	*column3
	xi: reg gI lL_calib ly_calib lI_obs i.ICC if Ngoods==10, vce(cluster ICC)
	eststo table_app_2
	qui	estadd local dep "Investment"
	qui	estadd local ngoods "10"
	
	*column5
	xi: reg gI lL_calib ly_calib lI_obs i.ICC if Ngoods==15, vce(cluster ICC)
	eststo table_app_3
	qui	estadd local dep "Investment"
	qui	estadd local ngoods "15"
	
	*column2
	xi: reg gL lL_calib ly_calib lI_obs lc_calib gI diff_producer i.ICC if Ngoods==5, vce(cluster ICC)
	eststo table_app_4
	qui	estadd local dep "Pop. Growth"
	qui	estadd local ngoods "5"
	
	*column3
	xi: reg gL lL_calib ly_calib lI_obs lc_calib gI diff_producer i.ICC if Ngoods==10, vce(cluster ICC)
	eststo table_app_5
	qui	estadd local dep "Pop. Growth"
	qui	estadd local ngoods "10"
	
	*column6
	xi: reg gL lL_calib ly_calib lI_obs lc_calib gI diff_producer i.ICC if Ngoods==15, vce(cluster ICC)
	eststo table_app_6
	qui	estadd local dep "Pop. Growth"
	qui	estadd local ngoods "15"
	
	
	esttab table_app* using "$save_draft2/table_a4.tex", b(%9.3f)  keep(lL_calib ly_calib lc_calib lI_obs gI diff_producer) ///
	                                                     coeflabels( lL_calib "Population"  ly_calib "Tradeable Income per Capita" lc_calib "Consumption per Capita" ///
																     lI_obs "Infrastructure" gI "Infrastructure Growth" diff_producer "Differentiated Producer" ) ///
																	 mtitle( "Investment" "Investment" "Investment" "Pop. Growth" "Pop. Growth" "Pop. Growth") ///
		                                                             noleg nonotes /// 
														replace not ///
														s(N r2_a, label("Observations" "Adjusted R-squared" ) fmt(0 %9.2f))
	
restore

************* Paper tables: welfare effects

g y_country = 1
 
* by country
forvalues mobil=0/1{

	preserve

		keep if Mobility==`mobil'

		* one observation counterfactual
		keep if tag_cfac!=0
		
		* number of goods and city allocation
		keep if Ngoods == 10
		keep if city_allocation == "largest cities"
		
		* keep unique ID
		keep welfare_gain ICC Country convex y_country cfac_type
		isid ICC Country convex y_country cfac_type
		
		
		reshape wide welfare_gain,i(ICC Country convex y_country) j(cfac_type) string
		reshape wide welfare_gainexp_eng welfare_gainexp_til welfare_gainmis_eng,i(ICC Country y_country) j(convex)
		
		sort Country
		export excel using "$save_draft/tables.xls", sheet("walfare_mobil`mobil'_raw") sheetmodify first(var)	
		
	restore	

}

* correlation of welfare effects across countries between exp_eng and exp_FOC
preserve

	keep if tag_cfac!=0
	keep if Ngoods == 10
	keep if city_allocation == "largest cities"
	keep if cfac_type=="exp_eng" || cfac_type=="exp_til"
	
	g delta_X=(cfac_type=="exp_eng")
	
	keep welfare_gain ICC Country convex Mobility delta_X
	
	reshape wide welfare_gain,i(ICC Country convex Mobility) j(delta_X)
	
	bysort convex Mobility: spearman welfare_gain0 welfare_gain1
	
restore

* correlation of welfare effects across countries across pairwise comparisons
preserve

	keep if tag_cfac!=0
	keep if Ngoods == 10
	keep if city_allocation == "largest cities"
	
	g delta=(cfac_type=="exp_eng")
	g mis=(cfac_type=="mis_eng")
	
	keep welfare_gain ICC Country convex Mobility delta mis
	
	egen unique=group(convex Mobility delta mis)
	
restore

* correlation of welfare effects across countries between Ngoods=10 and Ngoods=15
preserve

	keep if tag_cfac!=0
	keep if Ngoods == 10 | Ngoods == 15
	keep if city_allocation == "largest cities"

	keep welfare_gain ICC Country convex cfac_type Mobility Ngoods
	
	bysort cfac_type Mobility convex: table Ngoods
	
	reshape wide welfare_gain,i(ICC Country convex cfac_type Mobility) j(Ngoods)
	
	bysort convex cfac_type Mobility: correl welfare_gain10 welfare_gain15
	
restore

* correlation of welfare effects between NUTS and largest cities
preserve

	keep if tag_cfac!=0
	keep if Ngoods == 10 | Ngoods == 15
	keep if city_allocation == "largest cities" || city_allocation == "nuts"

	g city_allocation_X=(city_allocation=="largest cities")
	
	
	keep welfare_gain ICC Country convex cfac_type Mobility Ngoods city_allocation_X
	
	
	reshape wide welfare_gain,i(ICC Country convex cfac_type Mobility Ngoods) j(city_allocation_X)
	
	bysort cfac_type convex Mobility Ngoods: correl welfare_gain0 welfare_gain1
	
restore


forvalues mobil=0/1{

	preserve

		keep if Mobility==`mobil'

		* one observation counterfactual
		keep if tag_cfac!=0
		
		* number of goods and city allocation
		keep if Ngoods == 10
		keep if city_allocation == "nuts"
		
		* keep unique ID
		keep welfare_gain ICC Country convex y_country cfac_type
		isid ICC Country convex y_country cfac_type
		
		reshape wide welfare_gain,i(ICC Country convex y_country) j(cfac_type) string
		reshape wide welfare_gainexp_eng welfare_gainexp_til welfare_gainmis_eng,i(ICC Country y_country) j(convex)
		
		export excel using "$save_draft/tables.xls", sheet("walfare_mobil`mobil'_raw_nuts") sheetmodify first(var)	
		
	restore	

}

* average welfare effects with N=10
	preserve

		* one observation counterfactual
		keep if tag_cfac!=0
		
		* number of goods and city allocation
		keep if Ngoods == 10
		keep if city_allocation == "largest cities"
		
		* keep unique ID
		keep welfare_gain ICC Country convex y_country cfac_type Mobility
		isid ICC Country convex y_country cfac_type Mobility
		g N=1
		collapse (mean) welfare_gain (sum) N, by(convex cfac_type Mobility)

		reshape wide welfare_gain N,i(cfac_type convex) j(Mobility)
		sort convex cfac_type 
		order cfac_type convex welfare_gain0 welfare_gain1 N0 N1
		
		forvalues convex=0/1{
			export excel using "$save_draft/tables.xls" if convex==`convex', sheet("walfare_convex`convex'_raw_average") sheetmodify first(var)	
		}
		
	restore	

* average welfare effects in convex case with different number of goods and largest cities
	preserve

		* one observation counterfactual
		keep if tag_cfac!=0
		
		* convex and city allocation
		keep if convex == 1
		keep if city_allocation == "largest cities"
		
		* keep unique ID
		keep welfare_gain ICC Country Ngoods y_country cfac_type Mobility
		isid ICC Country Ngoods y_country cfac_type Mobility Ngoods
		g N=1
		collapse (mean) welfare_gain (sum) N, by(Ngoods cfac_type Mobility)

		reshape wide welfare_gain N,i(cfac_type Ngoods) j(Mobility)
		order cfac_type Ngoods welfare_gain* N*
		reshape wide welfare_gain* N0 N1,i(cfac_type) j(Ngoods)
		
		order cfac_type welfare_gain05 welfare_gain010 welfare_gain015 ///
		                welfare_gain15 welfare_gain110 welfare_gain115 ///
						N05 N010 N015 ///
						N15 N110 N115 ///
						
		export excel using "$save_draft/tables.xls", sheet("walfare_convex_51015_largest") sheetmodify first(var)	
		
	restore		

* average welfare effects in convex case with different number of goods and NUTS
	preserve

		* one observation counterfactual
		keep if tag_cfac!=0
		
		* convex and city allocation
		keep if convex == 1
		keep if city_allocation == "nuts"
		
		* keep unique ID
		keep welfare_gain ICC Country Ngoods y_country cfac_type Mobility
		isid ICC Country Ngoods y_country cfac_type Mobility Ngoods
		g N=1
		collapse (mean) welfare_gain (sum) N, by(Ngoods cfac_type Mobility)

		reshape wide welfare_gain N,i(cfac_type Ngoods) j(Mobility)
		order cfac_type Ngoods welfare_gain* N*
		reshape wide welfare_gain* N0 N1,i(cfac_type) j(Ngoods)
		
		order cfac_type welfare_gain010 welfare_gain015 ///
		                welfare_gain110 welfare_gain115 ///
						N010 N015 ///
						N110 N115 ///
						
		export excel using "$save_draft/tables.xls", sheet("walfare_convex_51015_nuts") sheetmodify first(var)	
		
	restore	
	
************* Paper table: fundamental law
preserve

	* keep one observation for each country n benchmark only
	keep if benchmark & tag_cfac
	keep if cfac_type=="mis_eng"
	destring flaw_loc_elast, replace
	sum flaw_loc_elast
	table Country, c(mean flaw_loc_elast)
	
	keep Country flaw_loc_elast
	
	export excel using "$save_draft/tables.xls", sheet("flaw") sheetmodify first(var)	

restore	
	
************* Paper figure: welfare effects

* plot labor fit
			preserve
			
				* keep benchmark
				keep if Ngoods==10 & city_allocation == "largest cities"
				keep if cfac_type=="mis_eng"
				keep if Mobility ==1 
				keep if convex==1
				
				* labor: model vs. data
				
				local x = "L_data"
				local y = "L_calib"
				reg `y' `x', r
				local b = round(_b[`x']*1000)/1000
				local se = round( _se[`x']*1000)/1000
				
				twoway  (scatter `y' `x', msize(vsmall) mcolor(blue) ) ( lfit `x' `x', lwidth(medium) color(black) ), ///
						xtitle("Labor Shares (Data)") ytitle("Labor Shares (Model)") graphregion( color(white) ) ///
						legend( order( 1 2 ) lab( 1 "Mobile Labor" ) lab( 2 "45 Degree" ) rows(2) cols(1) ) ///
						note( "Linear regression slope (robust SE): Mobile Labor: `b' (`se')"  )
				graph export "$save_draft/fit_L_all_ngoods10_Mobil1Convex1.pdf", replace
				
			restore

* plot income fit
			preserve						
				
				* keep benchmark
				keep if Ngoods==10 & city_allocation == "largest cities"
				keep if cfac_type=="mis_eng"
				*keep if Mobility ==1 
				keep if convex==1
				
				sort Mobility ICC
				by Mobility ICC: egen totY_calib=total(Y_calib)
				by Mobility ICC: egen totY_data=total(Y_data)
				
				replace Y_calib = Y_calib/totY_calib
				replace Y_data = Y_data/totY_data
				
				* income: model vs. data
				local x = "Y_data"
				local y = "Y_calib"
				
				reg `y' `x' if Mobility==1,  r
				local b_mobile = round(_b[`x']*1000)/1000
				local se_mobile = round( _se[`x']*1000)/1000
				
				reg `y' `x' if Mobility==0,  r
				local b_fixed = round(_b[`x']*1000)/1000
				local se_fixed = round( _se[`x']*1000)/1000
				
				twoway  (scatter `y' `x' if Mobility==1, msize(vsmall) mcolor(blue) ) (lfit `x' `x' if Mobility==1, lwidth(medium) color(black) ) ///
				        (scatter `y' `x' if Mobility==0, msize(vsmall) mcolor(red) ) (lfit `x' `x' if Mobility==0, lwidth(medium) color(black) ), ///
						xtitle("Income Shares (Data)") ytitle("Income Shares (Model)") graphregion( color(white) ) ///
						legend( order( 1 3 2 ) lab( 1 "Mobile Labor" ) lab( 3 "Fixed Labor" ) lab( 2 "45 Degree" ) rows(2) cols(2) ) ///
						note( "Linear regression slope (robust SE): Mobile Labor: `b_mobile' (`se_mobile'); Fixed Labor: `b_fixed' (`se_fixed')"  )
				graph export "$save_draft/fit_Y_all_ngoods10_Convex1.pdf", replace
		
				
			restore
			
* plot fundamentals versus data
			preserve
			
				* keep benchmark
				keep if Ngoods==10 & city_allocation == "largest cities"
				keep if cfac_type=="mis_eng"
				keep if convex==1
				
				* log fundamentals
			    g lZ = log(Z_calib)
			    g lh = log(h_calib)
				g lH = log(H_calib)
			
				forvalues mobil=0/1{
					
					* productivity versus income
					local y = "lZ"
					local x = "lY_calib"
					reg `y' `x' if Mobility==`mobil',  r
					local b_zy`mobil' = round(_b[`x']*1000)/1000
					local se_zy`mobil' = round( _se[`x']*1000)/1000
					
					* h versus income
					local y = "lH"
					local x = "lY_calib"
					reg `y' `x' if Mobility==`mobil',  r
					local b_hy`mobil' = round(_b[`x']*1000)/1000
					local se_hy`mobil' = round( _se[`x']*1000)/1000
					
					* productivity versus population
					local y = "lZ"
					local x = "lL_calib"
					reg `y' `x' if Mobility==`mobil',  r
					local b_zl`mobil' = round(_b[`x']*1000)/1000
					local se_zl`mobil' = round( _se[`x']*1000)/1000
					
					* h versus income
					local y = "lH"
					local x = "lL_calib"
					reg `y' `x' if Mobility==`mobil',  r
					local b_hl`mobil' = round(_b[`x']*1000)/1000
					local se_hl`mobil' = round( _se[`x']*1000)/1000
					
				}
					
					
			    local y1 = "lZ"
				local y2 = "lH_calib"
				local x = "lY_calib"
				twoway  (scatter `y1' `x' if Mobility==1, msize(vsmall) mcolor(blue) ) (lfit `y1' `x' if Mobility==1, lwidth(medium) color(blue) ) ///
						(scatter `y2' `x' if Mobility==1, msize(vsmall) mcolor(red) ) (lfit `y2' `x' if Mobility==1, lwidth(medium) color(red) ), ///
						xtitle("Log Income Shares (Data)") ytitle("Fundamentals") graphregion( color(white) ) ///
						legend( order( 1 3 ) lab( 1 "Log Productivity (z)" ) lab( 3 "Log Non-Traded Endowment" ) rows(1) cols(2) ) ///
						note( "Linear regression slope (robust SE): Productivity: `b_zy1' (`se_zy1'); Endowment: `b_hy1' (`se_hy1')"  )
						
				graph export "$save_draft/ZH_Y_all_mobility1_ngoods10.pdf", replace
			
				local y1 = "lZ"
				local y2 = "lH_calib"
				local x = "lL_calib"
				twoway  (scatter `y1' `x' if Mobility==1, msize(vsmall) mcolor(blue) ) (lfit `y1' `x' if Mobility==1, lwidth(medium) color(blue) ) ///
						(scatter `y2' `x' if Mobility==1, msize(vsmall) mcolor(red) ) (lfit `y2' `x' if Mobility==1, lwidth(medium) color(red) ), ///
						xtitle("Log Population Shares (Data)") ytitle("Fundamentals") graphregion( color(white) ) ///
						legend( order( 1 3 ) lab( 1 "Log Productivity (z)" ) lab( 3 "Log Non-Traded Endowment" ) rows(1) cols(2) ) ///
						note( "Linear regression slope (robust SE): Productivity: `b_zl1' (`se_zl1'); Endowment: `b_hl1' (`se_hl1')"  )
				
				graph export "$save_draft/ZH_L_all_mobility1_ngoods10.pdf", replace				
				
				
			restore			
			
***** table A6
import excel using "$raw_calibration\all_cfactuals.xls",sheet("all_results") firstrow clear

********** keep cfac **********
keep if cfac_type=="exp_eng" || cfac_type=="exp_til"  || cfac_type=="mis_eng" 

* define mobility as 0/1
g temp=Mobility
drop Mobility
g Mobility=1 if temp=="on"
replace Mobility=0 if temp=="off"

* define parameter space
g low_beta = (beta==0.13)
g convex = (beta>gamma)

* express welfare gain in % points
replace welfare_gain = ( welfare_gain-1 )*100	

* unique identifier of type of counterfactual
egen cfactual_type_id = group( convex cfac_type Ngoods city_allocation Mobility cong )
isid cfactual_type_id Country location
bys cfactual_type_id Country location: g rep=_N
drop cfactual_type_id

* check number of counterfactuals for each case
egen tag_cfac = tag( cfactual_id ICC )

* keep benchmark
keep if Ngoods==10
keep if convex==1
keep if city_allocation=="largest cities"

* keep 1 counterfactual 
keep if tag_cfac==1

* verify unique ID
isid Country cong Mobil cfac_type

* table with average welfare gains
preserve
	
	collapse (mean) welfare_gain, by(cong Mobil cfac_type)
	reshape wide welfare, i(cfac_type Mobil) j(cong) string
	reshape wide welfare_gainoff welfare_gainon,i(cfac_type) j(Mobil)

	export excel using "$save_draft/tables.xls", sheet("welfare_cong_on_off_raw") sheetmodify first(var)	
	
restore
