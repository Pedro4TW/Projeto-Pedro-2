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

* Benchmark only

	import excel using "$load_gravity\gravity_Spain_nuts_Ngoods15.xls",sheet("trade_Spain") firstrow clear

	sort origin distance

	foreach var of varlist x_real x_calib{

		egen `var'_tot=total(`var')
		replace `var'=`var'/`var'_tot
		
		bys origin: egen `var'_tot_exports = total(`var')
		bys destination: egen `var'_tot_imports = total(`var')
		
		g l`var'=log(`var')
		g l`var'_tot_exports=log(`var'_tot_exports)
		g l`var'_tot_imports=log(`var'_tot_imports)
		
		g l`var'_import_share =log(`var'/`var'_tot_imports)
		g l`var'_export_share =log(`var'/`var'_tot_exports)
		
		sum l`var'_import_share    
		g l`var'_import_share_d = l`var'_import_share-`r(mean)' // demeaned
		
		sum l`var'_export_share
		g l`var'_export_share_d = l`var'_export_share-`r(mean)' // demeaned
	}	

	g ldist = log(distance)
	g ldist1 = log(distance+1)
	
	* dummy for own share
	g own = (origin==destination)
	
	* gravity equations with FE
	xi: reg lx_calib ldist i.origin i.destination, r	
	xi: reg lx_real ldist i.origin i.destination, r
	
	* import share on FE
	xi: reg lx_calib_import_share ldist i.origin, r
	xi: reg lx_real_import_share ldist i.origin, r
	
	* plot actual versus calibrated
	correl x_calib x_real  // correlation  = 0.79
	local correl = round( r(rho)*1000 )/1000
	twoway (lfit lx_calib lx_real, lwidth(medium) color(blue) ) ///
		   (lfit lx_real lx_real, lwidth(thin) lp(dash) color(black)) ///
		   (scatter lx_calib lx_real if origin!=destination, msymbol(o) color(blue) )  ///
		   (scatter lx_calib lx_real if origin==destination, msymbol(X) color(red) ), ///
		   xtitle(Log of Bilateral Trade (Data)) ytitle(Log of Bilateral Trade (Model)) graphregion( color(white) ) ///		   
		   legend( order( 3 4 2 ) lab(2 "45 degree") lab(3 "Between-NUTS Flow") lab(4 "Within-NUTS Flow") rows(1) cols(3) ) 

	local x="$save_draft"+"actual_predicted_trade_Spain_nuts_Ngoods15_paper.eps"
	graph export "`x'", as(eps) replace	

	
	* plot residualized import share versus distance, after taking out fixed effects
	xi: reg lx_calib_import_share ldist i.origin, r
	local b_calib = round(_b[ldist]*1000)/1000
	local se_calib = round( _se[ldist]*1000)/1000
					
	reg lx_calib_import_share i.origin, r
	predict lx_calib_import_share_res, residuals
		
	reg ldist i.origin, r
	predict ldist_res, residuals
		
	xi: reg lx_real_import_share ldist i.origin, r
	local b_real = round(_b[ldist]*1000)/1000
	local se_real = round( _se[ldist]*1000)/1000
					
	reg lx_real_import_share i.origin, r
	predict lx_real_import_share_res, residuals		

	twoway (lfit lx_calib_import_share_res ldist_res if lx_calib_import_share_res<2, lwidth(medium) color(blue) ysc(r(-2 2)) ) (scatter lx_calib_import_share_res ldist_res if lx_calib_import_share_res<2, msize(vsmall) color(blue) ysc(r(-2 2))) ///
		   (lfit lx_real_import_share_res ldist_res if lx_calib_import_share_res<2, lwidth(medium) color(red) ysc(r(-2 2))) (scatter lx_real_import_share_res ldist_res if lx_calib_import_share_res<2, msize(vsmall) color(red) ysc(r(-2 2))), ///
			xtitle(Log of Distance (KM)) ytitle(Log of Import Share) graphregion( color(white) ) ///
			legend( order( 1 3 ) lab(1 "Model") lab(3 "Data") rows(1))
			
    local x="$save_draft"+"gravity_Spain_nuts_Ngoods15_residualized_paper.eps"
	graph export "`x'", as(eps) replace	 			
