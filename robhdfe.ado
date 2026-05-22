*! version 1.2.0 20260522 David Veenman

/*
20260522: 1.2.0     Added option for Driscoll-Kraay standard errors for spatial dependence
20260429: 1.1.1     Fixed minor bug when using keepsin option (error because of missing Ndrop scalar)
20260416: 1.1.0     Added python option in IRWLS for faster execution with pyfixest
                    Fixed minor collinearity effects and aligned convergence criterian for julia option
20260408: 1.0.0     First version

Dependencies (other Stata packages):
   moremata
   reghdfe
   hdfe
   julia [for julia option only]
   reghdfejl [for julia option only]
*/

program define robhdfe, eclass sortpreserve
	version 15
	syntax [anything] [in] [if], absorb(varlist) eff(real) [cluster(varlist) dkraay(string) tol(real 0) weightvar(str) keepsin julia python]

	capture findfile mf_mm_aqreg.hlp
	if _rc {
		di as error "Program requires the {bf:moremata} package: type {stata ssc install moremata, replace}"
		error 499
	}

	capture findfile reghdfe.ado 
	if _rc {
		di as error "Program requires the {bf:reghdfe} package: type {stata ssc install reghdfe, replace}"
		error 499
	}
	
	capture findfile hdfe.ado 
	if _rc {
		di as error "Program requires the {bf:hdfe} package: type {stata ssc install hdfe, replace}"
		error 499
	}

	if ("`julia'" != "") {
		capture findfile reghdfejl.ado 
		if _rc {
			di as error "Program requires the {bf:reghdfejl} package: type {stata ssc install reghdfejl, replace}"
			error 499
		}
		capture findfile jl.ado 
		if _rc {
			di as error "Program requires the {bf:julia} package: type {stata ssc install julia, replace}"
			error 499
		}		
	}
	
	if ("`python'" != "") {
		capture python: import pyfixest
		if _rc {
			di as error "Program requires Python installation and the {bf:pyfixest} Python package: pip install pyfixest"
			error 499
		}
	}

	if ("`julia'" != "" & "`python'" != "") {
		di as err "ERROR: Cannot specify both julia and python options simultaneously"
		exit
	}

	if ("`keepsin'" != "" & "`python'" != "") {
		di as err "ERROR: Cannot combine python option with keepsin"
		exit
	}
	
	marksample touse
		
	tokenize `anything'
	local subcmd `"`1'"'

	local cmdlist "m"
	if !`:list subcmd in cmdlist' {
		di as err `"Invalid subcommand: `subcmd'"'
		exit 
	}
		
	macro shift 1
	local depv `"`1'"'
	unab depv: `depv'
	local varlist `"`*'"'

	// Ensure dv is not a factor variable:
	_fv_check_depvar `depv'
	macro shift 1
	local indepv "`*'"

	// Ensure dv is not an indicator variable:
	qui capture assert inlist(`depv', 0, 1)
	if _rc==0 {
        di as err "ERROR: Dependent variable should not be an indicator (0/1) variable"
        exit 		
	}	
	
	// Ensure iv list does not contain a factor variable:
    fvexpand `indepv'
    if "`r(fvops)'" == "true" {
        di as err "ERROR: Independent variable list may not contain factor variables"
        exit 
    }
	else {
		local indepv `r(varlist)'
	}
	
	// Mark out missing observations:
	markout `touse' `depv' `indepv'

	// Check number of independent variables:
	local varcount = 0
	foreach v of local indepv {
		local `varcount++'
	}
	scalar k0 = `varcount'
	
	// Check absorb variables:
	local nabs: word count `absorb'
	local j = 0
	foreach var of local absorb {
		local `j++'
		local absvar`j' = "`var'"
		markout `touse' `var'
	}
	
	local abs_unique
	foreach var of local absorb {
		if !`: list var in abs_unique' {
			local abs_unique `abs_unique' `var'
		}
	}
	local nabs_un: word count `abs_unique'
	if (`nabs' != `nabs_un') {
	    di as err "ERROR: Option absorb() must contain different variables"
		exit				
	}

	// Ensure absorb dimensions are not nested:
	foreach ai of local absorb {
		local others : list absorb - ai
		foreach aj of local others {
			capture bysort `ai': assert `aj' == `aj'[1] if !missing(`ai', `aj')
			if (_rc == 0) {
				di as err "ERROR: `ai' is nested within `aj'"
				exit
			}
		}
	}
	
	// Process information in dkraay() option:
	local dkn: word count `dkraay'
	if ("`dkraay'" != "" & `dkn' != 2) {
		di as err "ERROR: Option dkraay() incorrectly specified"
		exit 198
	}
	if ("`dkraay'"!="") {
		if ("`cluster'" != "") {
			di as err "ERROR: Options cluster() and dkraay() may not be combined"
			exit 
		}
		// Time dimension:
		local dk_time: word 1 of `dkraay'
		capture confirm variable `dk_time'
		if _rc {
			di as err "ERROR: Time variable `dk_time' in option dkraay() not found"
			exit
		}
		capture confirm numeric variable `dk_time'
		if _rc {
			di as err "ERROR: Time variable `dk_time' in option dkraay() not numeric"
			exit
		}
		markout `touse' `dk_time'
		tempvar dk_time_id
		qui egen double `dk_time_id' = group(`dk_time') if `touse'
		// Number of lags:
		local dk_lags: word 2 of `dkraay'
		capture confirm integer number `dk_lags'
		if _rc {
			di as err "ERROR: Lag length in option dkraay() must be integer"
			exit
		}
		if (`dk_lags' < 1) {
			di as err "ERROR: Lag length in option dkraay() must be positive"
			exit			
		}
		scalar dk_lags = `dk_lags'
		qui sum `dk_time_id'
		if r(max) < dk_lags {
			di as err "ERROR: Lag length in option dkraay() must be smaller than number of time periods"
			exit						
		}
	}
	
	// Process information in cluster() option:
	local nc: word count `cluster'
	
	if (`nc'>2){
	    di as err "ERROR: Maximum number of dimensions to cluster on is two"
		exit
	}
	local clusterdim1: word 1 of `cluster'
	local clusterdim2: word 2 of `cluster'
	
	// Check nesting of FE in clusters and create indicators for dof adjustments:
	local allnest=1
	local j = 0
	local all1 = 1
	if ("`cluster'" == "" & "`dkraay'" == "") {
		local nocluster = 1
		foreach abs of local absorb {
			local `j++'
			local nest`j' = 1		
			local nest`j'dof = 1
		}
		local nest1dof = 0
	}
	else {
		if ("`cluster'" != "") {
			foreach abs of local absorb {
				local `j++'
				local nest`j' = 1
				foreach cl of local cluster {
					capture bysort `abs': assert `cl'==`cl'[1] if !missing(`abs', `cl')
					if (_rc == 0) {
						local nest`j' = 0
						continue, break
					}
				}			
				local nest`j'dof = `nest`j''
				local all1 = `all1' * `nest`j'dof'
				local allnest = `allnest' * `nest`j''
			}
		}
		if ("`dkraay'" != "") { // Currently redundant as small-sample correction uses full K only as in fixest
			foreach abs of local absorb {
				local `j++'
				local nest`j' = 1
				capture bysort `abs': assert `dk_time'==`dk_time'[1] if !missing(`abs', `dk_time')
				if (_rc == 0) {
					local nest`j' = 0
				}
				local nest`j' = 1 // Override for consistency with fixest
				local nest`j'dof = `nest`j''
				local all1 = `all1' * `nest`j'dof'
				local allnest = `allnest' * `nest`j''
			}
		}
		if (`all1' == 1) {
			local nest1dof = 0
		}		
	}
		
	// Set tolerance:
	if (`tol' != 0){
		local tolerance = `tol'
	}
	else {
		local tolerance = 1e-10
	}	

	// Check efficiency:	
	if (`eff' < 63.7 | `eff' > 99.9) {
		di as err "ERROR: Normal efficiency must be between 63.7 and 99.9"
		exit
	}
		
	// Create temporary variables: 
	tempvar clus1
	if ("`cluster'" != "") {
		qui egen double `clus1' = group(`clusterdim1') if `touse'
		if (`nc'>1){
			tempvar clus2 clus12
			qui egen double `clus2' = group(`clusterdim2') if `touse'
			qui egen double `clus12' = group(`clusterdim1' `clusterdim2') if `touse'
		}
	}
	else {
		 qui gen double `clus1' = 1 if `touse'
	}
			
	// Checking collinearity (including fixed effects):
	qui hdfe `indepv' if `touse', absorb(`absorb') gen(_stub_) keepsin
	foreach v of local indepv {
		tempvar `v'_temp
		qui ren `v' `v'_temp
		qui ren _stub_`v' `v'
	}
	_rmcoll `indepv'
	local k_omitted = r(k_omitted)
	local indepv_temp = r(varlist)
	foreach v of local indepv {
		drop `v'
		qui ren `v'_temp `v'
	}
	local indepv `indepv_temp'

	// Prepare list for phi-residualized X matrix for VCE calculation in the presence of collinear variables:
	local indepv0 ""
	foreach v of local indepv {
		if strpos("`v'", "o.") {
			local indepv0 "`indepv0'"
		}
		else {
			local indepv0 "`indepv0' `v'"			
		}
	}
	
	// Adjust variable lists for Julia and Python options in the presence of collinear variables:
	if ("`julia'" != "") {
		local indepv_jl ""
		local j = 1
		foreach v of local indepv {
			if strpos("`v'", "o.") {
				tempvar _collin_temp_`j'
				gen `_collin_temp_`j'' = 1
				local indepv_jl "`indepv_jl' `_collin_temp_`j''"
			}
			else {
				local indepv_jl "`indepv_jl' `v'"			
			}
		}
	}		
	if ("`python'" != "") {
		local indepv_py ""
		foreach v of local indepv {
			if strpos("`v'", "o.") == 0 {
				local indepv_py "`indepv_py' `v'"
			}
		}
	}	
	
	di ""
    /////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	di as text "STEP 1: Estimating initial MM-QR and obtaining scale estimate"
    /////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////

	// Location stage MM-QR (Machado and Santos Silva 2019):
	tempvar e Ipos r_raw denom u resid_tau

	if ("`julia'" != "") {
		capture drop _reghdfejl_*
	}
	else {
		capture sum _reghdfe_resid
		if (_rc == 0) {
			ren _reghdfe_resid _temp_reghdfe_resid
		}
	}
	
	qui sum `touse' if `touse'>0
    local N0=r(N)
	if ("`julia'" != "") {
		qui capture reghdfejl `depv' `indepv' if `touse', absorb(`absorb') notable nofootnote noheader resid `keepsin'
		qui predict _reghdfejl_res, res
		markout `touse' _reghdfejl_res
		qui ren _reghdfejl_res `e'
		drop _reghdfejl_*
	}
	else {
		qui capture reghdfe `depv' `indepv' if `touse', absorb(`absorb') dof(none) notable nofootnote noheader resid `keepsin'
		markout `touse' _reghdfe_resid
		qui ren _reghdfe_resid `e'
	}
	
	qui sum `touse' if `touse'>0
    local N = r(N)
	
	if ("`keepsin'" != "") {
		qui replace `e'=0 if abs(`e')<1e-10
		local Ndrop = 0	
	}
	else {
		if (`N' < `N0') {
			local Ndrop = `N0'-`N'
			if (`Ndrop' ==1 ) {
				di "note: dropped 1 singleton observation."
			}
			else {
				di "note: dropped `Ndrop' singleton observations."			
			}
		}
		else {
			local Ndrop = 0			
		}
	} 
	
	// Scale stage:
	qui gen `Ipos' = (`e' > =0) if `touse'
	qui sum `Ipos' if `touse', meanonly
	scalar Ibar = r(mean)
	qui gen double `r_raw' = 2 * `e' * (`Ipos' - Ibar) if `touse'
	
	if "`julia'"!="" {
		qui capture reghdfejl `r_raw' `indepv' if `touse', absorb(`absorb') notable nofootnote noheader resid `keepsin'  
		qui predict double `denom' if `touse', xbd
		drop _reghdfejl_*
	}
	else {
		qui capture reghdfe `r_raw' `indepv' if `touse', absorb(`absorb') dof(none) notable nofootnote noheader resid `keepsin'  
		qui predict double `denom' if `touse', xbd
		drop _reghdfe_resid
	}
	
	// Standardized residuals and create qhat:
	qui gen double `u' = `e' / `denom' if `touse'
	qui sum `u' if `touse', d // Note: xtqreg and mmqreg use qreg on constant; I use percentile approach instead for consistency with robreg and Mata function mm_aqreg()
	scalar qhat = r(p50)
	
	// Residuals:
	qui gen double `resid_tau' = `e' - qhat*`denom' if `touse'
	
	// Get relevant information from the data before creating scale estimate:
	local j = 0
	foreach abs of local absorb {
		local `j++'
		tempvar absvar`j'id 
		qui egen double `absvar`j'id' = group(`abs') if `touse'
	}
	local j = 0
	local ntotal = 0
	local ntotal_est = 0
	foreach abs of local absorb {
		local `j++'
		qui sum `absvar`j'id'
		local n`j' = r(max)
		local ntotal = `ntotal'+r(max)
		local n`j'_red = (1-`nest`j'')*`n`j'' + `nest`j'dof'
		local n`j'_est = `n`j'' - (1-`nest`j'')*`n`j'' - `nest`j'dof'
		local ntotal_est = `ntotal_est' + `n`j'_est'
	}	
	local Kinit: word count `indepv' 
	local Kinit = `Kinit' - `k_omitted' 
	if (`nabs'>1) {
		scalar df_initial = `N' - `ntotal' - (`Kinit' - 1) 
	}
	else {
		scalar df_initial = `N' - `ntotal' - `Kinit' 
	}
	local K = `Kinit' + 1 + `ntotal_est'
	if ("`dkraay'" != "") {
		local K_dk_full = `Kinit' + `ntotal' - `nabs' + 1
	}
	
	// Get scale estimate and initial weights:
	tempvar w 
	scalar eff = `eff'
	mata: _scale_initial()
	
    /////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	di as text "STEP 2: Iterating IRWLS"
    /////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	tempvar _resid_temp phi
	qui gen double `phi' = .
    local diff = 100
	local maxiter = c(maxiter)
	
	// Julia setup - transfer variables and build formula:
	if ("`julia'" != "") {
		reghdfejl_load
		reghdfejl_parse_absorb `absorb' if `touse'
		local jl_feterms `"`r(feterms)'"'
		local jl_absorbvars `r(absorbvars)'
		local jl_putvars `depv' `indepv_jl' `w' `jl_absorbvars'
		local jl_putvars: list uniq jl_putvars
		unab jl_putvars: `jl_putvars'
		jl PutVarsToDF `jl_putvars' if `touse', nomissing doubleonly nolabel
		local jl_indepfmla: subinstr local indepv_jl " " " + ", all
		_jl: jl_f = @formula(`depv' ~ `jl_indepfmla' `jl_feterms')
	}
	
	// Python setup - transfer variables to Python and build pyfixest formula:
	if ("`python'" != "") {
		local py_indepfmla: subinstr local indepv_py " " " + ", all
		local py_feterms:   subinstr local absorb " " " + ", all
		local py_formula "`depv' ~ `py_indepfmla' | `py_feterms'"
		python: from sfi import Data, Scalar, Macro
		python: import numpy as np
		python: import pandas as pd
		python: import pyfixest as pf
		python: _touse = Macro.getLocal("touse")
		python: _depv = Macro.getLocal("depv")
		python: _indepv = Macro.getLocal("indepv_py").split()
		python: _absorb = Macro.getLocal("absorb").split()
		python: _w = Macro.getLocal("w")
		python: _phi = Macro.getLocal("phi")
		python: _formula = Macro.getLocal("py_formula")
		python: _all_vars = [_depv] + _indepv + [_w] + _absorb
		python: _seen = set()
		python: _all_vars = [x for x in _all_vars if x not in _seen and not _seen.add(x)]
		python: _py_mask = np.array(Data.get(_touse)).astype(bool)
		python: _py_obs_stata = np.where(_py_mask)[0].tolist()
		python: _py_df = pd.DataFrame({var: np.array(Data.get(var))[_py_mask] for var in _all_vars})
		python: _py_scale = Scalar.getValue("scale")
		python: _py_krob = Scalar.getValue("krob")
		python: _py_w = _py_df[_w].values.copy()
		python: _py_phi = np.ones(_py_mask.sum())
		python: _py_b0 = None
	}
	
    forvalues i=1(1)`maxiter'{
        if (`diff' > `tolerance') {
			qui capture drop `_resid_temp'
			if ("`julia'" != "") {
				_jl: df[!, :`w'] = vec(st_data("`w'", "`touse'"))
				_jl: m = reg(df, jl_f, weights = :`w', drop_singletons=false, save=:residuals, maxiter=16000, tol=1e-8)
				_jl: res = residuals(m)
				_jl: replace!(res, missing=>NaN)
				qui jl GetVarsFromMat `_resid_temp' if `touse', source(res)
				if `i'>1 {
					_jl: jl_diff = maximum(abs.(coef(m) .- jl_b0) ./ (abs.(jl_b0) .+ 1.0))
					_jl: st_numscalar("jl_diff", jl_diff)
					local diff = jl_diff
				}
				_jl: jl_b0 = copy(coef(m))				
			}
			if ("`python'" != "") {
				python: _py_df[_w] = _py_w
				python: _fit = pf.feols(_formula, data=_py_df, weights=_w, vcov="iid")
				python: _coefs = np.array(_fit.coef())
				python: _resid = np.array(_fit.resid())
				python: _py_diff = float(np.max(np.abs(_coefs - _py_b0) / (np.abs(_py_b0) + 1.0))) if _py_b0 is not None else 100.0
				python: Scalar.setValue("py_diff", _py_diff)
				python: _py_b0 = _coefs.copy()
				if `i'>1 {
					local diff = py_diff
				}
			}
			if ("`julia'" == "" & "`python'" == "") {
				if (_caller() < 19) {
					qui capture reghdfe `depv' `indepv' [aw = `w'] if `touse', absorb(`absorb') dof(none) notable nofootnote noheader resid keepsin
					qui ren _reghdfe_resid `_resid_temp' 
				}
				else {
					qui capture areg `depv' `indepv' [aw = `w'] if `touse', absorb(`absorb') noabs 
					qui predict `_resid_temp', res 
				}
				matrix b = e(b)            
				if (`i' > 1) {
					local diff = mreldif(b0, b)
				}
				matrix b0 = b
			}
			if (`i' == `maxiter' & `diff' > `tolerance'){
				di as err "ERROR: Convergence not achieved"
				exit
			}
			if (`diff' > `tolerance'){
				if ("`python'" != "") {
					python: _z = _resid / _py_scale
					python: _py_w = np.where(np.abs(_z) <= _py_krob, 1.0, _py_krob / np.abs(_z))
					python: _py_phi = (np.abs(_z) <= _py_krob).astype(float)
					python: print(".", end="", flush=True)
				}
				else {
					qui gen double _z_temp = `_resid_temp'/scale
					mata: _update_weights()
					drop _z_temp
				}
			}
        }
    }
	
	if "`julia'"!="" {
		_jl: jl_b0_mat = reshape(jl_b0, 1, length(jl_b0))
		jl GetMatFromMat b0, source(jl_b0_mat)
		matrix b=b0
	}

	if "`python'"!="" {
		qui gen double `_resid_temp' = .
		python: Data.store(Macro.getLocal("_resid_temp"), _py_obs_stata, _resid.tolist())
		python: Data.store(_phi, _py_obs_stata, _py_phi.tolist())
		python: Data.store(_w, _py_obs_stata, _py_w.tolist())
		python: from sfi import Matrix
		python: _indepv_full = Macro.getLocal("indepv").split()
		python: _py_b0_iter = iter(_py_b0.tolist())
		python: _py_b0_full = [0.0 if v.startswith("o.") else next(_py_b0_iter) for v in _indepv_full]
		python: Matrix.store("b0", [_py_b0_full])
		matrix b=b0
	}
	
	if ("`weightvar'" != "") {
		capture drop `weightvar'
		if (_rc == 0) {
			local replaceweightvar "yes"
		}
		gen double `weightvar' = `w' 
	}
	
	mata: ""
    /////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	di as text "STEP 3: Computing standard errors"
    /////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	qui replace `phi' = 1e-20 if `phi' == 0 // Ensure that residualized values are also created for phi=0 cases 
	if "`julia'"!="" {
		qui partialhdfejl `indepv0' if `touse' [aw = `phi'], absorb(`absorb') prefix(_stub_)
	}
	else {
		qui hdfe `indepv0' if `touse' [aw = `phi'], absorb(`absorb') gen(_stub_) keepsin
	}
	local indepvr ""
	foreach v of local indepv {
		if strpos("`v'", "o.") {
			local indepvr "`indepvr' `v'"
			local vclean = subinstr("`v'", "o.", "", 1)
		}
		else {
			tempvar _tilde_`v'
			qui gen `_tilde_`v'' = _stub_`v'
			drop _stub_`v'
			local indepvr "`indepvr' `_tilde_`v''"
		}
	}

	// Calculation of Pseudo R2:
	scalar maxiter = `maxiter'
	scalar tol = `tolerance'
	mata: _huber_location()
	mata: _pseudo_r2()
	
	// VCE:
	if ("`dkraay'" != "") {
		sort `dk_time_id'
		local tvar "`dk_time_id'"	
		mata: _vce_dkraay()
		matrix beta = b0[.,1..k0]
		matrix Vc = Vdk
		local e_df_r = mata_ntime-1
	}
	else {
		sort `clus1' 
		local cvar "`clus1'"	
		mata: _vce_cluster()    
		local nclusterdim = mata_nclusters
		if ("`cluster'" == "") {
			local e_df_r = df_initial
		}
		else {
			local e_df_r = mata_nclusters-1
		}
		
		matrix beta = b0[.,1..k0]
		matrix Vc = Vclust
				
		if (`nc' > 1) {
			// Second clustering dimension:
			matrix V1 = Vclust
			sort `clus2' 
			local cvar "`clus2'"	
			mata: _vce_cluster()
			local nclusterdim1 = `nclusterdim'
			local nclusterdim2 = mata_nclusters
			local e_df_r2 = mata_nclusters - 1
			if (`nclusterdim2' < `nclusterdim') {
				local nclusterdim = `nclusterdim2'
				local e_df_r = `e_df_r2'
			}
			matrix V2 = Vclust
			// Intersection of clustering dimensions:
			sort `clus12' 
			local cvar "`clus12'"	
			mata: _vce_cluster()    
			matrix V12 = Vclust
			matrix Vc = V1 + V2 - V12
			matrix drop V1 V2 V12
		}		
	}
				
	if ("`cluster'" == "") {
		if ("`dkraay'" == "") {
			local factor = (`N'/`e_df_r')
		}
		else {
			local factor = ((`N' - 1)/(`N' - `K_dk_full'))*(mata_ntime / (mata_ntime-1))
		}
	}
	else{
		local factor = (`nclusterdim' / (`nclusterdim' - 1))*((`N' - 1)/(`N' - `K'))
	}
	matrix Vc = `factor' * Vc
    
	ereturn clear
	tempname b V

	matrix colnames Vc = `indepv'
	matrix rownames Vc = `indepv'
    matrix colnames beta = `indepv'	
	matrix rownames beta = `depv'
	
	matrix `b' = beta
	matrix `V' = Vc
	
	ereturn post `b' `V'
	ereturn scalar N = `N'
	if (`Ndrop' > 0) {
		ereturn scalar N_singletons = `Ndrop'		
	}
	if "`cluster'"!="" {
		if (`nc' == 1) {
			ereturn scalar N_clust=`nclusterdim'
		}
		else {
			ereturn scalar N_clust1 = `nclusterdim1'
			ereturn scalar N_clust2 = `nclusterdim2'
		}
	}
	ereturn scalar df_r = `e_df_r'
	ereturn scalar r2_p = r2_p
	ereturn scalar scale = scale 
	ereturn scalar ssc = `factor'
	if ("`dkraay'" != "") {
		ereturn scalar df_k = `K_dk_full'
		ereturn local vcetype "Driscoll-Kraay"
		ereturn scalar dk_lags = `dk_lags'
	}
	else {
		ereturn scalar df_k = `K'
		if ("`cluster'" == "") {
			ereturn local vcetype "robust"
		}
		else {
			ereturn local vcetype "cluster-robust"
		}
	}
	
    ereturn local depvar "`depv'"
    ereturn local indepvars "`indepv'"
    ereturn local cmd "robhdfe"
    ereturn local subcmd "`subcmd'"
    ereturn local clustvar "`cluster'"
	
	di ""
	di in green "Huber M-estimation with `eff'% normal efficiency and fixed effects"
	if ("`cluster'" == "" & "`dkraay'" == "") {
		di in green "Heteroskedasticity-robust standard errors" 
	}
	if (`nc' == 1){
		di in green "Standard errors adjusted for clustering by `clusterdim1'" 
	}
	if (`nc' == 2){
		di in green "Standard errors adjusted for clustering by `clusterdim1' and `clusterdim2'"
	}
	if ("`dkraay'" != "") {
		di in green "Driscoll-Kraay standard errors (`dk_lags' lags)"		
	}
	di ""
	di _column(51) in green "Number of obs = " %12.0fc in yellow e(N)
	di _column(51) in green "Pseudo R2" _column(65) "= " %12.4f in yellow e(r2_p)
	
    ereturn display
    
	if ("`weightvar'" != "") {
		di in green "Robust regression weights stored in " in yellow "`weightvar'" 	
	}
	
	if ("`replaceweightvar'" != "") {
		di in green "Careful: " in yellow "`weightvar'" in green " already existed and now replaced with new data"
	}

	di ""
	di in green "Degrees of freedom used by FE:"
	di "{hline 17}{c TT}{hline 36}{c TRC}"
	di "FE dimension: {col 18}{c |}  Categories - Redundant: {col 55}{c |}"
	di "{hline 17}{c +}{hline 36}{c RT}"
	local j = 0
	foreach abs of local absorb {
		local `j++'
		local offset1`j' = 28 - strlen("`n`j''")
		local offset2`j' = 40 - strlen("`n`j'_red'")
		local offset3`j' = 52 - strlen("`n`j'_est'")
		if (`nest`j'' == 0) {
			local star`j' = "*"
		}
		else {
			local star`j' = " "		
		}
		di in green "`absvar`j'' {col 17} {c |}" _column(`offset1`j'') " `n`j''" "   - " _column(`offset2`j'') (1-`nest`j'')*`n`j'' + `nest`j'dof' "   = " _column(`offset3`j'') in yellow `n`j'' - (1-`nest`j'')*`n`j'' - `nest`j'dof' " `star`j'' {col 53}{c |}"
	}	
	di "{hline 17}{c BT}{hline 36}{c BRC}"
	if (`allnest' == 0) {
		di in green "* FE nested within cluster; treated as redundant for DoF calculation"
	}
	
	if "`dkraay'"=="" {
		matrix drop beta Vc Vclust b b0  
		scalar drop df_initial eff mata_nclusters scale krob r2_p mu maxiter qhat Ibar k0
	}
	else {
		matrix drop beta Vc Vdk b b0  
		scalar drop df_initial eff mata_ntime scale krob r2_p mu maxiter qhat Ibar k0 dk_lags
	}

	capture sum _temp_reghdfe_resid
	if (_rc == 0) {
		ren _temp_reghdfe_resid _reghdfe_resid
	}
	
end

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
// Mata programs
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

mata:
	void _vce_cluster() {

		real vector r, cvar
		real matrix Xr
		real scalar scale, krob
		real scalar k, n, nc, nocluster
		real scalar i
		real vector z, psi, phi, psii, psi2
		real matrix XphiXinv, info, M, xi, Vclust
 
		st_view(Xr = ., ., tokens(st_local("indepvr")), st_local("touse"))
		st_view(r = ., ., st_local("_resid_temp"), st_local("touse"))
		st_view(cvar = ., ., st_local("cvar"), st_local("touse"))
		scale = st_numscalar("scale")		
		krob = st_numscalar("krob")
		nocluster = (st_local("nocluster") != "")
		
		// Process input:
		k = cols(Xr)
		n = rows(r)
		z = r:/scale
		psi = mm_huber_psi(z,krob)
		phi = mm_huber_phi(z,krob)	
		
		// Compute VCE:
		XphiXinv = invsym(quadcross(Xr,phi,Xr))
		info = panelsetup(cvar, 1)
        nc = rows(info)
		if (nocluster == 1) {
			nc = n
		}
        M = J(k, k, 0)
		if (nc < n) { // Loop over clusters:
			for(i=1; i<=nc; i++) {
				xi = panelsubmatrix(Xr, i, info)
				psii = panelsubmatrix(psi, i, info)
				M = M + (xi' * psii) * (psii' * xi) 
			}			
		}
		else { //Else use heteroskedasticity-robust version:
			psi2 = psi :* psi
			M = quadcross(Xr, psi2, Xr)
			nc = rows(r)
		}
		
		// Combine:
		Vclust = makesymmetric(scale^2 * XphiXinv * M * XphiXinv)
		
		// Export to Stata:
		st_matrix("Vclust", Vclust)
		st_numscalar("mata_nclusters", nc)
	}
	
	void _vce_dkraay() {

		real vector r, tvar
		real matrix Xr
		real scalar scale, krob, lags
		real scalar k, nt
		real scalar t, l, wl
		real vector z, psi, phi
		real matrix XphiXinv, info, M, S, Vdk, xrt
		real vector psit, st, sl
 
		st_view(Xr = ., ., tokens(st_local("indepvr")), st_local("touse"))
		st_view(r = ., ., st_local("_resid_temp"), st_local("touse"))
		st_view(tvar = ., ., st_local("tvar"), st_local("touse"))
		scale = st_numscalar("scale")		
		krob = st_numscalar("krob")
		lags = st_numscalar("dk_lags")
		
		// Process input:
		k = cols(Xr)
		z = r:/scale
		psi = mm_huber_psi(z,krob)
		phi = mm_huber_phi(z,krob)	
		
		// Compute Driscoll-Kraay VCE:
		XphiXinv = invsym(quadcross(Xr,phi,Xr))
		
		info = panelsetup(tvar, 1)
        nt = rows(info)
		
		S = J(nt, k, 0)
		for (t=1; t<=nt; t++) {
			xrt = panelsubmatrix(Xr, t, info)
			psit = panelsubmatrix(psi, t, info)
			S[t,.] = (xrt' * psit)'
		}
		
        M = J(k, k, 0)
		for (t=1; t<=nt; t++) {
			st = S[t,.]'
			M = M + st * st'
		}
		
		for (l=1; l<=lags; l++) {
			wl = 1 - l/(lags+1)
			for (t=l+1; t<=nt; t++) {
				st = S[t,.]'
				sl = S[t-l,.]'
				M = M + wl * (st * sl' + sl * st')
			}
		}
		
		// Combine:
		Vdk = makesymmetric(scale^2 * XphiXinv * M * XphiXinv)
		
		// Export to Stata:
		st_matrix("Vdk", Vdk)
		st_numscalar("mata_ntime", nt)
	}
    
	void _scale_initial() {

		real vector e, z, w
		real scalar df, eff, n, p, scale, krob
	
		st_view(e=., ., tokens(st_local("resid_tau")), st_local("touse"))
        df=st_numscalar("df_initial")
		eff=st_numscalar("eff")
        n=rows(e)
        p = (2*n - df) / (2*n) 
        scale=mm_quantile(abs(e), 1, p) / invnormal(0.75) // For consistency with robreg
        z = e / scale
		krob=mm_huber_k(eff)
		w=mm_huber_w(z, krob)
		st_store(., st_addvar("double", st_local("w")), st_local("touse"), w)
        st_numscalar("scale", scale)
		st_numscalar("krob", krob)
	}
	
	void _update_weights() {

		real vector z, phi, w
		real scalar eff, krob
	
		z = st_data(., "_z_temp")
		eff = st_numscalar("eff")
		krob = mm_huber_k(eff)
		phi = mm_huber_phi(z,krob)			
		w = mm_huber_w(z, krob)
        st_store(., st_local("w"), w)
        st_store(., st_local("phi"), phi)
        printf(".")
    }

	void _huber_location() {

		real vector y, u, w
		real scalar eff, maxiter, tol, k, mu, mu_new, n, df, p, scale, i
	
		st_view(y = ., ., tokens(st_local("depv")), st_local("touse"))
		eff = st_numscalar("eff")
		maxiter = st_numscalar("maxiter")
		tol = st_numscalar("tol")
		k = mm_huber_k(eff)
		mu = mm_median(y)
        n = rows(y)
		df = n - 1
        p = (2*n - df) / (2*n) 
        scale = mm_quantile(abs(y :- mu), 1, p) / invnormal(0.75) 
		for (i=1; i<=maxiter; i++) {
			u = (y :- mu) :/ scale
			w = mm_huber_w(u, k)
			mu_new = sum(w:*y) / sum(w)
			if (abs(mu_new - mu) < tol) break
			mu = mu_new
		}
		st_numscalar("mu", mu)
	}

	void _pseudo_r2() {

		real vector y, r
		real scalar scale, krob, mu
		real vector z, z0, rho, rho0
		real scalar r2_p
 
		st_view(y = ., ., st_local("depv"), st_local("touse"))
		st_view(r = ., ., st_local("_resid_temp"), st_local("touse"))
		scale = st_numscalar("scale")		
		krob = st_numscalar("krob")
		mu = st_numscalar("mu")
		
		z = r:/scale
		z0 = (y :- mu) :/ scale
		rho = mm_huber_rho(z, krob)			
		rho0 = mm_huber_rho(z0, krob)
		
		r2_p = 1 - (colsum(rho) / colsum(rho0))
		st_numscalar("r2_p", r2_p)		
	}
    
end
	
	
	
