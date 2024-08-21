*! version 1.0.0   10mar2020
/*
	parse syntax for -msreg- and store result in mata object
*/
program _msreg_parse_syntax, sclass
	_on_colon_parse `0'
	local before `s(before)'
	local after `s(after)'

	local 0 `before'
	syntax , touse1(passthru) 	///
		touse2(passthru)	

	local 0 `after'
	syntax anything(equalok name=eq) 	///
		using			 	///
		[if] [in]			///
		[, vce(string)			///
		estimator(string)		///
		metric(string)			///
		NNeighbor(passthru)		///
		order(passthru)			///
		noCONstant			///
		*]
	
						// create msreg object	
	_msreg_obj create
						// append samples
	AppendSample `if' `in' `using', `touse1' `touse2' 
						// parse eq
	ParseEq `eq', `touse1' `touse2'  
						// parse estimator
	ParseEstimator, `estimator'
						// parse metric
	ParseMetric, `metric'
						// parse nneighbor
	ParseNN, `nneighbor'
						// parse order
	ParseOrder, `order'
						// parse constant
	ParseConst, `constant'
						// diopts
	_get_diopts diopts extra, `options'
						// parse vce
	ParseVCE, `vce' estimator(`estimator')
						// check extra
	CheckExtra, `extra' check

	sret local diopts `diopts'
end
					//----------------------------//
					// parse equation
					//----------------------------//
program ParseEq
	syntax anything(equalok name=eq)	///
		, touse1(string)		///
		touse2(string)

	_iv_parse `eq'
	local st_y `s(lhs)'
	local st_x1 `s(exog)'
	local st_x2 `s(endog)'
	local st_z `s(inst)'

	CheckEqSyntax, depvar(`st_y') x2(`st_x2') z(`st_z')

	_rmcoll `st_x1' if `touse1', expand
	local st_x1 `r(varlist)'

	_rmcoll `st_x2' if `touse2', expand
	local st_x2 `r(varlist)'

	_rmcoll `st_z', expand
	local st_z `r(varlist)'
						// depvar
	_fv_check_depvar `st_y'
						// x2
	CheckVarSource `st_x2', touse(`touse2') i(2) 
						// z
	CheckVarSource `st_z', touse(`touse1') i(1) 
	CheckVarSource `st_z', touse(`touse2') i(2) 
						// markout variables
	markout `touse1' `st_y' `st_x1' `st_z'
	markout `touse2' `st_z' `st_x2'

	sum `touse1' if `touse1', meanonly
	if (r(N) == 0)  error 2000

	sum `touse2' if `touse2', meanonly
	if (r(N) == 0)  error 2000

	_msreg_obj name
	mata: `OBJ'._st_y = `"`st_y'"'
	mata: `OBJ'._st_x1 = `"`st_x1'"'
	mata: `OBJ'._st_x2 = `"`st_x2'"'
	mata: `OBJ'._st_z = `"`st_z'"'
end

					//----------------------------//
					// mark sample
					//----------------------------//
program AppendSample 
	syntax [if] [in] using, ///
		touse1(string)	///
		touse2(string)	

	tempvar src
	qui append `using', generate(`src')

	marksample mytouse
	qui gen byte `touse1' = 0
	qui replace `touse1' = `mytouse' if `src' == 0

	qui gen byte `touse2' = 0
	qui replace `touse2' = `mytouse' if `src' == 1

	_msreg_obj name
	mata: `OBJ'._touse1 = `"`touse1'"'
	mata: `OBJ'._touse2 = `"`touse2'"'
end
					//----------------------------//
					// Check vars source
					//----------------------------//
program CheckVarSource
	syntax varlist(numeric fv),	///
		touse(string)		///
		i(string)		

	sum `varlist' if `touse', meanonly
	if (r(N) == 0) {
		di as err "variables {it:`varlist'} are missing "	///
			"in sample `i'"
		exit 2000
	}
end
					//----------------------------//
					// Parse estimator
					//----------------------------//
program ParseEstimator
	syntax [, ONEstep	///
		TWOstep		///
		msols]

	local method `onestep' `twostep' `msols'
	local n_method : list sizeof method

	if (`n_method' == 0) {
		local method twostep
	}
	else if (`n_method' > 1) {
		di as err "only one of {bf:onestep} or {bf:twostep} is allowed"
		exit 198
	}

	_msreg_obj name
	mata: `OBJ'._method = `"`method'"'
end
					//----------------------------//
					// parse metric
					//----------------------------//
program ParseMetric
	cap syntax [, MAHAlanobis 	///
		EUCLidean]
	local rc = _rc

	local metric `mahalanobis' `euclidean'
	local n_metric : list sizeof metric

	if (`n_metric' > 1 | `rc') {
		di as err "only one of {bf:mahalanobis} or {bf:euclidean} " ///
			"can be specified in option {bf:metric()}"
		exit 198
	}

	if (`n_metric' == 0) {
		local metric mahalanobis
	}

	_msreg_obj name
	mata: `OBJ'._metric = `"`metric'"'
end
					//----------------------------//
					// parse nneighbor
					//----------------------------//
program ParseNN
	syntax [, nneighbor(integer 1)]

	if (`nneighbor' < 1) {
		di as err "must specify a positive integer in option "	///
			"{bf:nneighbor()}"
		exit 198
	}

	if (`nneighbor' >= 11) {
		di as err "must specify a positive integer no greater "	///
			"than 10 in option {bf:nneighbor()}"
		exit 198
	}

	_msreg_obj name
	mata: `OBJ'._nn = `nneighbor'
end
					//----------------------------//
					// parse order
					//----------------------------//
program ParseOrder
	syntax [, order(integer 2)]

	if (`order' < 1) {
		di as err "must specify a positive integer in option "	///
			"{bf:order()}"
		exit 198
	}

	if (`order' >= 6) {
		di as err "must specify a positive integer no greater "	///
			"than 5 in option {bf:order()}"
		exit 198
	}

	_msreg_obj name
	mata: `OBJ'._order = `order'
end
					//----------------------------//
					//parse constant
					//----------------------------//
program ParseConst
	syntax [, noCONstant]

	if (`"`constant'"' == "noconstant") {
		local cons = 0	
	}
	else {
		local cons = 1
	}

	_msreg_obj name
	mata: `OBJ'._cons = `cons'
end
					//----------------------------//
					// check extra
					//----------------------------//
program CheckExtra
	syntax , check
end

					//----------------------------//
					// check equation syntax
					//----------------------------//
program CheckEqSyntax
	syntax [, depvar(string)	///
		x2(string)		///
		z(string)]
	
	if (`"`depvar'"' == "") {
		di as err "must specify the dependent variable"
		exit 198
	}
	
	if (`"`x2'"' == "") {
		di as err "must specify the matched variables {it:varlist_x2}"
		exit 198
	}

	if (`"`z'"' == "") {
		di as err "must specify the common variables {it:varlist_z}"
		exit 198
	}
end
					//----------------------------//
					// parse vce
					//----------------------------//
program ParseVCE
	syntax [, vi	///
		vii	///
		viii	///
		estimator(string)]
	
	if (`"`estimator'"' == "msols") {
		exit
		// NotReached
	}
	
	local vce `vi' `vii' `viii'
	
	local n_vce : list sizeof vce

	if (`n_vce' >1) {
		di as err "must specify one of {bf:vi}, {bf:vii} "	///
			"or {bf:viii} in option vce()"
		exit 198
	}

	if (`n_vce' == 0) {
		local vce vi
	}

	_msreg_obj name
	mata: `OBJ'._vce = `"`vce'"'
	mata: `OBJ'._vcetype = strupper(`"`vce'"')
end
