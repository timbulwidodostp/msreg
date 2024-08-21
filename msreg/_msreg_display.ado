*! version 1.0.0  21feb2020 
					//----------------------------//
					// msreg display
					//----------------------------//
program _msreg_display
	syntax [, * ]

	_get_diopts diopts extra , `options'

	if (`"`extra'"' != "") {
		di as err "option {bf:`extra'} not allowed"
		exit 198
	}
						//  header
	Head
						//  coef table
	_coef_table, `diopts'
						//  footnote
	Footnote
end
					//----------------------------//
					// head
					//----------------------------//
program Head
//--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----
						// title
	local title _n as txt `"`e(title)'"'

	local col1 = 40

						// Number of obs
	local n_obs1 as txt _col(`col1') `"Number of obs in sample 1"' 	///
		_col(67) "=" _col(69) as res %10.0fc e(N_1)

	local n_obs2 as txt _col(`col1') `"Number of obs in sample 2"' 	///
		_col(67) "=" _col(69) as res %10.0fc e(N_2)

						//  Distance metric
	local metric as txt "Distance metric" _col(17) ":" "{bf: `e(metric)'}"

						// Estimator
	local estimator as txt "Estimator" _col(17) ":" "{bf: `e(estimator)'}"

						//  nneighbor	
	local nneighbor as txt _col(`col1') "Number of matches " 	///
		_col(67) "=" _col(69) as res %10.0fc e(nneighbor)

						//  order	
	local order as txt _col(`col1') "Order of power-series " 	///
		_col(67) "=" _col(69) as res %10.0fc e(order)

						//  p-value and test
	if !missing(e(df_r)) {
		local model as txt _col(`col1') 		///
			"F(" 				///
			as res %4.0f e(df_m) 		///
			as txt "," 			///
			as res %7.0f e(df_r) as txt ")"	///
			_col(67) "= " as res %10.2f e(F)

		local pvalue _col(`col1') as txt "Prob > F"	///
			_col(67) "= " as res %10.4f e(p)
	}
	else {
		if "`e(chi2type)'" == "" {
			local chitype Wald
		}
		else	local chitype `e(chi2type)'
		local model _col(`col1') as txt `"`chitype' chi2("' ///
			as res e(df_m) as txt ")" 		///
			_col(67) "= " as res %10.2f e(chi2)
		local pvalue _col(`col1') as txt "Prob > chi2"	///
			_col(67) "= " as res %10.4f e(p)
	}

						// Display Header
	di `title' `n_obs1'
	di `n_obs2'
	di `nneighbor'
	if (`"`e(estimator)'"' == "twostep") {
		di `order'
	}
	di `estimator' `model'
	di `metric' `pvalue'
	
	di 
end

program Footnote
	di "{p 0 19 2}"
	di as txt "Common variables : " `"`e(common)'"' 
	di "{p_end}"

	di "{p 0 19 2}"
	di as txt "Matched variables: " `"`e(matched)'"'
	di "{p_end}"

	if (`"`e(footnote)'"' != "") {
		di "{p 0 4 2}"
		di `"`e(footnote)'"'
		di "{p_end}"
	}
end

