*! version 1.0.0  21feb2020 
program msreg
	if replay() {
		if (`"`e(cmd)'"' == "msreg") {
			_msreg_display `0'
		}
		else {
			error 301
		}
	}
	else {
		SetRngState
		local state `s(state)'

		cap noi Estimate `0'
		local rc = _rc

		nobreak _msreg_obj clean, rc(`rc') state(`state')
	}
end

					//----------------------------//
					// Estimate
					//----------------------------//
program Estimate
						// preserve start
	preserve
						// parse syntax	
	tempvar touse1 touse2 
	_msreg_obj parse_syntax, touse1(`touse1') touse2(`touse2'): `0'
	local diopts `s(diopts)'
						//  compute
	_msreg_obj compute	
						//  display
	_msreg_display, `diopts'
						//  preserve end
	restore
end

					//----------------------------//
					//  set random seed 
					//----------------------------//
program SetRngState, sclass
	sret local state = c(rngstate)
	set seed 12344671
end
