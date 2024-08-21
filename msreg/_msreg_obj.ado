*! version 1.0.0  21feb2020 
/*
	utility for msreg
*/
program _msreg_obj
	gettoken subcmd 0 : 0 , parse(" ,")
	local 0 = subinstr(`"`0'"', `","' , `" , "' , .)
	
 	local len = length(`"`subcmd'"')
        if `len'==0 {
               di as err "No subcommand specified"
               exit 198
        }

	if (`"`subcmd'"' == "create") {
		Create
	}
	else if (`"`subcmd'"' ==  "clean") {
		Clean `0'	
	}
	else if (`"`subcmd'"' ==  "name") {
		Name
	}
	else if (`"`subcmd'"' ==  "parse_syntax") {
		_msreg_parse_syntax `0'
	}
	else if (`"`subcmd'"' ==  "compute") {
		Compute
	}
	else {
		di as err "unknown subcommand {bf:`subcmd'}"
		exit 198
	}

	if (`"`s(set_local)'"' != "") {
		c_local `s(name)' `s(val)'
	}
end

					//----------------------------//
					// create object
					//----------------------------//
program Create
	mata : st_msreg_create_obj()
end
					//----------------------------//
					// clean object
					//----------------------------//
program Clean
	syntax , rc(string) state(string)
	cap mata : st_msreg_clean_obj()
	set rngstate `state'
	if (`rc') exit `rc' 
end
					//----------------------------//
					// object name
					//----------------------------//
program Name, sclass
	mata : st_msreg_objname()	

	sret local name OBJ
	sret local val `s(msreg_objname)'
	sret local set_local set_local
end

					//----------------------------//
					// compute
					//----------------------------//
program Compute
	_msreg_obj name
	mata: `OBJ'.compute()
end
