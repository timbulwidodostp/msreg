*! version 1.0.0  21feb2020

findfile "_msreg_macdefs.matah", path(".")
include `"`r(fn)'"'

mata :
mata set matastrict on

/*----------------------------------------------------------------------------
		Stata interface to MSREG class
----------------------------------------------------------------------------*/

					//----------------------------//
					// create msreg object
					//----------------------------//
void st_msreg_create_obj() 
{
	pointer scalar	p

	if ((p = findexternal("`TMP_MSREG'")) == NULL) {
		p = crexternal("`TMP_MSREG'")
	}
	else {
		errprintf("_msreg_obj already exists\n")
		exit(3001)
	}

	*p = `MSREG'()
}
					//----------------------------//
					// clean msreg object
					//----------------------------//
void st_msreg_clean_obj() 
{
	rmexternal("`TMP_MSREG'")
}
					//----------------------------//
					// get msreg object name
					//----------------------------//
void st_msreg_objname()
{
	pointer p

	p = findexternal("`TMP_MSREG'")
	if (p == NULL) {
		errprintf("%s object not found\n", "`MSREG'")
		exit(3001)
	}

	st_global("s(msreg_objname)", `"`TMP_MSREG'"')
}

end

