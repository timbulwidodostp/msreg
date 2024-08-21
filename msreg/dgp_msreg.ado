*!version 1.0.0  10jun2017

program dgp_msreg
	syntax , beta0(string)		///
		beta1(string)		///
		beta2(string)		///
		gamma(string)		///
		[sample1(passthru)	///
		sample2(passthru)	///
		sample0(passthru)	///
		n_obs(integer 1000)	///
		ratio(real 0.5) ]	

	tempvar z1 z2 z3	
	tempname C
	matrix `C' = (1, 1/sqrt(2), 1/sqrt(3) \	///
		1/sqrt(2), 1, sqrt(2)/sqrt(3) \	///
		1/sqrt(3), sqrt(2)/sqrt(3), 1)
	drawnorm `z1' `z2' `z3', double n(`n_obs') cov(`C')
	local zs `z1' `z2' `z3'

	mata : dgp_msreg(	///
		`"`beta0'"',	///
		`"`beta1'"',	///
		`"`beta2'"',	///
		`"`gamma'"',	///
		`n_obs',	///
		`"`zs'"',	///
		`ratio')
	
	makeTwoSample, `sample1' `sample2' `sample0'
end
					//-- make Two sample from one sample -//
program makeTwoSample
	syntax [,sample0(string)	///
		sample1(string)		///
		sample2(string)]

	if (`"`sample0'"' == "") {
		local sample0 sample0
	}
	
	if (`"`sample1'"' == "") {
		local sample1 sample1
	}

	if (`"`sample2'"' == "") {
		local sample2 sample2
	}

	preserve 
	drop if t == 1
	keep y x1* z* 
	save `sample1', replace
	restore

	preserve
	drop if t == 0
	keep x2* z* 
	save `sample2', replace
	restore

	keep y x1* x2* z* t
	save `sample0', replace
end

mata :
mata set matastrict on

void dgp_msreg(			///
	string scalar 	_b0,	///
	string scalar 	_b1,	///
	string scalar 	_b2,	///
	string scalar 	_r,	///
	real scalar	n_obs,	///
	string scalar	_zs,	///
	real scalar	ratio)
{
	real matrix	X1, X2, Z, Zs
	real colvector	y
	real rowvector	b1, b2, r, idx, t
	real scalar	b0, i

	Zs	= st_data(., _zs)	
	b0	= st_matrix(_b0)		
	b1	= st_matrix(_b1)
	b2	= st_matrix(_b2)
	r	= st_matrix(_r)

	Z	= J(n_obs, 0, .)
	for (i=1; i<=length(r); i++) {
		Z	= (Z, 4*normal(Zs[.,i]):-2)
	}

	X1	= J(n_obs, 0, .)
	for (i=1; i<=length(b1); i++) {
		X1	= (X1, rowsum(Z)+ rnormal(n_obs, 1, 0, 1))
	}


	X2	= J(n_obs, length(b2), .)
 	X2[.,1]	= __gzx2(Z, 0.25)+rnormal(n_obs, 1, 0, 1)	
//!!rm 	X2[.,2]	= __gzx2(Z, 0.25)+rnormal(n_obs, 1, 0, 1)	
	X2[.,2] = __gzx3(Z) + rnormal(n_obs, 1, 0, 1)

	y	= b0:+ X1*b1' + X2*b2' + Z*r' + rnormal(n_obs, 1, 0, 1)

	t	= (1..rows(y))'
	t	= t:<= ratio*rows(y)
	
	dgp_store("y", y)
	dgp_store("x1", X1)
	dgp_store("x2", X2)
	dgp_store("z", Z)
	dgp_store("t", t)
}

real colvector	__gzx3(		///
	real matrix	Z)
{
	real colvector	res
	real matrix 	absZ

	absZ	= abs(Z/2)

	res	= rowsum(4*sqrt(absZ:*(1:-absZ))	///
		:*sin(2*pi()*(1.05):/(absZ:+0.05)))
	return(res)
}
					//-- g(z) = z + 5/tau*phi(z/tau)--//
real colvector __gzx2(		///
	real matrix Z,		///
	real scalar tau)
{
	real colvector res

	res	= rowsum(Z :+ (5/tau):*normalden(Z:/tau))
	return(res)
}
					//-- store variable --//
void dgp_store(			///	
	string scalar prefix,	///
	real matrix  _X)
{
	real rowvector	idx
	real scalar	i
	string scalar	s

	if (cols(_X) == 1) {
		idx	= st_addvar("double", prefix)
		st_store(., idx, _X)
		return
		// NotReached
	}

	for (i=1; i<= cols(_X); i++) {
		s	= prefix+strofreal(i)
		idx	= st_addvar("double", s)
		st_store(., idx, _X[., i])
	}
}

end
