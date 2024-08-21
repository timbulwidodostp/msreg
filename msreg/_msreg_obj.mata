*! version 1.0.0  21feb2020

findfile "_msreg_macdefs.matah", path(".")
include `"`r(fn)'"'

mata :
mata set matastrict on

class `MSREG' {
	public :
		string scalar	_touse1		// touse in sample 1
		string scalar	_touse2		// touse in sample 2
		string scalar	_st_y		// depvar
		string scalar	_st_x1		// x1
		string scalar	_st_x2		// x2
		string scalar	_st_z		// z
		string scalar	_method		// estimator
		string scalar	_metric		// metric
		real scalar	_nn		// neighbor of neighbor match
		real scalar	_order		// order of polynomial approx 
		real scalar	_cons		// constant or not
		string scalar	_vce		// vce (vi, vii or viii)
		string scalar	_vcetype	// vce (VI, VII OR VIII)

		void	compute()

	private:
		real colvector	_y1		// y from sample 1
		real matrix	_X1		// x1 from sample 1
		real matrix	_X2		// x2 from sample 2
		real matrix	_Z1		// Z from sample 1
		real matrix	_Z2		// Z from sample 2
		real matrix	_X2_md		// x2 from sample 1 matched
		real matrix	_W		// matched sample
		real scalar	_n		// sample size in sample 1
		real scalar	_m		// sample size in sample 2
		real scalar	_p		// number of variables
		pointer(real vector) vector	_idx_md // matched idx

		real vector	_b		// e(b)
		real matrix	_V		// e(V)
		string matrix	_bcs		// e(b) colstripe
		string scalar	_footnote	// footnote for vce

		void	get_data()
		void	get_metric()
		void	match_X2()
		void	match_zi()
		void	msols()
		void	norm_xax()
		void	onestep()
		void	twostep()
		void	post_result()
		void	post_b_V()
		void	post_aux()
		void	load_data()
		void	get_W()
		void	adjust_y()
		void	power_series()
		void	rmcoll()
		void	get_footnote()
}

					//----------------------------//
					// compute
					//----------------------------//
void `MSREG'::compute()
{
	// get data
	get_data()

	// compute b and V based on method
	if (_method == "msols") {
		msols()
	}
	else if (_method == "onestep") {
		onestep()
	}
	else if (_method == "twostep") {
		twostep()
	}
	
	// post result
	post_result()
}
					//----------------------------//
					// compute msols
					//----------------------------//
void `MSREG'::msols()
{
	// msols
	_msreg_msols(_y1, _W, _b, _V)
}
					//----------------------------//
					// load data
					//----------------------------//
void `MSREG'::load_data()
{
	//  read data
	_y1 = st_data(., _st_y, _touse1)

	if (_st_x1 != "") {
		_X1 = st_data(., _st_x1, _touse1)
	}
	else {
		_X1 = J(rows(_y1), 0, 0)
	}
	_X2 = st_data(., _st_x2, _touse2)
	_Z1 = st_data(., _st_z, _touse1)
	_Z2 = st_data(., _st_z, _touse2)
	_n = rows(_y1)
	_m = rows(_Z2)
}
					//----------------------------//
					// get W matrix
					//----------------------------//
void `MSREG'::get_W()
{
	//  matched sample
	_W = (_X1, _X2_md, _Z1)

	_p = cols(_W)
	_bcs = J(_p, 2, "")
	_bcs[., 2] = tokens(invtokens((_st_x1, _st_x2, _st_z)))'

	// add constant if necessary
	if (_cons) {
		_W = _W, J(_n, 1, 1)
		_bcs = _bcs \ ("", "_cons")
		_p = _p + 1
	}

}
					//----------------------------//
					// get data
					//----------------------------//
void `MSREG'::get_data()
{
	// load data
	load_data()

	//  match X2
	match_X2()

	// get W matrix
	get_W()

	// get footnote
	get_footnote()
}

					//----------------------------//
					// match X2
					//----------------------------//
void `MSREG'::match_X2()
{
	real matrix	A
	real scalar	i
	real rowvector	zi
	real colvector	sidx
	pointer(vector) scalar	ptr_idi

	sidx = runiform(rows(_Z2), 1)

	get_metric(A)

	_X2_md = J(rows(_y1), cols(_X2), .)
	_idx_md = J(_n, 1, NULL)

	for (i=1; i<= rows(_Z1); i++) {
		zi = _Z1[i, .]	
		match_zi(sidx, zi, A, ptr_idi)
		_idx_md[i] = ptr_idi
		_X2_md[i, .] = mean(_X2[*ptr_idi, .])
	}
}
					//----------------------------//
					// get metric
					//----------------------------//
void `MSREG'::get_metric(real matrix	A)
{
	real matrix	Z
	real rowvector	Zbar
	real scalar	N

	Z = (_Z1 \ _Z2)
	N = rows(Z)
	Zbar = mean(Z)
	Z = Z :- Zbar

	A = cross(Z, Z)/N

	if (_metric == "mahalanobis") {
		A = invsym(A)
	}
	else if (_metric == "euclidean") {
		A = invsym(diag(A))
	}
}
					//----------------------------//
					// norm A = (x'Ax)^{1/2}
					//----------------------------//
void `MSREG'::match_zi(		///
	real colvector	sidx,	///
	real rowvector	zi,	///
	real matrix	A,	///
	pointer(vector) scalar	ptr_id)
{
	real colvector	idx
	real scalar	i, n_tie
	real matrix	res

	// compute norm xax
	norm_xax(zi, A, res)

	// sort norm
	idx = order((res, sidx), (1, 2))

	n_tie = sum(res:==0)

	// get the first K neareast neighbors

	if (n_tie <= 1) {
		ptr_id = &idx[1.._nn]
	}
	else {
		ptr_id = &idx[1..(_nn + n_tie - 1)]
	}
}
					//----------------------------//
					// compute norm
					//----------------------------//
void `MSREG'::norm_xax(		///
	real rowvector	zi,	///
	real matrix	A,	///
	real matrix 	res)
{
	res = _Z2 :- zi
	res = sqrt(rowsum((res*A):*res))
}
					//----------------------------//
					// post result
					//----------------------------//
void `MSREG'::post_result()
{
	post_b_V()
	post_aux()
}

					//----------------------------//
					// post b and V
					//----------------------------//
void `MSREG'::post_b_V()
{
	string scalar	s, bs, Vs

	bs = st_tempname()
	st_matrix(bs, _b')
	st_matrixcolstripe(bs, _bcs)

	Vs = st_tempname()
	st_matrix(Vs, _V)
	st_matrixcolstripe(Vs, _bcs)
	st_matrixrowstripe(Vs, _bcs)

	s = sprintf("eret post %s %s, esample(%s) buildfvinfo", bs, Vs, _touse1)
	stata(s)
}

					//----------------------------//
					// one step estimation msii
					//----------------------------//
void `MSREG'::onestep()
{
	// msii
	_msreg_msii(	///
		_vce,	///
		_y1, 	///
		_W, 	///
		_X1, 	///
		_X2, 	///
		_Z1, 	///
		_Z2, 	///
		_n, 	///
		_m, 	///
		_cons, 	///
		_nn, 	///
		_b,	///
		_V)
}
					//----------------------------//
					// two step estimation msii
					//----------------------------//
void `MSREG'::twostep()
{
	real colvector	y1_adj, b0
	real matrix	X2_orig, Z2_orig

	// note: msii will resort X2 and Z2, we want to keep _id_md the same,
	// so in the first step, we should not change the order of _X2 and _Z2
	X2_orig = _X2
	Z2_orig = _Z2

	// first step msii
	_msreg_msii(		///
		_vce,		///
		_y1,		///
		_W, 		///
		_X1, 		///
		X2_orig, 	///
		_Z1, 		///
		Z2_orig, 	///
		_n, 		///
		_m, 		///
		_cons, 		///
		_nn, 		///
		b0)

	// adjust y
	adjust_y(b0, y1_adj)

	// second msii
	_msreg_msii(	///
		_vce,	///
		y1_adj,	///
		_W, 	///
		_X1, 	///
		_X2, 	///
		_Z1, 	///
		_Z2, 	///
		_n, 	///
		_m, 	///
		_cons, 	///
		_nn, 	///
		_b,	///
		_V)
}



					//----------------------------//
					// post auxilary 
					//----------------------------//
void `MSREG'::post_aux()
{
	string scalar s

	st_global("e(title)", "Matched sample regression")	

	st_numscalar("e(N_1)", rows(_y1))
	st_numscalar("e(N_2)", rows(_Z2))
	st_numscalar("e(nneighbor)", _nn)
	st_global("e(estimator)", _method)
	st_global("e(metric)", _metric)
	st_global("e(matched)", _st_x2)
	st_global("e(common)", _st_z)
	if (_method == "msols") {
		st_global("e(vcetype)", "Robust")
		st_global("e(vce)", "robust")
	}
	else {
		st_global("e(vce)", _vce)
		st_global("e(vcetype)", _vcetype)
	}

	s = sprintf("qui test %s %s %s", _st_x1, _st_x2, _st_z)
	stata(s)

	st_numscalar("e(chi2)", st_numscalar("r(chi2)"))
	st_numscalar("e(p)", st_numscalar("r(p)"))
	st_numscalar("e(df_m)", st_numscalar("r(df)"))
	st_global("e(chi2type)", "Wald")

	st_numscalar("e(order)", _order)

	st_global("e(footnote)", _footnote)

	st_global("e(cmd)", "msreg")
}
					//----------------------------//
					//  adjust y to remove the second order
					//  bias
					//----------------------------//
void `MSREG'::adjust_y(		///
	real colvector	b0,	///
	real colvector	y_adj)
{
	real matrix	Zp_1, Zp_2, B_p, g_Z1, g_Z2
	real vector	p, b2, lbda
	real scalar	i

					// get power series of Z in sample 1
	Zp_1 = _Z1	
	power_series(_order, _Z1, Zp_1)
	rmcoll(Zp_1, p)
	Zp_1 = (J(_n, 1, 1), Zp_1)

					// get power series of Z in sample 2
	Zp_2 = _Z2
	power_series(_order, _Z2, Zp_2)
	Zp_2 = select(Zp_2, p)
	Zp_2 = (J(_m, 1, 1), Zp_2)

					// get B from power series regression 
	B_p = invsym(cross(Zp_2, Zp_2))*cross(Zp_2, _X2)

					// get estimated x2_hat = g(z)
	g_Z1 = Zp_1*B_p
	g_Z2 = Zp_2*B_p

					// get b_msii for X2
	b2 = b0[cols(_X1)+1..cols(_X1) + cols(_X2)]

					// get lambda	
	lbda = J(_n, 1, .)	
	for (i=1; i<=_n; i++) {
		lbda[i] = (g_Z1[i, .] - mean(g_Z2[*(_idx_md[i]), .]))*b2
	}
	y_adj = _y1 - lbda
}
					//----------------------------//
					//  construct power series
					//----------------------------//
void `MSREG'::power_series(	///
	real scalar	i,	///
	real matrix	Z_orig,	///
	real matrix	Zp)
{
	real matrix	tmp
	real scalar	j, n, k

	n	= rows(Z_orig)
	k	= cols(Z_orig)
	
	if (i == 1) {
		Zp	= Z_orig
	}
	else {
		power_series(i-1, Z_orig, Zp)
		tmp	= J(n, 0, .)
		for (j=1; j<=k; j++) {
			tmp	= (tmp,  Z_orig[., j]:*Zp)
		}
		Zp	= (Zp, tmp)
	}
}
					//----------------------------//
					//  remove collinearity
					//----------------------------//
void `MSREG'::rmcoll(			///
	real matrix 	Z,		///
	real rowvector	p)
{
	p	= (diagonal(invsym(cross(Z, Z))):!=0)'
	Z	= select(Z, p) 
}

					//----------------------------//
					// get footnote
					//----------------------------//
void `MSREG'::get_footnote()
{
	string scalar	 s, s1, s2, s3

	if (_method == "msols") {
		return
		// NotReached
	}

	s1 = "there is at most one continuous common variable."
	s2 = "there are at most two continuous common variables."
	s3 = "there are at most three continuous common variables."

	if (_vce == "vi") {
		if (_method == "onestep") {
			s = s1
		}
		else if (_method == "twostep") {
			s = s3
		}
		_footnote = "Note: VI Std. Err. assumes that " + ///
			"sample size ratio between sample 1 and sample " + ///
			"2 converges to a nonzero constant and " + s
	}
	else if (_vce == "vii") {
		if (_method == "onestep") {
			s = s2
		}
		else if (_method == "twostep") {
			s = s3
		}
		_footnote = "Note: VII Std. Err. assumes that " + ///
			"sample size ratio between sample 1 and sample " + ///
			"2 converges to zero and " + s
	}
	else if (_vce == "viii") {
		if (_method == "onestep") {
			s = s1
		}
		else if (_method == "twostep") {
			s = s3
		}
		_footnote = "Note: VIII Std. Err. assumes that " + ///
			"sample size ratio between sample 1 and sample " + ///
			"2 converges to infinity and there " + s
	}
}


/*----------------------------------------------------------------------------
	computation utilities
----------------------------------------------------------------------------*/

					//----------------------------//
					// matched sample OLS 
					//----------------------------//
void _msreg_msols(		///
	real colvector	_y1,	///
	real matrix	_W,	///
	real vector	_b,	///
	real matrix	_V)
{
	real matrix	invWpW, uhat
	real colvector	eps

	invWpW = invsym(cross(_W, _W))

	_b = invWpW*cross(_W, _y1)

	eps = _y1 - _W*_b
	uhat = _W:*eps
	_V = invWpW*(cross(uhat, uhat))*invWpW
}

					//----------------------------//
					// get b and V for onestep (msii)
					//----------------------------//
/*
	one-step bias-corrected estimator (matched sample indirect inference)

	b_{II} = Pw^{-1}*Rw

	where,	Rw = 1/n * W'y
		Pw = Qw - 1/K * S
		Qw = 1/n * W'W

		S = diag(0_{d1+1}, S2, 0_{d3})
*/

void _msreg_msii(		///
	string scalar	_vce,	///
	real colvector	_y1,	///
	real matrix	_W,	///
	real matrix	_X1,	///
	real matrix	_X2,	///
	real matrix	_Z1,	///
	real matrix	_Z2,	///
	real scalar	_n,	///
	real scalar	_m,	///
	real scalar	_cons,	///
	real scalar	_nn,	///
	real vector	_b,	///
	| real matrix	_V)
{
	real matrix	invPw, Sigma

	// reorder X2 w.r.t Z2
	_msreg_reorder_S2(_Z2, _X2, _m)

	// get Sigma 
	_msreg_get_Sigma(_X1, _X2, _Z1, _cons, _m, Sigma)

	// get point estimate
	_msreg_msii_b(_y1, _W, _n, _nn, Sigma, _b, invPw)

	// get V estimate
	_msreg_msii_V(_vce, _y1, _W, _X1, _X2, _Z1, _Z2,	///
		Sigma, _b, _n, _m, _cons, _nn, invPw, _V)
}
					//----------------------------//
					//  reorder sample 2 w.r.t. Z
					//----------------------------//
void _msreg_reorder_S2(		///
	real matrix	_Z2,	///
	real matrix	_X2,	///
	real scalar	_m)
{
	real scalar	zmin, j, found, k
	real vector	p, res, z_jm1, pres, sidx
	

	p = J(_m, 1, 0)
	sidx = runiform(_m, 1)

	zmin = min(_Z2[.,1])

	p[1] = min(selectindex(_Z2[.,1]:<=zmin))

	for(j = 2; j<= _m; j++) {
		z_jm1 = _Z2[p[j-1],.]
		_msreg_tr_norm(_Z2, z_jm1, res)
		pres = order((res, sidx), (1, 2) )

		found = 0
		k = 1
		while (!found) {
			if (sum(pres[k]:==p) == 0) {
				p[j] = pres[k]
				found = 1
			}
			else {
				k++	
			}
		}
	}

	_Z2 = _Z2[p, .]
	_X2 = _X2[p, .]
}

					//----------------------------//
					//  trace norm
					//----------------------------//
void _msreg_tr_norm(		///
	real matrix A,		///
	real rowvector	z,	///
	real vector res)
{
	res = sqrt(rowsum((A:-z):^2))
}
					//----------------------------//
					//  get Sigma
					//----------------------------//
void _msreg_get_Sigma(		///
	real matrix	_X1,	///
	real matrix	_X2,	///
	real matrix	_Z1,	///
	real scalar	_cons,	///
	real scalar	_m,	///
	real matrix	_Sigma)
{
	real matrix	S2, D_X2
	real scalar	j, d1, d2, dz


	S2 = J(cols(_X2), cols(_X2), 0)
	for (j=2; j<= _m; j++) {
		D_X2 = _X2[j, .] - _X2[j-1, .]
		S2 = S2 + cross(D_X2, D_X2)
	}
	S2 = S2/(2*(_m -1))

	d1 = cols(_X1)
	d2 = cols(_X2)
	dz = cols(_Z1) + _cons

	_Sigma = J(d1, d1, 0)
	_Sigma = blockdiag(_Sigma, S2)
	_Sigma = blockdiag(_Sigma, J(dz, dz, 0))
}
					//----------------------------//
					//  msii beta
					//----------------------------//
void _msreg_msii_b(		///
	real colvector	_y1,	///
	real matrix	_W,	///
	real scalar	_n,	///
	real scalar	_nn,	///
	real matrix	_Sigma,	///
	real vector	_b,	///
	real matrix	invPw)
{
	real matrix	Rw, Qw, Pw

	Rw = cross(_W, _y1)/_n
	Qw = cross(_W, _W)/_n
	Pw = Qw - _Sigma/_nn

	invPw = qrinv(Pw)
	_b = invPw*Rw
}
					//----------------------------//
					// msii V
					//----------------------------//
void _msreg_msii_V(		///
	string scalar	_vce,	///
	real colvector	_y1,	///
	real matrix	_W,	///
	real matrix	_X1,	///
	real matrix	_X2,	///
	real matrix	_Z1,	///
	real matrix	_Z2,	///
	real matrix	_Sigma,	///
	real vector	_b,	///
	real scalar	_n,	///
	real scalar	_m,	///
	real scalar	_cons,	///
	real scalar	_nn,	///
	real matrix	invPw,	///
	real matrix	_V)
{
	real colvector	ehat, uhat, b2, Wbar
	real matrix	Om, O_11A, O_22, G1, G0, Sigma2, bsb, Vg2
	real scalar	d1, d2, dz

	// Omega_11A
	ehat = _y1 - _W*_b
	uhat = _W:*ehat :+ (_Sigma*_b/_nn)'
	O_11A = cross(uhat, uhat)/_n

	// Omega_22
	d1 = cols(_X1)
	d2 = cols(_X2)
	dz = cols(_Z1) + _cons
	Sigma2 = _Sigma[d1+1::d1+d2, d1+1::d1+d2]
	b2 = _b[d1+1::d1+d2] 

	G1 = _msreg_msii_Gamma(1, _m, _X2, Sigma2, b2)
	G0 = _msreg_msii_Gamma(0, _m, _X2, Sigma2, b2)

	O_22 = J(d1, d1, 0)
	O_22 = blockdiag(O_22, G0+2*G1)
	O_22 = blockdiag(O_22, J(dz, dz, 0))
	O_22 = O_22/_nn^2
	O_22 = (O_22 + O_22')/2

	// Omega
	bsb = b2'*Sigma2*b2
	Wbar = mean(_X1), mean(_X2),  mean((_Z1 \ _Z2)), 1
	Wbar = Wbar'
	Vg2 = _X2:-mean(_X2)
	Vg2 = cross(Vg2, Vg2)/(_m - 1) - Sigma2

	Om = J(d1, d1, 0)
	Om = blockdiag(Om, bsb*Vg2 + G0 - 2*G1)
	Om = blockdiag(Om, J(dz, dz, 0))
	Om = Om/_nn^2
	Om = Om + bsb*Wbar*Wbar'
	Om = Om*_n/_m + O_11A
	
	Om = (Om + Om')/2

	// V
	if (_vce == "vi") {
		_V = invPw*Om*invPw/_n
	}
	else if (_vce == "vii") {
		_V = invPw*O_11A*invPw/_n
	}
	else if (_vce == "viii") {
		_V = invPw*O_22*invPw/_m
	}
	else {
		errprintf("misspecified vce()\n")
		exit(198)
	}
}
					//----------------------------//
					//  gamma function
					//----------------------------//
real matrix _msreg_msii_Gamma(		///
	real scalar	_l,		///
	real scalar	_m,		///
	real matrix	_X2,		///
	real matrix	_Sigma2,	///
	real colvector	_b2)
{
	real matrix	D_X2, G, Gj, Gjl
	real scalar	s0, e0, j, d2

	D_X2 = _X2[2::_m, .] - _X2[1::_m-1, .]

	s0 = max((2, 2+_l)') -1
	e0 = min((_m, _m+_l)')-1
	d2 = cols(_X2)
	G = J(d2, d2, 0)

	for (j=s0; j<=e0; j++) {
		Gj = cross(D_X2[j, .], D_X2[j, .])/2 - _Sigma2
		Gj = Gj*_b2

		Gjl = cross(D_X2[j-_l, .], D_X2[j-_l, .])/2 - _Sigma2
		Gjl = Gjl*_b2
		G = G + Gj * Gjl'
	}

	G = G/(_m - 1)

	return(G)
}

end
