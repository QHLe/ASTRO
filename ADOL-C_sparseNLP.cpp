
#define SPARSE_HESS

#include <cassert>
#include "ADOL-C_sparseNLP.hpp"
//#include "examples/lowthr/model.hpp"
#include "model.hpp"
using namespace Ipopt;

/* Constructor. */
MyADOLC_sparseNLP::MyADOLC_sparseNLP() {
	x_lam			= NULL;
	cind_g 			= NULL;
	rind_g			= NULL;
	jacval 			= NULL;

	NLP_x_lb		= NULL;
	NLP_x_ub		= NULL;
	NLP_g_lb		= NULL;
	NLP_g_ub		= NULL;
	OCP_structure 	= NULL;
	NLP_x_opt		= NULL;
}

MyADOLC_sparseNLP::~MyADOLC_sparseNLP()
{
	delete[] NLP_x_opt;
}

template<class T> bool  MyADOLC_sparseNLP::eval_obj(Index n, const T *x, T& obj_value) {
  // return the value of the objective function

//	Index n_phases		= OCP_structure[0];
	Index n_nodes 		= OCP_structure[1];
	Index n_states 		= OCP_structure[2];
	Index n_controls 	= OCP_structure[3];
	Index n_param 		= OCP_structure[4];
//	Index n_events		= OCP_structure[5];
//	Index n_path		= OCP_structure[6];
//	Index n_linkages	= OCP_structure[7];

	T *y0		= new T [n_states];
	T *yf		= new T [n_states];
	T **y 		= new T *[n_nodes];
	T **u 		= new T *[n_nodes];
	T *param	= new T [n_param];
	T *delta	= new T [n_nodes - 1];
	T *t 		= new T [n_nodes];

	T tf 		= x[n-1]*NLP_x_sf[n-1];;
	T t0		= x[n-2]*NLP_x_sf[n-2];;

	t[0]		= t0;
	for (Index i 	= 0; i < n_nodes - 1; i++) {
		delta[i]	= (tf-t0)*node_str[i];
		t[i+1]		= t[i] + delta[i];
	}

	for (Index i = 0; i < n_nodes; i += 1) {
		y[i]	= new T [n_states];
		u[i]	= new T [n_controls];
	}
	Index idx = 0;

	while (idx < n_nodes*(n_states + n_controls) + n_param) {

		for (Index i = 0; i < n_nodes; i += 1)	{
			for (Index j = 0; j < n_states; j += 1){
				y[i][j] 	= x[idx]*NLP_x_sf[idx];
				idx++;
			}
			for (Index j = 0; j < n_controls; j += 1){
				u[i][j] 	= x[idx]*NLP_x_sf[idx];
				idx++;
			}
		}
		for (Index i = 0; i < n_param; i += 1)
		{
			param[i]			= x[idx]*NLP_x_sf[idx];
			idx++;
		}
	}
	for (Index i = 0; i < n_states; i += 1)
	{
		y0[i] 		= y[0][i];
		yf[i]		= y[n_nodes - 1][i];
	}
	obj_value = 0.0;
	for (Index i = 0; i<n_nodes-1; i++)
		obj_value += (lagrange_cost(y[i], u[i], param, t[i], 1)
					+ lagrange_cost(y[i+1], u[i+1], param, t[i+1], 1))*delta[i]*0.5;

	obj_value += endpoint_cost (y0, yf, param, t0, tf, 1);

	for (Index i = 0; i < n_nodes; i += 1) {
		delete[] y[i];
		delete[] u[i];
	}
	delete[] y;
	delete[] u;
	delete[] y0;
	delete[] yf;
	delete[] param;
   	delete[] delta;
	delete[] t;

	return true;
}

template<class T> bool  MyADOLC_sparseNLP::eval_constraints(Index n, const T *x, Index m, T* g) {

//	Index n_phases		= OCP_structure[0];
	Index n_nodes 		= OCP_structure[1];
	Index n_states 		= OCP_structure[2];
	Index n_controls 	= OCP_structure[3];
	Index n_param 		= OCP_structure[4];
	Index n_events		= OCP_structure[5];
	Index n_path		= OCP_structure[6];
//	Index n_linkages	= OCP_structure[7];

	T *y_start		= new T [n_states];
	T *y_end		= new T [n_states];
	T **y 			= new T *[n_nodes];
	T **f		 	= new T *[n_nodes];
	T **u 			= new T *[n_nodes];
	T **path		= new T *[n_nodes];
	T *param		= new T [n_param];
	T *e			= new T [n_events];

	T **y_m			= new T *[n_nodes-1];
	T **f_m			= new T *[n_nodes-1];
	T **u_m			= new T *[n_nodes-1];
	T *t_m			= new T [n_nodes - 1];


	T *delta	 	= new T [n_nodes - 1];
	T *t	 		= new T [n_nodes];

	T tf 			= x[n-1]*NLP_x_sf[n-1];
	T t0			= x[n-2]*NLP_x_sf[n-2];

	t[0]			= t0;
	for (Index i = 0; i < n_nodes - 1; i++) {
		delta[i]	= (tf-t0)*node_str[i];
		t[i+1]		= t[i] + delta[i];
		t_m[i]		= (t[i] + t[i+1])/2;
	}

	for (Index i = 0; i < n_nodes; i += 1) {
		y[i]		= new T [n_states];
		path[i]		= new T [n_path];
		f[i] 		= new T [n_states];
		u[i]		= new T [n_controls];
	}

	for (Index i = 0; i < n_nodes - 1; i++) {
		y_m[i]		= new T [n_states];
		f_m[i]		= new T [n_states];
		u_m[i]		= new T [n_controls];
	}

	Index idx_n = 0;
	while (idx_n < n_nodes*(n_states + n_controls) + n_param) {
		for (Index i = 0; i < n_nodes; i += 1)	{
			for (Index j = 0; j < n_states; j += 1){
				y[i][j] 	= x[idx_n]*NLP_x_sf[idx_n];
				idx_n++;
			}
			for (Index j = 0; j < n_controls; j += 1){
				u[i][j] 	= x[idx_n]*NLP_x_sf[idx_n];
				idx_n++;
			}
			derivatives(f[i], 	y[i], 	u[i], 	param, 	t[i], 	1);
		}
		for (Index i = 0; i < n_param; i += 1)
		{
			param[i]			= x[idx_n]*NLP_x_sf[idx_n];
			idx_n++;
		}
	}

	for (Index i = 0; i < n_nodes - 1; i += 1)	{
		for (Index j = 0; j < n_states; j += 1){
			y_m[i][j] 	= (y[i][j]+y[i+1][j])/2 + delta[i]/8*(f[i][j]-f[i+1][j]);
		}
		for (Index j = 0; j < n_controls; j += 1){
			u_m[i][j] 	= (u[i][j]+u[i+1][j])/2;
		}
		derivatives(f_m[i], y_m[i], u_m[i], param, t_m[i], 1);
	}

	for (Index i = 0; i < n_states; i += 1)
	{
		y_start[i] 		= y[0][i];
		y_end[i]		= y[n_nodes - 1][i];
	}
	events(e, y_start, y_end, param, t0, tf, 1);

	Index idx_m = 0;
	while( idx_m < m) {
		for (Index i = 0; i < n_nodes; i += 1) {


			for (Index j = 0; j < n_path; j += 1) {
				cout<<"path idx_m = "<<idx_m<<"\n";
				g[idx_m]	= 	0.0/NLP_g_sf[idx_m];	//need to be implemented
				idx_m++;
			}
			for (Index j = 0; j < n_states; j += 1) {
				if(i < n_nodes - 1) {
					//trapezoidal
					//g[idx]	= 	y[i+1][j] - y[i][j] - delta[i]/2.0*(f[i][j] + f[i+1][j]);
					//hermite simpson
					g[idx_m] 		= 	(y[i+1][j] - y[i][j] - delta[i]/6*(f[i][j]+4*f_m[i][j]+f[i+1][j]))/NLP_g_sf[idx_m];
					idx_m++;
				}
			}
		}
		for (Index i = 0; i < n_events; i += 1)
		{
			g[idx_m]	= 	e[i]/NLP_g_sf[idx_m];
			idx_m++;
		}
	}

	for (Index i = 0; i < n_nodes; i += 1)
	{
		delete[] y[i];
		delete[] path[i];
		delete[] f[i];
		delete[] u[i];
	}

	for (Index i = 0;i< n_nodes - 1;i++) {
		delete[] y_m[i];
		delete[] f_m[i];
		delete[] u_m[i];
	}
	delete[] y_m;
	delete[] f_m;
	delete[] u_m;


	delete[] y;
	delete[] path;
	delete[] u;
	delete[] f;
   	delete[] y_start;
   	delete[] y_end;
   	delete[] param;
   	delete[] e;
   	delete[] delta;
  	delete[] t;
  	delete[] t_m;

	return true;
}

bool MyADOLC_sparseNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{
	n = NLP_n;
	m = NLP_m;
	cout<<"n = "<<n<<"\n";
	cout<<"m = "<<m<<"\n";
	cout<<"generate tapes\n";
	generate_tapes(n, m, nnz_jac_g, nnz_h_lag);
	cout<<"end ini\n";
  // use the C style indexing (0-based)
	index_style = C_STYLE;

	return true;
}

bool MyADOLC_sparseNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
										Index m, Number* g_l, Number* g_u)
{

	for (Index i = 0; i < NLP_n; i += 1) {
		x_l[i] = NLP_x_lb[i];
		x_u[i] = NLP_x_ub[i];
	}
	for (Index i = 0; i < NLP_m; i += 1) {
		g_l[i] = NLP_g_lb[i];
		g_u[i] = NLP_g_ub[i];
	}

	return true;
}

bool MyADOLC_sparseNLP::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{

	assert(init_x == true);
	assert(init_z == false);
	assert(init_lambda == false);
	cout<<"get starting point\n";

	for (Index i = 0; i < n; i += 1) {
		x[i]	= guess[i];
	}
	cout<<"end of getting starting point\n";
	return true;
}

void 	MyADOLC_sparseNLP::setNLP_structure(Index n, Index m, Index* structure) {
	NLP_n 			= n;
	NLP_m 			= m;
	OCP_structure 	= structure;
}

void 	MyADOLC_sparseNLP::setBounds (Number* x_lb, Number* x_ub, Number* g_lb, Number* g_ub){

	NLP_x_lb = x_lb;
	NLP_x_ub = x_ub;
	NLP_g_lb = g_lb;
	NLP_g_ub = g_ub;

}

void 	MyADOLC_sparseNLP::setSF (Number* x_sf, Number* g_sf){
	NLP_x_sf = x_sf;
	NLP_g_sf = g_sf;

}

bool MyADOLC_sparseNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
	eval_obj(n,x,obj_value);
	return true;
}

bool MyADOLC_sparseNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
	gradient(tag_f,n,x,grad_f);
	return true;
}

bool MyADOLC_sparseNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
	eval_constraints(n,x,m,g);
	return true;
}

bool MyADOLC_sparseNLP::eval_jac_g(Index n, const Number* x, bool new_x,
                       Index m, Index nele_jac, Index* iRow, Index *jCol,
                       Number* values)
{
	if (values == NULL) {
	// return the structure of the jacobian

		for(Index idx=0; idx<nnz_jac; idx++)
		{
			iRow[idx] = rind_g[idx];
			jCol[idx] = cind_g[idx];
		}
	}
	else {
	// return the values of the jacobian of the constraints

		sparse_jac(tag_g, m, n, 1, x, &nnz_jac, &rind_g, &cind_g, &jacval, options_g);

		for(Index idx=0; idx<nnz_jac; idx++) {
			values[idx] = jacval[idx];
		}
	}
  return true;
}

bool MyADOLC_sparseNLP::eval_h(Index n, const Number* x, bool new_x,
                   Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{

#ifdef SPARSE_HESS

	if (values == NULL) {
	    // return the structure. This is a symmetric matrix, fill the lower left
	    // triangle only.

		for(Index idx=0; idx<nnz_L; idx++) {
			iRow[idx] = rind_L[idx];
			jCol[idx] = cind_L[idx];
		}
	}
	else {
	    // return the values. This is a symmetric matrix, fill the lower left
	    // triangle only

		for(Index idx = 0; idx<n ; idx++)
			x_lam[idx] = x[idx];
		for(Index idx = 0; idx<m ; idx++)
	    	x_lam[n+idx] = lambda[idx];
		x_lam[n+m] = obj_factor;

		sparse_hess(tag_L, n+m+1, 1, x_lam, &nnz_L_total, &rind_L_total, &cind_L_total, &hessval, options_L);

		Index idx = 0;
		for(Index idx_total = 0; idx_total <nnz_L_total; idx_total++) {
			if((rind_L_total[idx_total] < (unsigned int) n) && (cind_L_total[idx_total] < (unsigned int) n)) {
				values[idx] = hessval[idx_total];
				idx++;
			}
		}
	}
	return true;
#endif
	return false;
}

void MyADOLC_sparseNLP::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number obj_value,
			      const IpoptData* ip_data,
			      IpoptCalculatedQuantities* ip_cq)
{
	if (NLP_x_opt != NULL) {
		delete[] NLP_x_opt;
	}
	NLP_x_opt = new Number[n];
	for (Index i = 0; i < n; i++) {
		NLP_x_opt[i] 	= x[i]*NLP_x_sf[i];
	}
// memory deallocation of ADOL-C variables

	delete[] x_lam;
	delete[] OCP_structure;
	delete[] NLP_x_lb;
	delete[] NLP_x_ub;
	delete[] NLP_x_sf;

	delete[] NLP_g_lb;
	delete[] NLP_g_ub;
	delete[] NLP_g_sf;

	delete[] guess;
	delete[] node_str;

#ifdef SPARSE_HESS
	delete[] (rind_L);
	delete[] (cind_L);
	free(cind_L_total);
	free(rind_L_total);
	free(hessval);
	free(rind_g);
	free(cind_g);
	free(jacval);
#endif
}

void MyADOLC_sparseNLP::generate_tapes(Index n, Index m, Index& nnz_jac_g, Index& nnz_h_lag)
{
	Number *xp    = new double[n];
	Number *lamp  = new double[m];
	Number *zl    = new double[m];
	Number *zu    = new double[m];

	adouble *xa   = new adouble[n];
	adouble *g    = new adouble[m];
	adouble *lam  = new adouble[m];
	adouble sig;
	adouble obj_value;

	double dummy;

	x_lam   = new double[n+m+1];

	get_starting_point(n, 1, xp, 0, zl, zu, m, 0, lamp);

	trace_on(tag_f);
    
    for(Index idx=0;idx<n;idx++)
      xa[idx] <<= xp[idx];

    eval_obj(n,xa,obj_value);

    obj_value >>= dummy;

    trace_off();
    
    trace_on(tag_g);

    for(Index idx=0;idx<n;idx++)
    	xa[idx] <<= xp[idx];

    eval_constraints(n,xa,m,g);


    for(Index idx=0;idx<m;idx++)
		g[idx] >>= dummy;
    trace_off();

    trace_on(tag_L);
    
    for(Index idx=0;idx<n;idx++)
    	xa[idx] <<= xp[idx];
    for(Index idx=0;idx<m;idx++)
    	lam[idx] <<= 1.0;
    sig <<= 1.0;

    eval_obj(n,xa,obj_value);

    obj_value = obj_value*sig;
    eval_constraints(n,xa,m,g);
 
    for(Index idx=0;idx<m;idx++)
		obj_value = obj_value + g[idx]*lam[idx];

    obj_value >>= dummy;

    trace_off();

	rind_g = NULL;
	cind_g = NULL;

	options_g[0] = 0;          /* sparsity pattern by index domains (default) */
	options_g[1] = 0;          /*                         safe mode (default) */
	options_g[2] = -1;         /*                     &jacval is not computed */
	options_g[3] = 0;          /*                column compression (default) */


	this->jacval=NULL;


	sparse_jac(tag_g, m, n, 0, xp, &nnz_jac, &rind_g, &cind_g, &jacval, options_g);

	options_g[2] = 0;
	nnz_jac_g = nnz_jac;


	sparse_jac(tag_g, m, n, 1, xp, &nnz_jac, &rind_g, &cind_g, &jacval, options_g);

	double *grad_f = new double[n];
	gradient(tag_f,n,xp,grad_f);

	delete[] grad_f;

/*
	cout<<"eval_grad_f\n";
	for(Index i = 0; i < n; i++) {
		printf("grad_f[%d] = %f\n",i,grad_f[i]);
	}
	for (int i = 0; i < nnz_jac; i += 1)
	{
		printf("sparse_j[%d] = %f\n",i,jacval[i]);
	}
*/
#ifdef SPARSE_HESS
	hessval			= NULL;
	cind_L_total	= NULL;
	rind_L_total	= NULL;
	options_L[0]= 0;
	options_L[1]= 0;
	sparse_hess(tag_L, n+m+1, 0, x_lam, &nnz_L_total, &rind_L_total, &cind_L_total, &hessval, options_L);
	nnz_L = 0;
	for(Index idx_total = 0; idx_total <nnz_L_total ; idx_total++) {
		if((rind_L_total[idx_total] < (unsigned int) n) && (cind_L_total[idx_total] < (unsigned int) n)) {
			nnz_L++;
		}
	}
	nnz_h_lag = nnz_L;
	cout<<"nnz_h_lag = "<<nnz_h_lag<<"\n";
	rind_L 		= new unsigned int [nnz_L];
	cind_L 		= new unsigned int [nnz_L];
	Index idx 	= 0;
	for(Index idx_total = 0; idx_total <nnz_L_total ; idx_total++) {
		if((rind_L_total[idx_total] < (unsigned int) n) && (cind_L_total[idx_total] < (unsigned int) n)) {
			rind_L[idx]		= rind_L_total[idx_total];
			cind_L[idx]		= cind_L_total[idx_total];
			idx++;
		}
	}
#endif

	delete[] lam;
	delete[] g;
	delete[] xa;
	delete[] zu;
	delete[] zl;
	delete[] lamp;
	delete[] xp;

}
