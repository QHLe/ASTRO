

#include <cassert>
#include "ADOL-C_sparseNLP.hpp"

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
}

MyADOLC_sparseNLP::~MyADOLC_sparseNLP()
{}

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

	T *ini_states			= new T [n_states];
	T *fin_states			= new T [n_states];
	T **states 		= new T *[n_nodes];
	T **controls 	= new T *[n_nodes];
	T *param		= new T [n_param];
//	double delta	= (double)1/(n_nodes-1);

	T tf 			= x[n-1];
	T t0			= 0;
	T delta			= (tf-t0)/((double)n_nodes - 1);
	for (Index i = 0; i < n_nodes; i += 1) {
		states[i]	= new T [n_states];
		controls[i]	= new T [n_controls];
	}
	Index idx = 0;

	while (idx < n_nodes*(n_states + n_controls) + n_param) {

		for (Index i = 0; i < n_nodes; i += 1)	{
			for (Index j = 0; j < n_states; j += 1){
				states[i][j] 	= x[idx];
				idx++;
			}
			for (Index j = 0; j < n_controls; j += 1){
				controls[i][j] 	= x[idx];
				idx++;
			}
		}
		for (Index i = 0; i < n_param; i += 1)
		{
			param[i]			= x[idx];
			idx++;
		}
	}
	for (Index i = 0; i < n_states; i += 1)
	{
		ini_states[i] 		= states[0][i];
		fin_states[i]		= states[n_nodes - 1][i];
	}
	obj_value = 0.0;
	for (Index i = 0; i<n_nodes-1; i++)
		obj_value += (Lagrange_Cost(states[i], controls[i], param) + Lagrange_Cost(states[i+1], controls[i+1], param))*delta*0.5;

	obj_value += endpoint_cost (ini_states, fin_states, param, t0, tf, 1);

	for (Index i = 0; i < n_nodes; i += 1) {
		delete[] states[i];
		delete[] controls[i];
	}
	delete[] states;
	delete[] controls;
	delete[] ini_states;
	delete[] fin_states;
	delete[] param;

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

	T *ini_states			= new T [n_states];
	T *fin_states			= new T [n_states];
	T **states 				= new T *[n_nodes];
	T **states_dot		 	= new T *[n_nodes];
	T **controls 			= new T *[n_nodes];
	T **path				= new T *[n_nodes];
	T *param				= new T [n_param];
	T *e					= new T [n_events];
//	double delta			= (double)1/(n_nodes-1);
	T tf 					= x[n - 1];
	T delta 				= x[n - 1]/((double)n_nodes-1);

	for (Index i = 0; i < n_nodes; i += 1) {
		states[i]			= new T [n_states];
		path[i]				= new T [n_path];
		states_dot[i] 		= new T [n_states];
		controls[i]			= new T [n_controls];
	}

	Index idx = 0;
	while (idx < n_nodes*(n_states + n_controls) + n_param) {
		for (Index i = 0; i < n_nodes; i += 1)	{
			for (Index j = 0; j < n_states; j += 1){
				states[i][j] 	= x[idx];
				idx++;
			}
			for (Index j = 0; j < n_controls; j += 1){
				controls[i][j] 	= x[idx];
				idx++;
			}
			derivatives(states[i], controls[i], param, states_dot[i]);
			//euler(states[i],states_dot[i],delta,states_approx[i]);
			//trapezoidal(states[i],states_dot[i],states_dot[i+1],delta,states_approx[i]);
		}
		for (Index i = 0; i < n_param; i += 1)
		{
			param[i]			= x[idx];
			idx++;
		}
		//tf 	= x[idx];
		//idx++;

	}

	for (Index i = 0; i < n_states; i += 1)
	{
		ini_states[i] 		= states[0][i];
		fin_states[i]		= states[n_nodes - 1][i];
	}
	events(ini_states,fin_states,param,e);

	idx = 0;
	while( idx < m) {
		for (Index i = 0; i < n_nodes; i += 1) {


			for (Index j = 0; j < n_path; j += 1) {
				cout<<"path idx = "<<idx<<"\n";
				g[idx]	= 	0.0;	//need to be implemented
				idx++;
			}
			for (Index j = 0; j < n_states; j += 1) {
				if(i < n_nodes - 1) {
					g[idx]	= 	states[i+1][j] - states[i][j] - delta/2.0*(states_dot[i][j] + states_dot[i+1][j]);
					idx++;
				}
			}
		}
		for (Index i = 0; i < n_events; i += 1)
		{
			g[idx]	= 	e[i];
			idx++;
		}
	}

	for (Index i = 0; i < n_nodes; i += 1)
	{
		delete[] states[i];
		delete[] path[i];
		delete[] states_dot[i];
		delete[] controls[i];
	}
	delete[] states;
	delete[] path;
	delete[] controls;
	delete[] states_dot;
   	delete[] ini_states;
   	delete[] fin_states;
   	delete[] param;
   	delete[] e;
#ifdef DISPLAY_ON
   	cout<<"end of eval constraints\n";
#endif
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
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the others if
  // you wish.
	assert(init_x == true);
	assert(init_z == false);
	assert(init_lambda == false);
	cout<<"get starting point\n";
  // set all y's to the perfect match with y_d
/*		x[0] 	= 0.0;
		x[1] 	= 1.0;
		x[2]	= -2.0;

		x[3]	= 0.5;
		x[4]	= 0.0;
		x[5]	= -2.0;

		x[6]	= 0.0;
		x[7]	= -1.0;
		x[8]	= -2.0;
*/
	for (Index i = 0; i < n; i += 1) {
		x[i]	= 1.0;
	}

//	x[n-1]		= 1.0;
/*	for(Index i = 0;i<n ; i++) {
		printf("x[%d] = %f\n",i,x[i]);
	}*/
	cout<<"end of getting starting point\n";
	return true;
}

template<class T> T Lagrange_Cost(const T *states, const T *controls, const T *param) {
	return 0;
//	return controls[0]*controls[0]*0.5;
}
template<class T> T endpoint_cost (	const T* ini_states, const T* fin_states, const T* param, const T& t0, const T& tf, uint phase) {

	return tf;
}

template<class T> void derivatives(const T *states, const T *controls, const T *param, T *states_dot) {
	
	states_dot[0] 	= states[1];
	states_dot[1]	= controls[0];

}

template<class T> void events(const T *ini_states, const T *fin_states, const T *param, T *events){ 
	events [0]	= ini_states[0];
	events [1]	= ini_states[1];
	events [2]	= fin_states[0];
	events [3]	= fin_states[1];
	
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

	cout<<"eval h \n";
	return false;
/*  if (values == NULL) {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    for(Index idx=0; idx<nnz_L; idx++)
      {
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
    for(Index idx_total = 0; idx_total <nnz_L_total ; idx_total++)
      {
      //printf("rind_L_total[%d] = %d\trind_L_total[%d] = %d\thessval = %.4e\n",idx_total,rind_L_total[idx_total],idx_total,rind_L_total[idx_total],hessval[idx_total]);
	if((rind_L_total[idx_total] < (unsigned int) n) && (cind_L_total[idx_total] < (unsigned int) n))
	  {
	  	//printf("rind_L_total[%d] = %d\trind_L_total[%d] = %d\thessval = %.4e\n",idx_total,rind_L_total[idx_total],idx_total,rind_L_total[idx_total],hessval[idx_total]);
	    values[idx] = hessval[idx_total];
	    idx++;
	  }
      }
  }

	cout<<"end of eval h \n";
	return true;
	*/
}

void MyADOLC_sparseNLP::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number obj_value,
			      const IpoptData* ip_data,
			      IpoptCalculatedQuantities* ip_cq)
{
	cout<<"finalizing solution\n";
// memory deallocation of ADOL-C variables
/*//	Index n_phases		= OCP_structure[0];
	Index n_nodes 		= OCP_structure[1];
	Index n_states 		= OCP_structure[2];
	Index n_controls 	= OCP_structure[3];
	Index n_param 		= OCP_structure[4];
	Index n_events		= OCP_structure[5];
	Index n_path		= OCP_structure[6];
//	Index n_linkages	= OCP_structure[7];

	double *ini_states			= new double [n_states];
	double *fin_states			= new double [n_states];
	double **states 			= new double *[n_nodes];
	double **states_dot		 	= new double *[n_nodes];
	double **controls 			= new double *[n_nodes];
	double **path				= new double *[n_nodes];
	double *param				= new double [n_param];
	double *e					= new double [n_events];
//	double delta			= (double)1/(n_nodes-1);
	double tf 					= x[n - 1];
	double delta 				= x[n - 1]/((double)n_nodes-1);

	for (Index i = 0; i < n_nodes; i += 1) {
		states[i]			= new double [n_states];
		path[i]				= new double [n_path];
		states_dot[i] 		= new double [n_states];
		controls[i]			= new double [n_controls];
	}

	Index idx = 0;
	while (idx < n_nodes*(n_states + n_controls) + n_param) {
		for (Index i = 0; i < n_nodes; i += 1)	{
			for (Index j = 0; j < n_states; j += 1){
				states[i][j] 	= x[idx];
				idx++;
			}
			for (Index j = 0; j < n_controls; j += 1){
				controls[i][j] 	= x[idx];
				printf("controls[%d] = %f\n",i,controls[i][j]);
				idx++;
			}
			derivatives(states[i], controls[i], param, states_dot[i]);
			//euler(states[i],states_dot[i],delta,states_approx[i]);
			//trapezoidal(states[i],states_dot[i],states_dot[i+1],delta,states_approx[i]);
		}
		for (Index i = 0; i < n_param; i += 1)
		{
			param[i]			= x[idx];
			idx++;
		}
		//tf 	= x[idx];
		//idx++;

	}

	for (Index i = 0; i < n_states; i += 1)
	{
		ini_states[i] 		= states[0][i];
		fin_states[i]		= states[n_nodes - 1][i];
	}
	events(ini_states,fin_states,param,e);

	delete[] ini_states;
	delete[] fin_states;

	delete[] param;
	delete[] e;

	for (Index i = 0; i < n_nodes; i += 1) {
		delete[] states[i];
		delete[] path[i];
		delete[] states_dot[i];
		delete[] controls[i];
	}

	delete[] states;
	delete[] states_dot;
	delete[] controls;
	delete[] path;
*/
	delete[] x_lam;
	delete[] OCP_structure;
	delete[] NLP_x_lb;
	delete[] NLP_x_ub;
	delete[] NLP_g_lb;
	delete[] NLP_g_ub;
	//  free(rind_L);
	//  free(cind_L);
	//  free(hessval);
	free(rind_g);
	free(cind_g);
	free(jacval);



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
//	double *jacval;

//	unsigned int i,j,k,l,ii;

	x_lam   = new double[n+m+1];
//	unsigned int  **hesspat=NULL; // compressed row storage
//	int           options_h=0; // options for the hessian patterns
//	double        *x;
//	int retv_h = -1; // return value

//	hesspat = new unsigned int* [(n+m+1)];
//	x = new double[(n+m+1)];


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

    obj_value *= sig;
    eval_constraints(n,xa,m,g);
 
    for(Index idx=0;idx<m;idx++)
		obj_value += g[idx]*lam[idx];

    obj_value >>= dummy;

    trace_off();

	rind_g = NULL;
	cind_g = NULL;

	options_g[0] = 0;          /* sparsity pattern by index domains (default) */
	options_g[1] = 0;          /*                         safe mode (default) */
	options_g[2] = -1;         /*                     &jacval is not computed */
	options_g[3] = 0;          /*                column compression (default) */


	this->jacval=NULL;
//	this->hessval=NULL;

	sparse_jac(tag_g, m, n, 0, xp, &nnz_jac, &rind_g, &cind_g, &jacval, options_g);

	options_g[2] = 0;
	nnz_jac_g = nnz_jac;
/*
	sparse_jac(tag_g, m, n, 1, xp, &nnz_jac, &rind_g, &cind_g, &jacval, options_g);

	double *grad_f = new double[n];
	gradient(tag_f,n,xp,grad_f);

	cout<<"eval_grad_f\n";
	for(Index i = 0; i < n; i++) {
		printf("grad_f[%d] = %f\n",i,grad_f[i]);
	}
	for (int i = 0; i < nnz_jac; i += 1)
	{
		printf("sparse_j[%d] = %f\n",i,jacval[i]);
	}
*/

/*	hess_pat(tag_L,n+m+1,x,hesspat,options_h);

	for (Index i=0;i<n;i++) {
		for (Index j=1;j<=hesspat[i][0];j++){
			if(hesspat[i][j]<=i){ // hess_pat returns all the non-zeros, not just the lower-left corner ones
				nnz_L++;
			}
		}
	}
	rind_L = new unsigned int [nnz_L];
	cind_L = new unsigned int [nnz_L];
	hessval = new double [nnz_L];
	Index inx = 0;
	for (Index i=0;i<n;i++) {
		for (Index j=1;j<=hesspat[i][0];j++){
			if(hesspat[i][j]<=i){ // hess_pat returns all the non-zeros, not just the lower-left corner ones
				rind_L[inx] 	= i;
				cind_L[inx]		= j;
			}
		}
	}
*/
	delete[] lam;
	delete[] g;
	delete[] xa;
	delete[] zu;
	delete[] zl;
	delete[] lamp;
	delete[] xp;

}
