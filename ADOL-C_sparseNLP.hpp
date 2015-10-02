
#ifndef __MYADOLCNLP_HPP__
#define __MYADOLCNLP_HPP__

enum NLP_SOLVER 	{ma27=0, ma57, ma86, ma97, mumps};
enum OPT_ORDER		{first_order, second_order};
enum APPROX			{Hermite_Simpson=0, trapezoidal};

#include "IpTNLP.hpp"
#include <adolc/adolc.h>
#include <adolc/adolc_sparse.h>
#ifdef _OPENMP
#include <omp.h>
#include <adolc/adolc_openmp.h>
#endif
#include "SMatrix.hpp"

using namespace Ipopt;

class OCP;
class MyADOLC_sparseNLP : public TNLP
{
public:
	/** default constructor */
	MyADOLC_sparseNLP();

	/** default destructor */
	virtual ~MyADOLC_sparseNLP();

	/**@name Overloaded from TNLP */
	//@{
	/** Method to return some info about the nlp */
	virtual bool get_nlp_info(	Index& n, Index& m, Index& nnz_jac_g,
								Index& nnz_h_lag, IndexStyleEnum& index_style);

	/** Method to return the bounds for my problem */
	virtual bool get_bounds_info(	Index n, Number* x_l, Number* x_u,
									Index m, Number* g_l, Number* g_u);

	/** Method to return the starting point for the algorithm */
	virtual bool get_starting_point(Index n, bool init_x, Number* x,
		  	  	  	  	  	  		bool init_z, Number* z_L, Number* z_U,
		  	  	  	  	  	  		Index m, bool init_lambda,
		  	  	  	  	  	  		Number* lambda);

	/** Template to return the objective value */
	bool eval_obj(Index n, const double *x, double& obj_value);
	bool ad_eval_obj(Index n, const adouble *x, adouble& obj_value);

  
	/** Template to compute contraints */
	bool eval_constraints(Index n, const double *x, Index m, double *g);
	bool ad_eval_constraints(Index n, const adouble *x, Index m, adouble *g);
	/** Original method from Ipopt to return the objective value */
	/** remains unchanged */
	virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

	/** Original method from Ipopt to return the gradient of the objective */
	/** remains unchanged */
	virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

	/**  Original method from Ipopt to return the constraint residuals */
	/** remains unchanged */
	virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

	/** Original method from Ipopt to return:
	*   1) The structure of the jacobian (if "values" is NULL)
	*   2) The values of the jacobian (if "values" is not NULL)
	*/
	/** remains unchanged */
	virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
							Index m, Index nele_jac, Index* iRow, Index *jCol,
							Number* values);

  /** Original method from Ipopt to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  /** remains unchanged */
	virtual bool eval_h(Index n, const Number* x, bool new_x,
						Number obj_factor, Index m, const Number* lambda,
						bool new_lambda, Index nele_hess, Index* iRow,
						Index* jCol, Number* values);

  //@}

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
	virtual void finalize_solution(	SolverReturn status,
									Index n, const Number* x, const Number* z_L, const Number* z_U,
									Index m, const Number* g, const Number* lambda,
									Number obj_value,
									const IpoptData* ip_data,
									IpoptCalculatedQuantities* ip_cq);
  //@}

//***************    start ADOL-C part ***********************************

	/** Method to generate the required tapes */
	virtual void generate_tapes(Index n, Index m, Index& nnz_jac_g, Index& nnz_h_lag);

//	template<class T>	void euler(const T *states_0, const T *states_dot, const double delta, T *states_1);
//	template<class T>	void trapezoidal(const T *states_0, const T *states_dot_0, const T *states_dot_1, const double delta, T *states_1);

	//***************    end   ADOL-C part ***********************************
	void 	setNLP_structure(Index n, Index m, SMatrix<uint> structure, APPROX method);
	template<class T>
	void 	OCP_var_2_NLP_x(T*const* states, T*const* controls, const T* param, const T& t0, const T& tf, T* x, const T* sf);
	template<class T>
	void 	OCP_var_2_NLP_g( T*const* path, T*const* defects, const T* events, T* g, const T* g_sf);
	template<class T>
	void 	NLP_x_2_OCP_var(const T* x, const T* sf, T** states, T** controls, T* param, T& t0, T& tf);
	template<class T>
	void 	NLP_g_2_OCP_var(const T* g, const T* sf, T** path, T** defects, T* events);
	Index 	getNLP_n	() 				{ return NLP_n;}
	Index 	getNLP_m	() 				{ return NLP_m;}
	SMatrix<double> get_x_opt	()				{ return NLP_x_opt;}
	SMatrix<double> get_lam_opt	()				{ return NLP_lam_opt;}
	void 	setBounds	( double* x_lb, double* x_ub, double* g_lb, double* g_ub);
	void 	setSF		( double* x_sf, double* g_sf);
	void 	setguess	( double* x_guess)	{NLP_x_guess = x_guess;}
	void	setnodestr	(	SMatrix<double> str) 		{node_str = str;}
	SMatrix<double> getnode_str() 		{ return node_str;}
	double (*d_e_cost) 	(const  double* ini_states, const  double* fin_states, const  double* param, const  double& t0, const  double& tf, uint phase);
	adouble (*ad_e_cost)(const adouble* ini_states, const adouble* fin_states, const adouble* param, const adouble& t0, const adouble& tf, uint phase);
	double  (*d_l_cost)	(const  double *states, const  double *controls, const  double *param, const  double &time,	uint phase);
	adouble	(*ad_l_cost)(const adouble *states, const adouble *controls, const adouble *param, const adouble &time,	uint phase);
	void 	(*d_derv)	( double *states_dot,  double *path, const  double *states, const  double *controls, const  double *param, const  double &time, uint phase);
	void 	(*ad_derv)	(adouble *states_dot, adouble *path, const adouble *states, const adouble *controls, const adouble *param, const adouble &time, uint phase);
	void 	(*d_events)	( double *events, const  double *ini_states, const  double *fin_states, const  double *param, const  double &t0, const  double &tf, uint phase);
	void 	(*ad_events)(adouble *events, const adouble *ini_states, const adouble *fin_states, const adouble *param, const adouble &t0, const adouble &tf, uint phase);

private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually 
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *  
   */
  //@{
  //  MyADOLC_sparseNLP();
	MyADOLC_sparseNLP(const MyADOLC_sparseNLP&);
	MyADOLC_sparseNLP& operator=(const MyADOLC_sparseNLP&);
  //@}

	//@{
	double *x_lam;

	unsigned int *rind_g;        /* row indices    */
	unsigned int *cind_g;        /* column indices */
	double *jacval;              /* values         */
	unsigned int *rind_L;        /* row indices    */
	unsigned int *cind_L;        /* column indices */
	unsigned int *rind_L_total;  /* row indices    */
	unsigned int *cind_L_total;  /* column indices */
	double *hessval;             /* values */

	int nnz_jac;
	int nnz_L;//
	int nnz_L_total;
	int options_g[4];
	int options_L[2];

  //@}

	Index NLP_n, NLP_m;
	Index n_nodes, n_states, n_controls, n_param, n_events, n_path, n_phases, n_linkages;
	APPROX disc_method;
	SMatrix<double> NLP_lam_guess;
	double *y0, *yf, **y, **u, *param, tf, t0, **f, **path, **defects, *e, *t, *delta;
	double  *NLP_x_lb, *NLP_x_ub, *NLP_x_sf, *NLP_x_guess, *NLP_g_lb, *NLP_g_ub, *NLP_g_sf;


	SMatrix<double> node_str;
	SMatrix<uint> OCP_structure;
	SMatrix<double> NLP_x_opt, NLP_lam_opt;
};

template<class T>
void 	MyADOLC_sparseNLP::OCP_var_2_NLP_x(T*const* states, T*const* controls, const T* param, const T& t0, const T& tf, T* x, const T* sf) {

	Index idx = 0;

	if (disc_method == Hermite_Simpson) {
		for (Index i = 0; i < n_nodes; i += 1) {
			for (Index j = 0; j < n_states; j += 1) {
				x[idx]		= states[i][j]/sf[idx];
				idx++;
			}
			for (Index j = 0; j < n_controls; j += 1) {
				x[idx]		= controls[2*i][j]/sf[idx];
				idx++;
				if (i < n_nodes - 1) {
					x[idx]		= controls[2*i+1][j]/sf[idx];
					idx++;
				}
			}
		}

		for (Index i = 0; i < n_param; i += 1)	{
			x[idx]		= param[i]/sf[idx];
			idx++;
		}

		x[idx]			= t0/sf[idx];
		idx++;
		x[idx]			= tf/sf[idx];
		idx++;

		if (idx != NLP_n)
			printf("something went wrong in OCP_var_2_NLP_x\n");
	}
	else {
		for (Index i = 0; i < n_nodes; i += 1) {
			for (Index j = 0; j < n_states; j += 1) {
				x[idx]		= states[i][j]/sf[idx];
				idx++;
			}
			for (Index j = 0; j < n_controls; j += 1) {
				x[idx]		= controls[i][j]/sf[idx];
				idx++;
			}
		}
		for (Index i = 0; i < n_param; i += 1)	{
			x[idx]		= param[i]/sf[idx];
			idx++;
		}

		x[idx]			= t0/sf[idx];
		idx++;

		x[idx]			= tf/sf[idx];
		idx++;

		if (idx != NLP_n)
			printf("something went wrong in OCP_var_2_NLP_x\n");

	}
}

template<class T>
void 	MyADOLC_sparseNLP::OCP_var_2_NLP_g(T*const* path, T*const* defects, const T* events, T* g, const T* g_sf) {

	Index idx_m = 0;
	if (disc_method == Hermite_Simpson){
		for (Index i = 0; i < n_nodes; i += 1) {
			for (Index j = 0; j < n_path; j += 1) {
				g[idx_m]	= 	path[2*i][j]/g_sf[idx_m];
				idx_m++;
				if ( i < n_nodes - 1) {
					g[idx_m]	= 	path[2*i+1][j]/g_sf[idx_m];
					idx_m++;
				}
			}
			for (Index j = 0; j < n_states; j += 1) {
				if(i < n_nodes - 1) {
					g[idx_m] 		= 	defects[i][j]/g_sf[idx_m];
					idx_m++;
				}
			}
		}
		for (Index i = 0; i < n_events; i += 1)
		{
			g[idx_m]	= 	events[i]/g_sf[idx_m];
			idx_m++;
		}

		if (idx_m !=NLP_m)
			printf("something went wrong in OCP_var_2_NLP_g\n");
	}
	else {
		for (Index i = 0; i < n_nodes; i += 1) {
			for (Index j = 0; j < n_path; j += 1) {
				g[idx_m]	= 	path[i][j]/g_sf[idx_m];
				idx_m++;
			}
			for (Index j = 0; j < n_states; j += 1) {
				if(i < n_nodes - 1) {
					g[idx_m] 		= 	defects[i][j]/g_sf[idx_m];
					idx_m++;
				}
			}
		}
		for (Index i = 0; i < n_events; i += 1)
		{
			g[idx_m]	= 	events[i]/g_sf[idx_m];
			idx_m++;
		}

		if (idx_m !=NLP_m)
			printf("something went wrong in OCP_var_2_NLP_g\n");
	}
}

template<class T>
void 	MyADOLC_sparseNLP::NLP_x_2_OCP_var(const T* x, const T* sf, T** states, T** controls, T* param, T& t0, T& tf) {

	Index idx_n = 0;

	if (disc_method == Hermite_Simpson){
		for (Index i = 0; i < n_nodes; i += 1)	{
			for (Index j = 0; j < n_states; j += 1){
				states[i][j] 	= x[idx_n]*sf[idx_n];
				idx_n++;
			}
			for (Index j = 0; j < n_controls; j += 1){
				controls[2*i][j] 	= x[idx_n]*sf[idx_n];
				idx_n++;
				if (i < n_nodes - 1) {
					controls[2*i+1][j] 	= x[idx_n]*sf[idx_n];
					idx_n++;
				}
			}
		}
		for (Index i = 0; i < n_param; i += 1)
		{
			param[i]			= x[idx_n]*sf[idx_n];
			idx_n++;
		}

		t0 					= x[idx_n]*sf[idx_n];
		idx_n++;
		tf					= x[idx_n]*sf[idx_n];
		idx_n++;

		if (idx_n != NLP_n)
			printf("something went wrong in NLP_x_2_OCP_var\n");
	}
	else {
		for (Index i = 0; i < n_nodes; i += 1)	{
			for (Index j = 0; j < n_states; j += 1){
				states[i][j] 	= x[idx_n]*sf[idx_n];
				idx_n++;
			}
			for (Index j = 0; j < n_controls; j += 1){
				controls[i][j] 	= x[idx_n]*sf[idx_n];
				idx_n++;
			}
		}
		for (Index i = 0; i < n_param; i += 1)
		{
			param[i]			= x[idx_n]*sf[idx_n];
			idx_n++;
		}
		t0 						= x[idx_n]*sf[idx_n];
		idx_n++;
		tf		 				= x[idx_n]*sf[idx_n];
		idx_n++;

		if (idx_n != NLP_n)
			printf("something went wrong in NLP_x_2_OCP_var\n");
	}
}

template<class T>
void 	MyADOLC_sparseNLP::NLP_g_2_OCP_var(const T* g, const T* sf, T** path, T** defects, T* events) {

	Index idx_m = 0;
	if (disc_method == Hermite_Simpson){
		for (Index i = 0; i < n_nodes; i += 1) {
			for (Index j = 0; j < n_path; j += 1) {
				path[2*i][j]	= g[idx_m]*sf[idx_m];
				idx_m++;
				if ( i < n_nodes - 1) {
					path[2*i+1][j]	= g[idx_m]*sf[idx_m];
					idx_m++;
				}
			}
			for (Index j = 0; j < n_states; j += 1) {
				if(i < n_nodes - 1) {
					defects[i][j]	= g[idx_m]*sf[idx_m];
					idx_m++;
				}
			}
		}

		for (Index i = 0; i < n_events; i += 1)	{
			events[i]		= g[idx_m]*sf[idx_m];
			idx_m++;
		}
		if (idx_m != NLP_m)
			printf("something went wrong in NLP_g_2_OCP_var\n");
	}
	else {
		for (Index i = 0; i < n_nodes; i += 1) {
			for (Index j = 0; j < n_path; j += 1) {
				path[i][j]	= g[idx_m]*sf[idx_m];
				idx_m++;
			}
			for (Index j = 0; j < n_states; j += 1) {
				if(i < n_nodes - 1) {
					defects[i][j]	= g[idx_m]*sf[idx_m];
					idx_m++;
				}
			}
		}

		for (Index i = 0; i < n_events; i += 1)	{
			events[i]		= g[idx_m]*sf[idx_m];
			idx_m++;
		}
		if (idx_m != NLP_m)
			printf("something went wrong in NLP_g_2_OCP_var\n");
	}
}

#endif
