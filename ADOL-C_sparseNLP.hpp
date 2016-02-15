
#ifndef __MYADOLCNLP_HPP__
#define __MYADOLCNLP_HPP__

enum NLP_SOLVER 	{ma27=0, ma57, ma86, ma97, mumps};
enum OPT_ORDER		{first_order, second_order};
enum APPROX			{Hermite_Simpson=0, trapezoidal};

#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include <adolc/adolc.h>
#include <adolc/adolc_sparse.h>
#ifdef _OPENMP
#include <omp.h>
#include <adolc/adolc_openmp.h>
#endif
#include "SMatrix.hpp"

using namespace Ipopt;

class Phase;
class Guess {
public:
	Guess();
	~Guess();
	SMatrix <double> nodes;
	SMatrix <double> z;
	SMatrix <double> x;
	SMatrix <double> u;
	SMatrix <double> u_full;
	SMatrix <double> parameters;
	SMatrix <double> lam_z;
	SMatrix <double> lam_x;
	SMatrix <double> lam_path;
	SMatrix <double> lam_events;
	SMatrix <double> lam_defects;
};

class Config {
public:
	Config();
	~Config();
	Index max_iter;
	NLP_SOLVER NLP_solver;
	bool warmstart;
	bool with_mgl;
	double NLP_tol;
	APPROX disc_method;
	Index print_level;
	bool H_approximation;

};
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

	Index n_nodes, n_states, n_controls, n_parameters, n_events, n_path_constraints, n_phases, n_linkages;
	SMatrix<double> lb_states, ub_states, lb_controls, ub_controls, lb_parameters, ub_parameters, lb_path, ub_path, lb_events, ub_events, lb_t0, ub_t0, lb_tf, ub_tf;

	Guess 	guess;
	Guess 	results;
	Config 	config;
	double 	*constants;

	template<class T>
	void 	OCP_var_2_NLP_x(T*const* states, T*const* controls, const T* param, const T& t0, const T& tf, T* x, const T* sf);
	template<class T>
	void 	OCP_var_2_NLP_g( T*const* path, T*const* defects, const T* events, T* g);
	template<class T>
	void 	NLP_x_2_OCP_var(const T* x, T** states, T** controls, T* param, T& t0, T& tf);
	template<class T>
	void 	NLP_g_2_OCP_var(const T* g, T** path, T** defects, T* events);

	void 	mem_allocation();
	ApplicationReturnStatus 	initialization();
	ApplicationReturnStatus 	solve();

	void set_endpoint_cost(double (*)	(const  double* ini_states, const  double* fin_states, const double* param, const double& t0, const double& tf, Index phase, const double* constants),
						  adouble (*)	(const adouble* ini_states, const adouble* fin_states, const adouble* param, const adouble& t0, const adouble& tf, Index phase, const double* constants));
	void set_derivatives(void (*)( double *states_dot,  double *path, const  double *states, const  double *controls, const  double *param, const  double &time, Index phase, const double* constants),
						 void (*)(adouble *states_dot, adouble *path, const adouble *states, const adouble *controls, const adouble *param, const adouble &time, Index phase, const double* constants));
	void set_events(void (*)( double *events, const  double *ini_states, const  double *fin_states, const  double *param, const  double &t0, const  double &tf, Index phase, const double* constants),
						 void (*)(adouble *events, const adouble *ini_states, const adouble *fin_states, const adouble *param, const adouble &t0, const adouble &tf, Index phase, const double* constants));
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

	//MyADOLC_sparseNLP();
	MyADOLC_sparseNLP(const MyADOLC_sparseNLP&);
	MyADOLC_sparseNLP& operator=(const MyADOLC_sparseNLP&);

	SmartPtr<IpoptApplication> app;
	ApplicationReturnStatus status;

	double  (*d_e_cost) (const  double* ini_states, const  double* fin_states, const  double* param, const  double& t0, const  double& tf, Index phase, const double* constants);
	adouble (*ad_e_cost)(const adouble* ini_states, const adouble* fin_states, const adouble* param, const adouble& t0, const adouble& tf, Index phase, const double* constants);
	void 	(*d_derv)	( double *states_dot,  double *path, const  double *states, const  double *controls, const  double *param, const  double &time, Index phase, const double* constants);
	void 	(*ad_derv)	(adouble *states_dot, adouble *path, const adouble *states, const adouble *controls, const adouble *param, const adouble &time, Index phase, const double* constants);
	void 	(*d_events)	( double *events, const  double *ini_states, const  double *fin_states, const  double *param, const  double &t0, const  double &tf, Index phase, const double* constants);
	void 	(*ad_events)(adouble *events, const adouble *ini_states, const adouble *fin_states, const adouble *param, const adouble &t0, const adouble &tf, Index phase, const double* constants);

	Index				nlp_n, nlp_m;
	double 				*nlp_sf_x, *nlp_sf_g, *nlp_lb_x, *nlp_ub_x, *nlp_lb_g, *nlp_ub_g, *nlp_guess_x;
	Index 				**state_idx, **control_idx, *parameter_idx, t0_idx, tf_idx,
						**defect_idx, **path_constraint_idx, *event_idx;

	double 				**states, **controls, *parameters, t0, tf,
						**defects, **path_constraints, *events;

	SMatrix<double> 	node_str;

	double 				*x_lam;

	unsigned int 		*rind_g;        /* row indices    */
	unsigned int 		*cind_g;        /* column indices */
	double 				*jacval;              /* values         */
	unsigned int	 	*rind_L;        /* row indices    */
	unsigned int	 	*cind_L;        /* column indices */
	unsigned int	 	*rind_L_total;  /* row indices    */
	unsigned int	 	*cind_L_total;  /* column indices */
	double 				*hessval;             /* values */

	int nnz_jac;
	int nnz_L;//
	int nnz_L_total;
	int options_g[4];
	int options_L[2];

	void 	set_indexes();
	void 	set_sf();
	void 	set_bounces();
	void 	set_guess();
	void 	guess_gen();
	//	double *y0, *yf, **y, **u, *param, tf, t0, **f, **path, **defects, *e, *t, *delta;
	//	double  *NLP_x_lb, *NLP_x_ub, *NLP_x_sf, *NLP_x_guess, *NLP_g_lb, *NLP_g_ub, *NLP_g_sf;
};
/*
template<class T>
void 	MyADOLC_sparseNLP::OCP_var_2_NLP_x(T*const* states, T*const* controls, const T* param, const T& t0, const T& tf, T* x, const T* sf) {
	if (n_nodes != 0 && n_states != 0) { // no nodes or states => static optimization problem
		x[t0_idx]		= t0/sf[t0_idx];
		x[tf_idx]		= tf/sf[tf_idx];
	}
	for (Index i = 0; i < n_parameters; i++) {
		x[parameter_idx[i]]		= param[i]/sf[parameter_idx[i]];
	}

	if (config.disc_method == Hermite_Simpson) {
		for (Index i = 0; i < n_nodes; i += 1) {
			for (Index j = 0; j < n_states; j += 1) {
				x[state_idx[i][j]]		= states[i][j]/sf[state_idx[i][j]];
			}
			for (Index j = 0; j < n_controls; j += 1) {
				x[control_idx[2*i][j]]		= controls[2*i][j]/sf[control_idx[2*i][j]];
				if (i < n_nodes - 1) {
					x[control_idx[2*i+1][j]]		= controls[2*i+1][j]/sf[control_idx[2*i+1][j]];
				}
			}
		}
	}
	else {
		for (Index i = 0; i < n_nodes; i += 1) {
			for (Index j = 0; j < n_states; j += 1) {
				x[state_idx[i][j]]		= states[i][j]/sf[state_idx[i][j]];
			}
			for (Index j = 0; j < n_controls; j += 1) {
				x[control_idx[i][j]]		= controls[i][j]/sf[control_idx[i][j]];
			}
		}
	}
}
*/
template<class T>
void 	MyADOLC_sparseNLP::OCP_var_2_NLP_g(T*const* path, T*const* defects, const T* events, T* g) {

	if (config.disc_method == Hermite_Simpson){
			for (Index i = 0; i < n_nodes; i += 1) {
				for (Index j = 0; j < n_path_constraints; j += 1) {
					g[path_constraint_idx[2*i][j]]	= 	path[2*i][j]/nlp_sf_g[path_constraint_idx[2*i][j]];
					if ( i < n_nodes - 1) {
						g[path_constraint_idx[2*i+1][j]]	= 	path[2*i+1][j]/nlp_sf_g[path_constraint_idx[2*i+1][j]];
					}
				}
				for (Index j = 0; j < n_states; j += 1) {
					if(i < n_nodes - 1) {
						g[defect_idx[i][j]] 		= 	defects[i][j]/nlp_sf_g[defect_idx[i][j]];
					}
				}
			}
		}
	else {
		for (Index i = 0; i < n_nodes; i += 1) {
			for (Index j = 0; j < n_path_constraints; j += 1) {
				g[path_constraint_idx[i][j]]	= 	path[i][j]/nlp_sf_g[path_constraint_idx[i][j]];
			}
			for (Index j = 0; j < n_states; j += 1) {
				if(i < n_nodes - 1) {
					g[defect_idx[i][j]] 		= 	defects[i][j]/nlp_sf_g[defect_idx[i][j]];
				}
			}
		}
	}
	for (Index i = 0; i < n_events; i++) {
		g[event_idx[i]]		= events[i]/nlp_sf_g[event_idx[i]];
	}
}

template<class T>
void 	MyADOLC_sparseNLP::NLP_x_2_OCP_var(const T* x, T** states, T** controls, T* param, T& t0, T& tf) {
	if (n_nodes  && n_states ) { // no nodes or states => static optimization problem
		t0 		= x[t0_idx]*nlp_sf_x[t0_idx];
		tf 		= x[tf_idx]*nlp_sf_x[tf_idx];
	}
	for (Index i = 0; i < n_parameters; i++) {
		param[i]		= x[parameter_idx[i]]*nlp_sf_x[parameter_idx[i]];
	}

	if (config.disc_method == Hermite_Simpson){
		for (Index i = 0; i < n_nodes; i += 1)	{
			for (Index j = 0; j < n_states; j += 1){
				states[i][j] 	= x[state_idx[i][j]]*nlp_sf_x[state_idx[i][j]];
			}
			for (Index j = 0; j < n_controls; j += 1){
				controls[2*i][j] 	= x[control_idx[2*i][j]]*nlp_sf_x[control_idx[2*i][j]];
				if (i < n_nodes - 1) {
					controls[2*i+1][j] 	= x[control_idx[2*i+1][j]]*nlp_sf_x[control_idx[2*i+1][j]];
				}
			}
		}
	}
	else {
		for (Index i = 0; i < n_nodes; i += 1)	{
			for (Index j = 0; j < n_states; j += 1){
				states[i][j] 	= x[state_idx[i][j]]*nlp_sf_x[state_idx[i][j]];
			}
			for (Index j = 0; j < n_controls; j += 1){
				controls[i][j] 	= x[control_idx[i][j]]*nlp_sf_x[control_idx[i][j]];
			}
		}
	}
}


template<class T>
void 	MyADOLC_sparseNLP::NLP_g_2_OCP_var(const T* g, T** path, T** defects, T* events) {
	cout<<"nlp2ocp\n";
	for (Index i = 0; i < n_events; i++) {
		events[i]		= g[event_idx[i]]*nlp_sf_g[event_idx[i]];
	}

	if (config.disc_method == Hermite_Simpson){
		for (Index i = 0; i < n_nodes; i += 1) {
			for (Index j = 0; j < n_path_constraints; j += 1) {
				path[2*i][j]		= g[path_constraint_idx[2*i][j]]*nlp_sf_g[path_constraint_idx[2*i][j]];
				if ( i < n_nodes - 1) {
					path[2*i+1][j]	= g[path_constraint_idx[2*i][j]]*nlp_sf_g[path_constraint_idx[2*i][j]];
				}
			}
			for (Index j = 0; j < n_states; j += 1) {
				if(i < n_nodes - 1) {
					defects[i][j]	= g[defect_idx[i][j]]*nlp_sf_g[defect_idx[i][j]];
				}
			}
		}
	}
	else {
		for (Index i = 0; i < n_nodes; i += 1) {
			for (Index j = 0; j < n_path_constraints; j += 1) {
				path[i][j]	= g[path_constraint_idx[i][j]]*nlp_sf_g[path_constraint_idx[i][j]];
			}
			for (Index j = 0; j < n_states; j += 1) {
				if(i < n_nodes - 1) {
					defects[i][j]	= g[defect_idx[i][j]]*nlp_sf_g[defect_idx[i][j]];
				}
			}
		}
	}
	cout<<"end nlp2ocp\n";
}

template<class T>
SMatrix<T>	lin_interpol(const SMatrix<T> x, const SMatrix<T> y, const SMatrix<T> x_new) {
	SMatrix <T> y_new(x_new.getRowDim(),1);
	for (uint i = 1; i <= x_new.getRowDim(); i++) {
		for(uint j = 1; j < x.getRowDim(); j++) {
/*			if (x_new(i) < x(1) || x(x.getsize()) < x_new(i)) {
				printf("x_new[%d] = %.10e\t	x[%d] = %.10e\n",i,x_new(i),x.getsize(),x(x.getsize()));
				printf("lin_interpot out of range\n");
				exit(1);
			}*/

			if (x(j) <= x_new(i) && x(j+1) >= x_new(i)) {
				y_new(i)	= y(j) + (y(j+1) - y(j))/(x(j+1) - x(j))*(x_new(i) - x(j));
				break;
			}
			if (j == 1 || j == x.getRowDim() - 1) {
				if (x_new(i) < x(1) || x(x.getRowDim()) < x_new(i)) {
	//				cout<<i<<"\t"<<j<<endl;
					printf("x_new(%d) = %.10e\t	x(%d) = %.10e\t	x(%d) = %.10e\n",i,x_new(i),j,x(j),x.getRowDim(),x(x.getRowDim()));
					printf("lin_interpot out of range\n");
					//exit(1);
				}
			}
		}
	}
	return y_new;
}
#endif
