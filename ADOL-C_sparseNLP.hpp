#ifndef __MYADOLCNLP_HPP__
#define __MYADOLCNLP_HPP__

#include "IpTNLP.hpp"
#include <adolc/adolc.h>
#include <adolc/adolc_sparse.h>
#include "Matrix.hpp"

#define tag_f 1
#define tag_g 2
#define tag_L 3
#define HPOFF 30

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
	void 	setNLP_structure(Index n, Index m, Matrix<uint> structure);

	Index 	getNLP_n	() 				{ return NLP_n;}
	Index 	getNLP_m	() 				{ return NLP_m;}
	Matrix<double> get_x_opt	()				{ return NLP_x_opt;}
	Matrix<double> get_lam_opt	()				{ return NLP_lam_opt;}
	void 	setBounds	(	Matrix<double> x_lb, Matrix<double> x_ub, Matrix<double> g_lb, Matrix<double> g_ub);
	void 	setSF		(	Matrix<double> x_sf, Matrix<double> g_sf);
	void 	setguess	(	Matrix<double> x_guess)	{NLP_x_guess = x_guess;}
	void	setnodestr	(	Matrix<double> str) 		{node_str = str;}
	Matrix<double> getnode_str() 		{ return node_str;}
	double 	(*d_e_cost) (const   double* ini_states, const   double* fin_states, const   double* param, const   double& t0, const   double& tf, uint phase);
	adouble (*ad_e_cost)(const  adouble* ini_states, const  adouble* fin_states, const  adouble* param, const  adouble& t0, const  adouble& tf, uint phase);
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
	int nnz_L;
	int nnz_L_total;
	int options_g[4];
	int options_L[2];

  //@}

	Index NLP_n;
	Index NLP_m;
	Matrix<double> NLP_x_lb;
	Matrix<double> NLP_x_ub;
	Matrix<double> NLP_x_sf;
	Matrix<double> NLP_x_guess;
	Matrix<double> NLP_g_lb;
	Matrix<double> NLP_g_ub;
	Matrix<double> NLP_g_sf;
	Matrix<double> NLP_lam_guess;

	Matrix<double> node_str;
	Matrix<uint> OCP_structure;
	Matrix<double> NLP_x_opt;
	Matrix<double> NLP_lam_opt;
};

#endif


