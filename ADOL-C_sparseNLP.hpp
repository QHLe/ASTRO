
#ifndef __MYADOLCNLP_HPP__
#define __MYADOLCNLP_HPP__

#include "IpTNLP.hpp"
#include <adolc/adolc.h>
#include <adolc/adolc_sparse.h>
//#include "OCP.hpp"

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
	template<class T> bool eval_obj(Index n, const T *x, T& obj_value);

  
	/** Template to compute contraints */
	template<class T> bool eval_constraints(Index n, const T *x, Index m, T *g);

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
	void 	setNLP_structure(Index n, Index m, Index* structure);

	Index 	getNLP_n	() 				{ return NLP_n;}
	Index 	getNLP_m	() 				{ return NLP_m;}
	Number* getx_opt	()				{ return NLP_x_opt;}
	void 	setBounds	(	Number* x_lb, Number* x_ub, Number* g_lb, Number* g_ub);
	void 	setguess	(	Number* x_guess)	{guess = x_guess;}
	void	setnodestr	(	Number* str) 		{node_str = str;}

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


//	unsigned int **HP_t;         /* compressed block row storage */
	unsigned int *rind_g;        /* row indices    */
	unsigned int *cind_g;        /* column indices */
	double *jacval;              /* values         */
//	unsigned int *rind_L;        /* row indices    */
//	unsigned int *cind_L;        /* column indices */
//	unsigned int *rind_L_total;  /* row indices    */
//	unsigned int *cind_L_total;  /* column indices */
//	double *hessval;             /* values */
//	double **Hess;
//	double **hesspat;
	int nnz_jac;
	int nnz_L;
	int options_g[4];
	int options_L[4];

  //@}

	Index NLP_n;
	Index NLP_m;
	Number* NLP_x_lb;
	Number* NLP_x_ub;
	Number* NLP_g_lb;
	Number* NLP_g_ub;
	Number* guess;
	Number* node_str;
	Index* OCP_structure;
	Number* NLP_x_opt;
};

#endif
template<class T> T Lagrange_Cost(const T *states, const T *controls, const T *param);
template<class T> T endpoint_cost (	const T* ini_states, const T* fin_states, const T* param, const T& t0, const T& tf, uint phase);
template<class T> void derivatives(const T *states, const T *controls, const T *param, T *states_dot);
template<class T> void events(const T *ini_states, const T *fin_states, const T *param, T *events);
