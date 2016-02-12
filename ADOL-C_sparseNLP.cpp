#define SPARSE_HESS
//#define RETAPE

//#define DCOPT_DEBUG

#include <cassert>
#include "ADOL-C_sparseNLP.hpp"
using namespace Ipopt;

#define tag_f 1
#define tag_g 2
#define tag_L 3
#define HPOFF 30

const char* nlp_solver(int enumval) {
	const char* solver_string[] = {"ma27", "ma57", "ma86", "ma97", "mumps"};
	return solver_string[enumval];
}

/* Constructor. */
MyADOLC_sparseNLP::MyADOLC_sparseNLP() {

	x_lam				= NULL;
	cind_g 				= NULL;
	rind_g				= NULL;
	jacval 				= NULL;

	rind_L				= NULL;
	cind_L				= NULL;
	rind_L_total		= NULL;
	cind_L_total		= NULL;
	hessval				= NULL;


}

MyADOLC_sparseNLP::~MyADOLC_sparseNLP() {}


void MyADOLC_sparseNLP::set_endpoint_cost(double (*d_e_cost)	(const  double* ini_states, const  double* fin_states, const  double* param, const  double& t0, const  double& tf, Index phase, const double* constants),
																adouble (*ad_e_cost)	(const adouble* ini_states, const adouble* fin_states, const adouble* param, const adouble& t0, const adouble& tf, Index phase, const double* constants)) {
	this->d_e_cost 	= d_e_cost;
	this->ad_e_cost = ad_e_cost;
}

void MyADOLC_sparseNLP::set_derivatives(void (*d_derv)(  double *states_dot,  double *path, const  double *states, const  double *controls, const  double *param, const  double &time, Index phase, const double* constants),
										void (*ad_derv)(adouble *states_dot, adouble *path, const adouble *states, const adouble *controls, const adouble *param, const adouble &time, Index phase, const double* constants)){
	this->d_derv 	= d_derv;
	this->ad_derv 	= ad_derv;

}
void MyADOLC_sparseNLP::set_events(	void (*d_events)(  double *events, const  double *ini_states, const  double *fin_states, const  double *param, const  double &t0, const  double &tf, Index phase, const double* constants),
									void (*ad_events)(adouble *events, const adouble *ini_states, const adouble *fin_states, const adouble *param, const adouble &t0, const adouble &tf, Index phase, const double* constants)) {
	this->d_events 	= d_events;
	this->ad_events = ad_events;
}

void MyADOLC_sparseNLP::mem_allocation(){
	cout<<"start mem_alloc\n";
	if (n_states == 0 || n_nodes == 0) {
		nlp_n 	= n_parameters;
		nlp_m 	= n_events;
	}
	else {
		if (config.disc_method == Hermite_Simpson){
			nlp_n = (n_states + 2*n_controls)*n_nodes - n_controls + n_parameters + 2;
			nlp_m = ((n_nodes - 1)*n_states + n_events + n_path_constraints*(2*n_nodes-1));
		}
		else {
			nlp_n = (n_states + n_controls)*n_nodes + n_parameters + 2;
			nlp_m = ((n_nodes - 1)*n_states + n_events + n_path_constraints*n_nodes);
		}
	}
	lb_states.resize(n_states,1);
	ub_states.resize(n_states,1);

	lb_controls.resize(n_controls,1);
	ub_controls.resize(n_controls,1);

	lb_parameters.resize(n_parameters,1);
	ub_parameters.resize(n_parameters,1);

	lb_events.resize(n_events,1);
	ub_events.resize(n_events,1);

	lb_path.resize(n_path_constraints,1);
	ub_path.resize(n_path_constraints,1);

	app 				= new IpoptApplication();

	nlp_sf_x 			= new double [nlp_n];
	nlp_sf_g	 		= new double [nlp_m];

	nlp_lb_x		 	= new double[nlp_n];
	nlp_ub_x			= new double[nlp_n];
	nlp_lb_g			= new double[nlp_m];
	nlp_ub_g			= new double[nlp_m];

	state_idx 			= new Index *[n_nodes];
	control_idx			= new Index *[n_nodes];
	parameter_idx		= new Index [n_parameters];

	for (Index i = 0; i < n_nodes; i++) {
		state_idx[i]	= new Index[n_states];
		control_idx[i]	= new Index[n_controls];
	}

	path_constraint_idx	= new Index *[n_nodes];
	defect_idx			= new Index *[n_nodes-1];
	event_idx			= new Index [n_events];

	for (Index i = 0; i < n_nodes; i++) {
		path_constraint_idx[i]	= new Index [n_path_constraints];
		if (i < n_nodes - 1)
			defect_idx[i]		= new Index [n_states];
	}

	cout<<"end mem_alloc\n";
}

ApplicationReturnStatus MyADOLC_sparseNLP::initialization(){
	cout<<"start ini\n";
	this->set_indexes();
	this->set_sf();
	this->set_bounces();

	app->Options()->SetNumericValue("tol", 1e-6);
	app->Options()->SetStringValue("mu_strategy", "adaptive");
//	app->Options()->SetStringValue("output_file", "ipopt.out");
//	app->Options()->SetStringValue("nlp_scaling_method","gradient-based");
	app->Options()->SetStringValue("linear_solver", nlp_solver(config.NLP_solver));//	ma86 & ma57 with memmory leakage
	app->Options()->SetIntegerValue("max_iter", config.max_iter);
	app->Options()->SetIntegerValue("print_level", config.print_level);
	if (config.warmstart) {
		app->Options()->SetStringValue("warm_start_init_point", "yes");
	}

	if (config.H_approximation){
		app->Options()->SetStringValue("hessian_approximation", "limited-memory");
	}
	else {
		app->Options()->SetStringValue("hessian_approximation", "exact");
	}

	status = app->Initialize();
  	if (status != Solve_Succeeded) {
  		printf("\n\n*** Error during initialization!\n");
  	}
  	return status;
	cout<<"end ini\n";
}

void MyADOLC_sparseNLP::set_indexes() {
	cout<<"start indexing\n";
	Index idx_n = 0;

	if (n_nodes == 0 || n_states == 0) {
		for (Index i = 0; i < n_parameters; i++) {
			parameter_idx[i]	= idx_n;
			idx_n++;
		}
		if (idx_n != nlp_n)
			printf("something went wrong in setting indexes (static case)\n");
	}
	else {
		for (Index i = 0; i < n_nodes; i++) {
			for (Index j = 0; j < n_states; j++) {
				state_idx[i][j] 		= idx_n;
				idx_n++;
			}
			for (Index j = 0; j < n_controls; j++) {
				control_idx[i][j] 		= idx_n;
				idx_n++;
			}
		}
		for (Index i = 0; i < n_parameters; i++) {
			parameter_idx[i]			= idx_n;
			idx_n++;
		}
		t0_idx							= idx_n;
		idx_n++;
		tf_idx							= idx_n;
		idx_n++;
		if (idx_n != nlp_n)
			printf("something went wrong in setting indexes\n");
	}

	Index idx_m = 0;
	if (n_nodes == 0 || n_states == 0) {
		for (Index i = 1; i < n_events; i++) {
			event_idx[i]		= idx_m;
			idx_m++;
		}
		if (idx_m != nlp_m)
			printf("something went wrong in setting indexes (g static case)\n");
	}
	else {
		for (Index i = 0; i < n_nodes; i++) {
			for (Index j = 0; j < n_path_constraints; j++) {
				path_constraint_idx[i][j]	= idx_m;
				idx_m++;
			}
			if ( i < n_nodes - 1) {
				for (Index j = 0; j < n_states; j++) {
					defect_idx[i][j] 			= idx_m;
					idx_m++;
				}
			}
		}
		for (Index i = 0; i < n_events; i++) {
			event_idx[i]					= idx_m;
			idx_m++;
		}
		if (idx_m != nlp_m)
			printf("something went wrong in NLP_g_2_OCP_var\n");
	}
	cout<<"end indexing\n";
}

void MyADOLC_sparseNLP::set_sf(){
	cout<<"start set_sf\n";
	for (Index j = 1; j <= n_states; j++) {
		if(lb_states(j) != 0 || ub_states(j) != 0) {
			for (Index i = 0; i < n_nodes; i++)
				nlp_sf_x[state_idx[i][j-1]] 	= max(fabs(lb_states(j)), fabs(ub_states(j)));
		}
		else {
			for (Index i = 0; i < n_nodes; i++)
				nlp_sf_x[state_idx[i][j-1]] 	= 1;
		}
	}
	for (Index j = 1; j <= n_controls; j++) {
		if(lb_controls(j) != 0 || ub_controls(j) != 0) {
			for (Index i = 0; i < n_nodes; i++) {
				nlp_sf_x[control_idx[i][j-1]] = max(fabs(lb_controls(j)), fabs(ub_controls(j)));
			}
		}
		else {
			for (Index i = 0; i < n_nodes; i++) {
				nlp_sf_x[control_idx[i][j-1]] = 1;
			}
		}
	}
	for (Index i = 1; i <= n_parameters; i++) {
		if(lb_parameters(i) != 0 || ub_parameters(i) != 0) {
			nlp_sf_x[parameter_idx[i-1]] = max(fabs(lb_parameters(i)), fabs(ub_parameters(i)));
		}
		else {
			nlp_sf_x[parameter_idx[i-1]] = 1;
		}
	}
	if(lb_t0 != 0 || ub_tf != 0) {
		nlp_sf_x[t0_idx] = max(fabs(lb_t0(1)), fabs(ub_tf(1)));
		nlp_sf_x[tf_idx] = max(fabs(lb_t0(1)), fabs(ub_tf(1)));
	}
	else {
		nlp_sf_x[t0_idx] = 1;
		nlp_sf_x[tf_idx] = 1;
	}

	for (Index j = 1; j <= n_path_constraints; j++) {
		if(lb_path(j) != 0 || ub_path(j) != 0) {
			for (Index i = 0; i < n_nodes; i++)
				nlp_sf_g[path_constraint_idx[i][j-1]] 	= max(fabs(lb_path(j)), fabs(ub_path(j)));
		}
		else {
			for (Index i = 0; i < n_nodes; i++)
				nlp_sf_g[path_constraint_idx[i][j-1]] 	= 1;
		}
	}

	for (Index i = 0; i < n_nodes-1; i++)
		for (Index j = 0; j < n_states; j++)
			nlp_sf_g[defect_idx[i][j]]	= 1;

	for (Index j = 1; j <= n_events; j++) {
		if(lb_events(j) != 0 || ub_events(j) != 0) {
			nlp_sf_g[event_idx[j-1]] = max(fabs(lb_events(j)), fabs(ub_events(j)));
		}
		else {
			nlp_sf_g[event_idx[j-1]] = 1;
		}
	}

	cout<<"end set_sf\n";
}

void MyADOLC_sparseNLP::set_bounces() {

	for (Index i = 0; i < n_parameters; i++) {
		nlp_lb_x[parameter_idx[i]] = lb_parameters(i+1)/nlp_sf_x[parameter_idx[i]];
		nlp_ub_x[parameter_idx[i]] = ub_parameters(i+1)/nlp_sf_x[parameter_idx[i]];
	}

	if (n_nodes != 0 && n_states != 0) {
		nlp_lb_x[t0_idx] 		= lb_t0/nlp_sf_x[t0_idx];
		nlp_lb_x[tf_idx] 		= lb_tf/nlp_sf_x[tf_idx];
		nlp_ub_x[t0_idx] 		= ub_t0/nlp_sf_x[t0_idx];
		nlp_ub_x[tf_idx] 		= ub_tf/nlp_sf_x[tf_idx];
	}

	for (Index i = 0; i < n_nodes; i++) {
		for (Index j = 0; j < n_states; j++) {
			nlp_lb_x[state_idx[i][j]] = lb_states(j+1)/nlp_sf_x[state_idx[i][j]];
			nlp_ub_x[state_idx[i][j]] = ub_states(j+1)/nlp_sf_x[state_idx[i][j]];
		}
		if (config.disc_method == Hermite_Simpson) {
			for (Index j = 0; j < n_controls; j++) {
				nlp_lb_x[control_idx[2*i][j]] = lb_controls(j+1)/nlp_sf_x[control_idx[2*i][j]];
				nlp_ub_x[control_idx[2*i][j]] = ub_controls(j+1)/nlp_sf_x[control_idx[2*i][j]];
				if (i < n_nodes - 1) {
					nlp_lb_x[control_idx[2*i+1][j]] = lb_controls(j+1)/nlp_sf_x[control_idx[2*i+1][j]];
					nlp_ub_x[control_idx[2*i+1][j]] = ub_controls(j+1)/nlp_sf_x[control_idx[2*i+1][j]];
				}
			}
		}
		else {
			for (Index j = 0; j < n_controls; j++) {
				nlp_lb_x[control_idx[i][j]] = lb_controls(j+1)/nlp_sf_x[control_idx[i][j]];
				nlp_ub_x[control_idx[i][j]] = ub_controls(j+1)/nlp_sf_x[control_idx[i][j]];
			}
		}
	}

	for (Index i = 0; i < n_nodes; i++) {
		if (config.disc_method == Hermite_Simpson) {
			for (Index j = 0; j < n_path_constraints; j++) {
				nlp_lb_g[path_constraint_idx[2*i][j]] = lb_path(j+1)/nlp_sf_g[path_constraint_idx[2*i][j]];
				nlp_ub_g[path_constraint_idx[2*i][j]] = ub_path(j+1)/nlp_sf_g[path_constraint_idx[2*i][j]];
				if (i < n_nodes - 1) {
					nlp_lb_g[path_constraint_idx[2*i+1][j]] = lb_path(j+1)/nlp_sf_g[path_constraint_idx[2*i+1][j]];
					nlp_ub_g[path_constraint_idx[2*i+1][j]] = ub_path(j+1)/nlp_sf_g[path_constraint_idx[2*i+1][j]];
				}
			}
		}
		else {
			for (Index j = 0; j < n_path_constraints; j++) {
				nlp_lb_g[path_constraint_idx[i][j]] = lb_path(j+1)/nlp_sf_g[path_constraint_idx[i][j]];
				nlp_ub_g[path_constraint_idx[i][j]] = ub_path(j+1)/nlp_sf_g[path_constraint_idx[i][j]];
			}
		}
		if ( i < n_nodes - 1) {
			for (Index j = 0; j < n_states; j++) {
				nlp_lb_g[defect_idx[i][j]] = 0.0;
				nlp_ub_g[defect_idx[i][j]] = 0.0;
			}
		}
	}

	for (Index i = 0; i < n_events; i++) {
		nlp_lb_g[event_idx[i]] = lb_events(i+1)/nlp_sf_g[event_idx[i]];
		nlp_ub_g[event_idx[i]] = ub_events(i+1)/nlp_sf_g[event_idx[i]];
	}
}

ApplicationReturnStatus MyADOLC_sparseNLP::solve(){

	status = app->OptimizeTNLP(this);

	if (status == Solve_Succeeded) {
	// Retrieve some statistics about the solve
		Index iter_count = app->Statistics()->IterationCount();
		printf("\n\n*** The problem solved in %d iterations!\n", iter_count);

		Number final_obj = app->Statistics()->FinalObjective();
		printf("\n\n*** The final value of the objective function is %e.\n", final_obj);
	}

	return status;
}

bool  MyADOLC_sparseNLP::ad_eval_obj(Index n, const adouble *z, adouble& obj_value) {

#ifdef DCOPT_DEBUG
	printf("ad_eval_obj()\n");
#endif
	if ( n_states == 0 || n_nodes == 0){
		adouble *y0, *yf, **y, **u, t0, tf;
		adouble *param	= new adouble [n_parameters];

		this->NLP_x_2_OCP_var(z, y, u, param, t0, tf);

		obj_value = ad_e_cost (y0, yf, param, t0, tf, 1, constants);

		delete[] param;
	}
	else {

		adouble *y0		= new adouble [n_states];
		adouble *yf		= new adouble [n_states];
		adouble **y 	= new adouble *[n_nodes];
		adouble **u;
		adouble *param	= new adouble [n_parameters];

		adouble tf, t0;

		if (config.disc_method == Hermite_Simpson) {
			u 	= new adouble *[2*n_nodes - 1];
			for (Index i = 0; i < 2*n_nodes - 1; i += 1) {
				u[i]	= new adouble [n_controls];
			}
		}
		else {
			u 	= new adouble *[n_nodes];
			for (Index i = 0; i < n_nodes; i += 1) {
				u[i]	= new adouble [n_controls];
			}
		}

		for (Index i = 0; i < n_nodes; i += 1) {
			y[i]	= new adouble [n_states];
		}
		NLP_x_2_OCP_var(z, y, u, param, t0, tf);

		for (Index i = 0; i < n_states; i += 1)
		{
			y0[i] 		= y[0][i];
			yf[i]		= y[n_nodes - 1][i];
		}

		obj_value = ad_e_cost (y0, yf, param, t0, tf, 1, constants);

		for (Index i = 0; i < n_nodes; i += 1) {
			delete[] y[i];
		}

		if (config.disc_method == Hermite_Simpson) {
			for (Index i = 0; i < 2*n_nodes - 1; i += 1) {
				delete[] u[i];
			}
		}
		else {
			for (Index i = 0; i < n_nodes; i += 1) {
				delete[] u[i];
			}
		}

		delete[] y;
		delete[] u;
		delete[] y0;
		delete[] yf;
		delete[] param;
	}

#ifdef DCOPT_DEBUG
	printf("end of ad_eval_obj()\n");
#endif

	return true;

}

bool  MyADOLC_sparseNLP::eval_obj(Index n, const double *z, double& obj_value) {
  // return the value of the objective function
#ifdef DCOPT_DEBUG
	printf("eval_obj()\n");
#endif

	if ( n_states == 0 || n_nodes == 0){
		double *y0, *yf, **y, **u, t0, tf;
		double *param	= new double [n_parameters];

		this->NLP_x_2_OCP_var(z, y, u, param, t0, tf);

		obj_value = d_e_cost (y0, yf, param, t0, tf, 1, constants);

		delete[] param;
	}
	else {

		double *y0		= new double [n_states];
		double *yf		= new double [n_states];
		double **y 		= new double *[n_nodes];
		double **u;
		double *param	= new double [n_parameters];

		double tf, t0;

		if (config.disc_method == Hermite_Simpson) {
			u 	= new double *[2*n_nodes - 1];
			for (Index i = 0; i < 2*n_nodes - 1; i += 1) {
				u[i]	= new double [n_controls];
			}
		}
		else {
			u 	= new double *[n_nodes];
			for (Index i = 0; i < n_nodes; i += 1) {
				u[i]	= new double [n_controls];
			}
		}

		for (Index i = 0; i < n_nodes; i += 1) {
			y[i]	= new double [n_states];
		}

		NLP_x_2_OCP_var(z, y, u, param, t0, tf);

		for (Index i = 0; i < n_states; i += 1)
		{
			y0[i] 		= y[0][i];
			yf[i]		= y[n_nodes - 1][i];
		}
		this->NLP_x_2_OCP_var(z,y,u,param,t0,tf);

		for (Index i = 0; i < n_states; i += 1)
		{
			y0[i] 		= y[0][i];
			yf[i]		= y[n_nodes - 1][i];
		}

		obj_value = this->d_e_cost (y0, yf, param, t0, tf, 1, constants);

		for (Index i = 0; i < n_nodes; i += 1) {
			delete[] y[i];
		}

		if (config.disc_method == Hermite_Simpson) {
			for (Index i = 0; i < 2*n_nodes - 1; i += 1) {
				delete[] u[i];
			}
		}
		else {
			for (Index i = 0; i < n_nodes; i += 1) {
				delete[] u[i];
			}
		}

		delete[] y;
		delete[] u;
		delete[] y0;
		delete[] yf;
		delete[] param;
	}

#ifdef DCOPT_DEBUG
	printf("end of eval_obj()\n");
#endif
	return true;
}

bool  MyADOLC_sparseNLP::ad_eval_constraints(Index n, const adouble *z, Index m, adouble* g) {

#ifdef DCOPT_DEBUG
	printf("ad_eval_constraints()\n");
#endif
	if(n_states == 0 || n_nodes == 0) {
		adouble *y_start, *y_end, **path, **defects, t0, tf;

		adouble *param		= new adouble [n_parameters];
		adouble *e			= new adouble [n_events];

		ad_events(e, y_start, y_end, param, t0, tf, 1, constants);
		OCP_var_2_NLP_g(path, defects, e, g);

		delete[] param;
		delete[] e;
	}
	else {
		adouble *y_start	= new adouble [n_states];
		adouble *y_end		= new adouble [n_states];
		adouble **y 		= new adouble *[n_nodes];
		adouble **f		 	= new adouble *[n_nodes];
		adouble **u;
		adouble **path;
		adouble *param		= new adouble [n_parameters];
		adouble *e			= new adouble [n_events];
		adouble **defects	= new adouble *[n_nodes-1];

		adouble *y_m		= new adouble [n_states];
		adouble *f_m		= new adouble [n_states];

		adouble *delta	 	= new adouble [n_nodes - 1];
		adouble *t	 		= new adouble [n_nodes];

		adouble tf, t0;

		if (config.disc_method == Hermite_Simpson) {
			u 		= new adouble* [2*n_nodes - 1];
			path	= new adouble* [2*n_nodes - 1];
			for (Index i = 0; i < 2*n_nodes - 1; i += 1) {
				u[i]		= new adouble [n_controls];
				path[i]		= new adouble [n_path_constraints];
			}
		}
		else {
			u 		= new adouble *[n_nodes];
			path	= new adouble *[n_nodes];
			for (Index i = 0; i < n_nodes; i += 1) {
				u[i]		= new adouble [n_controls];
				path[i]		= new adouble [n_path_constraints];
			}
		}

		for (Index i = 0; i < n_nodes; i += 1) {
			y[i]		= new adouble [n_states];
			f[i] 		= new adouble [n_states];
		}

		for (Index i = 0; i < n_nodes - 1; i++) {
			defects[i]	= new adouble [n_states];
		}
		NLP_x_2_OCP_var(z,y,u,param,t0,tf);

		t[0]				= t0;
		for (Index i = 0; i < n_nodes - 1; i++) {
			delta[i]	= (tf-t0)*node_str(i+1);
			t[i+1]		= t[i] + delta[i];
		}
		t[n_nodes - 1]	= tf;

		for (Index i = 0; i < n_states; i += 1) {
				y_start[i] 		= y[0][i];
				y_end[i]		= y[n_nodes - 1][i];
		}


		ad_events(e, y_start, y_end, param, t0, tf, 1, constants);

		if (config.disc_method == Hermite_Simpson) {
			ad_derv(f[0], path[0], y[0], u[0], param, t[0], 1, constants);
			for (Index i = 0; i < n_nodes - 1; i += 1)	{
				ad_derv(f[i+1], path[2*(i+1)], y[i+1], u[2*(i+1)], param, t[i+1], 1, constants);
				for (Index j = 0; j < n_states; j += 1){
					y_m[j] 	= (y[i][j]+y[i+1][j])/2. + delta[i]/8.*(f[i][j]-f[i+1][j]);
				}
				ad_derv(f_m, path[2*i+1], y_m, u[2*i+1], param, (t[i] + t[i+1])/2., 1, constants);
				for (Index j = 0; j < n_states; j++)
					defects[i][j]	= (y[i+1][j] - y[i][j] - delta[i]/6.*(f[i][j]+4.*f_m[j]+f[i+1][j]));
			}
		}
		else if (config.disc_method == trapezoidal) {
			for (Index i = 0; i < n_nodes; i += 1)	{
				ad_derv(f[i], path[i], y[i], u[i], param, t[i], 1, constants);
			}
			for (Index i = 0; i < n_nodes - 1; i++)
				for (Index j = 0; j < n_states; j++)
					defects[i][j] 	= y[i+1][j] - y[i][j] - delta[i]/2.0*(f[i][j] + f[i+1][j]);
		}

		OCP_var_2_NLP_g(path, defects, e, g);

		if (config.disc_method == Hermite_Simpson) {
			for (Index i = 0; i < 2*n_nodes - 1; i += 1) {
				delete[] u[i];
				delete[] path[i];
			}
		}
		else {
			for (Index i = 0; i < n_nodes; i += 1) {
				delete[] u[i];
				delete[] path[i];
			}
		}

		for (Index i = 0; i < n_nodes; i += 1)
		{
			delete[] y[i];
			delete[] f[i];
			if (i < n_nodes - 1) {
				delete[] defects[i];
			}
		}

		delete[] y_m;
		delete[] f_m;

		delete[] y;
		delete[] path;
		delete[] defects;

		delete[] u;
		delete[] f;
		delete[] y_start;
		delete[] y_end;
		delete[] param;
		delete[] e;
		delete[] delta;
		delete[] t;
	//  	delete[] t_m;
	}

#ifdef DCOPT_DEBUG
	printf("end of ad_eval_constraints()\n");
#endif
	return true;
}

bool  MyADOLC_sparseNLP::eval_constraints(Index n, const double *x, Index m, double* g) {

#ifdef DCOPT_DEBUG
	printf("eval_constraints()\n");
#endif
	if(n_states == 0 || n_nodes == 0) {
		double *y0, *yf, **path, **defects, t0, tf;

		double *param		= new double [n_parameters];
		double *e			= new double [n_events];

		d_events(e, y0, yf, param, t0, tf, 1, constants);
		OCP_var_2_NLP_g(path, defects, e, g);

	}
	else {
		double *y0		= new double [n_states];
		double *yf		= new double [n_states];
		double **y 			= new double *[n_nodes];
		double **f		 	= new double *[n_nodes];
		double **u;
		double **path;
		double *param		= new double [n_parameters];
		double *e			= new double [n_events];
		double **defects	= new double *[n_nodes-1];

		double *delta	 	= new double [n_nodes - 1];
		double *t	 		= new double [n_nodes];

		double tf, t0;

		if (config.disc_method == Hermite_Simpson) {
			u 		= new double* [2*n_nodes - 1];
			path	= new double* [2*n_nodes - 1];
			for (Index i = 0; i < 2*n_nodes - 1; i += 1) {
				u[i]		= new double [n_controls];
				path[i]		= new double [n_path_constraints];
			}
		}
		else {
			u 		= new double *[n_nodes];
			path	= new double *[n_nodes];
			for (Index i = 0; i < n_nodes; i += 1) {
				u[i]		= new double [n_controls];
				path[i]		= new double [n_path_constraints];
			}
		}

		for (Index i = 0; i < n_nodes; i += 1) {
			y[i]		= new double [n_states];
			f[i] 		= new double [n_states];
		}

		for (Index i = 0; i < n_nodes - 1; i++) {
			defects[i]	= new double [n_states];
		}
		NLP_x_2_OCP_var(x,y,u,param,t0,tf);

		for (Index i = 0; i < n_states; i += 1)
		{
			y0[i] 		= y[0][i];
			yf[i]		= y[n_nodes - 1][i];
		}

		d_events(e, y0, yf, param, t0, tf, 1, constants);

		t[0]				= t0;
		for (Index i = 0; i < n_nodes - 1; i++) {
			delta[i]	= (tf-t0)*node_str(i+1);
			t[i+1]		= t[i] + delta[i];
		}
		t[n_nodes - 1]	= tf;

		if (config.disc_method == Hermite_Simpson) {
			double* y_m			= new double [n_states];
			double* f_m			= new double [n_states];

			d_derv(f[0], path[0], y[0], u[0], param, t[0], 1, constants);

			for (Index i = 0; i < n_nodes - 1; i += 1)	{
				d_derv(f[i+1], path[2*(i+1)], y[i+1], u[2*(i+1)], param, t[i+1], 1, constants);
				for (Index j = 0; j < n_states; j += 1){
					y_m[j] 	= (y[i][j]+y[i+1][j])/2. + delta[i]/8.*(f[i][j]-f[i+1][j]);
				}
				d_derv(f_m, path[2*i+1], y_m, u[2*i+1], param, (t[i] + t[i+1])/2., 1, constants);
				for (Index j = 0; j < n_states; j++)
					defects[i][j]	= (y[i+1][j] - y[i][j] - delta[i]/6.*(f[i][j]+4.*f_m[j]+f[i+1][j]));
			}

			delete[] y_m;
			delete[] f_m;

		}
		else if (config.disc_method == trapezoidal) {
			for (Index i = 0; i < n_nodes; i += 1)	{
				d_derv(f[i], path[i], y[i], u[i], param, t[i], 1, constants);
			}
			for (Index i = 0; i < n_nodes - 1; i++)
				for (Index j = 0; j < n_states; j++)
					defects[i][j] 	= y[i+1][j] - y[i][j] - delta[i]/2.0*(f[i][j] + f[i+1][j]);
		}

		OCP_var_2_NLP_g(path, defects, e, g);

	}

#ifdef DCOPT_DEBUG
	printf("end of eval_constraints()\n");
#endif
	return true;
}

bool MyADOLC_sparseNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                         Index& nnz_h_lag, IndexStyleEnum& index_style)
{

#ifdef DCOPT_DEBUG
	printf("get_nlp_info()\n");
#endif
	n = nlp_n;
	m = nlp_m;
	generate_tapes(n, m, nnz_jac_g, nnz_h_lag);
  // use the C style indexing (0-based)
	index_style = C_STYLE;

#ifdef DCOPT_DEBUG
	printf("exit get_nlp_info()\n");
#endif
	return true;
}

bool MyADOLC_sparseNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
										Index m, Number* g_l, Number* g_u)
{


	for (Index i = 0; i < nlp_n; i += 1) {
		x_l[i] = nlp_lb_x[i];
		x_u[i] = nlp_ub_x[i];
	}

	for (Index i = 0; i < nlp_m; i += 1) {
		g_l[i] = nlp_lb_g[i];
		g_u[i] = nlp_ub_g[i];
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

	for (Index i = 0; i < n; i += 1) {
		x[i]	= 0;//NLP_x_guess[i];
	}
	return true;
}

bool MyADOLC_sparseNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
	eval_obj(n,x,obj_value);
	return true;
}

bool MyADOLC_sparseNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
	gradient(tag_f + omp_get_thread_num()*10,n,x,grad_f);

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

#ifdef DCOPT_DEBUG
	cout<<"call eval_jac_g\n";
#endif
	if (values == NULL) {
	// return the structure of the jacobian
		for(Index idx=0; idx<nnz_jac; idx++) {
			iRow[idx] = rind_g[idx];
			jCol[idx] = cind_g[idx];
		}
	}
	else {

	#ifdef RETAPE
/*		cout<<"retape\n";
		// return the values of the jacobian of the constraints
		adouble *xa   = new adouble[n];
		adouble *ga    = new adouble[m];
		double dummy;

		trace_on(tag_g);

	    for(Index idx=0;idx<n;idx++)
	    	xa[idx] <<= x[idx];

	    ad_eval_constraints(n,xa,m,ga);

	    for(Index idx=0;idx<m;idx++)
			ga[idx] >>= dummy;
	    trace_off();

	    delete[] xa;
	    delete[] ga;
*/
		free(rind_g);
		free(cind_g);
		free(jacval);

		rind_g = NULL;
		cind_g = NULL;
		jacval = NULL;

		sparse_jac(tag_g, m, n, 0, x, &nnz_jac, &rind_g, &cind_g, &jacval, options_g);
		for(Index idx=0; idx<nnz_jac; idx++) {
			values[idx] = jacval[idx];
		}
	}
#else

		sparse_jac(tag_g + omp_get_thread_num()*10, m, n, 1, x, &nnz_jac, &rind_g, &cind_g, &jacval, options_g);
		for(Index idx=0; idx<nnz_jac; idx++)
			values[idx] = jacval[idx];
	}
#endif

#ifdef DCOPT_DEBUG
	cout<<"return eval_jac_g\n";
#endif
	return true;
}

bool MyADOLC_sparseNLP::eval_h(Index n, const Number* x, bool new_x,
                   Number obj_factor, Index m, const Number* lambda,
                   bool new_lambda, Index nele_hess, Index* iRow,
                   Index* jCol, Number* values)
{

#ifdef DCOPT_DEBUG
	cout<<"call eval_h\n";
#endif
	if (config.H_approximation) {
		#ifdef DCOPT_DEBUG
				cout<<"end eval_h\n";
		#endif
		return false;
	}
	else{
		if (values == NULL) {
			// return the structure. This is a symmetric matrix, fill the lower left
			// triangle only.

			for(Index idx=0; idx<nnz_L; idx++) {
				iRow[idx] = rind_L[idx];
				jCol[idx] = cind_L[idx];
			}

			//cout<<"nele_hess = "<<nele_hess<<"\n";

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
		//			cout<<"rind_L_total["<<rind_L_total[idx_total]<<"]\t";
		//			cout<<"cind_L_total["<<cind_L_total[idx_total]<<"]\t";
		//			cout<<hessval[idx_total]<<endl;
				}
			}
	//		cout<<"nnz_h_lag = "<<idx<<"\n";

		}

		#ifdef DCOPT_DEBUG
			cout<<"end eval_h\n";
		#endif
		return true;
	}
}

void MyADOLC_sparseNLP::finalize_solution(SolverReturn status,
                              Index n, const Number* x, const Number* z_L, const Number* z_U,
                              Index m, const Number* g, const Number* lambda,
                              Number obj_value,
			      const IpoptData* ip_data,
			      IpoptCalculatedQuantities* ip_cq)
{

	/*
	NLP_x_opt.resize(n,1);
	for (Index i = 0; i < n; i++) {
		NLP_x_opt(i+1) 	= x[i]*NLP_x_sf[i];
	}
	NLP_lam_opt.resize(m,1);
	for (Index i = 0; i < m; i++) {
		NLP_lam_opt(i+1) 	= lambda[i]*NLP_g_sf[i];
	}
	*/
//	cout<<"memory deallocation\n";

//	delete[] app;

	delete[] nlp_sf_x;
	delete[] nlp_sf_g;
	delete[] nlp_lb_x;
	delete[] nlp_ub_x;
	delete[] nlp_lb_g;
	delete[] nlp_ub_g;

	delete[] parameter_idx;
	delete[] event_idx;


	for (Index i = 0; i < n_nodes; i++) {
		delete[] state_idx[i];
		delete[] control_idx[i];
	}

	for (Index i = 0; i < n_nodes; i++) {
		delete[] path_constraint_idx[i];
		if (i < n_nodes - 1)
			delete[] defect_idx[i];
	}


	delete[] state_idx;
	delete[] control_idx;
	delete[] path_constraint_idx;
	delete[] defect_idx;

	free(rind_g);
	free(cind_g);
	free(jacval);

	if (!config.H_approximation) {
		delete[] (rind_L);
		delete[] (cind_L);
		free(cind_L_total);
		free(rind_L_total);
		free(hessval);
	}

//	cout<<"memory deallocation done \n";
}

void MyADOLC_sparseNLP::generate_tapes(Index n, Index m, Index& nnz_jac_g, Index& nnz_h_lag)
{
#ifdef DCOPT_DEBUG
	printf("generate_tapes()\n");
#endif
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

	trace_on(tag_f + omp_get_thread_num()*10);
    
    for(Index idx=0;idx<n;idx++)
      xa[idx] <<= xp[idx];

    ad_eval_obj(n,xa,obj_value);

    obj_value >>= dummy;

    trace_off();

    trace_on(tag_g + omp_get_thread_num()*10);

    for(Index idx=0;idx<n;idx++)
    	xa[idx] <<= xp[idx];

    ad_eval_constraints(n,xa,m,g);

    for(Index idx=0;idx<m;idx++)
		g[idx] >>= dummy;
    trace_off();
    
    rind_g = NULL;
	cind_g = NULL;

	options_g[0] = 0;          /* sparsity pattern by index domains (default) */
	options_g[1] = 0;          /*                         safe mode (default) */
	options_g[2] = -1;         /*                     &jacval is not computed */
	options_g[3] = 0;          /*                column compression (default) */


	this->jacval=NULL;

	sparse_jac(tag_g + omp_get_thread_num()*10, m, n, 0, xp, &nnz_jac, &rind_g, &cind_g, &jacval, options_g);

	options_g[2] = 0;
	nnz_jac_g = nnz_jac;


	/*
	sparse_jac(tag_g, m, n, 1, xp, &nnz_jac, &rind_g, &cind_g, &jacval, options_g);

	Number* NLP_constraint_sf = new Number[m];
	for (Index i=0; i<m; i++) {
		NLP_constraint_sf[i] = 0.0;
	}


	for (Index i=0;i<nnz_jac;i++) {
		NLP_constraint_sf[rind_g[i]] += jacval[i]*jacval[i];
	}
	cout<<endl;

	for(Index i = 0; i < m; i++) {
		printf("constraint_sf = %e\n",NLP_constraint_sf[i]);
	}

	delete[] NLP_constraint_sf;
*/
/*
	double *grad_f = new double[n];
	gradient(tag_f + omp_get_thread_num()*10,n,xp,grad_f);

	double enorm_grad_f = 0;
	for (Index i = 0; i < n ; i++)
		enorm_grad_f = enorm_grad_f + grad_f[i]*grad_f[i];
	enorm_grad_f = sqrt(enorm_grad_f);

	if(enorm_grad_f != 0 )
		NLP_obj_sf = enorm_grad_f;
*/
/*
	trace_on(tag_f);

    for(Index idx=0;idx<n;idx++)
      xa[idx] <<= xp[idx];

    ad_eval_obj(n,xa,obj_value);

    obj_value >>= dummy;

    trace_off();
*/
/*

	cout<<"eval_grad_f\n";
	for(Index i = 0; i < n; i++) {
		printf("grad_f[%d] = %f\n",i,grad_f[i]);
	}
	for (int i = 0; i < nnz_jac; i += 1)
	{
		printf("sparse_j[%d] = %f\n",i,jacval[i]);
	}

	delete[] grad_f;
*/
	if (!config.H_approximation) {
		trace_on(tag_L);

		for(Index idx=0;idx<n;idx++)
			xa[idx] <<= xp[idx];
		for(Index idx=0;idx<m;idx++)
			lam[idx] <<= 1.0;

		sig <<= 1.0;

		ad_eval_obj(n,xa,obj_value);

		obj_value = obj_value*sig;
		ad_eval_constraints(n,xa,m,g);

		for(Index idx=0;idx<m;idx++)
			obj_value = obj_value + g[idx]*lam[idx];

		obj_value >>= dummy;

		trace_off();

		hessval			= NULL;
		cind_L_total	= NULL;
		rind_L_total	= NULL;
		options_L[0]= 0;
		options_L[1]= 0;
		sparse_hess(tag_L, n+m+1, 0, x_lam, &nnz_L_total, &rind_L_total, &cind_L_total, &hessval, options_L);
		nnz_L = 0;

		for(Index idx_total = 0; idx_total < nnz_L_total ; idx_total++) {
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
	}


	delete[] lam;
	delete[] g;
	delete[] xa;
	delete[] zu;
	delete[] zl;
	delete[] lamp;
	delete[] xp;

#ifdef DCOPT_DEBUG
	printf("end of generate_tapes()\n");
#endif
}

Guess::Guess() {

}

Guess::~Guess() {

}

Config::Config() {
	max_iter 		= 5000;
	NLP_solver	 	= mumps;
	warmstart 		= false;
	NLP_tol			= 1e-6;
	with_mgl		= false;
	disc_method		= trapezoidal;
	print_level		= 5;
	H_approximation = true;
}

Config::~Config() {

}
