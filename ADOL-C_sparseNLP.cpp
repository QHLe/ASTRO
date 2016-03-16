#define SPARSE_HESS
//#define RETAPE

//#define DCOPT_DEBUG

#include <cassert>
#include "ADOL-C_sparseNLP.hpp"

double step = std::numeric_limits<double>::epsilon();
time_t rand_time;

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

	cind_g 				= NULL;
	rind_g				= NULL;
	jacval 				= NULL;

	rind_L				= NULL;
	cind_L				= NULL;
	hessval				= NULL;

	nnz_jac				= 0;
	nnz_L				= 0;


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
#ifdef DCOPT_DEBUG
	cout<<"start mem_alloc\n";
#endif
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

	nlp_guess_x 		= new double [nlp_n];
	nlp_sf_x 			= new double [nlp_n];
	nlp_sf_g	 		= new double [nlp_m];

	nlp_lb_x		 	= new double[nlp_n];
	nlp_ub_x			= new double[nlp_n];
	nlp_lb_g			= new double[nlp_m];
	nlp_ub_g			= new double[nlp_m];

	if (n_parameters) {
		parameter_idx		= new Index [n_parameters];
	}

	if (n_events) {
		event_idx			= new Index [n_events];
	}

	if (n_nodes && n_states) {
		state_idx 			= new Index *[n_nodes];

		defect_idx			= new Index *[n_nodes - 1];

		if (n_controls && config.disc_method == Hermite_Simpson) {
			control_idx			= new Index *[2*n_nodes-1];
			path_constraint_idx	= new Index *[2*n_nodes-1];
		}
		else {
			control_idx			= new Index *[n_nodes];
			path_constraint_idx	= new Index *[n_nodes];
		}
	}

	for (Index i = 0; i < n_nodes; i++) {
		if (n_states) {
			state_idx[i]	= new Index[n_states];
			if (i < n_nodes - 1) {
				defect_idx[i]		= new Index [n_states];
			}
		}
		if (n_controls && config.disc_method == Hermite_Simpson) {
			control_idx[2*i]				= new Index[n_controls];
			if (n_path_constraints) {
				path_constraint_idx[2*i]	= new Index [n_path_constraints];
			}
			else {
				path_constraint_idx[2*i]	= NULL;
			}
			if (i < n_nodes - 1) {
				control_idx[2*i+1]			= new Index[n_controls];
				if (n_path_constraints) {
					path_constraint_idx[2*i+1]	= new Index [n_path_constraints];
				}
				else {
					path_constraint_idx[2*i+1]	= NULL;
				}
			}
		}
		else if (n_controls) {
			control_idx[i]	= new Index[n_controls];
			if (n_path_constraints) {
				path_constraint_idx[i]	= new Index [n_path_constraints];
			}
			else {
				path_constraint_idx[i]	= NULL;
			}
		}

	}
#ifdef DCOPT_DEBUG
	cout<<"end mem_alloc\n";
#endif
}

ApplicationReturnStatus MyADOLC_sparseNLP::initialization(SmartPtr<IpoptApplication> app){
#ifdef DCOPT_DEBUG
	cout<<"start ini\n";
#endif
	this->set_indexes();
	this->set_sf();
	this->set_bounces();
	this->set_guess();

	ApplicationReturnStatus status;

	app->Options()->SetNumericValue("tol", config.NLP_tol);
	app->Options()->SetNumericValue("constr_viol_tol",config.constr_viol_tol);
	app->Options()->SetStringValue("mu_strategy", "adaptive");
//	app->Options()->SetStringValue("output_file", "ipopt.out");
//	app->Options()->SetStringValue("nlp_scaling_method","gradient-based");
	app->Options()->SetStringValue("linear_solver", nlp_solver(config.NLP_solver));//	ma86 & ma57 with memmory leakage
	app->Options()->SetIntegerValue("max_iter", config.max_iter);
	app->Options()->SetIntegerValue("print_level", config.print_level);
	app->Options()->SetIntegerValue("print_frequency_iter", config.print_freq);
	if (config.warmstart) {
		app->Options()->SetStringValue("warm_start_init_point", "yes");
	}
	app->Options()->SetNumericValue("mu_max",config.mu_max);
	if (config.H_approximation){
		app->Options()->SetStringValue("hessian_approximation", "limited-memory");
	}
	else {
		app->Options()->SetStringValue("hessian_approximation", "exact");
	}
	if (config.warmstart) {
		app->Options()->SetStringValue("warm_start_init_point","yes");
		app->Options()->SetNumericValue("warm_start_bound_push",1e-4);
		app->Options()->SetNumericValue("warm_start_mult_bound_push",1e-4);
		app->Options()->SetNumericValue("mu_init",config.ipopt_mu_ini);
	}

	status = app->Initialize();
  	if (status != Solve_Succeeded) {
  		printf("\n\n*** Error during initialization!\n");
  	}
  	return status;
#ifdef DCOPT_DEBUG
	cout<<"end ini\n";
#endif
}

void MyADOLC_sparseNLP::set_indexes() {
#ifdef DCOPT_DEBUG
	cout<<"start indexing\n";
#endif
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
			if (config.disc_method == Hermite_Simpson) {
				for (Index j = 0; j < n_controls; j++) {
					control_idx[2*i][j] 		= idx_n;
					idx_n++;
					if (i < n_nodes-1) {
						control_idx[2*i+1][j] 		= idx_n;
						idx_n++;
					}
				}
			}
			else {
				for (Index j = 0; j < n_controls; j++) {
					control_idx[i][j] 		= idx_n;
					idx_n++;
				}
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
		for (Index i = 0; i < n_events; i++) {
			event_idx[i]		= idx_m;
			idx_m++;
		}
		if (idx_m != nlp_m)
			printf("something went wrong in setting indexes (g static case)\n");
	}
	else {
		for (Index i = 0; i < n_nodes; i++) {
			if ( config.disc_method == Hermite_Simpson) {
				for (Index j = 0; j < n_path_constraints; j++) {
					path_constraint_idx[2*i][j]		= idx_m;
					idx_m++;
					if (i < n_nodes - 1) {
						path_constraint_idx[2*i+1][j]	= idx_m;
						idx_m++;
					}
				}
			}
			else {
				for (Index j = 0; j < n_path_constraints; j++) {
						path_constraint_idx[i][j]	= idx_m;
						idx_m++;
				}
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


#ifdef DCOPT_DEBUG
	cout<<"end indexing\n";
#endif
}

void MyADOLC_sparseNLP::set_sf(){
#ifdef DCOPT_DEBUG
	cout<<"start set_sf\n";
#endif
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

	if (config.disc_method == Hermite_Simpson) {
		for (Index j = 1; j <= n_controls; j++) {
			if(lb_controls(j) != 0 || ub_controls(j) != 0) {
				for (Index i = 0; i < 2*n_nodes - 1; i++) {
					nlp_sf_x[control_idx[i][j-1]] = max(fabs(lb_controls(j)), fabs(ub_controls(j)));
				}
			}
			else {
				for (Index i = 0; i < n_nodes; i++) {
					nlp_sf_x[control_idx[i][j-1]] = 1;
				}
			}
		}
	}
	else {
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
	}

	for (Index i = 1; i <= n_parameters; i++) {
		if(lb_parameters(i) != 0 || ub_parameters(i) != 0) {
			nlp_sf_x[parameter_idx[i-1]] = max(fabs(lb_parameters(i)), fabs(ub_parameters(i)));
		}
		else {
			nlp_sf_x[parameter_idx[i-1]] = 1;
		}
	}
	if (n_nodes && n_states) {
		if(lb_t0 != 0 || ub_tf != 0) {
			nlp_sf_x[t0_idx] = max(fabs(lb_t0(1)), fabs(ub_tf(1)));
			nlp_sf_x[tf_idx] = max(fabs(lb_t0(1)), fabs(ub_tf(1)));
		}
		else {
			nlp_sf_x[t0_idx] = 1;
			nlp_sf_x[tf_idx] = 1;
		}
	}

	if (config.disc_method == Hermite_Simpson) {
		for (Index j = 1; j <= n_path_constraints; j++) {
			if(lb_path(j) != 0 || ub_path(j) != 0) {
				for (Index i = 0; i < 2*n_nodes-1; i++)
					nlp_sf_g[path_constraint_idx[i][j-1]] 	= max(fabs(lb_path(j)), fabs(ub_path(j)));
			}
			else {
				for (Index i = 0; i < n_nodes; i++)
					nlp_sf_g[path_constraint_idx[i][j-1]] 	= 1;
			}
		}
	}
	else {
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

#ifdef DCOPT_DEBUG
	cout<<"end set_sf\n";
#endif
}

void MyADOLC_sparseNLP::set_bounces() {

	for (Index i = 0; i < n_parameters; i++) {
		nlp_lb_x[parameter_idx[i]] = lb_parameters(i+1)/nlp_sf_x[parameter_idx[i]];
		nlp_ub_x[parameter_idx[i]] = ub_parameters(i+1)/nlp_sf_x[parameter_idx[i]];
	}

	if (n_nodes && n_states) {
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

void MyADOLC_sparseNLP::set_guess() {
#ifdef DCOPT_DEBUG
	cout<<"set_guess\n";
#endif
	if (config.warmstart) {
		if (guess.z.getRowDim() == nlp_n && guess.lam_g.getRowDim() == nlp_m && guess.nodes.getRowDim() == n_nodes) {

			double delta_t 	= guess.nodes(n_nodes,1) - guess.nodes(1,1);
			node_str.resize(n_nodes-1,1);
			for (Index i = 1; i <= n_nodes - 1; i++) {
				node_str(i)	= (guess.nodes(i+1,1) - guess.nodes(i,1))/delta_t;
			}

			nlp_guess_z_L 	= new double [nlp_n];
			nlp_guess_z_U 	= new double [nlp_n];
			nlp_guess_lam	= new double [nlp_m];

			for (Index i = 0; i < nlp_n; i++) {
				nlp_guess_x[i]		= guess.z(i+1,1);
				nlp_guess_z_L[i]	= guess.z_L(i+1,1);
				nlp_guess_z_U[i]	= guess.z_U(i+1,1);
			}
			for (Index i = 0; i < nlp_m; i++){
				nlp_guess_lam[i]	= guess.lam_g(i+1,1);
			}
		}
		else{
			cout<< "\n\nNo lambda provided\n\n";
			exit(1);
		}

	}
	else {
		if (n_nodes == 0 || n_states == 0) {
			if (guess.parameters.getRowDim() != (uint)n_parameters)  {

				cout<<"===== No guess provided or invalid guess =====\n"
					  "=====  Generating guess based on Bounds  =====\n\n";
				guess_gen();
			}
			for (Index i = 0; i < n_parameters; i++) {
				nlp_guess_x[parameter_idx[i]] = guess.parameters(i+1)/nlp_sf_x[parameter_idx[i]];
			}
		}
		else {
			if ((guess.nodes.getRowDim() != (uint)n_nodes) 	||
				(guess.parameters.getRowDim() != (uint)n_parameters) 	||
				((guess.u.getColDim() != (uint)n_controls)	&&	guess.u_full.getColDim() != (uint)n_controls) ||
				(guess.x.getColDim() != (uint)n_states)	||
				((guess.u.getRowDim() != (uint)n_nodes)	&&	guess.u_full.getRowDim() != (uint)(n_nodes*2-1)) ||
				(guess.x.getRowDim() != (uint)n_nodes)	) {

				cout<<"===== No guess provided or invalid guess =====\n"
						"=====  Generating guess based on Bounds  =====\n\n";
				guess_gen();
			}

			double delta_t 	= guess.nodes(n_nodes,1) - guess.nodes(1,1);

			node_str.resize(n_nodes-1,1);
			for (Index i = 1; i <= n_nodes - 1; i++) {
				node_str(i)	= (guess.nodes(i+1,1) - guess.nodes(i,1))/delta_t;
			}

			for (Index i = 0; i < n_parameters; i++) {
				nlp_guess_x[parameter_idx[i]] = guess.parameters(i+1)/nlp_sf_x[parameter_idx[i]];
			}

			for (Index i = 0; i < n_nodes; i++) {
				for (Index j = 0; j < n_states; j++) {
					nlp_guess_x[state_idx[i][j]] = guess.x(i+1,j+1)/nlp_sf_x[state_idx[i][j]];
				}
				if (config.disc_method == Hermite_Simpson && guess.u_full.getColDim() == (2*n_nodes - 1)) {
					for (Index j = 0; j < n_controls; j++) {
						nlp_guess_x[control_idx[2*i][j]] = guess.u_full(2*i+1,j+1)/nlp_sf_x[control_idx[2*i][j]];
						if (i < n_nodes - 1) {
							nlp_guess_x[control_idx[2*i][j]] = guess.u_full(2*i+2,j+1)/nlp_sf_x[control_idx[2*i+1][j]];
						}
						nlp_guess_x[control_idx[2*i][j]] = guess.u(i+1,j+1)/nlp_sf_x[control_idx[2*i][j]];
					}
				}
				else if (config.disc_method == Hermite_Simpson){
					for (Index j = 0; j < n_controls; j++) {
						if (guess.u_full.getRowDim() == n_nodes && guess.u_full.getColDim() == n_controls) {
							nlp_guess_x[control_idx[2*i][j]] 	= guess.u_full(2*i+1,j+1)/nlp_sf_x[control_idx[2*i][j]];
							if (i < n_nodes - 1)
								nlp_guess_x[control_idx[2*i+1][j]] 	= guess.u_full((i+1)*2,j+1)/(nlp_sf_x[control_idx[2*i+1][j]]);
						}
						else {
						nlp_guess_x[control_idx[2*i][j]] 	= guess.u(i+1,j+1)/nlp_sf_x[control_idx[2*i][j]];
						if (i < n_nodes - 1)
							nlp_guess_x[control_idx[2*i+1][j]] 	= (guess.u(i+1,j+1)+guess.u(i+2,j+1))/(2*nlp_sf_x[control_idx[2*i+1][j]]);
						}
					}
				}
				else {
					for (Index j = 0; j < n_controls; j++) {
						nlp_guess_x[control_idx[i][j]] = guess.u(i+1,j+1)/nlp_sf_x[control_idx[i][j]];
					}
				}
			}

			nlp_guess_x[t0_idx]		= guess.nodes(1,1)/nlp_sf_x[t0_idx];
			nlp_guess_x[tf_idx]		= guess.nodes(n_nodes,1)/nlp_sf_x[tf_idx];
		}
	}
#ifdef DCOPT_DEBUG
	cout<<"end set_guess\n";
#endif
}

void MyADOLC_sparseNLP::guess_gen() {
#ifdef DCOPT_DEBUG
	cout<<"guess generation\n";
#endif
	if( n_nodes == 0 || n_states == 0) {
		if(n_parameters > 0) {
			guess.parameters.resize(n_parameters,1);
			for (Index i = 1; i <= n_parameters; i++) {
				guess.parameters(i) 	= (lb_parameters(i) + ub_parameters(i))/2;
			}
		}
	}
	else {
		guess.nodes = linspace<double>((lb_t0+ub_t0)/2.0,(lb_tf+ub_tf)/2.0,n_nodes);
		guess.x.resize(n_nodes, n_states);
		for (Index i = 1; i <= n_states; i++) {
			guess.x.setCol(i,linspace<double>(lb_states(i), ub_states(i),n_nodes));
		}
		if (n_controls) {
			guess.u.resize(n_nodes,n_controls);
			for (Index i = 1; i <= n_controls; i++) {
				guess.u.setCol(i,linspace<double>(lb_controls(i), ub_controls(i),n_nodes));
			}
			if (config.disc_method == Hermite_Simpson) {
				guess.u_full.resize(2*n_nodes-1,n_controls);
				for (Index i = 1; i <= n_controls; i++) {
					guess.u_full.setCol(i,linspace<double>(lb_controls(i), ub_controls(i),2*n_nodes-1));
				}
			}
		}
		else {
			SMatrix<double> Null_vector;
			guess.u = Null_vector;
		}

		if(n_parameters > 0) {
			guess.parameters.resize(n_parameters,1);
			for (Index i = 1; i <= n_parameters; i++) {
				guess.parameters(i) 	= (double)(lb_parameters(i) + ub_parameters(i))/2;
			}
		}
	}
#ifdef DCOPT_DEBUG
	cout<<"end guess generation\n";
#endif
}

ApplicationReturnStatus MyADOLC_sparseNLP::solve(SmartPtr<IpoptApplication> app){


	ApplicationReturnStatus status;

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

	adouble *d_z, *p_states_ini, *p_states_fin, *p_parameters, *p_t0, *p_tf;
	d_z 			= new adouble [n];

	for( Index i = 0; i < n; i++)
		d_z[i]		= z[i]*nlp_sf_x[i];

	if (n_nodes != 0 && n_states != 0) {
		p_states_ini	= &d_z[state_idx[0][0]];
		p_states_fin	= &d_z[state_idx[n_nodes - 1][0]];
		p_t0 			= &d_z[t0_idx];
		p_tf 			= &d_z[tf_idx];
	}
	else{
		p_states_ini	= NULL;
		p_states_fin	= NULL;
		p_t0 			= NULL;
		p_tf 			= NULL;
	}

	if (n_parameters)
		p_parameters 	= &d_z[parameter_idx[0]];
	else
		p_parameters	= NULL;

	obj_value = ad_e_cost (p_states_ini, p_states_fin, p_parameters, *p_t0, *p_tf, 1, constants);

	delete[] d_z;

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

	double 	*d_z, *p_states_ini, *p_states_fin, *p_parameters, *p_t0, *p_tf;
	d_z 			= new double [n];

	for( Index i = 0; i < n; i++)
		d_z[i]		= z[i]*nlp_sf_x[i];

	if (n_nodes != 0 && n_states != 0) {
		p_states_ini	= &d_z[state_idx[0][0]];
		p_states_fin	= &d_z[state_idx[n_nodes - 1][0]];
		p_t0 			= &d_z[t0_idx];
		p_tf 			= &d_z[tf_idx];
	}
	else{
		p_states_ini	= NULL;
		p_states_fin	= NULL;
		p_t0 			= NULL;
		p_tf 			= NULL;
	}

	if (n_parameters)
		p_parameters 	= &d_z[parameter_idx[0]];
	else
		p_parameters	= NULL;


	obj_value = d_e_cost (p_states_ini, p_states_fin, p_parameters, *p_t0, *p_tf, 1, constants);

	delete[] d_z;

#ifdef DCOPT_DEBUG
	printf("end of eval_obj()\n");
#endif
	return true;
}

bool  MyADOLC_sparseNLP::ad_eval_constraints(Index n, const adouble *x, Index m, adouble* g) {

#ifdef DCOPT_DEBUG
	printf("ad_eval_constraints()\n");
#endif

	adouble *d_z, *p_controls, *p_states0,*p_states1, *p_states_ini, *p_states_fin, *p_parameters, *p_t0, *p_tf;
	adouble *p_path = NULL, *p_events, *p_defects;

	d_z 			= new adouble [n];
	for( Index i = 0; i < n; i++)
		d_z[i]		= x[i]*nlp_sf_x[i];

	if (n_events)
		p_events		= &g[event_idx[0]];
	else
		p_events		= NULL;

	if (n_nodes != 0 && n_states != 0) {
		p_states_ini	= &d_z[state_idx[0][0]];
		p_states_fin	= &d_z[state_idx[n_nodes - 1][0]];
		p_t0 			= &d_z[t0_idx];
		p_tf 			= &d_z[tf_idx];
	}
	else{
		p_states_ini	= NULL;
		p_states_fin	= NULL;
		p_t0 			= NULL;
		p_tf 			= NULL;
	}

	if (n_parameters)
		p_parameters 	= &d_z[parameter_idx[0]];
	else
		p_parameters	= NULL;


	ad_events(p_events, p_states_ini, p_states_fin, p_parameters, *p_t0, *p_tf, 1, constants);

	if(n_states != 0 && n_nodes != 0) {
		adouble **f		 	= new adouble *[n_nodes];

		adouble *t	 		= new adouble [n_nodes];
		adouble *delta		= new adouble [n_nodes - 1];

		for (Index i = 0; i < n_nodes; i += 1) {
			f[i] 		= new adouble [n_states];
		}

		t[0]				= *p_t0;
		for (Index i = 0; i < n_nodes - 1; i++) {
			delta[i] 		= (*p_tf-*p_t0)*node_str(i+1);
			t[i+1]			= t[i] + (*p_tf-*p_t0)*node_str(i+1);
		}
		t[n_nodes - 1]	= *p_tf;

		if (n_path_constraints)
			p_path 			= &g[path_constraint_idx[0][0]];
		p_states0 			= &d_z[state_idx[0][0]];
		p_controls 			= &d_z[control_idx[0][0]];

		ad_derv(f[0], p_path, p_states0, p_controls, p_parameters, t[0], 1, constants);

		if (config.disc_method == Hermite_Simpson) {
			adouble* y_m			= new adouble [n_states];
			adouble* f_m			= new adouble [n_states];

			for (Index i = 0; i < n_nodes - 1; i += 1)	{
				if (n_path_constraints)
					p_path 		= &g[path_constraint_idx[2*(i+1)][0]];
				p_states1 			= &d_z[state_idx[i+1][0]];
				p_controls 			= &d_z[control_idx[2*(i+1)][0]];
				ad_derv(f[i+1], p_path, p_states1, p_controls, p_parameters, t[i+1], 1, constants);

				for (Index j = 0; j < n_states; j += 1){
					y_m[j] 	= (p_states0[j] + p_states1[j])/2. + delta[i]/8.*(f[i][j]-f[i+1][j]);
				}
				if (n_path_constraints)
					p_path 		= &g[path_constraint_idx[2*i+1][0]];
				p_controls 			= &d_z[control_idx[2*i+1][0]];
				ad_derv(f_m, p_path, y_m, p_controls, p_parameters, (t[i] + t[i+1])/2., 1, constants);

				p_defects 		= &g[defect_idx[i][0]];
				for (Index j = 0; j < n_states; j++) {
					p_defects[j]	= (p_states1[j] - p_states0[j] - delta[i]/6.*(f[i][j]+4.*f_m[j]+f[i+1][j]));
				}
				p_states0 = p_states1;
			}
			delete[] y_m;
			delete[] f_m;

		}
		else if (config.disc_method == trapezoidal) {

			for (Index i = 0; i < n_nodes - 1; i += 1)	{
				if (n_path_constraints)
					p_path 			= &g[path_constraint_idx[i+1][0]];
				p_states1 			= &d_z[state_idx[i+1][0]];
				p_controls 			= &d_z[control_idx[i+1][0]];

				ad_derv(f[i+1], p_path, p_states1, p_controls, p_parameters, t[i+1], 1, constants);

				p_defects 		= &g[defect_idx[i][0]];
				for (Index j = 0; j < n_states; j++) {
					p_defects[j] 	= p_states1[j] - p_states0[j] - delta[i]/2.0*(f[i][j] + f[i+1][j]);
				}
				p_states0 = p_states1;
			}
		}

		for (Index i = 0; i < nlp_m; i++) {
			g[i] = g[i]/nlp_sf_g[i]+1;
		}

		for (Index i = 0; i < n_nodes; i++) {
			delete[] f[i];
		}
		delete[] f;
		delete[] t;
		delete[] delta;
		delete[] d_z;
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

	double *d_z, *p_controls, *p_states0,*p_states1, *p_states_ini, *p_states_fin, *p_parameters, *p_t0, *p_tf;
	double *p_path = NULL, *p_events, *p_defects;

	d_z 			= new double [n];
	for( Index i = 0; i < n; i++)
		d_z[i]		= x[i]*nlp_sf_x[i];

	if (n_events)
		p_events		= &g[event_idx[0]];
	else
		p_events		= NULL;

	if (n_nodes != 0 && n_states != 0) {
		p_states_ini	= &d_z[state_idx[0][0]];
		p_states_fin	= &d_z[state_idx[n_nodes - 1][0]];
		p_t0 			= &d_z[t0_idx];
		p_tf 			= &d_z[tf_idx];
	}
	else{
		p_states_ini	= NULL;
		p_states_fin	= NULL;
		p_t0 			= NULL;
		p_tf 			= NULL;
	}

	if (n_parameters)
		p_parameters 	= &d_z[parameter_idx[0]];
	else
		p_parameters	= NULL;


	d_events(p_events, p_states_ini, p_states_fin, p_parameters, *p_t0, *p_tf, 1, constants);

	if(n_states != 0 && n_nodes != 0) {
		double **f		 	= new double *[n_nodes];

		double *t	 		= new double [n_nodes];
		double *delta		= new double [n_nodes - 1];

		for (Index i = 0; i < n_nodes; i += 1) {
			f[i] 		= new double [n_states];
		}

		t[0]				= *p_t0;
		for (Index i = 0; i < n_nodes - 1; i++) {
			delta[i] 		= (*p_tf-*p_t0)*node_str(i+1);
			t[i+1]			= t[i] + (*p_tf-*p_t0)*node_str(i+1);
		}
		t[n_nodes - 1]	= *p_tf;

		if (n_path_constraints)
			p_path 			= &g[path_constraint_idx[0][0]];
		p_states0 			= &d_z[state_idx[0][0]];
		p_controls 			= &d_z[control_idx[0][0]];

		d_derv(f[0], p_path, p_states0, p_controls, p_parameters, t[0], 1, constants);

		if (config.disc_method == Hermite_Simpson) {
			double* y_m			= new double [n_states];
			double* f_m			= new double [n_states];

			for (Index i = 0; i < n_nodes - 1; i += 1)	{
				if (n_path_constraints)
					p_path 		= &g[path_constraint_idx[2*(i+1)][0]];
				p_states1 			= &d_z[state_idx[i+1][0]];
				p_controls 			= &d_z[control_idx[2*(i+1)][0]];
				d_derv(f[i+1], p_path, p_states1, p_controls, p_parameters, t[i+1], 1, constants);

				for (Index j = 0; j < n_states; j += 1){
					y_m[j] 	= (p_states0[j] + p_states1[j])/2. + delta[i]/8.*(f[i][j]-f[i+1][j]);
				}
				if (n_path_constraints)
					p_path 		= &g[path_constraint_idx[2*i+1][0]];
				p_controls 			= &d_z[control_idx[2*i+1][0]];
				d_derv(f_m, p_path, y_m, p_controls, p_parameters, (t[i] + t[i+1])/2., 1, constants);

				p_defects 		= &g[defect_idx[i][0]];
				for (Index j = 0; j < n_states; j++) {
					p_defects[j]	= (p_states1[j] - p_states0[j] - delta[i]/6.*(f[i][j]+4.*f_m[j]+f[i+1][j]));
				}
				p_states0 = p_states1;
			}
			delete[] y_m;
			delete[] f_m;

		}
		else if (config.disc_method == trapezoidal) {

			for (Index i = 0; i < n_nodes - 1; i += 1)	{
				if (n_path_constraints)
					p_path 			= &g[path_constraint_idx[i+1][0]];
				p_states1 			= &d_z[state_idx[i+1][0]];
				p_controls 			= &d_z[control_idx[i+1][0]];

				d_derv(f[i+1], p_path, p_states1, p_controls, p_parameters, t[i+1], 1, constants);

				p_defects 		= &g[defect_idx[i][0]];
				for (Index j = 0; j < n_states; j++) {
					p_defects[j] 	= p_states1[j] - p_states0[j] - delta[i]/2.0*(f[i][j] + f[i+1][j]);
				}
				p_states0 = p_states1;
			}
		}

		for (Index i = 0; i < nlp_m; i++) {
			g[i] = g[i]/nlp_sf_g[i]+1;
		}

		for (Index i = 0; i < n_nodes; i++) {
			delete[] f[i];
		}
		delete[] f;
		delete[] t;
		delete[] delta;
		delete[] d_z;
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
	cout<<"\n";

	for (Index i = 0; i < nlp_m; i += 1) {
		g_l[i] = nlp_lb_g[i]+1;
		g_u[i] = nlp_ub_g[i]+1;
	}

	return true;
}

bool MyADOLC_sparseNLP::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{
//	assert(init_x == true);
//	assert(init_z == false);
//	assert(init_lambda == false);

	if (init_z == true) {
		for (Index i = 0; i < n; i++) {
			z_L[i]	= nlp_guess_z_L[i];
			z_U[i]	= nlp_guess_z_U[i];
		}
	}

	if (init_lambda == true) {
		for (Index i = 0; i < m; i++) {
			lambda[i]	= nlp_guess_lam[i];
		}
	}

	for (Index i = 0; i < n; i += 1) {
		x[i]	= nlp_guess_x[i];
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

		sparse_jac(tag_g + omp_get_thread_num()*10, m, n, 1, x, &nnz_jac, &rind_g, &cind_g, &jacval, options_g);
		for(Index idx=0; idx<nnz_jac; idx++) {
			values[idx] = jacval[idx];
		}
/*
 *
 *
 *
 *

		dcomp *p_controls1, *p_controls0, *p_path, *p_events, *p_states0,*p_states1, *p_states_ini, *p_states_fin, *p_parameters, *p_t0, *p_tf;
		dcomp *p_defects;

		dcomp *f0		= new dcomp [n_states-1];
		dcomp *f1		= new dcomp [n_states-1];
		dcomp *y_m		= new dcomp [n_states-1];
		dcomp *f_m		= new dcomp [n_states-1];

		dcomp *delta		= new dcomp [n_nodes - 1];
		dcomp *t	 	= new dcomp [n_nodes];

		dcomp *x_comp	= new dcomp[n];
		dcomp *g_comp	= new dcomp[m];

		for (int i = 0; i < n; i++)
			x_comp[i]	= x[i];


		t[0]				= t0;
			for (Index i = 0; i < n_nodes - 1; i++) {
				delta[i]	= (tf-t0)*node_str(i+1);
				t[i+1]		= t[i] + delta[i];
			}
			t[n_nodes - 1]	= tf;

//			NLP_x_2_OCP_var(x_comp,y,u,param,t0,tf);

			p_events		= &g_comp[event_idx[0]];
			p_states_ini	= &x_comp[state_idx[0][0]];
			p_states_fin	= &x_comp[state_idx[n_nodes - 1][0]];
			p_t0 			= &x_comp[t0_idx];
			p_t0 			= &x_comp[tf_idx];
			if (n_parameters)
				p_parameters= &x_comp[parameter_idx[0]];

			if (config.disc_method == Hermite_Simpson) {
				ad_derv(f[0], path[0], y[0], u[0], param, t[0], 1, constants);
				for (Index i = 0; i < n_nodes - 1; i += 1)	{
					ad_derv(f[i+1], path[2*(i+1)], y[i+1], u[2*(i+1)], param, t[i+1], 1, constants);

					for (Index j = 0; j < n_states; j += 1){
						y_m[j] 	= (y[i][j]+y[i+1][j])/2. + delta[i]/8.*(f[i][j]-f[i+1][j]);
					}

					ad_derv(f_m, path[2*i+1], y_m, u[2*i+1], param, (t[i] + t[i+1])/2., 1, constants);

					for (Index j = 0; j < n_states; j++) {
						p_defects	= (p_states1[j+1] - p_states0[j] - delta[i]/6.*(f[i][j]+4.*f_m[j]+f[i+1][j]));
						p_defects++;
						p_states0++;
						p_states1++;
					}
				}
			}
			else if (config.disc_method == trapezoidal) {
				for (Index i = 0; i < n_nodes; i += 1)	{
					ad_derv(f[i], p_path, y[i], u[i], param, t[i], 1, constants);
				}
				for (Index i = 0; i < n_nodes - 1; i++) {
					p_defects 	= g_comp[defect_idx[i][0]];
					p_states0	= x_comp[defect_idx[i][0]];
					p_states1	= x_comp[defect_idx[i+1][0]];

					for (Index j = 0; j < n_states; j++) {
						p_defects 	= p_states1[j+1] - p_states0[j] - delta[i]/2.0*(f1[j] + f0[j]);
						p_defects++;
						p_states0++;
						p_states1++;
					}
				}
			}

			d_events(p_events, p_states_ini, p_states_fin, p_parameters, p_t0, p_t0, 1, constants);

			for (Index i = 0; i < nlp_m; i++) {
				g_comp[i] = g_comp[i] + 1;
			}

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

			delete[] delta;
			delete[] t;


			int nnz = 0;
	//#pragma omp parallel
		{
			int this_thread = omp_get_thread_num(), num_threads = omp_get_num_threads();
			for (int col = this_thread; col < n; col+=num_threads) {
				x_comp[col] 	+= dcomp (0,step);



				eval_constraints(n,	x_comp, m, g_comp); // eval_constraint is not thread safe
				x_comp[col] 	= x[col];
				for (int row = 0; row < m; row++)
					if (rind_g[nnz] == row && cind_g[nnz] == col) {
						values[nnz]	= g_comp[row].imag()/step;
						nnz++;
					}
			}
			delete[] g_comp;
		}
			delete[] x_comp;
//*/
	}

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
		}
		else {
			// return the values. This is a symmetric matrix, fill the lower left
			// triangle only
			obj_lam[0] = obj_factor;
				for(Index idx = 0; idx<m ; idx++)
					obj_lam[1+idx] = lambda[idx];

			set_param_vec(tag_L,m+1,obj_lam);

		    sparse_hess(tag_L, n, 1, const_cast<double*>(x), &nnz_L, &rind_L, &cind_L, &hessval, options_L);
		    for(Index idx = 0; idx <nnz_L ; idx++) {
		    	values[idx] = hessval[idx];
			}


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
	setHamiltonian(x,lambda);
	results.z.resize(nlp_n,1);
	results.z_L.resize(nlp_n,1);
	results.z_U.resize(nlp_n,1);
	results.lam_g.resize(nlp_m,1);
	results.x.resize(n_nodes, n_states);
	results.u.resize(n_nodes, n_controls);
	if (config.disc_method == Hermite_Simpson)
		results.u_full.resize(2*n_nodes-1, n_controls);

	results.parameters.resize(n_parameters,1);
	results.nodes.resize(n_nodes, 1);

	for (Index i = 0; i < nlp_n; i++) {
		results.z(i+1) 		= x[i];
		results.z_L(i+1) 	= z_L[i];
		results.z_U(i+1) 	= z_U[i];
	}
	for (Index i = 0; i < nlp_m; i++){
		results.lam_g(i+1)	= lambda[i];
	}

	for (Index i = 0; i < n_nodes; i++)
		for (Index j = 0; j < n_states; j++)
			results.x(i+1,j+1) 		= x[state_idx[i][j]]*nlp_sf_x[state_idx[i][j]];

	for (Index i = 0; i < n_nodes; i++) {
		for (Index j = 0; j < n_controls; j++){
			if (config.disc_method == Hermite_Simpson){
				results.u(i+1,j+1) 		= x[control_idx[2*i][j]]*nlp_sf_x[control_idx[2*i][j]];
				results.u_full	(2*i+1,j+1) 		= x[control_idx[2*i][j]]*nlp_sf_x[control_idx[2*i][j]];
				if (i < n_nodes - 1)
					results.u_full(2*i+2,j+1) 		= x[control_idx[2*i+1][j]]*nlp_sf_x[control_idx[2*i+1][j]];
			}
			else
				results.u(i+1,j+1) 		= x[control_idx[i][j]]*nlp_sf_x[control_idx[i][j]];
		}
	}


	for (Index i = 0; i < n_parameters; i++)
		results.parameters(i+1)		= x[parameter_idx[i]]*nlp_sf_x[parameter_idx[i]];

	if (n_nodes && n_states) {
		for (Index i = 2; i < n_nodes; i++) {
			results.nodes(i,1) 	= results.nodes(i-1,1) + (x[tf_idx]*nlp_sf_x[tf_idx] - x[t0_idx]*nlp_sf_x[t0_idx])*node_str(i-1);
		}
		results.nodes(n_nodes,1) = x[tf_idx]*nlp_sf_x[tf_idx];
	}

	if(n_events) {
		results.lam_events.resize(n_events,1);
		for (Index i = 0; i < n_events; i++) {
			results.lam_events(i+1)		= -lambda[event_idx[i]]	;
		}
	}

	if (n_nodes && n_states) {

		results.lam_defects.resize(n_nodes - 1, n_states);
		for (Index i = 0; i < n_nodes-1; i++) {
			for (Index j = 0; j < n_states; j++) {
				results.lam_defects(i+1,j+1) 	= -lambda[defect_idx[i][j]];
			}
		}

		for (Index i = 0; i < n_nodes; i++) {
			delete[] state_idx[i];
			if (i< n_nodes - 1) {
				delete[] defect_idx[i];
			}
		}



		if (config.disc_method == Hermite_Simpson) {
			for (Index i = 0; i < 2*n_nodes - 1; i++) {
				delete[] control_idx[i];
				delete[] path_constraint_idx[i];
			}
		}
		else {
			for (Index i = 0; i < n_nodes; i++) {
				delete[] control_idx[i];
				delete[] path_constraint_idx[i];
			}
		}

		delete[] state_idx;

		delete[] defect_idx;

		delete[] control_idx;

		delete[] path_constraint_idx;
	}

	if (n_parameters) {
		delete[] parameter_idx;
	}

	if (n_events) {
		delete[] event_idx;
	}

	delete[] obj_lam;
	free(rind_g);
	free(cind_g);
	free(jacval);


		if (!config.H_approximation) {
			free(rind_L);
			free(cind_L);
			free(hessval);
		}

	delete[] nlp_sf_x;
	delete[] nlp_guess_z_L;
	delete[] nlp_guess_z_U;
	delete[] nlp_guess_lam;
	delete[] nlp_sf_g;
	delete[] nlp_guess_x;
	delete[] nlp_lb_x;
	delete[] nlp_ub_x;
	delete[] nlp_lb_g;
	delete[] nlp_ub_g;

}

void MyADOLC_sparseNLP::setHamiltonian(const double* x, const double* lam_g) {
	double *d_z, *p_controls, *p_states, *p_parameters, *p_t0, *p_tf, *p_path;

	d_z 			= new double [nlp_n];
	for( Index i = 0; i < nlp_n; i++)
		d_z[i]		= x[i]*nlp_sf_x[i];

	if (n_parameters)
		p_parameters 	= &d_z[parameter_idx[0]];
	else
		p_parameters	= NULL;

	if (n_path_constraints)
		p_path 			= new double [n_path_constraints];
	else
		p_path			= NULL;

	if(n_states != 0 && n_nodes != 0) {
		double *f		 	= new double [n_states];

		double *t	 		= new double [n_nodes];
		double *delta		= new double [n_nodes - 1];

		p_t0 				= &d_z[t0_idx];
		p_tf 				= &d_z[tf_idx];

		t[0]				= *p_t0;
		for (Index i = 0; i < n_nodes - 1; i++) {
			delta[i] 		= (*p_tf-*p_t0)*node_str(i+1);
			t[i+1]			= t[i] + (*p_tf-*p_t0)*node_str(i+1);
		}
		t[n_nodes - 1]	= *p_tf;


		results.Hamiltonian.resize(n_nodes - 1,1);

		for (Index i = 0; i < n_nodes - 1; i += 1)	{
			p_states 			= &d_z[state_idx[i+1][0]];
			p_controls 			= &d_z[control_idx[i+1][0]];

			d_derv(f, p_path, p_states, p_controls, p_parameters, t[i+1], 1, constants);

			for (Index j = 0; j < n_states; j++) {
				results.Hamiltonian(i+1) += f[j]*(-lam_g[defect_idx[i][j]]);
			}
		}
		delete[] f;
		delete[] t;
		delete[] delta;
		if (n_path_constraints)
			delete[] p_path;
	}


	delete[] d_z;
}

void MyADOLC_sparseNLP::generate_tapes(Index n, Index m, Index& nnz_jac_g, Index& nnz_h_lag)
{
#ifdef DCOPT_DEBUG
	printf("generate_tapes()\n");
#endif
	Number *xp   	= new double[n];
	Number *lamp 	= new double[m];
	Number *zl   	= new double[m];
	Number *zu   	= new double[m];

	adouble *xf  	= new adouble[n];
	adouble *xg  	= new adouble[n];
	adouble *xL   	= new adouble[n];

	adouble *g   	= new adouble[m];
	double *lam  	= new double[m];
	double sig;
	adouble obj_value;

	obj_lam   = new double[m+1];

	double dummy;

	get_starting_point(n, 1, xp, 0, zl, zu, m, 0, lamp);

	trace_on(tag_f + omp_get_thread_num()*10);
    
    for(Index idx=0;idx<n;idx++)
      xf[idx] <<= xp[idx];

    ad_eval_obj(n,xf,obj_value);

    obj_value >>= dummy;

    trace_off();

    trace_on(tag_g + omp_get_thread_num()*10);

    for(Index idx=0;idx<n;idx++)
    	xg[idx] <<= xp[idx];

    ad_eval_constraints(n,xg,m,g);

    for(Index idx=0;idx<m;idx++)
		g[idx] >>= dummy;
    trace_off();

	if (!config.H_approximation) {
		trace_on(tag_L);

		for(Index idx=0;idx<n;idx++)
			xL[idx] <<= xp[idx];
		for(Index idx=0;idx<m;idx++)
			lam[idx] = 1.0;

		sig = 1.0;

		ad_eval_obj(n,xL,obj_value);

		obj_value *= mkparam(sig);
		ad_eval_constraints(n,xL,m,g);

		for(Index idx=0;idx<m;idx++)
			obj_value += g[idx]*mkparam(lam[idx]);

		obj_value >>= dummy;

		trace_off();
	}
    rind_g = NULL;
	cind_g = NULL;

	options_g[0] = 0;          /* sparsity pattern by index domains (default) */
	options_g[1] = 0;          /*                         safe mode (default) */
	options_g[2] = 0;         /*                     &jacval is not computed */
	options_g[3] = 0;          /*                column compression (default) */


	this->jacval=NULL;

	sparse_jac(tag_g + omp_get_thread_num()*10, m, n, 0, xp, &nnz_jac, &rind_g, &cind_g, &jacval, options_g);
	nnz_jac_g = nnz_jac;


/*
	sparse_jac(tag_g, m, n, 1, xp, &nnz_jac, &rind_g, &cind_g, &jacval, options_g);

	double* NLP_constraint_sf = new double[m];
	for (Index i=0; i<m; i++) {
		NLP_constraint_sf[i] = 0.0;
	}


	for (Index i=0;i<nnz_jac;i++) {
		cout<<"rind_g"<<rind_g[i]<<endl;
		NLP_constraint_sf[rind_g[i]] += jacval[i]*jacval[i];

	}
	cout<<endl;

	for(Index i = 0; i < m; i++) {
		printf("constraint_sf[%d/%d] = %e\n",i+1,m,NLP_constraint_sf[i]);
	}
	cout<<"common\n";

	delete[] NLP_constraint_sf;
//*/
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
		hessval			= NULL;
		cind_L			= NULL;
		rind_L			= NULL;
		options_L[0]= 0;
		options_L[1]= 1;
		sparse_hess(tag_L, n, 0, xp, &nnz_L, &rind_L, &cind_L, &hessval, options_L);
	}
	nnz_h_lag = nnz_L;

	delete[] lam;
	delete[] g;
	delete[] xf;
	delete[] xg;
	delete[] xL;
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
	constr_viol_tol = 1e-4;
	with_mgl		= false;
	disc_method		= trapezoidal;
	print_level		= 5;
	print_freq		= 1;
	H_approximation = true;
	ipopt_mu_ini	= 0.1;
	mu_max			= 1e5;
}

Config::~Config() {

}
