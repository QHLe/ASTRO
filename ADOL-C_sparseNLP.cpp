#define SPARSE_HESS

//#define DCOPT_DEBUG

#include <cassert>
#include "ADOL-C_sparseNLP.hpp"
using namespace Ipopt;

double step = std::numeric_limits<double>::epsilon();

time_t rand_time;

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


void MyADOLC_sparseNLP::set_endpoint_cost(dcomp (*cost)(		const  dcomp* ini_states,
																const  dcomp* fin_states,
																const  dcomp* param,
																const  dcomp& t0,
																const  dcomp& tf,
																Index phase, const double* constants)) {
	this->d_e_cost 	= cost;
}

void MyADOLC_sparseNLP::set_derivatives(void (*derv)( 	dcomp *states_dot,
														dcomp *path,
														const  dcomp *states,
														const  dcomp *controls,
														const  dcomp *param,
														const  dcomp &time,
														Index phase, const double* constants)) {
	this->d_derv 	= derv;
}

void MyADOLC_sparseNLP::set_events(	void (*events)(	dcomp *events,
													const  dcomp *ini_states,
													const  dcomp *fin_states,
													const  dcomp *param,
													const  dcomp &t0,
													const  dcomp &tf,
													Index phase, const double* constants)) {
	this->d_events 	= events;
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
		parameters			= new dcomp[n_parameters];
	}

	if (n_events) {
		event_idx			= new Index [n_events];
		events				= new dcomp[n_events];
	}

	if (n_nodes && n_states) {
		state_idx 			= new Index *[n_nodes];
		states  			= new dcomp*[n_nodes];

		defect_idx			= new Index *[n_nodes - 1];
		defects				= new dcomp*[n_nodes - 1];

		if (n_controls && config.disc_method == Hermite_Simpson) {
			control_idx			= new Index *[2*n_nodes-1];
			controls			= new dcomp*[2*n_nodes-1];

			path_constraint_idx	= new Index *[2*n_nodes-1];
			path_constraints	= new dcomp*[2*n_nodes-1];
		}
		else {
			control_idx			= new Index *[n_nodes];
			controls			= new dcomp*[n_nodes];

			path_constraint_idx	= new Index *[n_nodes];
			path_constraints	= new dcomp*[n_nodes];
		}
	}

	for (Index i = 0; i < n_nodes; i++) {
		if (n_states) {
			state_idx[i]	= new Index[n_states];
			states[i]		= new dcomp[n_states];
			if (i < n_nodes - 1) {
				defect_idx[i]		= new Index [n_states];
				defects[i]			= new dcomp [n_states];
			}
		}
		if (n_controls && config.disc_method == Hermite_Simpson) {
			control_idx[2*i]				= new Index[n_controls];
			controls[2*i]					= new dcomp[n_controls];
			if (n_path_constraints) {
				path_constraint_idx[2*i]	= new Index [n_path_constraints];
				path_constraints[2*i]		= new dcomp [n_path_constraints];
			}
			else {
				path_constraint_idx[2*i]	= NULL;
				path_constraints[2*i]		= NULL;
			}
			if (i < n_nodes - 1) {
				control_idx[2*i+1]			= new Index[n_controls];
				controls[2*i+1]				= new dcomp[n_controls];
				if (n_path_constraints) {
					path_constraint_idx[2*i+1]	= new Index [n_path_constraints];
					path_constraints[2*i+1]		= new dcomp [n_path_constraints];
				}
				else {
					path_constraint_idx[2*i+1]	= NULL;
					path_constraints[2*i+1]		= NULL;
				}
			}
		}
		else if (n_controls) {
			control_idx[i]	= new Index[n_controls];
			controls[i]		= new dcomp[n_controls];
			if (n_path_constraints) {
				path_constraint_idx[i]	= new Index [n_path_constraints];
				path_constraints[i]		= new dcomp [n_path_constraints];
			}
			else {
				path_constraint_idx[i]	= NULL;
				path_constraints[i]		= NULL;
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
		if(lb_states(j) != 0. || ub_states(j) != 0.) {
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
			if(lb_controls(j) != 0. || ub_controls(j) != 0.) {
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
			if(lb_controls(j) != 0. || ub_controls(j) != 0.) {
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
		if(lb_parameters(i) != 0. || ub_parameters(i) != 0.) {
			nlp_sf_x[parameter_idx[i-1]] = max(fabs(lb_parameters(i)), fabs(ub_parameters(i)));
		}
		else {
			nlp_sf_x[parameter_idx[i-1]] = 1;
		}
	}
	if (n_nodes && n_states) {
		if(lb_t0 != 0. || ub_tf != 0.) {
			nlp_sf_x[t0_idx] = max(fabs(lb_t0), fabs(ub_tf));
			nlp_sf_x[tf_idx] = max(fabs(lb_t0), fabs(ub_tf));
		}
		else {
			nlp_sf_x[t0_idx] = 1;
			nlp_sf_x[tf_idx] = 1;
		}
	}

	if (config.disc_method == Hermite_Simpson) {
		for (Index j = 1; j <= n_path_constraints; j++) {
			if(lb_path(j) != 0. || ub_path(j) != 0.) {
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
			if(lb_path(j) != 0. || ub_path(j) != 0.) {
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
		if(lb_events(j) != 0. || ub_events(j) != 0.) {
			nlp_sf_g[event_idx[j-1]] = max(fabs(lb_events(j)), fabs(ub_events(j)));
		}
		else {
			nlp_sf_g[event_idx[j-1]] = 1;
		}
	}
/*
	for(int i = 0; i < nlp_n; i++)
		nlp_sf_x[i]		=	1;

	for(int i = 0; i < nlp_m; i++)
		nlp_sf_g[i]		=	1;
//*/
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
			(guess.u.getColDim() != (uint)n_controls)	||
			(guess.x.getColDim() != (uint)n_states)	||
			(guess.u.getRowDim() != (uint)n_nodes)	||
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
			if (config.disc_method == Hermite_Simpson && (Index)guess.u_full.getColDim() == (2*n_nodes - 1)) {
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
					nlp_guess_x[control_idx[2*i][j]] 	= guess.u(i+1,j+1)/nlp_sf_x[control_idx[2*i][j]];
					if (i < n_nodes - 1)
						nlp_guess_x[control_idx[2*i+1][j]] 	= (guess.u(i+1,j+1)+guess.u(i+2,j+1))/(2*nlp_sf_x[control_idx[2*i+1][j]]);
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

bool  MyADOLC_sparseNLP::eval_obj(Index n, const dcomp *z, dcomp& obj_value) {
  // return the value of the objective function
#ifdef DCOPT_DEBUG
	printf("eval_obj()\n");
#endif
	if(n_states == 0 || n_nodes == 0) {
		dcomp dummy;
		NLP_x_2_OCP_var(z, states, controls, parameters, t0, tf);
		obj_value = d_e_cost (&dummy, &dummy, parameters, t0, tf, 1, constants);
	}
	else {
		NLP_x_2_OCP_var(z, states, controls, parameters, t0, tf);
		obj_value = d_e_cost (states[0], states[n_nodes-1], parameters, t0, tf, 1, constants);
	}

#ifdef DCOPT_DEBUG
	printf("end of eval_obj()\n");
#endif
	return true;
}

bool  MyADOLC_sparseNLP::eval_constraints(Index n, const dcomp *x, Index m, dcomp* g) {

#ifdef DCOPT_DEBUG
	printf("eval_constraints()\n");
#endif
	if(n_states == 0 || n_nodes == 0) {
		dcomp *dummy;
		NLP_x_2_OCP_var(x,states,controls,parameters,t0,tf);
		d_events(events, dummy, dummy, parameters, t0, tf, 1, constants);
		OCP_var_2_NLP_g(path_constraints, defects, events, g);
		for (Index i = 0; i < nlp_m; i++) {
			g[i] = g[i] + dcomp (1,0);
		}
	}
	else {
		dcomp **f		 	= new dcomp *[n_nodes];

		dcomp *t	 		= new dcomp [n_nodes];
		dcomp *delta		= new dcomp [n_nodes - 1];

		for (Index i = 0; i < n_nodes; i += 1) {
			f[i] 		= new dcomp [n_states];
		}

/*		for (int i = 0; i < nlp_n; i++)
			printf("nlp_x[%d] = (%.2e,%.2e)\n",i,x[i].real(),x[i].imag());
//*/
		NLP_x_2_OCP_var(x,states,controls,parameters,t0,tf);
/*
		for (int i = 0; i < n_nodes; i++) {
			for (int j = 0; j < n_states; j++) {
				printf("states[%d][%d] = (%.2e,%.2e)\t",i,j,states[i][j].real(),states[i][j].imag());
			}
			printf("\n");
		}

		for (int i = 0; i < n_nodes; i++) {
			for (int j = 0; j < n_controls; j++) {
				printf("controls[%d][%d] = (%.2e,%.2e)\t",i,j,controls[i][j].real(),controls[i][j].imag());
			}
			printf("\n");
		}
//*/
		d_events(events, states[0], states[n_nodes-1], parameters, t0, tf, 1, constants);

		t[0]				= t0;
		for (Index i = 0; i < n_nodes - 1; i++) {
			delta[i] 		= (tf-t0)*node_str(i+1);
			t[i+1]			= t[i] + (tf-t0)*node_str(i+1);
		}
		t[n_nodes - 1]	= tf;

		if (config.disc_method == Hermite_Simpson) {
			dcomp* y_m			= new dcomp [n_states];
			dcomp* f_m			= new dcomp [n_states];

			d_derv(f[0], path_constraints[0], states[0], controls[0], parameters, t[0], 1, constants);

			for (Index i = 0; i < n_nodes - 1; i += 1)	{
				d_derv(f[i+1], path_constraints[2*(i+1)], states[i+1], controls[2*(i+1)], parameters, t[i+1], 1, constants);
				for (Index j = 0; j < n_states; j += 1){
					y_m[j] 	= (states[i][j]+states[i+1][j])/2. + delta[i]/8.*(f[i][j]-f[i+1][j]);
				}
				d_derv(f_m, path_constraints[2*i+1], y_m, controls[2*i+1], parameters, (t[i] + t[i+1])/2., 1, constants);
				for (Index j = 0; j < n_states; j++)
					defects[i][j]	= (states[i+1][j] - states[i][j] - delta[i]/6.*(f[i][j]+4.*f_m[j]+f[i+1][j]));
			}

			delete[] y_m;
			delete[] f_m;

		}
		else if (config.disc_method == trapezoidal) {
#pragma omp parallel for
			for (Index i = 0; i < n_nodes; i += 1)	{
				d_derv(f[i], path_constraints[i], states[i], controls[i], parameters, t[i], 1, constants);
			}

			for (Index i = 0; i < n_nodes - 1; i++)
				for (Index j = 0; j < n_states; j++)
					defects[i][j] 	= states[i+1][j] - states[i][j] - delta[i]/2.0*(f[i][j] + f[i+1][j]);
		}

/*
		for (int i = 0; i < n_nodes; i++) {
			for (int j = 0; j < n_path_constraints; j++) {
				printf("path[%d][%d] = (%.2e,%.2e)\t",i,j,path_constraints[i][j].real(),path_constraints[i][j].imag());
			}
			printf("\n");
		}

		for (int i = 0; i < n_nodes-1; i++) {
			for (int j = 0; j < n_states; j++) {
				printf("defects[%d][%d] = (%.2e,%.2e)\t",i,j,defects[i][j].real(),defects[i][j].imag());
			}
			printf("\n");
		}

		for (int i = 0; i < n_events; i++) {
			printf("events[%d] = (%.2e,%.2e)\n",i,events[i].real(),events[i].imag());

		}
//*/
		OCP_var_2_NLP_g(path_constraints, defects, events, g);

/*		for (int i = 0; i < nlp_m; i++)
			printf("nlp_g[%d] = (%.2e,%.2e)\n",i,g[i].real(),g[i].imag());
//*/
		for (Index i = 0; i < nlp_m; i++) {
			g[i] = g[i] + dcomp(1,0);
		}

		for (Index i = 0; i < n_nodes; i++) {
			delete[] f[i];
		}
		delete[] f;
		delete[] t;
		delete[] delta;
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
	cout<<"n = "<<n<<"\nm = "<<m<<endl;
	pattern_gen(n, m, nnz_jac_g, nnz_h_lag);
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

#ifdef DCOPT_DEBUG
	cout<<"get_bounds_info\n";
#endif

	for (Index i = 0; i < nlp_n; i += 1) {
		x_l[i] = nlp_lb_x[i];
		x_u[i] = nlp_ub_x[i];
	}

	for (Index i = 0; i < nlp_m; i += 1) {
		g_l[i] = nlp_lb_g[i]+1;
		g_u[i] = nlp_ub_g[i]+1;
	}

#ifdef DCOPT_DEBUG
	cout<<"end get_bounds_info\n";
#endif
	return true;
}

bool MyADOLC_sparseNLP::get_starting_point(Index n, bool init_x, Number* x,
                               bool init_z, Number* z_L, Number* z_U,
                               Index m, bool init_lambda,
                               Number* lambda)
{

#ifdef DCOPT_DEBUG
	cout<<"get_staring_point\n";
#endif
	assert(init_x == true);
	assert(init_z == false);
	assert(init_lambda == false);

	for (Index i = 0; i < n; i += 1) {
		x[i]	= nlp_guess_x[i];
	}

#ifdef DCOPT_DEBUG
	cout<<"end get_staring_point\n";
#endif
	return true;
}

bool MyADOLC_sparseNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
#ifdef DCOPT_DEBUG
	cout<<"eval_f\n";
#endif
	dcomp obj;
	dcomp *xp	= new dcomp[n];
	for (int i = 0; i < n; i++)
		xp[i]	= x[i];
	eval_obj(n,xp,obj);
	obj_value = obj.real();

	delete[] xp;
#ifdef DCOPT_DEBUG
	cout<<"end of eval_f\n";
#endif
	return true;
}

bool MyADOLC_sparseNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
#ifdef DCOPT_DEBUG
	cout<<"eval_grad_f\n";
#endif
	dcomp obj;
	dcomp *xp	= new dcomp[n];
	for (int i = 0; i < n; i++)
		xp[i]	= x[i];

	for (int i = 0; i < n; i++) {
		xp[i] 	+= dcomp (0,step);
		eval_obj(n,xp,obj);
		xp[i] 	= x[i];
		grad_f[i]	= obj.imag()/step;
	}
	delete[] xp;

#ifdef DCOPT_DEBUG
	cout<<"end of eval_grad_f\n";
#endif
	return true;
}

bool MyADOLC_sparseNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
#ifdef DCOPT_DEBUG
	cout<<"eval_g\n";
#endif
	dcomp *g_comp	= new dcomp[m];
	dcomp *x_comp	= new dcomp[n];
	for (int i = 0; i < n; i++)
		x_comp[i]	= x[i];
	eval_constraints(n,x_comp,m,g_comp);

	for (Index i = 0; i < m; i++)
		g[i]	= g_comp[i].real();

	delete[] g_comp;
	delete[] x_comp;
#ifdef DCOPT_DEBUG
	cout<<"end of eval_grad_g\n";
#endif
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
		for(Index idx = 0; idx < nnz_jac; idx++) {
			iRow[idx] = rind_g[idx];
			jCol[idx] = cind_g[idx];
		}
	}
	else {

		dcomp *x_comp	= new dcomp[n];

		for (int i = 0; i < n; i++)
			x_comp[i]	= x[i];

		int nnz = 0;
//#pragma omp parallel
	{
		dcomp *g_comp	= new dcomp[m];
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

	results.z.resize(nlp_n,1);
	results.x.resize(n_nodes, n_states);
	results.u.resize(n_nodes, n_controls);
	if (config.disc_method == Hermite_Simpson)
		results.u_full.resize(2*n_nodes-1, n_controls);

	results.parameters.resize(n_parameters,1);
	results.nodes.resize(n_nodes, 1);

	dcomp *x_comp		= new dcomp[n];

	for (Index i = 0; i < nlp_n; i++) {
		results.z(i+1) 	= x[i];
		x_comp[i]		= x[i];
	}

	NLP_x_2_OCP_var(x_comp, states, controls, parameters, t0, tf);

	for (Index i = 0; i < n_nodes; i++)
		for (Index j = 0; j < n_states; j++)
			results.x(i+1,j+1) 		= states[i][j].real();

	for (Index i = 0; i < n_nodes; i++)
		for (Index j = 0; j < n_controls; j++){
			if (config.disc_method == Hermite_Simpson){
				results.u(i+1,j+1) 		= controls[2*i][j].real();
				results.u_full	(2*i+1,j+1) 		= controls[2*i][j].real();
				if (i < n_nodes - 1)
					results.u_full(2*i+2,j+1) 		= controls[2*i+1][j].real();
			}
			else
				results.u(i+1,j+1) 		= controls[i][j].real();
		}

	for (Index i = 0; i < n_parameters; i++)
		results.parameters(i+1)		= parameters[i].real();

	if (n_nodes && n_states) {
		for (Index i = 2; i < n_nodes; i++) {
			results.nodes(i,1) 	= results.nodes(i-1,1) + (tf - t0).real()*node_str(i-1);
		}
		results.nodes(n_nodes,1) = tf.real();
	}

	if (n_nodes && n_states) {

		for (Index i = 0; i < n_nodes; i++) {
			delete[] state_idx[i];
			delete[] states[i];
			if (i< n_nodes - 1) {
				delete[] defect_idx[i];
				delete[] defects[i];
			}
		}

		if (config.disc_method == Hermite_Simpson) {
			for (Index i = 0; i < 2*n_nodes - 1; i++) {
				delete[] control_idx[i];
				delete[] controls[i];
				delete[] path_constraint_idx[i];
				delete[] path_constraints[i];
			}
		}
		else {
			for (Index i = 0; i < n_nodes; i++) {
				delete[] control_idx[i];
				delete[] controls[i];
				delete[] path_constraint_idx[i];
				delete[] path_constraints[i];
			}
		}

		delete[] state_idx;
		delete[] states;

		delete[] defect_idx;
		delete[] defects;

		delete[] control_idx;
		delete[] controls;

		delete[] path_constraint_idx;
		delete[] path_constraints;
	}

	if (n_parameters) {
		delete[] parameters;
		delete[] parameter_idx;
	}

	if (n_events) {
		delete[] events;
		delete[] event_idx;
	}

	delete[] x_lam;
	delete[] rind_g;
	delete[] cind_g;
	delete[] jacval;


		if (!config.H_approximation) {
			delete[] (rind_L);
			delete[] (cind_L);
			free(cind_L_total);
			free(rind_L_total);
			free(hessval);
		}

	delete[] nlp_sf_x;
	delete[] nlp_sf_g;
	delete[] nlp_guess_x;
	delete[] nlp_lb_x;
	delete[] nlp_ub_x;
	delete[] nlp_lb_g;
	delete[] nlp_ub_g;

}

void MyADOLC_sparseNLP::pattern_gen(Index n, Index m, Index& nnz_jac_g, Index& nnz_h_lag)
{
#ifdef DCOPT_DEBUG
	printf("generate_tapes()\n");
#endif
	double *xp		= new double[n];
	double *lamp	= new double[m];
	double *zl		= new double[m];
	double *zu		= new double[m];

	dcomp *lam  	= new dcomp[m];
	x_lam   		= new double[n+m+1];

	dcomp sig;
	dcomp obj_value;
	nnz_jac 		= 0;
	nnz_jac_g		= 0;

	get_starting_point(n, 1, xp, 0, zl, zu, m, 0, lamp);

	time(&rand_time);
	srand((unsigned int)rand_time);              /* Zufallsgenerator initialisieren */
	for(Index i = 0; i < n; i++) {
		xp[i] = (double) rand()/(RAND_MAX);//* nlp_sf_x[i];
//		printf("x_comp[%d] = %.2e\n",i, x_comp[i]);
	}
//#pragma omp parallel
{
	int nnz			= 0;
	dcomp *x_comp	= new dcomp[n];
	for(Index i = 0; i < n; i++) {
		x_comp[i]	= xp[i];
	}
	dcomp *g_comp	= new dcomp[m];
	int this_thread = omp_get_thread_num(), num_threads = omp_get_num_threads();
	for (int col = this_thread; col < n; col+=num_threads) {
		x_comp[col] 	+= dcomp (0,step);
		eval_constraints(n,	x_comp, m, g_comp);
		x_comp[col] 	= xp[col];
		for (int row = 0; row < m; row++) {
			if(g_comp[row].imag() != 0){
				nnz++;
			}
		}
	}

#pragma omp atomic
	nnz_jac	+= nnz;

	delete[] g_comp;
	delete[] x_comp;
}

	cind_g		= new uint [nnz_jac];
	rind_g		= new uint [nnz_jac];
	jacval		= new double [nnz_jac];
	cout<<"nnz_jac = "<<nnz_jac<<endl;
/*
	for (int row = 0; row < m; row++) {
		printf("row[%d]\t",row);
		for (int col = 0; col < n; col++) {
			printf("%.1e ",val[row][col]);
		}
		cout<<endl;
	}
//*/
	time(&rand_time);
	srand((unsigned int)rand_time);              /* Zufallsgenerator initialisieren */
	for(Index i = 0; i < n; i++) {
		xp[i] = (double) rand()/(RAND_MAX);//* nlp_sf_x[i];
//		printf("x_comp[%d] = %.2e\n",i, x_comp[i]);
	}
/*	exit(1);
	for(Index i = 0; i < n; i++) {
		x_comp[i] = xp[i];// = (double) rand()/(RAND_MAX);
	}
//*/
//#pragma omp parallel shared(val)
{
	dcomp *x_comp	= new dcomp[n];
	for(Index i = 0; i < n; i++) {
		x_comp[i]	= xp[i];
	}
	dcomp *g_comp	= new dcomp[m];
	int this_thread = omp_get_thread_num(), num_threads = omp_get_num_threads();
	for (int col = this_thread; col < n; col+=num_threads) {
		x_comp[col] 	+= dcomp (0,step);
		eval_constraints(n,	x_comp, m, g_comp);
		x_comp[col] 	= xp[col];
		for (int row = 0; row < m; row++) {
			if (g_comp[row].imag() != 0.) {
				cind_g[nnz_jac_g]		= col;
				rind_g[nnz_jac_g]		= row;
				#pragma omp atomic
				nnz_jac_g++;
			}
		}
	}
	delete[] g_comp;
	delete[] x_comp;
}
/*
	cout<<"nnz_jac_g = "<<nnz_jac_g<<endl;
	for (int row = 0; row < m; row++) {
		printf("row[%d]\t",row);
		for (int col = 0; col < n; col++) {
			printf("%.1e ",val[row][col]);
		}
		cout<<endl;
	}
	//*/
	if (nnz_jac == nnz_jac_g) {
		nnz_jac_g = nnz_jac;
		cout<<"yeah\n";
	}
	else{
		cout<<"no proper jac_g found\n";
		cout<<"1 = "<<nnz_jac_g<<endl<<"2 = "<<nnz_jac<<endl;
		exit(1);
	}
	nnz_h_lag = 0;

/*
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
*/

	delete[] lam;
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
