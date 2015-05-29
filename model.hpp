/*
 * cpp_example.hpp
 *
 *  Created on: May 17, 2015
 *      Author: zineus
 */

#define N_NODES		1000



#ifndef CPP_EXAMPLE_HPP_
#define CPP_EXAMPLE_HPP_

template<class T> T lagrange_cost(	const T *states,
									const T *controls,
									const T *param,
									const T &time,
									uint phase) {
	return 0;
//	return controls[0]*controls[0]*0.5;
}
template<class T> T endpoint_cost (	const T* ini_states,
									const T* fin_states,
									const T* param,
									const T& t0,
									const T& tf,
									uint phase) {

	return tf;//fin_states[2];
}

template<class T> void derivatives(	T *states_dot,
							 	 	const T *states,
							 	 	const T *controls,
							 	 	const T *param,
							 	 	const T &time,
							 	 	uint phase) {

	T x1 = states[0];
	T x2 = states[1];
//	T x3 = states[2];

	T u = controls[0];

	states_dot[0] 	= x2;
	states_dot[1]	= u;
//	states_dot[2]	= u*u/2;
}

template<class T> void events(	T *events,
								const T *ini_states,
								const T *fin_states,
								const T *param,
								const T &t0,
								const T &tf,
								uint phase) {

/*	events [0]	= ini_states[0];
	events [1]	= ini_states[1];
	events [2]	= ini_states[2];
	events [3]	= fin_states[0];
	events [4]	= fin_states[1];*/

	events [0]	= ini_states[0];
	events [1]	= ini_states[1];
	events [2]	= fin_states[0];
	events [3]	= fin_states[1];

}

#endif /* CPP_EXAMPLE_HPP_ */
