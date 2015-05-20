/*
 * cpp_example.hpp
 *
 *  Created on: May 17, 2015
 *      Author: zineus
 */

#ifndef CPP_EXAMPLE_HPP_
#define CPP_EXAMPLE_HPP_

template<class T> T lagrange_cost ( const double* states,
									const double* controls,
									const double* param,
									const double& time,
									uint phase);

template<class T> T endpoint_cost (	const double* ini_states,
									const double* fin_states,
									const double* param,
									const double& t0,
									const double& tf,
									uint phase);
/*
template <class T> void events(		T* e,
									const T *ini_states,
									const T *fin_states,
									const T *param);
*/


#endif /* CPP_EXAMPLE_HPP_ */
