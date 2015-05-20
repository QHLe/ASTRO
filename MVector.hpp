/*
 * MVector.hpp
 *
 *  Created on: May 14, 2015
 *      Author: zineus
 */

#ifndef MULTIVECTOR_HPP_
#define MULTIVECTOR_HPP_

#include<stdlib.h>
#include "SVector.hpp"

using namespace std;

class SVector;

class MVector {
	friend class SVector;
public:
	MVector();
	MVector(uint v_length, uint v_num);
	virtual ~MVector();
	void 	Print 		()											const;
	void 	Print		(string line) 								const;
	void 	resize		(uint v_length, uint v_num);
	SVector 	getSVector	(uint num)									const;
	SVector 	getRow		(uint num)									const;
	void 	load		(const char* filename);
	void 	save		(const char* filename)						const;
	const MVector& operator= (const MVector& rhs);
	const MVector& operator= (const SVector& rhs);
	MVector operator() 	(uint l_limit, uint u_limit)				const;
	MVector	operator,	(const SVector& rhs)							const;
	MVector	operator,	(const MVector& rhs)						const;
private:
	uint vector_num;
	uint vector_length;
	SVector* mvector;
};

#endif /* MULTIVECTOR_HPP_ */
