/*
 * SVector.h
 *
 *  Created on: May 13, 2015
 *      Author: zineus
 */

#ifndef SVECTOR_H_
#define SVECTOR_H_

#include<stdlib.h>
#include "MVector.hpp"

using namespace std;

class MVector;

class SVector {
	friend class MVector;
public:
	SVector();
	SVector(const SVector& cp_SVector);
	SVector(uint size);
	virtual ~SVector();
	void 	resize		(uint n_new);
	uint 	getsize		() 									const	{return n;}
	double 	getelement	(uint num) 							const;
	void 	setelement	(uint num, double value_new);
	void 	Print		() 									const;
	void 	Print		(string line) 						const;

	double 	enorm		();
	void 	load		(const char* filename);
	void 	save		(const char* filename)				const;

	const SVector& operator= (const SVector& rhs);
	double operator() (uint num)						const;
	SVector operator() (uint l_limit, uint u_limit)		const;

	/* elementary operations */
	SVector 		operator+ 	(const SVector& rhs) 			const;
	SVector 		operator- 	(const SVector& rhs)				const;
	SVector 		operator* 	(const SVector& rhs)				const;
	SVector 		operator/ 	(const SVector& rhs)				const;
	friend SVector sqrt		(const SVector& vector);
	friend SVector sin		(const SVector& vector);
	friend SVector cos		(const SVector& vector);
	friend SVector tan		(const SVector& vector);
	friend SVector atan2		(const SVector&, const SVector&);

	/* operations with a scalar */
	SVector 		operator+ 	(double sum)					const;
	SVector 		operator- 	(double sum)					const;
	SVector 		operator* 	(double factor)					const;
	SVector 		operator/ 	(double factor)					const;

	SVector 		operator- 	()								const;

	/* concatenation operation */
	MVector	operator,	(const SVector& rhs)					const;
	MVector	operator,	(const MVector& rhs)				const;

	friend SVector operator+ (double, const SVector&);
	friend SVector operator- (double, const SVector&);
	friend SVector operator* (double, const SVector&);

	/* friend functions */

	friend double dot 			(const SVector& v1, const SVector& v2);
	friend SVector cross		(const SVector& v1, const SVector& v2);
	friend SVector ones 		(uint num);
	friend SVector zeros		(uint num);
	friend SVector linspace		(double x0, double xn, uint n);


private:
	double* values;
	uint n;
};

double 		dot 		(const SVector& v1, const SVector& v2);
SVector 	cross 		(const SVector& v1, const SVector& v2);
SVector 	ones 		(uint num);
SVector 	zeros 		(uint num);
SVector 	linspace	(double x0, double xn, uint n);
#endif /* VECTOR_H_ */
