/*
 * MVector.hpp
 *
 *  Created on: Jun 5, 2015
 *      Author: quanghoale
 */



#ifndef MVECTOR_HPP_
#define MVECTOR_HPP_

//#include "SVector.hpp"
template<class T>
class SVector;

class MVector {
public:
	MVector();
	MVector(uint v_length, uint v_num);
	virtual ~MVector();
	void 	Print 		()											const;
	void 	Print		(const char* line) 								const;
	void 	resize		(uint v_length, uint v_num);
	SVector<double> 	getSVector	(uint num)							const;
	SVector<double> 	getRow		(uint num)							const;
	uint 	getRowDim	()	const	{return vector_length;}
	uint 	getColDim	()	const 	{return vector_num;}
	void 	load		(const char* filename);
	void 	save		(const char* filename)						const;

	const MVector& operator= (const MVector& rhs);
	const MVector& operator= (const SVector<double>& rhs);

	double operator() 	(uint row, uint col)						const;
	double& operator() 	(uint row, uint col);

	MVector truncate_col	(uint l_limit, uint u_limit);
	MVector truncate_row	(uint l_limit, uint u_limit);


	MVector	operator,	(const SVector<double>& rhs)						const;
	MVector	operator,	(const MVector& rhs)						const;
	friend double min	(const MVector& mv);
	friend double max	(const MVector& mv);

private:
	uint vector_num;
	uint vector_length;
	SVector<double>* mvector;
};

//void Print(const SVector<double>);
//void Print(const SVector<double>, const char*);


#endif /* MVECTOR_HPP_ */
