/*
 * Matrix.h
 *
 *  Created on: Aug 11, 2015
 *      Author: zineus
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <adolc/adouble.h>

using namespace std;

template <class T>
class Matrix {};

template <>
class Matrix<adouble>;

template <>
class Matrix<double> {
public:
	Matrix();
	Matrix(const Matrix<double>& cp_Matrix);
	Matrix(const Matrix<adouble>& cp_Matrix);
	Matrix(uint v_length, uint v_num);
	virtual ~Matrix();
	void 			Print 		()									const;
	void 			Print		(const char* line) 					const;
	void 			resize		(uint n_row, uint n_col);
	uint 			getRowDim	()	const	{return n_row;}
	uint 			getColDim	()	const 	{return n_col;}

	double 			operator()	(uint row, uint col)				const;
	double&			operator()	(uint row, uint col);
	double 			operator()	(uint idx) 							const;
	double&			operator()	(uint idx);

	Matrix<double>  getCol		(uint num)							const;
	Matrix<double>	getRow		(uint num)							const;

	void 			load		(const char* filename);
	void 			save		(const char* filename)				const;

	const Matrix<double>& 	operator= (const Matrix<double>& rhs);
	const Matrix<double>& 	operator= (const double scalar);
	Matrix<double>			operator,	(const Matrix<double>& rhs)						const;

	Matrix<double> 	truncate_col	(uint l_limit, uint u_limit);
	Matrix<double> 	truncate_row	(uint l_limit, uint u_limit);

	Matrix<double> 	operator+ 	(double sum)						const;
	Matrix<double> 	operator- 	(double sum)						const;
	Matrix<double> 	operator* 	(double factor)						const;
	Matrix<double> 	operator/ 	(double factor)						const;
	Matrix<adouble>	operator+ 	(adouble sum)						const;
	Matrix<adouble>	operator- 	(adouble sum)						const;
	Matrix<adouble>	operator* 	(adouble factor)					const;
	Matrix<adouble>	operator/ 	(adouble factor)					const;
	Matrix<double> 	operator- 	()									const;

	Matrix<double> 	operator+ 	(const Matrix<double>& rhs) 		const;
	Matrix<double> 	operator- 	(const Matrix<double>& rhs)			const;
	Matrix<double> 	operator* 	(const Matrix<double>& rhs)			const;
	Matrix<double> 	operator/ 	(const Matrix<double>& rhs)			const;
	Matrix<adouble> operator+ 	(const Matrix<adouble>& rhs) 		const;
	Matrix<adouble> operator- 	(const Matrix<adouble>& rhs)		const;
	Matrix<adouble> operator* 	(const Matrix<adouble>& rhs)		const;
	Matrix<adouble> operator/ 	(const Matrix<adouble>& rhs)		const;
	double 			enorm		();

private:
	uint n_col;
	uint n_row;
	double** data;
};

template<>
class Matrix <adouble>{
public:
	Matrix();
	Matrix(const Matrix<adouble>& cp_Matrix);
	Matrix(const Matrix<double>& cp_Matrix);
	Matrix(uint v_length, uint v_num);
	virtual ~Matrix();
	void 			Print 		()									const;
	void 			Print		(const char* line) 					const;
	void 			resize		(uint n_row, uint n_col);
	uint 			getRowDim	()	const	{return n_row;}
	uint 			getColDim	()	const 	{return n_col;}

	adouble 		operator()	(uint row, uint col)				const;
	adouble& 		operator()	(uint row, uint col);
	adouble 		operator()	(uint idx) 							const;
	adouble& 		operator()	(uint idx);

	Matrix<adouble> getCol		(uint num)							const;
	Matrix<adouble> getRow		(uint num)							const;

	void 			load		(const char* filename);
	void 			save		(const char* filename)				const;

	const Matrix& 	operator= (const Matrix& rhs);
	const Matrix& 	operator= (const adouble scalar);
	const Matrix& 	operator= (const Matrix<double>& rhs);
	const Matrix&	operator= (const double scalar);

	Matrix			operator,	(const Matrix& rhs)					const;

	Matrix	 		truncate_col	(uint l_limit, uint u_limit);
	Matrix 			truncate_row	(uint l_limit, uint u_limit);

	Matrix 			operator+ 	(adouble sum)				const;
	Matrix 			operator- 	(adouble sum)				const;
	Matrix 			operator* 	(adouble factor)			const;
	Matrix 			operator/ 	(adouble factor)			const;

	Matrix 			operator+ 	(double sum)				const;
	Matrix 			operator- 	(double sum)				const;
	Matrix 			operator* 	(double factor)				const;
	Matrix 			operator/ 	(double factor)				const;

	Matrix 			operator- 	()							const;

	Matrix 			operator+ 	(const Matrix& rhs) 		const;
	Matrix 			operator- 	(const Matrix& rhs)			const;
	Matrix 			operator* 	(const Matrix& rhs)			const;
	Matrix 			operator/ 	(const Matrix& rhs)			const;

	Matrix 			operator+ 	(const Matrix<double>& rhs) const;
	Matrix 			operator- 	(const Matrix<double>& rhs)	const;
	Matrix 			operator* 	(const Matrix<double>& rhs)	const;
	Matrix 			operator/ 	(const Matrix<double>& rhs)	const;

	adouble 		enorm		();

private:
	uint n_col;
	uint n_row;
	adouble** data;
};


inline Matrix<double>::Matrix() {
	n_col = 0;
	n_row = 0;
	data = NULL;
}

inline Matrix<double>::Matrix(uint m_length, uint m_num) {
	n_col = m_num;
	n_row = m_length;
	if (n_col == 0 || n_row == 0) {
		cout<<"WARNING: column or row size is zero!\n";
		data = NULL;
	}
	else {
		data 	= new double*[n_row];
		for (uint i = 0; i < n_row; i ++)
			data[i] 	= new double[n_col];

		for (uint i = 0; i < n_row; i ++)
			for (uint j = 0; j < n_col; j++)
				data[i][j] = 0.0;
	}
}

inline Matrix<double>::Matrix(const Matrix<double>& cp_Matrix){
	n_col = cp_Matrix.getColDim();
	n_row = cp_Matrix.getRowDim();

	data 	= new double*[n_row];
	for (uint i = 0; i < n_row; i ++)
		data[i] 	= new double[n_col];

	for (uint i = 0; i < n_row; i++)
		for (uint j = 0; j < n_col; j++)
			data[i][j] = cp_Matrix(i+1,j+1);
}

inline Matrix<double>::Matrix(const Matrix<adouble>& cp_Matrix){
	n_col = cp_Matrix.getColDim();
	n_row = cp_Matrix.getRowDim();

	data 	= new double*[n_row];
	for (uint i = 0; i < n_row; i ++)
		data[i] 	= new double[n_col];

	for (uint i = 0; i < n_row; i++)
		for (uint j = 0; j < n_col; j++)
			data[i][j] = cp_Matrix(i+1,j+1).getValue();
}

inline Matrix<double>::~Matrix() {
	if (this->data != NULL) {
		for ( uint i=0; i < n_row; i++) {
			delete[] data[i];
		}
		delete[] data;
	}
}

inline void Matrix<double>::Print() const {
	for (uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++)
			printf("%+.4e\t",data[i][j]);
		printf("\n");
	}
}

inline void Matrix<double>::Print(const char* printline) const {
	printf("\n*** %s ***\n\n",printline);
	for (uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++)
			printf("%+.4e\t",data[i][j]);
		printf("\n");

	}
}

inline void Matrix<double>::resize(uint row_num, uint col_num) {
	Matrix <double> temp_matrix(*this);
	if (this->data != NULL) {
		for (uint i = 0; i < n_row; i++)
			delete[] data[i];
		delete[] data;
	}
	n_row 	= row_num;
	n_col 	= col_num;
	data 	= new double* [row_num];
	for (uint i = 0; i<row_num; i++)
		data[i]	= new double[col_num];

	for (uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++) {
			if (i < temp_matrix.getRowDim() && j < temp_matrix.getColDim())
				data[i][j]	= temp_matrix.data[i][j];
			else
				data[i][j] 	= 0.0;
		}
	}
}

inline double Matrix<double>::operator ()(uint row, uint col) const {
	if (col > n_col || row > n_row) {
		cout<<"===== Overloading operator (row, col) =====\n"
			  "=====          invalid range          =====\n";
		exit(1);
	}
	if (col < 1 || row < 1) {
		cout<<"===== Overloading operator (row, col) =====\n"
			  "=====       Index starts from 1       =====\n";
		exit(1);
	}
	return data[row-1][col-1];
}

inline double& Matrix<double>::operator ()(uint row, uint col) {
	if (col > n_col || row > n_row) {
		cout<<"===== Overloading operator (row, col) =====\n"
			  "=====          invalid range          =====\n";
		exit(1);
	}
	if (col < 1 || row < 1) {
		cout<<"===== Overloading operator (row, col) =====\n"
			  "=====       Index starts from 1       =====\n";
		exit(1);
	}
	return data[row-1][col-1];
}

inline double Matrix<double>::operator()(uint idx) const{
	if (n_col > 1 && n_row > 1) {
		cout<<"=====    Overloading operator (idx)     =====\n"
			  "===== Matrix is not a row or col vector =====\n";
		exit(1);
	}
	if (n_col == 0 || n_row == 0) {
		cout<<"=====  Overloading operator (idx)   =====\n"
			  "===== Matrix is not initialized yet =====\n";
		exit(1);
	}
	if (idx < 1 ) {
		cout<<"===== Overloading operator (idx) =====\n"
			  "=====     Index starts from 1    =====\n";
		exit(1);
	}
	if (n_col == 1)
		return data[idx-1][1];
	else
		return data[1][idx-1];
}

inline double& Matrix<double>::operator()(uint idx) {
	if (n_col > 1 && n_row > 1) {
		cout<<"=====    Overloading operator (idx)     =====\n"
			  "===== Matrix is not a row or col vector =====\n";
		exit(1);
	}
	if (n_col == 0 || n_row == 0) {
		cout<<"=====  Overloading operator (idx)   =====\n"
			  "===== Matrix is not initialized yet =====\n";
		exit(1);
	}
	if (idx < 1 ) {
		cout<<"===== Overloading operator (idx) =====\n"
			  "=====     Index starts from 1    =====\n";
		exit(1);
	}
	if (n_col == 1)
		return data[idx-1][1];
	else
		return data[1][idx-1];
}

inline void Matrix<double>::load (const char* filename) {

	FILE *fp;
	if ( (fp = fopen(filename,"r")) == NULL ) {
		printf("===== Error opening file to load Matrix =====");
	}
	for (uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++) {
			fscanf(fp,"%le", &data[i][j]);
		}
	}
	fclose(fp);
}

inline void Matrix<double>::save (const char* filename) const {
	FILE *fp;
	if ( (fp = fopen(filename,"w")) == NULL ) {
		printf("===== Error opening file to save Matrix =====\n");
	}

	for (uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++) {
			fprintf(fp,"%.8e\t",data[i][j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

inline const Matrix<double>& Matrix<double>::operator= (const Matrix<double>& rhs) {
	this->resize(rhs.getRowDim(),rhs.getColDim());
	for(uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++)
			data[i][j] = rhs(i+1,j+1);
	}
	return *this;
}

inline const Matrix<double>& Matrix<double>::operator= (const double scalar) {
	this->resize(1,1);
	data[0][0] = scalar;
	return *this;
}

inline Matrix<double> Matrix<double>::operator, (const Matrix<double>& rhs) const{
	if (rhs.getRowDim() != n_row) {
		cout<<"===== Error in matrix concatenation  =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);
	}
	Matrix<double> temp(n_row, rhs.getColDim() + n_col);

	for (uint i = 1; i <= n_row; i++) {
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = this->operator ()(i,j);
	}

	for (uint i = 1; i <= n_row; i++) {
		for (uint j = 1; j <= rhs.getColDim(); j++)
			temp(i,j+n_col) = rhs(i,j);
	}
	return temp;
}

inline Matrix<double> Matrix<double>::truncate_col	(uint l_limit, uint u_limit) {
	if (l_limit > u_limit){
		cout<<"===== Truncating columns of matrix =====\n"
			  "=====  lower limit > upper limit   =====\n";
		exit(1);
	}
	if (u_limit > n_col) {
		cout<<"===== Truncating columns of matrix =====\n"
			  "=====        invalid range         =====\n";
		exit(1);
	}
	if (l_limit < 1) {
		cout<<"===== Truncating columns of matrix =====\n"
			  "=====     Index starts from 1      =====\n";
		exit(1);
	}
	Matrix<double> temp(n_row, u_limit - l_limit + 1);
	for (uint i = 1; i <= n_row; i++) {
		for (uint j = l_limit; j <= u_limit; j++) {
			temp(i,j - l_limit + 1) = this->operator()(i,j);
		}
	}
	return temp;
}

inline Matrix<double> Matrix<double>::truncate_row	(uint l_limit, uint u_limit) {
	if (l_limit > u_limit){
		cout<<"===== Truncating rows of matrix =====\n"
			  "===== lower limit > upper limit =====\n";
		exit(1);
	}
	if (u_limit > n_col) {
		cout<<"===== Truncating rows of matrix =====\n"
			  "=====       invalid range       =====\n";
		exit(1);
	}
	if (l_limit < 1) {
		cout<<"===== Truncating rows of matrix =====\n"
			  "=====    Index starts from 1    =====\n";
		exit(1);
	}
	Matrix<double> temp(u_limit - l_limit + 1,n_col);
	for (uint i = l_limit; i <= u_limit; i++) {
		for (uint j = 1; j <= n_col; j++) {
			temp(i - l_limit + 1,j) = this->operator()(i,j);
		}
	}
	return temp;
}

inline double Matrix<double>::enorm() {
	if (n_row > 1 && n_col > 1) {
		cout<<"=====			  !!! WARNING !!!			   =====\n"
			  "=====  matrix is neither column nor row vector  =====\n"
			  "===== enorm of a matrix does not make any sense =====\n";
	}
	double temp = 0;
	for(uint i = 0; i < n_row; i++)
		for (uint j = 0; j < n_col; j++)
			temp += data[i][j]*data[i][j];
	return sqrt(temp);
}

inline Matrix<double> Matrix<double>::operator+ (double sum) const{
	Matrix<double> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) += sum;
	return temp;
}

inline Matrix<double> Matrix<double>::operator- (double sum) const{
	Matrix<double> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) -= sum;
	return temp;
}

inline Matrix<double> Matrix<double>::operator* (double factor) const{
	Matrix<double> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) *= factor;
	return temp;
}

inline Matrix<double> Matrix<double>::operator/ (double factor) const{
	Matrix<double> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) /= factor;
	return temp;
}

inline Matrix<double> Matrix<double>::operator- () const{
	Matrix<double> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) *= -1;
	return temp;
}

inline Matrix<double> Matrix<double>::operator+ (const Matrix<double>& rhs) const{
	Matrix<double> temp(*this);
	if (rhs.getColDim() != n_col || rhs.getRowDim() != n_row) {
		cout<<"=====     Element-wise addition      =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);
	}
	for( uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j) + rhs(i,j);
	return temp;
}

inline Matrix<double> Matrix<double>::operator- (const Matrix<double>& rhs) const{
	Matrix<double> temp(*this);
	if (rhs.getColDim() != n_col || rhs.getRowDim() != n_row) {
		cout<<"=====   Element-wise substitution    =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);	}
	for( uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j) - rhs(i,j);
	return temp;
}

inline Matrix<double> Matrix<double>::operator* (const Matrix<double>& rhs) const{
	Matrix<double> temp(*this);
	if (rhs.getRowDim() != n_row || rhs.getColDim() != n_col) {
		cout<<"=====  Element-wise multiplication   =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);
	}
	if (n_row > 1 && n_col > 1) {
		cout<<"=====             !!! WARNING !!!             =====\n"
			  "===== matrix is neither column nor row vector =====\n"
			  "===== Element-wise multiplication of a matrix =====\n"
			  "=====    Do you know what you are doing???    =====\n";
	}
	for (uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j)*rhs(i,j);
	return temp;
}

inline Matrix<double> Matrix<double>::operator/ (const Matrix<double>& rhs) const{
	Matrix<double> temp(*this);
	if (rhs.getRowDim() != n_row || rhs.getColDim() != n_col) {
		cout<<"=====     Element-wise division      =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);
	}
	if (n_row > 1 && n_col > 1) {
		cout<<"=====             !!! WARNING !!!             =====\n"
			  "===== matrix is neither column nor row vector =====\n"
			  "=====    Element-wise division of a matrix    =====\n"
			  "=====    Do you know what you are doing???    =====\n";
	}
	for (uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j)/rhs(i,j);
	return temp;
}

inline Matrix<adouble> Matrix<double>::operator+ (const Matrix<adouble>& rhs) const{
	Matrix<adouble> temp(*this);
	if (rhs.getColDim() != n_col || rhs.getRowDim() != n_row) {
		cout<<"=====     Element-wise addition      =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);
	}
	for( uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j) + rhs(i,j);
	return temp;
}

inline Matrix<adouble> Matrix<double>::operator- (const Matrix<adouble>& rhs) const{
	Matrix<adouble> temp(*this);
	if (rhs.getColDim() != n_col || rhs.getRowDim() != n_row) {
		cout<<"=====   Element-wise substitution    =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);	}
	for( uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j) - rhs(i,j);
	return temp;
}

inline Matrix<adouble> Matrix<double>::operator* (const Matrix<adouble>& rhs) const{
	Matrix<adouble> temp(*this);
	if (rhs.getRowDim() != n_row || rhs.getColDim() != n_col) {
		cout<<"=====  Element-wise multiplication   =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);
	}
	if (n_row > 1 && n_col > 1) {
		cout<<"=====             !!! WARNING !!!             =====\n"
			  "===== matrix is neither column nor row vector =====\n"
			  "===== Element-wise multiplication of a matrix =====\n"
			  "=====    Do you know what you are doing???    =====\n";
	}
	for (uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j)*rhs(i,j);
	return temp;
}

inline Matrix<adouble> Matrix<double>::operator/ (const Matrix<adouble>& rhs) const{
	Matrix<adouble> temp(*this);
	if (rhs.getRowDim() != n_row || rhs.getColDim() != n_col) {
		cout<<"=====     Element-wise division      =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);
	}
	if (n_row > 1 && n_col > 1) {
		cout<<"=====             !!! WARNING !!!             =====\n"
			  "===== matrix is neither column nor row vector =====\n"
			  "=====    Element-wise division of a matrix    =====\n"
			  "=====    Do you know what you are doing???    =====\n";
	}
	for (uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j)/rhs(i,j);
	return temp;
}

inline Matrix<double> Matrix<double>::getCol(uint num) const{
	Matrix<double> temp(n_row,1);
	for (uint i = 1; i <= n_row; i++)
		temp(i,1)	= data[i-1][num-1];
	return temp;
}

inline Matrix<double> Matrix<double>::getRow(uint num) const{
	Matrix<double> temp(1,n_col);
	for (uint i = 1; i <= n_col; i++)
		temp(1,i)	= data[num-1][i-1];
	return temp;
}

/*
 * friends operators
 */

template<class T>
T min(const Matrix<T>& matrix) {
	double val = 0;

	for (uint i = 1;i <= matrix.getRowDim(); i++) {
		for (uint j = 1; j <= matrix.getColDim(); j++)
			if (val < matrix(i,j))
				val = matrix(i,j);
	}
	return val;
}

template<class T>
T max(const Matrix<T>& matrix) {
	double val = 0;

	for (uint i = 1;i <= matrix.getRowDim(); i++) {
		for (uint j = 1; j <= matrix.getColDim(); j++)
			if (val > matrix(i,j))
				val = matrix(i,j);
	}
	return val;
}

template<class T>
Matrix<T> sin		 		 (const Matrix<T>& mat) {
	Matrix<T> temp(mat);
	for( uint i = 1 ; i <= mat.getRowDim();i++) {
		for (uint j = 1; j <= mat.getColDim(); j++)
			temp(i,j) 	= sin(mat(i,j));
	}
	return temp;
}

template<class T>
Matrix<T> cos		 		 (const Matrix<T>& mat) {
	Matrix<T> temp(mat);
	for( uint i = 1 ; i <= mat.getRowDim();i++) {
		for (uint j = 1; j <= mat.getColDim(); j++)
			temp(i,j) 	= cos(mat(i,j));
	}
	return temp;
}

template<class T>
Matrix<T> tan		 		 (const Matrix<T>& mat) {
	Matrix<T> temp(mat);
	for( uint i = 1 ; i <= mat.getRowDim();i++) {
		for (uint j = 1; j <= mat.getColDim(); j++)
			temp(i,j) 	= tan(mat(i,j));
	}
	return temp;
}

template<class T>
Matrix<T> sqrt		 		 (const Matrix<T>& mat) {
	Matrix<T> temp(mat);
	for( uint i = 1 ; i <= mat.getRowDim();i++) {
		for (uint j = 1; j <= mat.getColDim(); j++)
			temp(i,j) 	= sqrt(mat(i,j));
	}
	return temp;
}

template<class T>
Matrix<T> operator+ (double sum, const Matrix<T>& rhs) {
	return rhs + sum;
}

template<class T>
Matrix<T> operator- (double sum, const Matrix<T>& rhs){
	return (-rhs) + sum;
}

template<class T>
Matrix<T> operator* (double factor, const Matrix<T>& rhs){
	return rhs * factor;
}

template<class T>
Matrix<T> operator/ (double factor, const Matrix<T>& rhs){
	Matrix<T> temp(rhs.getRowDim(),rhs.getColDim());
	for(uint i = 1; i <= rhs.getRowDim(); i++)
		for (uint j = 1; j <= rhs.getColDim(); j++)
			temp(i,j) = factor/rhs(i,j);
	return temp;
}

template<class T>
Matrix<T> ones(uint row, uint col) {
	Matrix<T> temp(row, col);
	return temp + 1;
}

Matrix<double> linspace (double x1, double xn, uint n) {
	Matrix<double> temp(n,1);
	temp(1,1) 	= x1;
	double delta = (xn - x1)/(n-1);
	for (uint i = 2; i < n; i++) {
		temp(i,1)	= temp(i-1,1) + delta;
	}
	temp(n,1) = xn;
	return temp;
}

template<class T>
Matrix<T> cross(const Matrix<T>& v1, const Matrix<T>& v2) {
	if (v1.getRowDim() != v2.getRowDim() || v1.getColDim() != v2.getColDim()) {
		cout<<"===== Cross product of 2 col/row matrixes =====\n"
			  "=====      Matrix sizes do not match      =====\n\n";
	exit(1);
	}
	else if (v1.getColDim() != 1 && 1 != v1.getRowDim()) {
		cout<<"=====  Cross product of 2 col/row matrixes   =====\n"
			  "===== Matrixes must be column or row vectors =====\n\n";
		exit(1);
	}
	else if (v1.getRowDim() == 3) {
		Matrix<T> temp(3,1);

		temp(1,1) = v1(2,1) * v2(3,1) - v1(3,1) * v2(2,1);
		temp(2,1) = v1(3,1) * v2(1,1) - v1(1,1) * v2(3,1);
		temp(3,1) = v1(1,1) * v2(2,1) - v1(2,1) * v2(1,1);

		return temp;
	}
	else if (v1.getColDim() == 3) {
		Matrix<T> temp(1,3);

		temp(1,1) = v1(1,2) * v2(1,3) - v1(1,3) * v2(1,2);
		temp(1,2) = v1(1,3) * v2(1,1) - v1(1,1) * v2(1,3);
		temp(1,3) = v1(1,1) * v2(1,2) - v1(1,2) * v2(1,1);

		return temp;
	}
	else {
		cout<<"===== Cross product of 2 col/row matrixes =====\n"
			  "=====  Matrixes must have a length of 3   =====\n\n";
		exit(1);

	}
}

template<class T>
T dot(const Matrix<T>& v1, const Matrix<T>& v2) {
	if (v1.getRowDim() != v2.getRowDim() || v1.getColDim() != v2.getColDim()) {
		cout<<"===== Dot product of 2 col/row matrixes =====\n"
			  "=====     Matrix sizes do not match     =====\n\n";
	exit(1);
	}
	else if (v1.getColDim() != 1 && 1 != v1.getRowDim()) {
		cout<<"=====   Dot product of 2 col/row matrixes    =====\n"
			  "===== Matrixes must be column or row vectors =====\n\n";
		exit(1);
	}
	else {
		T temp_value = 0;
		for(uint i = 1; i <= v1.getRowDim(); i++)
			for (uint j = 1; j <= v2.getColDim(); j++)
				temp_value += v1(i,j) * v2(i,j);
		return temp_value;
	}
}

template<class T>
Matrix<T> transpose (const Matrix<T>& mat) {
	Matrix<T> temp(mat.getColDim(),mat.getRowDim());
	for (uint i = 1; i <= mat.getRowDim(); i++)
		for (uint j = 1; j <= mat.getColDim(); j++)
			temp(j,i) 	= mat(i,j);
	return temp;
}

inline Matrix<adouble>::Matrix() {
	n_col = 0;
	n_row = 0;
	data = NULL;
}

inline Matrix<adouble>::~Matrix() {
	if (this->data != NULL) {
		for ( uint i=0; i < n_row; i++) {
			delete[] data[i];
		}
		delete[] data;
	}
}

inline Matrix<adouble>::Matrix(uint m_length, uint m_num) {
	n_col = m_num;
	n_row = m_length;
	if (n_col == 0 || n_row == 0) {
		cout<<"WARNING: column or row size is zero!\n";
		data = NULL;
	}
	else {
		data 	= new adouble*[n_row];
		for (uint i = 0; i < n_row; i ++)
			data[i] 	= new adouble[n_col];

		for (uint i = 0; i < n_row; i ++)
			for (uint j = 0; j < n_col; j++)
				data[i][j] = 0.0;
	}
}

inline Matrix<adouble>::Matrix(const Matrix<adouble>& cp_Matrix){
	n_col = cp_Matrix.n_col;
	n_row = cp_Matrix.n_row;

	data 	= new adouble*[n_row];
	for (uint i = 0; i < n_row; i ++)
		data[i] 	= new adouble[n_col];

	for (uint i = 0; i < n_row; i++)
		for (uint j = 0; j < n_col; j++)
			data[i][j] = cp_Matrix.data[i][j];
}

inline Matrix<adouble>::Matrix(const Matrix<double>& cp_Matrix){
	n_col = cp_Matrix.getColDim();
	n_row = cp_Matrix.getRowDim();

	data 	= new adouble*[n_row];
	for (uint i = 0; i < n_row; i ++)
		data[i] 	= new adouble[n_col];

	for (uint i = 0; i < n_row; i++)
		for (uint j = 0; j < n_col; j++)
			data[i][j] = cp_Matrix(i+1,j+1);
}

inline void Matrix<adouble>::Print() const {
	for (uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++)
			printf("%+.4e\t",data[i][j].getValue());
		printf("\n");
	}
}

inline void Matrix<adouble>::Print(const char* printline) const {
	printf("\n*** %s ***\n\n",printline);
	for (uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++)
			printf("%+.4e\t",data[i][j].getValue());
		printf("\n");

	}
}

inline void Matrix<adouble>::resize(uint row_num, uint col_num) {
	Matrix <adouble> temp_matrix(*this);
	if (this->data != NULL) {
		for (uint i = 0; i < n_row; i++)
			delete[] data[i];
		delete[] data;
	}
	n_row 	= row_num;
	n_col 	= col_num;
	data 	= new adouble* [row_num];
	for (uint i = 0; i<row_num; i++)
		data[i]	= new adouble[col_num];

	for (uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++) {
			if (i < temp_matrix.getRowDim() && j < temp_matrix.getColDim())
				data[i][j]	= temp_matrix.data[i][j];
			else
				data[i][j] 	= 0.0;
		}
	}
}

inline adouble Matrix<adouble>::operator ()(uint row, uint col) const {
	if (col > n_col || row > n_row) {
		cout<<"===== Overloading operator (row, col) =====\n"
			  "=====          invalid range          =====\n";
		exit(1);
	}
	if (col < 1 || row < 1) {
		cout<<"===== Overloading operator (row, col) =====\n"
			  "=====       Index starts from 1       =====\n";
		exit(1);
	}
	return data[row-1][col-1];
}

inline adouble& Matrix<adouble>::operator ()(uint row, uint col) {
	if (col > n_col || row > n_row) {
		cout<<"===== Overloading operator (row, col) =====\n"
			  "=====          invalid range          =====\n";
		exit(1);
	}
	if (col < 1 || row < 1) {
		cout<<"===== Overloading operator (row, col) =====\n"
			  "=====       Index starts from 1       =====\n";
		exit(1);
	}
	return data[row-1][col-1];
}

inline adouble Matrix<adouble>::operator()(uint idx) const{
	if (n_col > 1 && n_row > 1) {
		cout<<"=====    Overloading operator (idx)     =====\n"
			  "===== Matrix is not a row or col vector =====\n";
		exit(1);
	}
	if (n_col == 0 || n_row == 0) {
		cout<<"=====  Overloading operator (idx)   =====\n"
			  "===== Matrix is not initialized yet =====\n";
		exit(1);
	}
	if (idx < 1 ) {
		cout<<"===== Overloading operator (idx) =====\n"
			  "=====     Index starts from 1    =====\n";
		exit(1);
	}
	if (n_col == 1)
		return data[idx-1][1];
	else
		return data[1][idx-1];
}

inline adouble& Matrix<adouble>::operator()(uint idx) {
	if (n_col > 1 && n_row > 1) {
		cout<<"=====    Overloading operator (idx)     =====\n"
			  "===== Matrix is not a row or col vector =====\n";
		exit(1);
	}
	if (n_col == 0 || n_row == 0) {
		cout<<"=====  Overloading operator (idx)   =====\n"
			  "===== Matrix is not initialized yet =====\n";
		exit(1);
	}
	if (idx < 1 ) {
		cout<<"===== Overloading operator (idx) =====\n"
			  "=====     Index starts from 1    =====\n";
		exit(1);
	}
	if (n_col == 1)
		return data[idx-1][1];
	else
		return data[1][idx-1];
}

inline void Matrix<adouble>::load (const char* filename) {
	double val;
	FILE *fp;
	if ( (fp = fopen(filename,"r")) == NULL ) {
		printf("===== Error opening file to load Matrix =====");
	}
	for (uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++) {
			fscanf(fp,"%le", &val);
			data[i][j] = val;
		}
	}
	fclose(fp);
}

inline void Matrix<adouble>::save (const char* filename) const {
	FILE *fp;
	if ( (fp = fopen(filename,"w")) == NULL ) {
		printf("===== Error opening file to save Matrix =====\n");
	}

	for (uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++) {
			fprintf(fp,"%.8e\t",data[i][j].getValue());
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

inline const Matrix<adouble>& Matrix<adouble>::operator= (const Matrix<adouble>& rhs) {
	this->resize(rhs.getRowDim(),rhs.getColDim());
	for(uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++)
			data[i][j] = rhs(i+1,j+1);
	}
	return *this;
}

inline const Matrix<adouble>& Matrix<adouble>::operator= (const Matrix<double>& rhs) {
	this->resize(rhs.getRowDim(),rhs.getColDim());
	for(uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++)
			data[i][j] = rhs(i+1,j+1);
	}
	return *this;
}

inline const Matrix<adouble>& Matrix<adouble>::operator= (const adouble scalar) {
	this->resize(1,1);
	data[0][0] = scalar;
	return *this;
}

inline const Matrix<adouble>& Matrix<adouble>::operator= (const double scalar) {
	this->resize(1,1);
	data[0][0] = scalar;
	return *this;
}

inline Matrix<adouble> Matrix<adouble>::operator, (const Matrix<adouble>& rhs) const{
	if (rhs.getRowDim() != n_row) {
		cout<<"===== Error in matrix concatenation  =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);
	}
	Matrix<adouble> temp(n_row, rhs.getColDim() + n_col);

	for (uint i = 1; i <= n_row; i++) {
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = this->operator ()(i,j);
	}

	for (uint i = 1; i <= n_row; i++) {
		for (uint j = 1; j <= rhs.getColDim(); j++)
			temp(i,j+n_col) = rhs(i,j);
	}
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::truncate_col	(uint l_limit, uint u_limit) {
	if (l_limit > u_limit){
		cout<<"===== Truncating columns of matrix =====\n"
			  "=====  lower limit > upper limit   =====\n";
		exit(1);
	}
	if (u_limit > n_col) {
		cout<<"===== Truncating columns of matrix =====\n"
			  "=====        invalid range         =====\n";
		exit(1);
	}
	if (l_limit < 1) {
		cout<<"===== Truncating columns of matrix =====\n"
			  "=====     Index starts from 1      =====\n";
		exit(1);
	}
	Matrix<adouble> temp(n_row, u_limit - l_limit + 1);
	for (uint i = 1; i <= n_row; i++) {
		for (uint j = l_limit; j <= u_limit; j++) {
			temp(i,j - l_limit + 1) = this->operator()(i,j);
		}
	}
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::truncate_row	(uint l_limit, uint u_limit) {
	if (l_limit > u_limit){
		cout<<"===== Truncating rows of matrix =====\n"
			  "===== lower limit > upper limit =====\n";
		exit(1);
	}
	if (u_limit > n_col) {
		cout<<"===== Truncating rows of matrix =====\n"
			  "=====       invalid range       =====\n";
		exit(1);
	}
	if (l_limit < 1) {
		cout<<"===== Truncating rows of matrix =====\n"
			  "=====    Index starts from 1    =====\n";
		exit(1);
	}
	Matrix<adouble> temp(u_limit - l_limit + 1,n_col);
	for (uint i = l_limit; i <= u_limit; i++) {
		for (uint j = 1; j <= n_col; j++) {
			temp(i - l_limit + 1,j) = this->operator()(i,j);
		}
	}
	return temp;
}

inline adouble Matrix<adouble>::enorm() {
	if (n_row > 1 && n_col > 1) {
		cout<<"=====			  !!! WARNING !!!			   =====\n"
			  "=====  matrix is neither column nor row vector  =====\n"
			  "===== enorm of a matrix does not make any sense =====\n";
	}
	adouble temp = 0;
	for(uint i = 0; i < n_row; i++)
		for (uint j = 0; j < n_col; j++)
			temp += data[i][j]*data[i][j];
	return sqrt(temp);
}

inline Matrix<adouble> Matrix<adouble>::operator+ (adouble sum) const{
	Matrix<adouble> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) += sum;
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::operator+ (double sum) const{
	Matrix<adouble> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) += sum;
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::operator- (adouble sum) const{
	Matrix<adouble> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) -= sum;
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::operator- (double sum) const{
	Matrix<adouble> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) -= sum;
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::operator* (adouble factor) const{
	Matrix<adouble> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) *= factor;
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::operator* (double factor) const{
	Matrix<adouble> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) *= factor;
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::operator/ (adouble factor) const{
	Matrix<adouble> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) /= factor;
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::operator/ (double factor) const{
	Matrix<adouble> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) /= factor;
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::operator- () const{
	Matrix<adouble> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) *= -1;
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::operator+ (const Matrix<adouble>& rhs) const{
	Matrix<adouble> temp(*this);
	if (rhs.getColDim() != n_col || rhs.getRowDim() != n_row) {
		cout<<"=====     Element-wise addition      =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);
	}
	for( uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j) + rhs(i,j);
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::operator- (const Matrix<adouble>& rhs) const{
	Matrix<adouble> temp(*this);
	if (rhs.getColDim() != n_col || rhs.getRowDim() != n_row) {
		cout<<"=====   Element-wise substitution    =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);	}
	for( uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j) - rhs(i,j);
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::operator* (const Matrix<adouble>& rhs) const{
	Matrix<adouble> temp(*this);
	if (rhs.getRowDim() != n_row || rhs.getColDim() != n_col) {
		cout<<"=====  Element-wise multiplication   =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);
	}
	if (n_row > 1 && n_col > 1) {
		cout<<"=====             !!! WARNING !!!             =====\n"
			  "===== matrix is neither column nor row vector =====\n"
			  "===== Element-wise multiplication of a matrix =====\n"
			  "=====    Do you know what you are doing???    =====\n";
	}
	for (uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j)*rhs(i,j);
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::operator/ (const Matrix<adouble>& rhs) const{
	Matrix<adouble> temp(*this);
	if (rhs.getRowDim() != n_row || rhs.getColDim() != n_col) {
		cout<<"=====     Element-wise division      =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);
	}
	if (n_row > 1 && n_col > 1) {
		cout<<"=====             !!! WARNING !!!             =====\n"
			  "===== matrix is neither column nor row vector =====\n"
			  "=====    Element-wise division of a matrix    =====\n"
			  "=====    Do you know what you are doing???    =====\n";
	}
	for (uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j)/rhs(i,j);
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::operator+ (const Matrix<double>& rhs) const{
	Matrix<adouble> temp(*this);
	if (rhs.getColDim() != n_col || rhs.getRowDim() != n_row) {
		cout<<"=====     Element-wise addition      =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);
	}
	for( uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j) + rhs(i,j);
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::operator- (const Matrix<double>& rhs) const{
	Matrix<adouble> temp(*this);
	if (rhs.getColDim() != n_col || rhs.getRowDim() != n_row) {
		cout<<"=====   Element-wise substitution    =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);	}
	for( uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j) - rhs(i,j);
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::operator* (const Matrix<double>& rhs) const{
	Matrix<adouble> temp(*this);
	if (rhs.getRowDim() != n_row || rhs.getColDim() != n_col) {
		cout<<"=====  Element-wise multiplication   =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);
	}
	if (n_row > 1 && n_col > 1) {
		cout<<"=====             !!! WARNING !!!             =====\n"
			  "===== matrix is neither column nor row vector =====\n"
			  "===== Element-wise multiplication of a matrix =====\n"
			  "=====    Do you know what you are doing???    =====\n";
	}
	for (uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j)*rhs(i,j);
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::operator/ (const Matrix<double>& rhs) const{
	Matrix<adouble> temp(*this);
	if (rhs.getRowDim() != n_row || rhs.getColDim() != n_col) {
		cout<<"=====     Element-wise division      =====\n"
			  "===== Matrix dimensions do not match =====\n";
		exit(1);
	}
	if (n_row > 1 && n_col > 1) {
		cout<<"=====             !!! WARNING !!!             =====\n"
			  "===== matrix is neither column nor row vector =====\n"
			  "=====    Element-wise division of a matrix    =====\n"
			  "=====    Do you know what you are doing???    =====\n";
	}
	for (uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j)/rhs(i,j);
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::getCol(uint num) const{
	Matrix<adouble> temp(n_row,1);
	for (uint i = 1; i <= n_row; i++)
		temp(i,1)	= data[i-1][num-1];
	return temp;
}

inline Matrix<adouble> Matrix<adouble>::getRow(uint num) const{
	Matrix<adouble> temp(1,n_col);
	for (uint i = 1; i <= n_col; i++)
		temp(1,i)	= data[num-1][i-1];
	return temp;
}

#endif /* MATRIX_H_ */
