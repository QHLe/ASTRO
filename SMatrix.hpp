/*
 * SMatrix.hpp
 *
 *  Created on: Aug 11, 2015
 *      Author: zineus
 */

#ifndef SMATRIX_H_
#define SMATRIX_H_

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

template <class T>
class SMatrix {
public:
	SMatrix();
	SMatrix(const SMatrix<T>& cp_SMatrix);
	SMatrix(uint v_length, uint v_num);
	virtual ~SMatrix();
	void 			Print 		()									const;
	void 			Print		(const char* line) 					const;
	void 			resize		(uint n_row, uint n_col);
	uint 			getRowDim	()	const	{return n_row;}
	uint 			getColDim	()	const 	{return n_col;}

	T 				operator()	(uint row, uint col)				const;
	T&				operator()	(uint row, uint col);
	T 				operator()	(uint idx) 							const;
	T&				operator()	(uint idx);

	SMatrix  		getCol		(uint num)							const;
	SMatrix			getRow		(uint num)							const;

	void 			setRow		(uint row_idx, SMatrix<T> row);
	void 			setCol		(uint col_idx, SMatrix<T> col);
	void			addRow		(SMatrix<T> row){
			if (row.getColDim() != this->getColDim()) {
				cout<<"  Error occurred by adding new row to SMatrix  \n"
					  "Column number of the two Matrixes doesn't match\n";
			}
			else {
				this->resize(n_row+1, n_col);
				this->setRow(n_row, row);
			}
		}

	SMatrix<T> 		operator+ 	(T sum)								const;
	SMatrix<T>	 	operator- 	(T sum)								const;
	SMatrix<T> 		operator* 	(T factor)							const;
	SMatrix<T> 		operator/ 	(T factor)							const;

	SMatrix<T> 		operator+ 	(const SMatrix<T>& rhs) 		const;
	SMatrix<T> 		operator- 	(const SMatrix<T>& rhs)			const;
	SMatrix<T> 		operator* 	(const SMatrix<T>& rhs)			const;
	SMatrix<T> 		operator/ 	(const SMatrix<T>& rhs)			const;
	SMatrix<T> operator- 	()								const;

	void 			load		(const char* filename);
	void 			save		(const char* filename)				const;

	const SMatrix& 	operator= 	(const SMatrix<T>& rhs);
	const SMatrix& 	operator= 	(const T scalar);
	SMatrix			operator,	(const SMatrix<T>& rhs)						const;

	SMatrix 			truncate_col	(uint l_limit, uint u_limit);
	SMatrix 			truncate_row	(uint l_limit, uint u_limit);

private:
	uint n_col;
	uint n_row;
	T** data;
};


template<class T>
SMatrix<T>::SMatrix() {
	n_col = 0;
	n_row = 0;
	data = NULL;
}

template<class T>
SMatrix<T>::SMatrix(uint m_length, uint m_num) {
	n_col = m_num;
	n_row = m_length;
	if (n_col == 0 || n_row == 0) {
		cout<<"WARNING in SMatrix constructor (row_num, col_num)\n"
			  "           column or row size is zero!           \n";
		data = NULL;
	}
	else {
		data 	= new T*[n_row];
		for (uint i = 0; i < n_row; i ++)
			data[i] 	= new T[n_col];

		for (uint i = 0; i < n_row; i ++)
			for (uint j = 0; j < n_col; j++)
				data[i][j] = 0.0;
	}
}

template<class T>
SMatrix<T>::SMatrix(const SMatrix<T>& cp_SMatrix){
	n_col = cp_SMatrix.getColDim();
	n_row = cp_SMatrix.getRowDim();

	data 	= new T*[n_row];
	for (uint i = 0; i < n_row; i ++)
		data[i] 	= new T[n_col];

	for (uint i = 0; i < n_row; i++)
		for (uint j = 0; j < n_col; j++)
			data[i][j] = cp_SMatrix(i+1,j+1);
}

template<class T>
SMatrix<T>::~SMatrix() {
	if (this->data != NULL) {
		for ( uint i=0; i < n_row; i++) {
			delete[] data[i];
		}
		delete[] data;
	}
}

template<class T>
void SMatrix<T>::Print() const {
	for (uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++)
			cout<<data[i][j]<<"\t";
		printf("\n");
	}
}

template<class T>
void SMatrix<T>::Print(const char* printline) const {
	printf("\n*** %s ***\n\n",printline);
	for (uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++)
			cout<<data[i][j]<<"\t";
		printf("\n");

	}
}

template<class T>
void SMatrix<T>::resize(uint row_num, uint col_num) {
	SMatrix <T> temp_matrix(*this);
	if (this->data != NULL) {
		for (uint i = 0; i < n_row; i++)
			delete[] data[i];
		delete[] data;
	}
	n_row 	= row_num;
	n_col 	= col_num;
	data 	= new T* [row_num];
	for (uint i = 0; i<row_num; i++)
		data[i]	= new T[col_num];

	for (uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++) {
			if (i < temp_matrix.getRowDim() && j < temp_matrix.getColDim())
				data[i][j]	= temp_matrix.data[i][j];
			else
				data[i][j] 	= 0.0;
		}
	}
}

template<class T>
T SMatrix<T>::operator ()(uint row, uint col) const {
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

template<class T>
T& SMatrix<T>::operator ()(uint row, uint col) {
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

template<class T>
T SMatrix<T>::operator()(uint idx) const{
	if (n_col > 1 && n_row > 1) {
		cout<<"=====    Overloading operator (idx)     =====\n"
			  "===== SMatrix is not a row or col vector =====\n";
		exit(1);
	}
	if (n_col == 0 || n_row == 0) {
		cout<<"=====  Overloading operator (idx)   =====\n"
			  "===== SMatrix is not initialized yet =====\n";
		exit(1);
	}
	if (idx < 1 ) {
		cout<<"===== Overloading operator (idx) =====\n"
			  "=====     Index starts from 1    =====\n";
		exit(1);
	}
	if (n_col == 1)
		return data[idx-1][0];
	else
		return data[0][idx-1];
}

template<class T>
T& SMatrix<T>::operator()(uint idx) {
	if (n_col > 1 && n_row > 1) {
		cout<<"=====    Overloading operator (idx)     =====\n"
			  "===== SMatrix is not a row or col vector =====\n";
		exit(1);
	}
	if (n_col == 0 || n_row == 0) {
		cout<<"=====  Overloading operator (idx)   =====\n"
			  "===== SMatrix is not initialized yet =====\n";
		exit(1);
	}
	if (idx < 1 ) {
		cout<<"===== Overloading operator (idx) =====\n"
			  "=====     Index starts from 1    =====\n";
		exit(1);
	}
	if (n_col == 1)
		return data[idx-1][0];
	else
		return data[0][idx-1];
}

template<class T>
void SMatrix<T>::load (const char* filename) {

	FILE *fp;
	if ( (fp = fopen(filename,"r")) == NULL ) {
		printf("===== Error opening file to load SMatrix =====");
	}
	for (uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++) {
			fscanf(fp,"%lf", &data[i][j]);
		}
	}
	fclose(fp);
}

template<class T>
void SMatrix<T>::save (const char* filename) const {
	FILE *fp;
	if ( (fp = fopen(filename,"w")) == NULL ) {
		printf("===== Error opening file to save SMatrix =====\n");
	}

	for (uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++) {
			fprintf(fp,"%.10e\t",data[i][j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

template<class T>
const SMatrix<T>& SMatrix<T>::operator= (const SMatrix<T>& rhs) {
	this->resize(rhs.getRowDim(),rhs.getColDim());
	for(uint i = 0; i < n_row; i++) {
		for (uint j = 0; j < n_col; j++)
			data[i][j] = rhs(i+1,j+1);
	}
	return *this;
}

template<class T>
const SMatrix<T>& SMatrix<T>::operator= (const T scalar) {
	this->resize(1,1);
	data[0][0] = scalar;
	return *this;
}

template<class T>
SMatrix<T> SMatrix<T>::operator, (const SMatrix<T>& rhs) const{
	if (rhs.getRowDim() != n_row) {
		cout<<"===== Error in matrix concatenation  =====\n"
			  "===== SMatrix dimensions do not match =====\n";
		exit(1);
	}
	SMatrix<T> temp(n_row, rhs.getColDim() + n_col);

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

template<class T>
SMatrix<T> SMatrix<T>::truncate_col	(uint l_limit, uint u_limit) {
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
	SMatrix<T> temp(n_row, u_limit - l_limit + 1);
	for (uint i = 1; i <= n_row; i++) {
		for (uint j = l_limit; j <= u_limit; j++) {
			temp(i,j - l_limit + 1) = this->operator()(i,j);
		}
	}
	return temp;
}

template<class T>
SMatrix<T> SMatrix<T>::truncate_row	(uint l_limit, uint u_limit) {
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
	SMatrix<T> temp(u_limit - l_limit + 1,n_col);
	for (uint i = l_limit; i <= u_limit; i++) {
		for (uint j = 1; j <= n_col; j++) {
			temp(i - l_limit + 1,j) = this->operator()(i,j);
		}
	}
	return temp;
}

template<class T>
SMatrix<T> SMatrix<T>::getCol(uint num) const{
	SMatrix<T> temp(n_row,1);
	for (uint i = 1; i <= n_row; i++)
		temp(i,1)	= data[i-1][num-1];
	return temp;
}

template<class T>
SMatrix<T> SMatrix<T>::getRow(uint num) const{
	SMatrix<T> temp(1,n_col);
	for (uint i = 1; i <= n_col; i++)
		temp(1,i)	= data[num-1][i-1];
	return temp;
}

template<class T>
void SMatrix<T>::setRow	(uint row_idx, SMatrix<T> row){
	if (row.getColDim() != n_col) {
		cout<<"=====    Error in function setRow    =====\n"
			  "===== Column dimensions do not match =====\n";
	}

	else if (row.getRowDim() != 1) {
		cout<<"=====       Error in function setRow       =====\n"
			  "===== Parameter is not a single row matrix =====\n";
	}
	else if (row_idx < 1 || row_idx > n_row) {
		cout<<"===== Error in function setRow =====\n"
			  "=====    Row_idx is invalid    =====\n";
	}

	for (uint i = 0; i < n_col; i++)
		data[row_idx-1][i] 	= row(i+1);
}


template<class T>
void SMatrix<T>::setCol	(uint col_idx, SMatrix<T> col){
	if (col.getRowDim() != n_row) {
		cout<<"=====  Error in function setCol   =====\n"
			  "===== Row dimensions do not match =====\n";
		exit(1);
	}

	else if (col.getColDim() != 1) {
		cout<<"=====        Error in function setCol         =====\n"
			  "===== Parameter is not a single column matrix =====\n";
		exit(1);
	}
	else if (col_idx < 1 || col_idx > n_col) {
		cout<<"===== Error in function setCol =====\n"
			  "=====    col_idx is invalid    =====\n";
		exit(1);
	}

	for (uint i = 0; i < n_row; i++)
		data[i][col_idx-1] 	= col(i+1);
}

template<class T>
SMatrix<T> SMatrix<T>::operator+ (T sum) const{
	SMatrix<T> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) += sum;
	return temp;
}

template<class T>
SMatrix<T> SMatrix<T>::operator- (T sum) const{
	SMatrix<T> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) -= sum;
	return temp;
}

template<class T>
SMatrix<T> SMatrix<T>::operator* (T factor) const{
	SMatrix<T> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) *= factor;
	return temp;
}

template<class T>
SMatrix<T> SMatrix<T>::operator/ (T factor) const{
	SMatrix<T> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) /= factor;
	return temp;
}

template<class T>
SMatrix<T> SMatrix<T>::operator+ (const SMatrix<T>& rhs) const{
	SMatrix<T> temp(*this);
	if (rhs.getRowDim() == 1 && rhs.getColDim() == 1) {
			for (uint i = 1; i <= n_row; i++)
				for (uint j = 1; j <= n_col; j++)
					temp(i,j) = temp(i,j)+rhs(1,1);
			return temp;
		}
		else if (n_row == 1 && n_col == 1) {
			double val = temp(1,1);
			temp.resize(rhs.getRowDim(), rhs.getColDim());
			for (uint i = 1; i <= rhs.getRowDim(); i++)
				for (uint j = 1; j <= rhs.getColDim(); j++)
					temp(i,j) = rhs(i,j)+val;
			return temp;
		}
		else if (rhs.getColDim() != n_col || rhs.getRowDim() != n_row) {
		cout<<"=====     Element-wise addition      =====\n"
			  "===== SMatrix dimensions do not match =====\n";
		exit(1);
	}
	for( uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j) + rhs(i,j);
	return temp;
}

template<class T>
SMatrix<T> SMatrix<T>::operator- (const SMatrix<T>& rhs) const{
	SMatrix<T> temp(*this);
	if (rhs.getRowDim() == 1 && rhs.getColDim() == 1) {
			for (uint i = 1; i <= n_row; i++)
				for (uint j = 1; j <= n_col; j++)
					temp(i,j) = temp(i,j)-rhs(1,1);
			return temp;
		}
		else if (n_row == 1 && n_col == 1) {
			double val = temp(1,1);
			temp.resize(rhs.getRowDim(), rhs.getColDim());
			for (uint i = 1; i <= rhs.getRowDim(); i++)
				for (uint j = 1; j <= rhs.getColDim(); j++)
					temp(i,j) = val-rhs(i,j);
			return temp;
		}
		else if (rhs.getColDim() != n_col || rhs.getRowDim() != n_row) {
		cout<<"=====   Element-wise substitution    =====\n"
			  "===== SMatrix dimensions do not match =====\n";
		exit(1);	}
	for( uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) = temp(i,j) - rhs(i,j);
	return temp;
}

template<class T>
SMatrix<T> SMatrix<T>::operator* (const SMatrix<T>& rhs) const{
	SMatrix<T> temp(*this);
	if (rhs.getRowDim() == 1 && rhs.getColDim() == 1) {
		for (uint i = 1; i <= n_row; i++)
			for (uint j = 1; j <= n_col; j++)
				temp(i,j) = temp(i,j)*rhs(1,1);
		return temp;
	}
	else if (n_row == 1 && n_col == 1) {
		double val = temp(1,1);
		temp.resize(rhs.getRowDim(), rhs.getColDim());
		for (uint i = 1; i <= rhs.getRowDim(); i++)
			for (uint j = 1; j <= rhs.getColDim(); j++)
				temp(i,j) = rhs(i,j)*val;
		return temp;
	}
	else if (rhs.getRowDim() != n_row || rhs.getColDim() != n_col) {
		cout<<"=====  Element-wise multiplication   =====\n"
			  "===== SMatrix dimensions do not match =====\n";
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

template<class T>
SMatrix<T> SMatrix<T>::operator/ (const SMatrix<T>& rhs) const{
	SMatrix<T> temp(*this);
	if (rhs.getRowDim() == 1 && rhs.getColDim() == 1) {
			for (uint i = 1; i <= n_row; i++)
				for (uint j = 1; j <= n_col; j++)
					temp(i,j) = temp(i,j)/rhs(1,1);
			return temp;
		}
		else if (n_row == 1 && n_col == 1) {
			double val = temp(1,1);
			temp.resize(rhs.getRowDim(), rhs.getColDim());
			for (uint i = 1; i <= rhs.getRowDim(); i++)
				for (uint j = 1; j <= rhs.getColDim(); j++)
					temp(i,j) = val/rhs(i,j);
			return temp;
		}
		else if (rhs.getRowDim() != n_row || rhs.getColDim() != n_col) {
		cout<<"=====     Element-wise division      =====\n"
			  "===== SMatrix dimensions do not match =====\n";
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

template<class T>
SMatrix<T> SMatrix<T>::operator- () const{
	SMatrix<T> temp(*this);
	for(uint i = 1; i <= n_row; i++)
		for (uint j = 1; j <= n_col; j++)
			temp(i,j) *= -1;
	return temp;
}

/*
 * friends operators
 */
template<class T>
T enorm(SMatrix<T> vec) {
	if (vec.getRowDim() > 1 && vec.getColDim() > 1) {
		cout<<"=====			  !!! WARNING !!!			   =====\n"
			  "=====  matrix is neither column nor row vector  =====\n"
			  "===== enorm of a matrix does not make any sense =====\n";
	}
	T temp = 0;
	for(uint i = 1; i <= vec.getRowDim(); i++)
		for (uint j = 1; j <= vec.getColDim(); j++)
			temp += vec(i,j)*vec(i,j);
	return sqrt(temp);
}

template<class T>
T min(const SMatrix<T>& matrix) {
	double val = matrix(1,1);

	for (uint i = 1;i <= matrix.getRowDim(); i++) {
		for (uint j = 1; j <= matrix.getColDim(); j++)
			if (val < matrix(i,j))
				val = matrix(i,j);
	}
	return val;
}

template<class T>
T max(const SMatrix<T>& matrix) {
	double val = matrix(1,1);

	for (uint i = 1;i <= matrix.getRowDim(); i++) {
		for (uint j = 1; j <= matrix.getColDim(); j++)
			if (val > matrix(i,j))
				val = matrix(i,j);
	}
	return val;
}

template<class T>
SMatrix<T> sin		 		 (const SMatrix<T>& mat) {
	SMatrix<T> temp(mat);
	for( uint i = 1 ; i <= mat.getRowDim();i++) {
		for (uint j = 1; j <= mat.getColDim(); j++)
			temp(i,j) 	= sin(mat(i,j));
	}
	return temp;
}

template<class T>
SMatrix<T> cos		 		 (const SMatrix<T>& mat) {
	SMatrix<T> temp(mat);
	for( uint i = 1 ; i <= mat.getRowDim();i++) {
		for (uint j = 1; j <= mat.getColDim(); j++)
			temp(i,j) 	= cos(mat(i,j));
	}
	return temp;
}

template<class T>
SMatrix<T> tan		 		 (const SMatrix<T>& mat) {
	SMatrix<T> temp(mat);
	for( uint i = 1 ; i <= mat.getRowDim();i++) {
		for (uint j = 1; j <= mat.getColDim(); j++)
			temp(i,j) 	= tan(mat(i,j));
	}
	return temp;
}

template<class T>
SMatrix<T> sqrt		 		 (const SMatrix<T>& mat) {
	SMatrix<T> temp(mat);
	for( uint i = 1 ; i <= mat.getRowDim();i++) {
		for (uint j = 1; j <= mat.getColDim(); j++)
			temp(i,j) 	= sqrt(mat(i,j));
	}
	return temp;
}

template<class T>
SMatrix<T> operator+ (double sum, const SMatrix<T>& rhs) {
	return rhs + sum;
}

template<class T>
SMatrix<T> operator- (double sum, const SMatrix<T>& rhs){
	return (-rhs) + sum;
}

template<class T>
SMatrix<T> operator* (double factor, const SMatrix<T>& rhs){
	return rhs * factor;
}


template<class T>
SMatrix<T> operator/ (double factor, const SMatrix<T>& rhs){
	SMatrix<T> temp(rhs.getRowDim(),rhs.getColDim());
	for(uint i = 1; i <= rhs.getRowDim(); i++)
		for (uint j = 1; j <= rhs.getColDim(); j++)
			temp(i,j) = factor/rhs(i,j);
	return temp;
}


template<class T>
SMatrix<T> ones(uint row, uint col) {
	SMatrix<T> temp(row, col);
	return temp + 1.0;
}

template<class T>
SMatrix<T> linspace (double x1, double xn, uint n) {
	SMatrix<T> temp(n,1);
	temp(1,1) 	= x1;
	double delta = (xn - x1)/(n-1);
	for (uint i = 2; i < n; i++) {
		temp(i,1)	= temp(i-1,1) + delta;
	}
	temp(n,1) = xn;
	return temp;
}

template<class T>
SMatrix<T> cross(const SMatrix<T>& v1, const SMatrix<T>& v2) {
	if (v1.getRowDim() != v2.getRowDim() || v1.getColDim() != v2.getColDim()) {
		cout<<"===== Cross product of 2 col/row matrixes =====\n"
			  "=====      SMatrix sizes do not match      =====\n\n";
	exit(1);
	}
	else if (v1.getColDim() != 1 && 1 != v1.getRowDim()) {
		cout<<"=====  Cross product of 2 col/row matrixes   =====\n"
			  "===== SMatrixes must be column or row vectors =====\n\n";
		exit(1);
	}
	else if (v1.getRowDim() == 3) {
		SMatrix<T> temp(3,1);

		temp(1,1) = v1(2,1) * v2(3,1) - v1(3,1) * v2(2,1);
		temp(2,1) = v1(3,1) * v2(1,1) - v1(1,1) * v2(3,1);
		temp(3,1) = v1(1,1) * v2(2,1) - v1(2,1) * v2(1,1);

		return temp;
	}
	else if (v1.getColDim() == 3) {
		SMatrix<T> temp(1,3);

		temp(1,1) = v1(1,2) * v2(1,3) - v1(1,3) * v2(1,2);
		temp(1,2) = v1(1,3) * v2(1,1) - v1(1,1) * v2(1,3);
		temp(1,3) = v1(1,1) * v2(1,2) - v1(1,2) * v2(1,1);

		return temp;
	}
	else {
		cout<<"===== Cross product of 2 col/row matrixes =====\n"
			  "=====  SMatrixes must have a length of 3   =====\n\n";
		exit(1);

	}
}

template<class T>
T dot(const SMatrix<T>& v1, const SMatrix<T>& v2) {
	if (v1.getRowDim() != v2.getRowDim() || v1.getColDim() != v2.getColDim()) {
		cout<<"===== Dot product of 2 col/row matrixes =====\n"
			  "=====     SMatrix sizes do not match     =====\n\n";
	exit(1);
	}
	else if (v1.getColDim() != 1 && 1 != v1.getRowDim()) {
		cout<<"=====   Dot product of 2 col/row matrixes    =====\n"
			  "===== SMatrixes must be column or row vectors =====\n\n";
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
SMatrix<T> transpose (const SMatrix<T>& mat) {
	SMatrix<T> temp(mat.getColDim(),mat.getRowDim());
	for (uint i = 1; i <= mat.getRowDim(); i++)
		for (uint j = 1; j <= mat.getColDim(); j++)
			temp(j,i) 	= mat(i,j);
	return temp;
}
template<class T>
SMatrix<T> abs (const SMatrix<T>& mat) {
	SMatrix<T> temp(mat.getColDim(),mat.getRowDim());
	for (uint i = 1; i <= mat.getRowDim(); i++)
		for (uint j = 1; j <= mat.getColDim(); j++)
			temp(j,i) 	= fabs(mat(i,j));
	return temp;
}

#endif /* SMATRIX_H_ */
