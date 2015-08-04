/*
 * Vector.h
 *
 *  Created on: Jun 3, 2015
 *      Author: quanghoale
 */


using namespace std;

#ifndef SVECTOR_HPP_
#define SVECTOR_HPP_

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <adolc/adouble.h>

#include "MVector.hpp"

template<class T>
class SVector {
public:
	SVector();
	SVector(const SVector<T>& cp_SVector);
	SVector(uint size);
	virtual ~SVector();

	SVector<double> 	double_conv	();

	void 	resize		(uint n_new);
	uint 	getsize		() 									const	{return n;}

	void 	Print		() 									const;
	void 	Print		(const char* filename) 				const;

	SVector truncate	(uint l_limit, uint u_limit)		const;
	T 		enorm		();

	void 	load		(const char* filename);
	void 	save		(const char* filename)				const;

	const SVector<T>& operator= (const SVector<T>& rhs);
	T operator() (uint num)						const;
	T& operator() (uint num);

	/* element wise operations */
	SVector<T> 		operator+ 	(const SVector<T>& rhs) 			const;
	SVector<T> 		operator- 	(const SVector<T>& rhs)				const;
	SVector<T> 		operator* 	(const SVector<T>& rhs)				const;
	SVector<T> 		operator/ 	(const SVector<T>& rhs)				const;

	/* operations with a scalar */

	SVector<T> 		operator+ 	(T sum)					const;
	SVector<T> 		operator- 	(T sum)					const;
	SVector<T> 		operator* 	(T factor)					const;
	SVector<T> 		operator/ 	(T factor)					const;

	SVector<T> 		operator- 	()								const;

	/* concatenation operation */
	MVector	operator,	(const SVector<T>& rhs)			const;
	MVector	operator,	(const MVector& rhs)					const;

private:
	T* values;
	uint n;
};

template<class T>
SVector<T>::SVector() {
	n = 0;
	values = NULL;
}

template<class T>
SVector<T>::SVector(const SVector<T>& cp_SVector){
	n = cp_SVector.n;
	values = new T [n];
	for(uint i = 0; i<n; i++) {
		values[i] = cp_SVector.values[i];
	}
}

template<class T>
SVector<T>::SVector(uint size) {
	n = size;
	if (size <= 0) {
		values = NULL;
	}
	else {
		values = new T [n];
		for(uint i = 0; i<n; i++) {
			values[i] = 0.0;
		}
	}
}

template<class T>
SVector<T>::~SVector() {
	if (values != NULL)
		delete[] values;
}

template<class T>
SVector<double> SVector<T>::double_conv() {
	return *this;
}


template<class T>
void SVector<T>::resize(uint n_new) {
	if (values != NULL) {
			delete[] values;
	}
	n = n_new;
	if (n_new > 0) {
		values = new T [n];
		for(uint i = 0; i<n; i++) {
			values[i] 	= 0.0;
		}
	}
	else {
		values = NULL;
	}
}

template<class T>
void SVector<T>::Print() const {
	for (uint i = 0; i < this->getsize(); i++) {
		printf("%+.4e\n",values[i]);
	}
}

template<class T>
void SVector<T>::Print(const char* print_line) const {
	printf("\n*** %s ***\n\n", print_line);
	for (uint i = 0; i < this->getsize(); i++) {
		printf("%+.4e\n",values[i]);
	}
}

template<class T>
void SVector<T>::load (const char* filename) {

	FILE *fp;
	if ( (fp = fopen(filename,"r")) == NULL ) {
		printf("===== Error opening file to load SVector =====");
	}
//	double *values = new double[vec.getsize()];
	for (uint i = 0; i < this->getsize(); i++) {
		fscanf(fp,"%le", &values[i]);
//		vec(i+1)	= values;
	}
	fclose(fp);
}

template<class T>
void SVector<T>::save (const char* filename) const {
	FILE *fp;
	if ( (fp = fopen(filename,"w")) == NULL ) {
		printf("===== Error opening file to save MVector =====\n");
	}
	for (uint i = 0; i < this->getsize(); i++) {
		fprintf(fp,"%.8e\n",values[i]);
	}
	fclose(fp);
}

template<class T>
SVector<T> SVector<T>::truncate (uint l_limit, uint u_limit) const {
	if (l_limit > u_limit){
		cout<<"===== Truncating SVector elements =====\n"
			  "=====       invalid range        =====\n";
		exit(1);
	}
	if (u_limit > n) {
		cout<<"===== Truncating SVector elements =====\n"
			  "=====       invalid range        =====\n";
		exit(1);
	}
	if (l_limit < 1) {
		cout<<"===== Truncating SVector elements =====\n"
			  "=====    Index starts from 1     =====\n";
		exit(1);
	}
	SVector temp(u_limit - l_limit + 1);
	for (uint i = l_limit; i <= u_limit; i++) {
		temp.values[i - l_limit] = values[i-1];
	}
	return temp;
}

template<class T>
T SVector<T>::enorm() {
	T temp = 0;
	for(uint i = 0;i<n;i++) {
		temp += values[i]*values[i];
	}
	return sqrt(temp);
}

template<class T>
const SVector<T>& SVector<T>::operator= (const SVector<T>& rhs) {
	n = rhs.n;
	this->resize(n);
	for(uint i = 0; i<n; i++) {
		values[i] = rhs.values[i];
	}
	return *this;
}

template<class T>
T SVector<T>::operator() (uint num) const {
	if (num > n) {
		cout<<"===== getting SVector element =====\n===== invalid element number =====\n\n";
		exit(1);
	}
	if (num < 1) {
		cout<<"===== Index starts from 1 =====\n";
		exit(1);
	}
	return values[num - 1];
}

template<class T>
T& SVector<T>::operator() (uint num) {
	if (num > n) {
		cout<<"===== getting SVector element =====\n===== invalid element number =====\n\n";
		exit(1);
	}
	if (num < 1) {
		cout<<"===== Index starts from 1 =====\n";
		exit(1);
	}
	return values[num - 1];
}

/*
 * element wise operations
 */

template<class T>
SVector<T> SVector<T>::operator+ (const SVector<T>& rhs) const{
	SVector<T> temp(rhs);
	if (rhs.n != n) {
		cout<<"SVector dimensions do not match\n";
		exit(1);
	}
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]+rhs.values[i];
	}
	return temp;
}

template<class T>
SVector<T> SVector<T>::operator- (const SVector<T>& rhs) const{
	SVector<T> temp(rhs);
	if (rhs.n != n) {
		cout<<"SVector dimensions do not match\n";
		exit(1);
	}
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]-rhs.values[i];
	}
	return temp;
}

template<class T>
SVector<T> SVector<T>::operator* (const SVector<T>& rhs) const{
	SVector<T> temp(rhs);
	if (rhs.n != n) {
		cout<<"SVector dimensions do not match\n";
		exit(1);
	}
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]*rhs.values[i];
	}
	return temp;
}

template<class T>
SVector<T> SVector<T>::operator/ (const SVector<T>& rhs) const{
	SVector<T> temp(rhs);
	if (rhs.n != n) {
		cout<<"SVector dimensions do not match\n";
		exit(1);
	}
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]/rhs.values[i];
	}
	return temp;
}

template<class T>
SVector<T> SVector<T>::operator- ()						 const{
	SVector<T> temp(n);
	for(uint i = 0;i<n;i++) {
		temp.values[i] = -values[i];
	}
	return temp;
}

template<class T>
SVector<T> sin		 		 (const SVector<T>& vec) {
	SVector<T> temp(vec);
	for(uint i = 0;i<vec.getsize();i++) {
		temp(i+1) = sin(vec(i+1));
	}
	return temp;
}

template<class T>
SVector<T> cos		 		 (const SVector<T>& vec) {
	SVector<T> temp(vec);
	for(uint i = 0;i<vec.getsize();i++) {
		temp(i+1) = cos(vec(i+1));
	}
	return temp;
}

template<class T>
SVector<T> tan		 		 (const SVector<T>& vec) {
	SVector<T> temp(vec);
	for(uint i = 0;i<vec.getsize();i++) {
		temp(i+1) = tan(vec(i+1));
	}
	return temp;
}

template<class T>
SVector<T> atan2		 	 (const SVector<T>& v1, const SVector<T>& v2) {
	if (v1.getsize() != v2.getsize()) {
		cout<<"===== Evaluating atan2(SVector1, SVector2 =====\n"
			  "=====  SVector dimensions do not match   =====\n";
		exit(1);
	}
	SVector<T> temp(v1);
	for(uint i = 0;i<v1.getsiz();i++) {
		temp(i+1) = atan2(v1(i+1),v2(i+1));
	}
	return temp;
}

template<class T>
SVector<T> sqrt		 		 (const SVector<T>& vec) {
	SVector<T> temp(vec);
	for(uint i = 0;i<vec.getsize();i++) {
		temp(i+1) = sqrt(vec(i+1));
	}
	return temp;
}

/*
 * operations with a scalar
 */

template<class T>
SVector<T> SVector<T>::operator+ (T sum) const{
	SVector<T> temp(n);
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]+sum;
	}
	return temp;
}

template<class T>
SVector<T> SVector<T>::operator- (T sum) const{
	SVector<T> temp(n);
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]-sum;
	}
	return temp;
}

template<class T>
SVector<T> SVector<T>::operator* (T factor) const{
	SVector<T> temp;
	temp.resize(n);
	for(uint i = 0; i<n; i++) {
		temp.values[i] = factor*values[i];
	}
	return temp;
}

template<class T>
SVector<T> SVector<T>::operator/ (T factor) const{
	SVector<T> temp;
	temp.resize(n);
	for(uint i = 0; i<n; i++) {
		temp.values[i] = values[i]/factor;
	}
	return temp;
}

/*
 * Concatenation operations
 */

template<class T>
MVector SVector<T>::operator, (const SVector<T>& rhs) const{
	if (rhs.n != n) {
		cout<<"===== Error in vector concatenation  =====\n"
			  "===== SVector dimensions do not match =====\n";
		exit(1);
	}

	MVector temp(n,2);
	for (uint i = 1; i <= n; i++) {
		temp(i,1) 	= values[i-1];
		temp(i,2)	= rhs(i);
	}
	return temp;
}

template<class T>
MVector SVector<T>::operator, (const MVector& rhs) const{
	if (rhs.getRowDim() != n) {
		cout<<"===== Error in vector concatenation  =====\n"
			  "===== SVector dimensions do not match =====\n";
		exit(1);
	}
	MVector temp(n, rhs.getColDim() + 1);
	temp.getSVector(1) = *this;
	for (uint i = 1; i <= n; i++) {
		temp(i,1) 	= values[i-1];
	}

	for (uint i = 2; i <= rhs.getColDim(); i++) {
		for (uint j = 1; j<=n; j++) {
			temp(j,i) = rhs.getSVector(i-1);
		}
	}
	return temp;
}

/*
 * Specialization of SVector class
 */

template<>
class SVector <adouble>{
public:
	SVector();
	SVector(const SVector<adouble>& cp_SVector);
	SVector(uint size);
	virtual ~SVector();

	SVector<double>		double_conv();

	void 	resize		(uint n_new);
	uint 	getsize		() 									const	{return n;}

	void 	Print		() 									const;
	void 	Print		(const char* filename) 				const;

	SVector truncate	(uint l_limit, uint u_limit)		const;
	adouble 		enorm		();

	void 	load		(const char* filename);
	void 	save		(const char* filename)				const;

	const SVector<adouble>& operator= (const SVector<adouble>& rhs);
	adouble operator() (uint num)						const;
	adouble& operator() (uint num);

	/* element wise operations */
	SVector<adouble> 		operator+ 	(const SVector<double>& rhs) 			const;
	SVector<adouble> 		operator- 	(const SVector<double>& rhs)			const;
	SVector<adouble> 		operator* 	(const SVector<double>& rhs)			const;
	SVector<adouble> 		operator/ 	(const SVector<double>& rhs)			const;
	SVector<adouble> 		operator+ 	(const SVector<adouble>& rhs) 			const;
	SVector<adouble> 		operator- 	(const SVector<adouble>& rhs)			const;
	SVector<adouble> 		operator* 	(const SVector<adouble>& rhs)			const;
	SVector<adouble> 		operator/ 	(const SVector<adouble>& rhs)			const;

	/* operations with a scalar */

	SVector<adouble> 		operator+ 	(adouble sum)					const;
	SVector<adouble> 		operator- 	(adouble sum)					const;
	SVector<adouble> 		operator* 	(adouble factor)				const;
	SVector<adouble> 		operator/ 	(adouble factor)				const;

	SVector<adouble> 		operator+ 	(double sum)					const;
	SVector<adouble> 		operator- 	(double sum)					const;
	SVector<adouble> 		operator* 	(double factor)					const;
	SVector<adouble> 		operator/ 	(double factor)					const;

	SVector<adouble> 		operator- 	()								const;


private:
	adouble* values;
	uint n;
};

inline SVector<adouble>::SVector(){
	values = NULL;
	n = 0;
}

inline SVector<adouble>::SVector(const SVector<adouble>& cp_SVector){
	n = cp_SVector.n;
	values = new adouble [n];
	for(uint i = 0; i<n; i++) {
		values[i] = cp_SVector.values[i];
	}
}

inline SVector<adouble>::SVector(uint size) {
	n = size;
	if (size <= 0) {
		values = NULL;
	}
	else {
		values = new adouble [n];
		for(uint i = 0; i<n; i++) {
			values[i] = 0.0;
		}
	}
}

inline SVector<adouble>::~SVector() {
	if (values != NULL)
		delete[] values;
}

inline SVector<double> SVector<adouble>::double_conv() {
	SVector<double> temp(n);
	for (uint i = 0; i<n; i++) {
		temp(i+1)	= this->values[i].getValue();
	}
	return temp;
}

inline void SVector<adouble>::resize(uint n_new) {
	if (values != NULL) {
			delete[] values;
	}
	n = n_new;
	if (n_new > 0) {
		values = new adouble [n];
		for(uint i = 0; i<n; i++) {
			values[i] 	= 0.0;
		}
	}
	else {
		values = NULL;
	}
}

inline void SVector<adouble>::Print() const {
	for (uint i = 0; i < this->getsize(); i++) {
		printf("%+.4e\n",values[i].getValue());
	}
}

inline void SVector<adouble>::Print(const char* print_line) const{
	printf("\n*** %s ***\n\n", print_line);
	for (uint i = 0; i < this->getsize(); i++) {
		printf("%+.4e\n",values[i].getValue());
	}
}

inline SVector<adouble> SVector<adouble>::truncate (uint l_limit, uint u_limit) const {
	if (l_limit > u_limit){
		cout<<"===== Truncating SVector elements =====\n"
			  "=====       invalid range        =====\n";
		exit(1);
	}
	if (u_limit > n) {
		cout<<"===== Truncating SVector elements =====\n"
			  "=====       invalid range        =====\n";
		exit(1);
	}
	if (l_limit < 1) {
		cout<<"===== Truncating SVector elements =====\n"
			  "=====    Index starts from 1     =====\n";
		exit(1);
	}
	SVector temp(u_limit - l_limit + 1);
	for (uint i = l_limit; i <= u_limit; i++) {
		temp.values[i - l_limit] = values[i-1];
	}
	return temp;
}

inline adouble SVector<adouble>::enorm() {
	adouble temp = 0;
	for(uint i = 0;i<n;i++) {
		temp += values[i]*values[i];
	}
	return sqrt(temp);
}

inline const SVector<adouble>& SVector<adouble>::operator= (const SVector<adouble>& rhs) {
	n = rhs.n;
	this->resize(n);
	for(uint i = 0; i<n; i++) {
		values[i] = rhs.values[i];
	}
	return *this;
}

inline adouble SVector<adouble>::operator() (uint num) const {
	if (num > n) {
		cout<<"===== getting SVector element =====\n===== invalid element number =====\n\n";
		exit(1);
	}
	if (num < 1) {
		cout<<"===== Index starts from 1 =====\n";
		exit(1);
	}
	return values[num - 1];
}

inline adouble& SVector<adouble>::operator() (uint num) {
	if (num > n) {
		cout<<"===== getting SVector element =====\n===== invalid element number =====\n\n";
		exit(1);
	}
	if (num < 1) {
		cout<<"===== Index starts from 1 =====\n";
		exit(1);
	}
	return values[num - 1];
}

/*
 * element wise operations
 */

inline SVector<adouble> SVector<adouble>::operator+ (const SVector<adouble>& rhs) const{
	SVector<adouble> temp(rhs);
	if (rhs.n != n) {
		cout<<"SVector dimensions do not match\n";
		exit(1);
	}
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]+rhs.values[i];
	}
	return temp;
}

inline SVector<adouble> SVector<adouble>::operator- (const SVector<adouble>& rhs) const{
	SVector<adouble> temp(rhs);
	if (rhs.n != n) {
		cout<<"SVector dimensions do not match\n";
		exit(1);
	}
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]-rhs.values[i];
	}
	return temp;
}

inline SVector<adouble> SVector<adouble>::operator* (const SVector<adouble>& rhs) const{
	SVector<adouble> temp(rhs);
	if (rhs.n != n) {
		cout<<"SVector dimensions do not match\n";
		exit(1);
	}
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]*rhs.values[i];
	}
	return temp;
}

inline SVector<adouble> SVector<adouble>::operator/ (const SVector<adouble>& rhs) const{
	SVector<adouble> temp(rhs);
	if (rhs.n != n) {
		cout<<"SVector dimensions do not match\n";
		exit(1);
	}
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]/rhs.values[i];
	}
	return temp;
}

inline SVector<adouble> SVector<adouble>::operator+ (const SVector<double>& rhs) const{
	SVector<adouble> temp(rhs.getsize());
	if (rhs.getsize() != n) {
		cout<<"SVector dimensions do not match\n";
		exit(1);
	}
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]+rhs(i+1);
	}
	return temp;
}

inline SVector<adouble> SVector<adouble>::operator- (const SVector<double>& rhs) const{
	SVector<adouble> temp(rhs.getsize());
	if (rhs.getsize() != n) {
		cout<<"SVector dimensions do not match\n";
		exit(1);
	}
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]-rhs(i+1);
	}
	return temp;
}

inline SVector<adouble> SVector<adouble>::operator* (const SVector<double>& rhs) const{
	SVector<adouble> temp(rhs.getsize());
	if (rhs.getsize() != n) {
		cout<<"SVector dimensions do not match\n";
		exit(1);
	}
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]*rhs(i+1);
	}
	return temp;
}

inline SVector<adouble> SVector<adouble>::operator/ (const SVector<double>& rhs) const{
	SVector<adouble> temp(rhs.getsize());
	if (rhs.getsize() != n) {
		cout<<"SVector dimensions do not match\n";
		exit(1);
	}
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]/rhs(i+1);
	}
	return temp;
}

inline SVector<adouble> SVector<adouble>::operator- ()						 const{
	SVector<adouble> temp(n);
	for(uint i = 0;i<n;i++) {
		temp.values[i] = -values[i];
	}
	return temp;
}


/*
 * operations with a scalar
 */

inline SVector<adouble> SVector<adouble>::operator+ (double sum) const{
	SVector<adouble> temp(n);
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]+sum;
	}
	return temp;
}

inline SVector<adouble> SVector<adouble>::operator- (double sum) const{
	SVector<adouble> temp(n);
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]-sum;
	}
	return temp;
}

inline SVector<adouble> SVector<adouble>::operator* (double factor) const{
	SVector<adouble> temp;
	temp.resize(n);
	for(uint i = 0; i<n; i++) {
		temp.values[i] = factor*values[i];
	}
	return temp;
}

inline SVector<adouble> SVector<adouble>::operator/ (double factor) const{
	SVector<adouble> temp;
	temp.resize(n);
	for(uint i = 0; i<n; i++) {
		temp.values[i] = values[i]/factor;
	}
	return temp;
}


inline SVector<adouble> SVector<adouble>::operator+ (adouble sum) const{
	SVector<adouble> temp(n);
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]+sum;
	}
	return temp;
}

inline SVector<adouble> SVector<adouble>::operator- (adouble sum) const{
	SVector<adouble> temp(n);
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]-sum;
	}
	return temp;
}

inline SVector<adouble> SVector<adouble>::operator* (adouble factor) const{
	SVector<adouble> temp;
	temp.resize(n);
	for(uint i = 0; i<n; i++) {
		temp.values[i] = factor*values[i];
	}
	return temp;
}

inline SVector<adouble> SVector<adouble>::operator/ (adouble factor) const{
	SVector<adouble> temp;
	temp.resize(n);
	for(uint i = 0; i<n; i++) {
		temp.values[i] = values[i]/factor;
	}
	return temp;
}


/*
 * friend functions
 */

template<class T>
SVector<T> operator+ (double sum, const SVector<T>& rhs) {
	SVector<T> temp(rhs.getsize());
	for(uint i = 0;i<rhs.getsize();i++) {
		temp(i+1) = rhs(i+1)+sum;
	}
	return temp;
}

template<class T>
SVector<T> operator- (double sum, const SVector<T>& rhs){
	SVector<T> temp(rhs.getsize());

	for(uint i = 0;i<rhs.getsize();i++) {
		temp(i+1) = sum - rhs(i+1);
	}
	return temp;
}

template<class T>
SVector<T> operator* (double factor, const SVector<T>& rhs){
	SVector<T> temp(rhs.getsize());

	for(uint i = 0;i<rhs.getsize();i++) {
		temp(i+1) = factor*rhs(i+1);
	}
	return temp;
}

template<class T>
SVector<T> operator/ (double factor, const SVector<T>& rhs){
	SVector<T> temp(rhs.getsize());

	for(uint i = 0;i<rhs.getsize();i++) {
		temp(i+1) = factor/rhs(i+1);
	}
	return temp;
}

template<class T>
SVector<T> operator+ (adouble sum, const SVector<T>& rhs) {
	SVector<T> temp(rhs.getsize());
	for(uint i = 0;i<rhs.getsize();i++) {
		temp(i+1) = rhs(i+1)+sum;
	}
	return temp;
}

template<class T>
SVector<T> operator- (adouble sum, const SVector<T>& rhs){
	SVector<T> temp(rhs.getsize());

	for(uint i = 0;i<rhs.getsize();i++) {
		temp(i+1) = sum - rhs(i+1);
	}
	return temp;
}

template<class T>
SVector<T> operator* (adouble factor, const SVector<T>& rhs){
	SVector<T> temp(rhs.getsize());

	for(uint i = 0;i<rhs.getsize();i++) {
		temp(i+1) = factor*rhs(i+1);
	}
	return temp;
}

template<class T>
SVector<T> operator/ (adouble factor, const SVector<T>& rhs){
	SVector<T> temp(rhs.getsize());

	for(uint i = 0;i<rhs.getsize();i++) {
		temp(i+1) = factor/rhs(i+1);
	}
	return temp;
}

template<class T>
SVector<T> linspace (double x0, double xn, uint n) {
	SVector<double> temp(n);
	temp(1) 	= x0;
	double delta = (xn - x0)/(n-1);
	for (uint i = 2; i < n; i++) {
		temp(i)	= temp(i-1) + delta;
	}
	temp(n) = xn;
	return temp;
}

template<class T>
T dot(const SVector<T>& v1, const SVector<T>& v2) {
	if (v1.getsize() != v2.getsize()) {
		cout<<"===== calculating dot product  =====\n===== SVector size do not match =====\n\n";
	exit(1);
	}
	else {
		T temp_value = 0;
		for(uint i = 0; i < v1.getsize(); i++) {
			temp_value += v1(i+1) * v2(i+1);
		}
		return temp_value;
	}
}

template<class T>
SVector<T> cross(const SVector<T>& v1, const SVector<T>& v2) {
	if (v1.getsize() != v2.getsize()) {
		cout<<"===== calculating cross product  =====\n===== SVector size do not match =====\n\n";
		exit(1);
	}
	else if (v1.getsize() != 3) {
		cout<<"===== calculating cross product  =====\n=====    SVector size is not 3    =====\n\n";
		exit(1);
	}
	else {
		SVector<T> temp(3);
		T temp_value;

		temp_value = v1(2) * v2(3) - v1(3) * v2(2);
		temp(1) = temp_value;
		temp_value = v1(3) * v2(1) - v1(1) * v2(3);
		temp(2)	= temp_value;
		temp_value = v1(1) * v2(2) - v1(2) * v2(1);
		temp(3)	= temp_value;
		return temp;
	}
}

template<class T>
SVector<T> ones(uint num) {
	SVector<T> temp(num);
	return temp + 1;
}

template<class T>
SVector<T> zeros(uint num) {
	SVector<T> temp(num);
	return temp;
}

template<class T>
T min(const SVector<T>& v) {
	T value;
	value = v(1);
	for (uint i = 1; i < v.getsize(); i++ ) {
		if (value > v(i+1)) {
			value = v(i+1);
		}
	}
	return value;
}

template<class T>
T max(const SVector<T>& v) {
	T value;
	value = v(1);
	for (uint i = 1; i < v.getsize(); i++ ) {
		if (value < v(i+1)) {
			value = v(i+1);
		}
	}
	return value;
}

adouble dot(const SVector<adouble>& v1, const SVector<double>& v2);
adouble dot(const SVector<double>& v1, const SVector<adouble>& v2);

SVector<adouble> cross(const SVector<adouble>& v1, const SVector<double>& v2);
SVector<adouble> cross(const SVector<double>& v1, const SVector<adouble>& v2);


#endif
