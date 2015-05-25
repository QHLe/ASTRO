/*
 * SVector.cpp
 *
 *  Created on: May 13, 2015
 *      Author: zineus
 */


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <math.h>

#include "SVector.hpp"
#include "MVector.hpp"

using namespace std;

SVector::SVector() {
	// TODO Auto-generated constructor stub
	values = NULL;
	n = 0;
}

SVector::SVector(const SVector& cp_SVector){
	n = cp_SVector.n;
	values = new double [n];
	for(uint i = 0; i<n; i++) {
		values[i] = cp_SVector.values[i];
	}
}
SVector::SVector(uint size) {
	n = size;
	if (size <= 0) {
		values = NULL;
	}
	else {
		values = new double [n];
		for(uint i = 0; i<n; i++) {
			values[i] = 0;
		}
	}
}
SVector::~SVector() {
	// TODO Auto-generated destructor stub
	if (values != NULL) {
		delete[] values;
	}
}
void SVector::Print() const{
	for (uint i = 0; i < n; i++) {
		printf("%+.4e\n",values[i]);
	}
}
void SVector::Print(string print_line) const{
	cout<<endl<<"*** "<<print_line<<" ***"<<endl<<endl;

	for (uint i = 0; i < n; i++) {
		printf("%+.4e\n",values[i]);
	}
}

void SVector::resize(uint n_new) {
	if (values != NULL) {
			delete[] values;
	}
	n = n_new;
	if (n_new > 0) {
		values = new double [n];
		for(uint i = 0; i<n; i++) {
			values[i] 	= 0.0;
		}
	}
	else {
		values = NULL;
	}
}

void 	SVector::setelement(uint num, double value_new) {
	if(num > n) {
		cout<<"===== setting SVector element =====\n"
			  "===== invalid element number =====\n";
		exit(1);
	}
	if (num < 1) {
		cout<<"===== setting SVector element =====\n"
			  "=====  Index starts from 1   =====\n";
		exit(1);
	}
	values[num - 1] = value_new;
}

double	SVector::getelement(uint num) const{
	if (num > n) {
		cout<<"===== getting SVector element =====\n===== invalid element number =====\n\n";
		exit(1);
	}
	if (num < 1) {
		cout<<"===== setting SVector element =====\n"
			  "=====  Index starts from 1   =====\n";
		exit(1);
	}
	return values[num - 1];
}

double SVector::enorm() {
	double temp = 0;
	for(uint i = 0;i<n;i++) {
		temp += values[i]*values[i];
	}
	return sqrt(temp);
}
double dot(const SVector& v1, const SVector& v2) {
	if (v1.getsize() != v2.getsize()) {
		cout<<"===== calculating dot product  =====\n===== SVector size do not match =====\n\n";
	exit(1);
	}
	else {
		double temp_value = 0;
		for(uint i = 0; i < v1.n; i++) {
			temp_value += v1.values[i] * v2.values[i];
		}
		return temp_value;
	}
}
SVector ones(uint num) {
	SVector temp(num);
	return temp + 1;
}
SVector zeros(uint num) {
	SVector temp(num);
	return temp;
}

SVector cross(const SVector& v1, const SVector& v2) {
	if (v1.getsize() != v2.getsize()) {
		cout<<"===== calculating cross product  =====\n===== SVector size do not match =====\n\n";
		exit(1);
	}
	else if (v1.getsize() != 3) {
		cout<<"===== calculating cross product  =====\n=====    SVector size is not 3    =====\n\n";
		exit(1);
	}
	else {
		SVector temp(3);
		double temp_value;

		temp_value = v1.values[1] * v2.values[2] - v1.values[2] * v2.values[1];
		temp.values[0] = temp_value;
		temp_value = v1.values[2] * v2.values[0] - v1.values[0] * v2.values[2];
		temp.values[1]	= temp_value;
		temp_value = v1.values[0] * v2.values[1] - v1.values[1] * v2.values[0];
		temp.values[2]	= temp_value;
		return temp;
	}
}

SVector linspace (double x0, double xn, uint n) {
	SVector temp(n);
	temp.values[0] 	= x0;
	double delta = (xn - x0)/(n-1);
	for (uint i = 1; i < n; i++) {
		temp.values[i]	= temp.values[i-1] + delta;
	}
	return temp;
}

/* Elementary operations */

SVector SVector::operator+ (const SVector& rhs) const{
	SVector temp(rhs);
	if (rhs.n != n) {
		cout<<"SVector dimensions do not match\n";
		exit(1);
	}
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]+rhs.values[i];
	}
	return temp;
}

SVector SVector::operator- (const SVector& rhs) const{
	SVector temp(rhs);
	if (rhs.n != n) {
		cout<<"SVector dimensions do not match\n";
		exit(1);
	}
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]-rhs.values[i];
	}
	return temp;
}

SVector SVector::operator* (const SVector& rhs) const{
	SVector temp(rhs);
	if (rhs.n != n) {
		cout<<"SVector dimensions do not match\n";
		exit(1);
	}
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]*rhs.values[i];
	}
	return temp;
}

SVector SVector::operator/ (const SVector& rhs) const{
	SVector temp(rhs);
	if (rhs.n != n) {
		cout<<"SVector dimensions do not match\n";
		exit(1);
	}
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]/rhs.values[i];
	}
	return temp;
}

SVector sqrt		 		 (const SVector& vec) {
	SVector temp(vec);
	for(uint i = 0;i<vec.n;i++) {
		temp.values[i] = sqrt(vec.values[i]);
	}
	return temp;
}

SVector sin		 		 (const SVector& vec) {
	SVector temp(vec);
	for(uint i = 0;i<vec.n;i++) {
		temp.values[i] = sin(vec.values[i]);
	}
	return temp;
}

SVector cos		 		 (const SVector& vec) {
	SVector temp(vec);
	for(uint i = 0;i<vec.n;i++) {
		temp.values[i] = cos(vec.values[i]);
	}
	return temp;
}

SVector tan		 		 (const SVector& vec) {
	SVector temp(vec);
	for(uint i = 0;i<vec.n;i++) {
		temp.values[i] = tan(vec.values[i]);
	}
	return temp;
}

SVector atan2		 	 (const SVector& v1, const SVector& v2) {
	if (v1.n != v2.n) {
		cout<<"===== Evaluating atan2(SVector1, SVector2 =====\n"
			  "=====  SVector dimensions do not match   =====\n";
		exit(1);
	}
	SVector temp(v1);
	for(uint i = 0;i<v1.n;i++) {
		temp.values[i] = atan2(v1.values[i],v2.values[i]);
	}
	return temp;
}

/* Operations with a scalar */

SVector SVector::operator+ (double sum) const{
	SVector temp(n);
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]+sum;
	}
	return temp;
}

SVector SVector::operator- (double sum) const{
	SVector temp(n);
	for(uint i = 0;i<n;i++) {
		temp.values[i] = values[i]-sum;
	}
	return temp;
}

SVector SVector::operator* (double factor) const{
	SVector temp;
	temp.resize(n);
	for(uint i = 0; i<n; i++) {
		temp.values[i] = factor*values[i];
	}
	return temp;
}

SVector SVector::operator/ (double factor) const{
	SVector temp;
	temp.resize(n);
	for(uint i = 0; i<n; i++) {
		temp.values[i] = values[i]/factor;
	}
	return temp;
}


SVector SVector::operator- () const{
	SVector temp(n);
	for(uint i = 0;i<n;i++) {
		temp.values[i] = -values[i];
	}
	return temp;
}

double SVector::operator() (uint num) const {
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

double& SVector::operator() (uint num) {
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

SVector SVector::truncate (uint l_limit, uint u_limit) const {
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

/* friend operations */

SVector operator+ (double sum, const SVector& rhs) {
	SVector temp(rhs.n);
	for(uint i = 0;i<rhs.n;i++) {
		temp.values[i] = rhs.values[i]+sum;
	}
	return temp;
}

SVector operator- (double sum, const SVector& rhs){
	SVector temp(rhs.n);

	for(uint i = 0;i<rhs.n;i++) {
		temp.values[i] = sum - rhs.values[i];
	}
	return temp;
}

SVector operator* (double factor, const SVector& rhs){
	SVector temp(rhs.n);

	for(uint i = 0;i<rhs.n;i++) {
		temp.values[i] = factor*rhs.values[i];
	}
	return temp;
}

/* Concatenation operations */

MVector SVector::operator, (const SVector& rhs) const{
	if (rhs.n != n) {
		cout<<"===== Error in vector concatenation  =====\n"
			  "===== SVector dimensions do not match =====\n";
		exit(1);
	}

	MVector temp(n,2);
	temp.mvector[0] = *this;
	temp.mvector[1] = rhs;
	return temp;
}

MVector SVector::operator, (const MVector& rhs) const{
	if (rhs.vector_length != n) {
		cout<<"===== Error in vector concatenation  =====\n"
			  "===== SVector dimensions do not match =====\n";
		exit(1);
	}
	MVector temp(n, rhs.vector_num + 1);
	temp.mvector[0] = *this;
	for (uint i = 0; i < rhs.vector_num; i++) {
		temp.mvector[i+1] = rhs.mvector[i];
	}
	return temp;
}

const SVector& SVector::operator= (const SVector& rhs) {
	n = rhs.n;
	this->resize(n);
	for(uint i = 0; i<n; i++) {
		values[i] = rhs.values[i];
	}
	return *this;
}

void SVector::load (const char* filename) {

	FILE *fp;
	if ( (fp = fopen(filename,"r")) == NULL ) {
		printf("===== Error opening file to load MVector =====");
	}
	for (uint i = 0; i < n; i++) {
		fscanf(fp,"%le", &values[i]);
	}
	fclose(fp);
}

void SVector::save (const char* filename) const {
	FILE *fp;
	if ( (fp = fopen(filename,"w")) == NULL ) {
		printf("===== Error opening file to save MVector =====\n");
	}
	for (uint i = 0; i < n; i++) {
		fprintf(fp,"%.8e\n",values[i]);
	}
	fclose(fp);
}
double min(const SVector& v) {
	double value;
	value = v.values[0];
	for (uint i = 1; i < v.n; i++ ) {
		if (value > v.values[i]) {
			value = v.values[i];
		}
	}
	return value;
}
double max(const SVector& v) {
	double value;
	value = v.values[0];
	for (uint i = 1; i < v.n; i++ ) {
		if (value < v.values[i]) {
			value = v.values[i];
		}
	}
	return value;
}
