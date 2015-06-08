/*
 * MVector.cpp
 *
 *  Created on: Jun 5, 2015
 *      Author: quanghoale
 */


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <adolc/adolc.h>

#include "SVector.hpp"

using namespace std;


MVector::MVector() {
	vector_num = 0;
	vector_length = 0;
	mvector = NULL;
}

MVector::MVector(uint v_length, uint v_num) {
	vector_num = v_num;
	vector_length = v_length;
	if (v_num == 0) {
		mvector = NULL;
	}
	else {
		SVector<double> zero_vec(v_length);
		SVector<double>* pvec = new SVector<double>[v_num];
		for (uint i = 0; i < v_num; i++) {
			pvec[i] 	= zero_vec;
		}
		mvector = pvec;
	}
}

MVector::~MVector() {
	if (mvector != NULL) {
		delete[] mvector;
	}
}

void MVector::Print() const {
	for (uint i = 0; i < vector_length; i++) {
		for (uint j = 0; j < vector_num; j++) {
			printf("%+.4e\t",mvector[j](i+1));
		}
		cout<<endl;

	}
}

void MVector::Print(const char* printline) const {
	printf("\n*** %s ***\n\n",printline);
	for (uint i = 0; i < vector_length; i++) {
		for (uint j = 0; j < vector_num; j++) {
			printf("%+.4e\t",mvector[j](i+1));
		}
		cout<<endl;

	}
}

void MVector::resize(uint v_length, uint v_num) {
	if (mvector != NULL) {
		delete[] mvector;
	}
	vector_num = v_num;
	vector_length = v_length;
	SVector<double> zero_vec(v_length);
	SVector<double>* pvec = new SVector<double>[v_num];
	for (uint i = 0; i < v_num; i++) {
		pvec[i] 	= zero_vec;
	}
	mvector = pvec;
}

SVector<double> MVector::getSVector(uint num) const{
	if (num > vector_num) {
		cout<<"===== Getting vector from mvector =====\n"
			  "==== requested number out of range ====\n";
		exit(1);
	}
	if (num > vector_num) {
			cout<<"===== Getting vector from mvector =====\n"
				  "====      Index start from 1      ====\n";
			exit(1);
	}
	return mvector[num-1];
}

SVector<double> MVector::getRow(uint num) const{
	if (num > vector_length) {
		cout<<"===== Getting vector from mvector =====\n"
			  "==== requested number out of range ====\n";
		exit(1);
	}
	if (num < 1) {
		cout<<"===== Getting vector from mvector =====\n"
			  "=====     Index starts from 1     ====\n";
		exit(1);
	}
	SVector<double> temp(vector_num);
	for (uint i = 0; i < vector_num; i++) {
		temp(i+1)	= mvector[i](num);
	}
	return temp;
}

void MVector::load (const char* filename) {

	FILE *fp;
	if ( (fp = fopen(filename,"r")) == NULL ) {
		printf("===== Error opening file to load MVector =====");
	}
	for (uint i = 0; i < vector_length; i++) {
		for (uint j = 0; j < vector_num; j++) {
			fscanf(fp,"%le", &mvector[j](i+1));
		}
	}
	fclose(fp);
}

void MVector::save (const char* filename) const {
	FILE *fp;
	if ( (fp = fopen(filename,"w")) == NULL ) {
		printf("===== Error opening file to save MVector =====\n");
	}
//	fprintf(fp,"comment\n");
	for (uint i = 0; i < vector_length; i++) {
		for (uint j = 0; j < vector_num; j++) {
			fprintf(fp,"%.8e\t",mvector[j](i+1));
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

const MVector& MVector::operator= (const MVector& rhs) {
	vector_num = rhs.vector_num;
	vector_length = rhs.vector_length;

	this->resize(vector_length,vector_num);
	for(uint i = 0; i<vector_num; i++) {
		mvector[i] = rhs.mvector[i];
	}
	return *this;
}

const MVector& MVector::operator= (const SVector<double>& rhs) {
	vector_num = 1;
	vector_length = rhs.getsize();

	this->resize(vector_length,vector_num);
	for(uint i = 0; i<vector_length; i++) {
		mvector[0](i+1) = rhs(i+1);
	}
	return *this;
}

double MVector::operator ()(uint row, uint col) const {
	if (col > vector_num || row > vector_length) {
		cout<<"===== Overloading operator (row, col) =====\n"
			  "=====          invalid range          =====\n";
		exit(1);
	}
	if (col < 1 || row < 1) {
		cout<<"===== Overloading operator (row, col) =====\n"
			  "=====       Index starts from 1       =====\n";
		exit(1);
	}
	return mvector[col-1](row);
}

double& MVector::operator ()(uint row, uint col) {
	if (col > vector_num || row > vector_length) {
		cout<<"===== Overloading operator (row, col) =====\n"
			  "=====          invalid range          =====\n";
		exit(1);
	}
	if (col < 1 || row < 1) {
		cout<<"===== Overloading operator (row, col) =====\n"
			  "=====       Index starts from 1       =====\n";
		exit(1);
	}
	return mvector[col-1](row);
}

MVector MVector::truncate_col	(uint l_limit, uint u_limit) {
	if (l_limit > u_limit){
		cout<<"===== Truncating MVector =====\n"
			  "=====   invalid range    =====\n";
		exit(1);
	}
	if (u_limit > vector_num) {
		cout<<"===== Truncating MVector =====\n"
			  "=====   invalid range    =====\n";
		exit(1);
	}
	if (l_limit < 1) {
		cout<<"===== Truncating MVector =====\n"
			  "===== Index starts from 1 =====\n";
		exit(1);
	}
	MVector temp(vector_length, u_limit - l_limit + 1);
	for (uint i = l_limit; i <= u_limit; i++) {
		for (uint j = 0; j < vector_length; j++) {
			temp.mvector[i - l_limit](j+1) = mvector[i-1](j+1);
		}
	}
	return temp;
}

MVector MVector::truncate_row	(uint l_limit, uint u_limit) {
	if (l_limit > u_limit){
		cout<<"===== Truncating MVector =====\n"
			  "=====   invalid range    =====\n";
		exit(1);
	}
	if (u_limit > vector_length) {
		cout<<"===== Truncating MVector =====\n"
			  "=====   invalid range    =====\n";
		exit(1);
	}
	if (l_limit < 1) {
		cout<<"===== Truncating MVector =====\n"
			  "===== Index starts from 1 =====\n";
		exit(1);
	}
	MVector temp(u_limit - l_limit + 1,vector_num);
	for (uint i = l_limit; i <= u_limit; i++) {
		for (uint j = 0; j < vector_num; j++) {
			temp.mvector[j](i - l_limit +1) = mvector[j](i);
		}
	}
	return temp;
}




MVector MVector::operator, (const SVector<double>& rhs) const{
	if (rhs.getsize() != vector_length) {
		cout<<"===== Error in vector concatenation  =====\n"
			  "===== SVector dimensions do not match =====\n";
		exit(1);
	}
	MVector temp(vector_length,vector_num + 1);
	temp.mvector[vector_num] = rhs;
	for (uint i = 0; i < vector_num; i++) {
		temp.mvector [i] = mvector[i];
	}
	return temp;
}

MVector MVector::operator, (const MVector& rhs) const{
	if (rhs.vector_length != vector_length) {
		cout<<"===== Error in vector concatenation  =====\n"
			  "===== SVector dimensions do not match =====\n";
		exit(1);
	}
	MVector temp(vector_length, rhs.vector_num + vector_num);

	for (uint i = 0; i < vector_num; i++) {
		temp.mvector [i] = rhs.mvector[i];
	}

	for (uint i = 0; i < rhs.vector_num; i++) {
		temp.mvector [i+vector_num] = rhs.mvector[i];
	}
	return temp;
}

double min(const MVector& mv) {
	SVector<double> temp(mv.vector_num);
	for (uint i = 1;i <= mv.vector_num; i++) {
		temp(i) = min(mv.getSVector(i));
	}
	double retval = min(temp);
	return retval;
}

double max(const MVector& mv) {
	SVector<double> temp(mv.vector_num);
		for (uint i = 1;i <= mv.vector_num; i++) {
			temp(i) = max(mv.getSVector(i));
		}
		double retval = max(temp);
		return retval;
}

