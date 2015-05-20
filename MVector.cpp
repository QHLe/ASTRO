/*
 * MVector.cpp
 *
 *  Created on: May 14, 2015
 *      Author: zineus
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include "MVector.hpp"
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
		SVector zero_vec(v_length);
		SVector* pvec = new SVector[v_num];
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
			printf("%+.4e\t",mvector[j].values[i]);
		}
		cout<<endl;

	}
}
void MVector::Print(string printline) const {
	cout<<endl<<"*** "<<printline<<" ***"<<endl<<endl;
	for (uint i = 0; i < vector_length; i++) {
		for (uint j = 0; j < vector_num; j++) {
			printf("%+.4e\t",mvector[j].values[i]);
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
	SVector zero_vec(v_length);
	SVector* pvec = new SVector[v_num];
	for (uint i = 0; i < v_num; i++) {
		pvec[i] 	= zero_vec;
	}
	mvector = pvec;
}

SVector MVector::getSVector(uint num) const{
	if (num >= vector_num) {
		cout<<"===== Getting vector from mvector =====\n"
			  "==== requested number out of range ====\n";
		exit(1);
	}
	return mvector[num];
}

SVector MVector::getRow(uint num) const{
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
	SVector temp(vector_num);
	for (uint i = 0; i < vector_num; i++) {
		temp.values[i]	= mvector[i].values[num - 1];
	}
	return temp;
}

MVector MVector::operator, (const SVector& rhs) const{
	if (rhs.n != vector_length) {
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


MVector MVector::truncate	(uint l_limit, uint u_limit) {
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
	MVector temp(u_limit - l_limit + 1,vector_length);
	for (uint i = l_limit; i <= u_limit; i++) {
		for (uint j = 0; j < vector_length; j++) {
			temp.mvector[i - l_limit].values[j] = mvector[i-1].values[j];
		}
	}
	return temp;
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
	return mvector[col-1].values[row-1];
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


const MVector& MVector::operator= (const SVector& rhs) {
	vector_num = 1;
	vector_length = rhs.n;

	this->resize(vector_length,vector_num);
	for(uint i = 0; i<vector_length; i++) {
		mvector[0].values[i] = rhs.values[i];
	}
	return *this;
}

void MVector::load (const char* filename) {

	FILE *fp;
	if ( (fp = fopen(filename,"r")) == NULL ) {
		printf("===== Error opening file to load MVector =====");
	}
	for (uint i = 0; i < vector_length; i++) {
		for (uint j = 0; j < vector_num; j++) {
			fscanf(fp,"%le", &mvector[j].values[i]);
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
			fprintf(fp,"%.8e\t",mvector[j].values[i]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}
