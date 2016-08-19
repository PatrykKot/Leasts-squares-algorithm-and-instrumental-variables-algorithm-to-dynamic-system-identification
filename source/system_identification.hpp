/*
 * object_identyfication_embedded.h
 *
 *  Created on: Mar 30, 2016
 *      Author: Patryk
 */

#ifndef OBJECT_IDENTYFICATION_EMBEDDED_H_
#define OBJECT_IDENTYFICATION_EMBEDDED_H_

#include "matrix_operations.hpp"

#define SYSTEM_IDENTYFICATION_METHOD_RLS 0
#define SYSTEM_IDENTYFICATION_METHOD_RIV 1

class SystemIdentification {
public:
	MatrixIdent objectParam;
private:
	USER_DATATYPE* inputSamples;
	USER_DATATYPE* outputSamples;
	int samplesIterator;
	int inputSamplesIterator;
	MatrixIdent covMatrix;
	USER_DATATYPE forgettingFactor;
	int rank;
	int identMethod;

public:
	SystemIdentification();
	SystemIdentification(const SystemIdentification & mat);
	SystemIdentification(int rank, USER_DATATYPE forgettingFactor = 1,
			USER_DATATYPE initialCovarianceMatrix = 1, int identMethod =
					SYSTEM_IDENTYFICATION_METHOD_RLS);
	~SystemIdentification();

	int getRank();
	
	SystemIdentification& operator=(const SystemIdentification& sys);

	void identStep(USER_DATATYPE inputSample, USER_DATATYPE outputSample);
	void print();
	bool isCorrect();

private:
	void insertNewSamples(USER_DATATYPE inSample,
			USER_DATATYPE outSammple);
	USER_DATATYPE getInputSample(int index); //
	USER_DATATYPE getOutputSample(int index);

	void RLSstep();
	void RIVstep();
};

#endif /* OBJECT_IDENTYFICATION_EMBEDDED_H_ */

