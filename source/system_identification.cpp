/*
 * object_identyfication_embedded.c
 *
 *  Created on: Mar 30, 2016
 *      Author: Patryk
 */

#include "system_identification.hpp"

SystemIdentification::SystemIdentification() {
	rank = forgettingFactor = 1;
	samplesIterator = 0;
	identMethod = SYSTEM_IDENTYFICATION_METHOD_RLS;
	inputSamples = new USER_DATATYPE[rank + 1];
	outputSamples = new USER_DATATYPE[rank + 1];

	for (int i = 0; i < (rank + 1); i++)
	{
		inputSamples[i] = outputSamples[i] = 0;
	}

	this->objectParam = MatrixIdent(rank * 2, 1);

	this->covMatrix = MatrixIdent(2 * rank, 2 * rank);
	for (int i = 0; i < covMatrix.getRows(); i++)
		for (int j = 0; j < covMatrix.getCols(); j++)
			if (i != j)
				covMatrix.insert(0, i, j);
			else
				covMatrix.insert(1, i, j);
}

SystemIdentification::SystemIdentification(const SystemIdentification & mat) {
	rank = mat.rank;
	forgettingFactor = mat.forgettingFactor;
	samplesIterator = mat.samplesIterator;
	identMethod = mat.identMethod;
	inputSamplesIterator = mat.inputSamplesIterator;

	if (identMethod == SYSTEM_IDENTYFICATION_METHOD_RIV)
	{
		inputSamples = new USER_DATATYPE[2 * rank + 1];
		for (int i = 0; i < (2* rank + 1); i++)
			inputSamples[i] = mat.inputSamples[i];
	}
	else
	{
		inputSamples = new USER_DATATYPE[rank + 1];
		for (int i = 0; i < (rank + 1); i++)
			inputSamples[i] = mat.inputSamples[i];
	}

	outputSamples = new USER_DATATYPE[rank + 1];
	for (int i = 0; i < (rank + 1); i++)
		outputSamples[i] = mat.outputSamples[i];

	covMatrix = mat.covMatrix;
	objectParam = mat.objectParam;
}

SystemIdentification::SystemIdentification(int rank,
USER_DATATYPE forgettingFactor, USER_DATATYPE initialCovarianceMatrix, int identMethod) {
	this->rank = rank;
	this->forgettingFactor = forgettingFactor;
	samplesIterator = 0;
	if (identMethod == SYSTEM_IDENTYFICATION_METHOD_RIV
			|| identMethod == SYSTEM_IDENTYFICATION_METHOD_RLS)
		this->identMethod = identMethod;
	else
		this->identMethod = SYSTEM_IDENTYFICATION_METHOD_RLS;

	this->objectParam = MatrixIdent(rank * 2, 1);

	this->covMatrix = MatrixIdent(2 * rank, 2 * rank);
	for (int i = 0; i < covMatrix.getRows(); i++)
		for (int j = 0; j < covMatrix.getCols(); j++)
			if (i != j)
				covMatrix.insert(0, i, j);
			else
				covMatrix.insert(initialCovarianceMatrix, i, j);

	if (identMethod == SYSTEM_IDENTYFICATION_METHOD_RIV) {
		inputSamples = new USER_DATATYPE[2 * rank + 1];
		inputSamplesIterator = 0;
	}
	else
	{
		inputSamples = new USER_DATATYPE[rank + 1];
	}

	outputSamples = new USER_DATATYPE[rank + 1];

	for (int i = 0; i < (rank + 1); i++) outputSamples[i] = 0;
	if (identMethod == SYSTEM_IDENTYFICATION_METHOD_RIV)
	{
		for (int i = 0; i < (2 * rank + 1); i++) inputSamples[i] = 0;
	}
	else
	{
		for (int i = 0; i < (rank + 1); i++) inputSamples[i] = 0;
	}
}

SystemIdentification::~SystemIdentification() {
	delete[] outputSamples;
	delete[] inputSamples;
}

SystemIdentification& SystemIdentification::operator=(const SystemIdentification& sys)
{
	if(inputSamples!=NULL)	delete[] inputSamples;
	if(outputSamples!=NULL) delete[] outputSamples;
	
	rank = sys.rank;
	forgettingFactor = sys.forgettingFactor;
	samplesIterator = sys.samplesIterator;
	identMethod = sys.identMethod;
	inputSamplesIterator = sys.inputSamplesIterator;

	if (identMethod == SYSTEM_IDENTYFICATION_METHOD_RIV)
	{
		inputSamples = new USER_DATATYPE[2*rank + 1];
		for (int i = 0; i < (2*rank + 1); i++)
			inputSamples[i] = sys.inputSamples[i];
	}
	else
	{
		inputSamples = new USER_DATATYPE[rank + 1];
		for (int i = 0; i < (rank + 1); i++)
			inputSamples[i] = sys.inputSamples[i];
	}

	outputSamples = new USER_DATATYPE[rank + 1];
	for (int i = 0; i < (rank + 1); i++)
		outputSamples[i] = sys.outputSamples[i];

	covMatrix = sys.covMatrix;
	objectParam = sys.objectParam;

	return *this;
}

void SystemIdentification::insertNewSamples(USER_DATATYPE inSample,
USER_DATATYPE outSample) {
	samplesIterator = (samplesIterator + 1) % (rank + 1);
	if (identMethod != SYSTEM_IDENTYFICATION_METHOD_RIV)
		inputSamples[samplesIterator] = inSample;
	else {
		inputSamplesIterator = (inputSamplesIterator + 1) % (2 * rank + 1);
		inputSamples[inputSamplesIterator] = inSample;
	}
	outputSamples[samplesIterator] = outSample;
}

USER_DATATYPE SystemIdentification::getInputSample(int index) {
	int tempIndex;

	if (identMethod != SYSTEM_IDENTYFICATION_METHOD_RIV) {
		tempIndex = samplesIterator + index;
	} else
		tempIndex = inputSamplesIterator + index;

	while (tempIndex < 0)
		if (identMethod != SYSTEM_IDENTYFICATION_METHOD_RIV)
			tempIndex += rank + 1;
		else
			tempIndex += 2 * rank + 1;

	return inputSamples[tempIndex];
}

USER_DATATYPE SystemIdentification::getOutputSample(int index) {
	int tempIndex = samplesIterator + index;
	while (tempIndex < 0)
		tempIndex += rank + 1;
	return outputSamples[tempIndex];
}

void SystemIdentification::identStep(USER_DATATYPE inputSample,
USER_DATATYPE outputSample) {
	insertNewSamples(inputSample, outputSample);
	if (identMethod == SYSTEM_IDENTYFICATION_METHOD_RLS)
		RLSstep();
	else
		RIVstep();
}

void SystemIdentification::print() {
	printf("\nSystem parameters:\n");
	this->objectParam.print();
	printf("\nCovariance matrix:\n");
	this->covMatrix.print();
	printf("\nRank:%d\nForgetting factor:%4.3f\n", this->rank,
			this->forgettingFactor);

	printf("\nLast input samples:");
	int tempIndex = 0;
	if (identMethod != SYSTEM_IDENTYFICATION_METHOD_RIV) {
		for (int i = 0; i < (this->rank + 1); i++) {
			printf("%4.2f, ", this->getInputSample(tempIndex));
			tempIndex--;
		}
	} else {
		for (int i = 0; i < (2 * this->rank + 1); i++) {
			printf("%4.2f, ", this->getInputSample(tempIndex));
			tempIndex--;
		}
	}

	printf("\nLast output samples:");
	tempIndex = 0;
	for (int i = 0; i < (this->rank + 1); i++) {
		printf("%4.2f, ", this->getOutputSample(tempIndex));
		tempIndex--;
	}
}

int SystemIdentification::getRank()
{
	return rank;
}

bool SystemIdentification::isCorrect()
{
	return (covMatrix.isCorrect() && objectParam.isCorrect());
}

void SystemIdentification::RLSstep() {
	MatrixIdent fiVec(2 * rank, 1);

	int tempIterator = -1;
	for (int i = 0; i < rank; i++) {
		fiVec.insert(-1 * getOutputSample(tempIterator), i, 0);
		fiVec.insert(getInputSample(tempIterator), rank + i, 0);
		tempIterator--;
	}

	MatrixIdent fiVecT = fiVec.trans();

	//Error
	USER_DATATYPE error = getOutputSample(0)
		- (fiVecT * (objectParam)).toScalar();

	//Dominator
	USER_DATATYPE dominator = forgettingFactor
		+ (fiVecT * (covMatrix * fiVec)).toScalar();

	if (dominator == 0) dominator = 0.0001;

	//New system parameters
	objectParam = objectParam
		+ (covMatrix * fiVec * (error / dominator));

	//New covariance matrix
	covMatrix = (covMatrix)
		-((covMatrix)) * fiVec * fiVecT * (covMatrix)
		/ dominator;
}

void SystemIdentification::RIVstep() {
	MatrixIdent fiVec(2 * rank, 1);

	int tempIterator = -1;
	for (int i = 0; i < rank; i++) {
		fiVec.insert(-1 * getOutputSample(tempIterator), i, 0);
		fiVec.insert(getInputSample(tempIterator), rank + i, 0);
		tempIterator--;
	}

	MatrixIdent zVec(2 * rank, 1);

	tempIterator = -1;
	for (int i = 0; i < 2 * rank; i++) {
		zVec.insert(getInputSample(tempIterator), i, 0);
		tempIterator--;
	}

	MatrixIdent fiVecT = fiVec.trans();

	//Error
	USER_DATATYPE error = getOutputSample(0)
			- (fiVecT * (objectParam)).toScalar();

	//Dominator
	USER_DATATYPE dominator = forgettingFactor
			+ (fiVecT * (covMatrix * zVec)).toScalar();

	if(dominator==0) dominator=0.0001;

	//New system parameters
	objectParam = objectParam
			+ (covMatrix * zVec * (error / dominator));

	//New covariance matrix
	covMatrix = (covMatrix)
			- ((covMatrix)) * zVec * fiVecT * (covMatrix)
					/ dominator;
}

