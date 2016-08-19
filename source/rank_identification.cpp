/*
 * rank_identification.cpp
 *
 *  Created on: Apr 29, 2016
 *      Author: Patryk
 */

#include "rank_identification.hpp"
#include <stdio.h>

RankIdent::RankIdent() {
	this->currentIndex = 0;
	this->finished = false;
	this->identMethod = SYSTEM_IDENTYFICATION_METHOD_RLS;
	this->maxRank = 5;
	this->maxSamplesNumber = 1024;
	this->forgettingFactor = 1;
	this->initialCovarianceMatrix = 1;

	this->systemTab = new SystemIdentification[this->maxRank];
	for (int i = 0; i < this->maxRank; i++) {
		SystemIdentification tempSystem(i + 1, this->forgettingFactor,
				this->initialCovarianceMatrix, this->identMethod);
		this->systemTab[i] = tempSystem;
	}

	this->inputSamples = new USER_DATATYPE[this->maxSamplesNumber];
	this->outputSamples = new USER_DATATYPE[this->maxSamplesNumber];
}

RankIdent::RankIdent(const RankIdent& sys) {
	this->currentIndex = sys.currentIndex;
	this->finished = sys.finished;
	this->identMethod = sys.identMethod;
	this->maxRank = sys.maxRank;
	this->maxSamplesNumber = sys.maxSamplesNumber;
	this->forgettingFactor = sys.forgettingFactor;

	this->systemTab = new SystemIdentification[this->maxRank];
	for (int i = 0; i < this->maxRank; i++) {
		this->systemTab[i] = sys.systemTab[i];
	}

	this->inputSamples = new USER_DATATYPE[this->maxSamplesNumber];
	for (int i = 0; i < this->maxSamplesNumber; i++) {
		this->inputSamples[i] = sys.inputSamples[i];
	}
	this->outputSamples = new USER_DATATYPE[this->maxSamplesNumber];
	for (int i = 0; i < this->maxSamplesNumber; i++) {
		this->outputSamples[i] = sys.outputSamples[i];
	}
}

RankIdent::RankIdent(int maxRank, USER_DATATYPE forgettingFactor,
USER_DATATYPE initialCovarianceMatrix, int identMethod, int maxSamples) {
	this->currentIndex = 0;
	this->finished = false;
	this->identMethod = identMethod;
	this->maxRank = maxRank;
	this->maxSamplesNumber = maxSamples;
	this->forgettingFactor = forgettingFactor;
	this->initialCovarianceMatrix = initialCovarianceMatrix;

	this->systemTab = new SystemIdentification[this->maxRank];
	
	for (int i = 0; i < this->maxRank; i++) {
		SystemIdentification tempSystem(i + 1, this->forgettingFactor,
				this->initialCovarianceMatrix, this->identMethod);
				
		this->systemTab[i] = tempSystem;
	}

	this->inputSamples = new USER_DATATYPE[this->maxSamplesNumber];
	this->outputSamples = new USER_DATATYPE[this->maxSamplesNumber];
}

RankIdent::~RankIdent() {
	delete[] this->inputSamples;
	delete[] this->outputSamples;	
	delete[] systemTab;
}

bool RankIdent::isFinished() {
	return finished;
}

SystemIdentification RankIdent::getBestSystem() {
	int minIndex=0;
	USER_DATATYPE lastAIC=this->getAIC(1);
	for (int i = 1; i < maxRank; i++)
		{
			USER_DATATYPE newAIC=this->getAIC(i+1);
			if (lastAIC > newAIC)
			{
				minIndex = i;
				lastAIC=newAIC;
			}
			else
			{
				return systemTab[minIndex];
			}
		}

	return systemTab[minIndex];
}

bool RankIdent::identStep(USER_DATATYPE inputSample,
USER_DATATYPE outputSample) {
	if (!finished) {
		inputSamples[currentIndex] = inputSample;
		outputSamples[currentIndex++] = outputSample;

		for (int i = 0; i < maxRank; i++)
			{
				systemTab[i].identStep(inputSample, outputSample);
			}

		if (currentIndex > maxSamplesNumber)

		{
			finished = true;
		}
	}

	return finished;
}

USER_DATATYPE RankIdent::getAIC(int rank) {
	if (rank < 0 || rank > this->maxRank)
		return 0;

	//Simulation
	int outLength = this->currentIndex;
	USER_DATATYPE* simulatedOutput = new USER_DATATYPE[outLength];
	for (int i = 0; i < outLength; i++)
		simulatedOutput[i] = 0;

	for (int outIterator = 0; outIterator < outLength; outIterator++) {
		for (int currOut = -1; currOut >= (-1 * rank); currOut--) {
			if ((outIterator + currOut) >= 0) {
				simulatedOutput[outIterator] += -1
					* simulatedOutput[outIterator + currOut]
					* this->systemTab[rank - 1].objectParam.elementAt(
					-1 * currOut - 1, 0);
			}
		}

		for (int currIn = -1; currIn >= (-1 * rank); currIn--) {
			simulatedOutput[outIterator] += this->getInputSample(
				outIterator + currIn)
				* this->systemTab[rank - 1].objectParam.elementAt(
				-1 * currIn - 1 + rank, 0);
		}

	}

	//Calculating Mean Square Error
	USER_DATATYPE mse = 0;
	for (int i = 0; i < outLength; i++) {
		simulatedOutput[i] -= this->getOutputSample(i);
		mse += simulatedOutput[i] * simulatedOutput[i];
	}
	delete[] simulatedOutput;

	//Calculating Akaike's Information Criterion
	return log((double)mse) + 6 * rank / outLength;
}

SystemIdentification RankIdent::getSystemByRank(int rank)
{
	return systemTab[rank - 1];
}

USER_DATATYPE RankIdent::getInputSample(int index) {
	if (index < 0)
		return 0;
	if (index > (this->currentIndex))
		return 0;

	return this->inputSamples[index];
}

USER_DATATYPE RankIdent::getOutputSample(int index) {
	if (index < 0)
		return 0;
	if (index > (this->currentIndex))
		return 0;

	return this->outputSamples[index];
}

void RankIdent::print()
{
	printf("\nRANK IDENTIFICATION PRINT\n");
	for(int i=0;i<this->maxRank;i++)
{
	printf("\n\nRANK:%d\n",i+1);
	this->systemTab[i].print();
}
printf("\nEND OF PRINTING RANK IDENTIFICATION\n");
}

void RankIdent::printAllAIC()
{
	printf("\nAkaike's Information Criterion vector printing:\n");
	for(int i=0;i<this->maxRank;i++)
	{
		printf("RANK:%d\tAIC=%f\n",i+1,this->getAIC(i+1));
	}
}

bool RankIdent::isCorrect()
{
	for (int i = 0; i < maxRank; i++)
	{
		if (!systemTab[i].isCorrect()) return false;
	}

	return true;
}

int RankIdent::getMaxRank()
{
	return maxRank;
}