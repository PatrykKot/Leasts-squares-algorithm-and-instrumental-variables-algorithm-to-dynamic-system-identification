/*
 * rank_identification.hpp
 *
 *  Created on: Apr 29, 2016
 *      Author: Patryk
 */

#ifndef RANK_IDENTIFICATION_HPP_
#define RANK_IDENTIFICATION_HPP_

#include "system_identification.hpp"
#include <math.h>

class RankIdent {
private:
	SystemIdentification* systemTab;
	USER_DATATYPE* outputSamples;
	USER_DATATYPE* inputSamples;
	int maxRank;
	int identMethod;
	bool finished;
	int currentIndex;
	int maxSamplesNumber;
	USER_DATATYPE forgettingFactor;
	int initialCovarianceMatrix;

public:
	RankIdent();
	RankIdent(const RankIdent&);
	RankIdent(int maxRank = 5, USER_DATATYPE forgettingFactor = 1,
		USER_DATATYPE initialCovarianceMatrix = 1,
		int identMethod = SYSTEM_IDENTYFICATION_METHOD_RLS, int maxSamples = 1000);
	~RankIdent();

	bool isFinished();
	SystemIdentification getBestSystem();
	bool identStep(USER_DATATYPE inputSample, USER_DATATYPE outputSample);
	void print();
	void printAllAIC();
	bool isCorrect();
	USER_DATATYPE getAIC(int rank);
	SystemIdentification getSystemByRank(int rank);
	int getMaxRank();

private:
	USER_DATATYPE getInputSample(int index);
	USER_DATATYPE getOutputSample(int index);
};

#endif /* RANK_IDENTIFICATION_HPP_ */
