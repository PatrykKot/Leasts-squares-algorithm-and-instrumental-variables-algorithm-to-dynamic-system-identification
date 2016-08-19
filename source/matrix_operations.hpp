/*
 * matrix_operations.h
 *
 *  Created on: Mar 30, 2016
 *      Author: Patryk
 */

#ifndef MATRIX_OPERATIONS_H_
#define MATRIX_OPERATIONS_H_

#include <stdio.h>

#define USER_DATATYPE long double

class MatrixIdent {
private:
	int rows;
	int cols;
	USER_DATATYPE* data;
public:
	MatrixIdent();
	MatrixIdent(int rows, int cols);
	MatrixIdent(const MatrixIdent & mat);
	~MatrixIdent();

	USER_DATATYPE getRows();
	USER_DATATYPE getCols();
	USER_DATATYPE elementAt(int row, int col) const;

	MatrixIdent operator+(MatrixIdent const& mat) const;
	MatrixIdent operator-(MatrixIdent const& mat) const;
	MatrixIdent operator*(MatrixIdent const& mat) const;
	MatrixIdent operator*(USER_DATATYPE val) const;
	MatrixIdent & operator=(const MatrixIdent & mat);
	friend MatrixIdent operator*(USER_DATATYPE val,const MatrixIdent& mat);
	MatrixIdent operator/(USER_DATATYPE val) const;
	MatrixIdent trans() const;

	USER_DATATYPE toScalar ();
	void insert(USER_DATATYPE val, int row, int col);
	bool isCorrect();

	void print();
};


#endif
