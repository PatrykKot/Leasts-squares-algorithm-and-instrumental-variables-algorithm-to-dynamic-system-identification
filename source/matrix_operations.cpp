/*
 * matrix_operations.c
 *
 *  Created on: Mar 30, 2016
 *      Author: Patryk
 */

#include "matrix_operations.hpp"

MatrixIdent::MatrixIdent() :
		rows(0), cols(0), data(NULL) {
}

MatrixIdent::MatrixIdent(int rows, int cols) {
	this->rows = rows;
	this->cols = cols;
	data = new USER_DATATYPE[rows * cols];
	for (int i = 0; i < (rows * cols); i++)
		data[i] = 0;
}

MatrixIdent::~MatrixIdent() {
	delete[] data;
}

MatrixIdent::MatrixIdent(const MatrixIdent & mat) {
	this->cols = mat.cols;
	this->rows = mat.rows;
	this->data = new USER_DATATYPE[this->cols * this->rows];
	for (int i = 0; i < this->rows; i++)
		for (int j = 0; j < this->cols; j++)
			this->insert(mat.elementAt(i, j), i, j);
}

USER_DATATYPE MatrixIdent::elementAt(int row, int col) const {
	return data[row * cols + col];
}

void MatrixIdent::insert(USER_DATATYPE val, int row, int col) {
	data[row * cols + col] = val;
}

MatrixIdent MatrixIdent::operator+(MatrixIdent const& mat) const {
	if (this->cols != mat.cols || this->rows != mat.rows)
		return mat;

	MatrixIdent tempMat(this->rows, this->cols);
	for (int i = 0; i < this->rows; i++)
		for (int j = 0; j < this->cols; j++)
			tempMat.insert(this->elementAt(i, j) + mat.elementAt(i, j), i, j);
	return tempMat;
}

MatrixIdent MatrixIdent::operator-(MatrixIdent const& mat) const {
	if (this->cols != mat.cols || this->rows != mat.rows)
		return mat;

	MatrixIdent tempMat(this->rows, this->cols);
	for (int i = 0; i < this->rows; i++)
		for (int j = 0; j < this->cols; j++)
			tempMat.insert(this->elementAt(i, j) - mat.elementAt(i, j), i, j);
	return tempMat;
}

MatrixIdent MatrixIdent::operator*(MatrixIdent const& mat) const {
	if (this->cols != mat.rows) {
		return mat;
	}
	MatrixIdent tempMat(this->rows, mat.cols);

	for (int col = 0; col < mat.cols; col++)
		for (int row = 0; row < this->rows; row++) {
			USER_DATATYPE tempVal = 0;
			for (int iterator = 0; iterator < this->cols; iterator++)
				tempVal += this->elementAt(row, iterator)
						* mat.elementAt(iterator, col);
			tempMat.insert(tempVal, row, col);
		}

	return tempMat;

}

USER_DATATYPE MatrixIdent::getRows()
{
	return rows;
}

USER_DATATYPE MatrixIdent::getCols()
{
	return cols;
}

MatrixIdent MatrixIdent::operator*(USER_DATATYPE val) const {
	MatrixIdent tempMat(this->rows, this->cols);
	for (int i = 0; i < this->rows; i++)
		for (int j = 0; j < this->cols; j++)
			tempMat.insert(this->elementAt(i, j) * val, i, j);
	return tempMat;
}

MatrixIdent& MatrixIdent::operator=(const MatrixIdent & mat)
{
	if(this->data!=NULL) delete[] data;
	
	this->cols=mat.cols;
	this->rows=mat.rows;
	this->data=new USER_DATATYPE[this->cols*this->rows];
	for (int i = 0; i < this->rows; i++)
		for (int j = 0; j < this->cols; j++)
			this->insert(mat.elementAt(i, j), i, j);

	return *this;
}

MatrixIdent operator*(USER_DATATYPE val, const MatrixIdent& mat) {
	return mat * val;
}

MatrixIdent MatrixIdent::trans() const {
	MatrixIdent tempMat(this->cols, this->rows);
	for (int i = 0; i < this->rows; i++)
		for (int j = 0; j < this->cols; j++)
			tempMat.insert(this->elementAt(i, j), j, i);
	return tempMat;
}

USER_DATATYPE MatrixIdent::toScalar() {
	return data[0];
}

MatrixIdent MatrixIdent::operator/(USER_DATATYPE val) const {
	return (*this) * (1 / val);
}

void MatrixIdent::print() {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++)
			printf("%4.7f\t", this->elementAt(i, j));
		printf("\n");
	}
}

bool MatrixIdent::isCorrect()
{
	for (int r = 0; r < rows; r++)
		for (int c = 0; c < cols; c++)
			if (elementAt(r, c) != elementAt(r, c))
				return false;
	return true;
}

