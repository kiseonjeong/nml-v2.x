#include "stdafx.h"
#include "Interpolation.h"

namespace nml
{
	void interpolation::checkEmpty(const algmat& X, const algmat& Y) const
	{
		// Check the matrices are empty or not
		assert(X.empty() == false && Y.empty() == false);
	}

	void interpolation::checkLength(const algmat& X, const algmat& Y) const
	{
		// Check the length of the matrices
		assert(X.dm.length() > 1 && Y.dm.length() > 1);
	}

	void interpolation::compareLength(const algmat& X, const algmat& Y) const
	{
		// Compare the length between the matrices
		assert(X.dm.length() == Y.dm.length());
	}

	const algmat interpolation::coeff() const
	{
		// Get the estimated coefficients
		return C;
	}

	vandermonde::vandermonde()
	{
		// Set an object
		setObject();
	}

	vandermonde::vandermonde(const algmat& X, const algmat& Y)
	{
		// Set an object
		setObject();

		// Estimate the coefficients
		estimate(X, Y);
	}

	vandermonde::vandermonde(const vandermonde& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	vandermonde::~vandermonde()
	{
		// Clear the object
		clearObject();
	}

	vandermonde& vandermonde::operator=(const vandermonde& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void vandermonde::setObject()
	{
		// Set the parameters
		setType(*this);
		deg = 0;

		// Set the memories
		X.release();
		Y.release();
		C.release();
	}

	void vandermonde::copyObject(const object& obj)
	{
		// Do down casting
		vandermonde* _obj = (vandermonde*)&obj;

		// Copy the parameters
		deg = _obj->deg;

		// Copy the memories
		X = _obj->X;
		Y = _obj->Y;
		C = _obj->C;
	}

	void vandermonde::clearObject()
	{
		// Clear the memories
		X.release();
		Y.release();
		C.release();
	}

	void vandermonde::estimate(const algmat& X, const algmat& Y)
	{
		// Check the matrices is empty or not
		checkEmpty(X, Y);

		// Check the length of the matrices
		checkLength(X, Y);

		// Compare the length between the matrices
		compareLength(X, Y);

		// Clear the old memories
		this->X.release();
		this->Y.release();
		C.release();

		// Backup the input vectors
		this->X = X.reshape(X.dm.length());
		this->Y = Y.reshape(Y.dm.length());

		// Set the degree
		deg = X.dm.length() - 1;

		// Calculate the coefficients
		calculate();
	}

	void vandermonde::calculate()
	{
		// Create the Vandermonde matrix X
		algmat vX(msize(deg + 1, deg + 1));
		for (int i = 0; i < vX.rows; i++)
		{
			for (int j = 0; j < vX.cols; j++)
			{
				vX[i][j] = pow(X[i], deg - j);
			}
		}

		// Create the Vandermonde matrix Y
		algmat vY(msize(deg + 1, 1));
		for (int i = 0; i < vY.dm.length(); i++)
		{
			vY[i] = Y[i];
		}

		// Calculate the coefficients
		C = vX.inv().dot(vY);
	}

	const double vandermonde::interpolate(const double x) const
	{
		// Check the coefficients
		assert(C.empty() == false);

		// Interpolate a value on the input X
		double y = 0.0;
		for (int i = 0; i < C.dm.length(); i++)
		{
			y += C[i] * pow(x, deg - i);
		}

		return y;
	}

	newton::newton()
	{
		// Set an object
		setObject();
	}

	newton::newton(const algmat& X, const algmat& Y)
	{
		// Set an object
		setObject();

		// Estimate the coefficients
		estimate(X, Y);
	}

	newton::newton(const newton& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	newton::~newton()
	{
		// Clear the object
		clearObject();
	}

	newton& newton::operator=(const newton& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void newton::setObject()
	{
		// Set the parameters
		setType(*this);
		deg = 0;
		idx = 0;

		// Set the memories
		X.release();
		Y.release();
		C.release();
	}

	void newton::copyObject(const object& obj)
	{
		// Do down casting
		newton* _obj = (newton*)&obj;

		// Copy the parameters
		deg = _obj->deg;
		idx = _obj->idx;

		// Copy the memories
		X = _obj->X;
		Y = _obj->Y;
		C = _obj->C;
	}

	void newton::clearObject()
	{
		// Clear the memories
		X.release();
		Y.release();
		C.release();
	}

	void newton::estimate(const algmat& X, const algmat& Y)
	{
		// Check the matrices is empty or not
		checkEmpty(X, Y);

		// Check the length of the matrices
		checkLength(X, Y);

		// Compare the length between the matrices
		compareLength(X, Y);

		// Clear the old memories
		this->X.release();
		this->Y.release();
		C.release();

		// Backup the input vectors
		this->X = X.reshape(X.dm.length());
		this->Y = Y.reshape(Y.dm.length());

		// Set the degree
		deg = X.dm.length() - 1;

		// Calculate the coefficients
		idx = 0;
		C = algmat::zeros(msize(deg + 1, 1));
		calculate(this->X, this->Y);
	}

	void newton::calculate(const algmat& X, const algmat& Y)
	{
		// Calculate a divided difference value
		algmat diff(msize(1, Y.dm.length() - 1));
		const int offset = X.dm.length() - Y.dm.length();
		for (int i = 0; i < diff.dm.length(); i++)
		{
			diff[i] = slope(X[i], Y[i], X[i + 1 + offset], Y[i + 1]);
		}

		// Check the index
		if (idx == 0)
		{
			C[idx] = Y[idx];
			idx++;
		}
		C[idx] = diff[0];
		idx++;

		// Calculate the others
		if (diff.dm.length() != 1)
		{
			calculate(X, diff);
		}
	}

	const double newton::slope(const double x0, const double y0, const double x1, const double y1) const
	{
		// Calculate the slope
		return (y1 - y0) / (x1 - x0);
	}

	const double newton::interpolate(const double x) const
	{
		// Check the coefficients
		assert(C.empty() == false);

		// Interpolate a value on the input X
		double y = C[0];
		for (int i = 1; i < C.dm.length(); i++)
		{
			double mp = x - X[0];
			for (int j = 1; j < i; j++)
			{
				mp *= x - X[j];
			}
			y += C[i] * mp;
		}

		return y;
	}

	lagrange::lagrange()
	{
		// Set an object
		setObject();
	}

	lagrange::lagrange(const algmat& X, const algmat& Y)
	{
		// Set an object
		setObject();

		// Estimate the coefficients
		estimate(X, Y);
	}

	lagrange::lagrange(const lagrange& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	lagrange::~lagrange()
	{
		// Clear the object
		clearObject();
	}

	lagrange& lagrange::operator=(const lagrange& obj)
	{
		// Clear the object
		clearObject();

		// Set an object
		setObject();

		// Copy the object
		copyObject(obj);

		return *this;
	}

	void lagrange::setObject()
	{
		// Set the parameters
		setType(*this);
		deg = 0;

		// Set the memories
		X.release();
		Y.release();
		C.release();
	}

	void lagrange::copyObject(const object& obj)
	{
		// Do down casting
		lagrange* _obj = (lagrange*)&obj;

		// Copy the parameters
		deg = _obj->deg;

		// Copy the memories
		X = _obj->X;
		Y = _obj->Y;
		C = _obj->C;
	}

	void lagrange::clearObject()
	{
		// Clear the memories
		X.release();
		Y.release();
		C.release();
	}

	void lagrange::estimate(const algmat& X, const algmat& Y)
	{
		// Check the matrices is empty or not
		checkEmpty(X, Y);

		// Check the length of the matrices
		checkLength(X, Y);

		// Compare the length between the matrices
		compareLength(X, Y);

		// Clear the old memories
		this->X.release();
		this->Y.release();
		C.release();

		// Backup the vectors
		this->X = X.reshape(X.dm.length());
		this->Y = Y.reshape(Y.dm.length());

		// Set the degree
		deg = X.dm.length() - 1;
	}

	const double lagrange::interpolate(const double x) const
	{
		// Check the backup vectors are empty or not
		assert(X.empty() == false && Y.empty() == false);

		// Interpolate a value on the input X
		double y = 0.0;
		for (int i = 0; i < deg + 1; i++)
		{
			y += calculate(x, i) * Y[i];
		}

		return y;
	}

	const double lagrange::calculate(const double x, const int idx) const
	{
		// Calculate an weighted value
		double weighted = 1.0;
		for (int i = 0; i < deg + 1; i++)
		{
			if (i != idx)
			{
				weighted *= (x - X[i]) / (X[idx] - X[i]);
			}
		}

		return weighted;
	}

	const algmat lagrange::coeff() const
	{
		// Do nothing
		return algmat();
	}
}