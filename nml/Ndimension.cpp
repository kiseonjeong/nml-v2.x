#include "stdafx.h"
#include "Ndimension.h"

namespace nml
{
	dim::dim() : d(_d), N(_N)
	{
		// Initialize information
		_N = 0;
		len = 0;
	}

	dim::dim(const int N, ...) : d(_d), N(_N)
	{
		// Create the new information
		_d.create(N);

		// Set dimension information
		_N = N;

		// Parse dimension values
		va_list vl;
		va_start(vl, N);
		for (int i = 0; i < N; i++)
		{
			_d[i] = va_arg(vl, int);
		}
		va_end(vl);

		// Calculate a total length
		len = 1;
		for (int i = 0; i < _d.size; i++)
		{
			len *= _d[i];
		}
	}

	dim::dim(const dim& obj) : d(_d), N(_N)
	{
		// Copy dimension information
		_d = obj._d;
		_N = obj._N;
		len = obj.len;
	}

	dim::~dim()
	{
		// Clear the dimension information
		_d.release();
	}

	dim& dim::operator=(const dim& obj)
	{
		// Copy from the input object
		_d = obj._d;
		_N = obj._N;
		len = obj.len;

		return *this;
	}

	const int& dim::operator[](const int idx) const
	{
		// Check the index
		assert(idx >= 0 && idx < _d.size);

		// Get a dimension length
		return _d[idx];
	}

	void dim::create(const int N)
	{
		// Set dimension information
		_N = N;

		// Create dimension information
		_d.create(N);
	}

	void dim::release()
	{
		// Initialize information
		_d.release();
		_N = 0;
	}

	const int dim::length() const
	{
		// Get a total length
		return len;
	}

	const int dim::get(const int idx) const
	{
		// Check the index
		assert(idx >= 0 && idx < _d.size);

		// Get a dimension length
		return _d[idx];
	}

	void dim::set(const int idx, const int len)
	{
		// Check the index
		assert(idx >= 0 && idx < _d.size);

		// Set a dimension length
		_d[idx] = len;

		// Update the total length
		this->len = 1;
		for (int i = 0; i < _d.size; i++)
		{
			this->len *= _d[i];
		}
	}

	void dim::reset(const int N, ...)
	{
		// Create the new information
		_d.create(N);

		// Set dimension information
		_N = N;

		// Parse dimension values
		va_list vl;
		va_start(vl, N);
		for (int i = 0; i < N; i++)
		{
			_d[i] = va_arg(vl, int);
		}
		va_end(vl);

		// Calculate a total length
		len = 1;
		for (int i = 0; i < _d.size; i++)
		{
			len *= _d[i];
		}
	}

	void dim::cout() const
	{
		// Check dimension information
		std::cout << "(";
		for (int i = 0; i < _N; i++)
		{
			std::cout << _d[i];
			if (i < _N - 1)
			{
				std::cout << ", ";
			}
		}
		std::cout << ")" << endl;
	}
}