#include "stdafx.h"
#include "Msize.h"

namespace nml
{
	msize::msize()
	{
		d0 = -1;
		d1 = -1;
	}

	msize::msize(const int length)
	{
		// Set a 1d size information
		set(length);
	}

	msize::msize(const int rows, const int cols)
	{
		// Set a 2d size information
		set(rows, cols);
	}

	msize::~msize()
	{

	}

	msize& msize::operator=(const msize& obj)
	{
		// Copy from the input object
		d0 = obj.d0;
		d1 = obj.d1;

		return *this;
	}

	void msize::set(const int length)
	{
		// Set the dimension parameters
		d0 = length;
		d1 = 1;
	}

	void msize::set(const int rows, const int cols)
	{
		// Set the dimension parameters
		d0 = rows;
		d1 = cols;
	}
}