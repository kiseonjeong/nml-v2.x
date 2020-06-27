#include "stdafx.h"
#include "Algmat.h"

namespace nml
{
	algmat::algmat() : rows(sz.d0), cols(sz.d1)
	{
		// Set size information
		sz.set(0, 0);
	}

	algmat::algmat(const msize& sz) : rows(this->sz.d0), cols(this->sz.d1)
	{
		// Set dimension information
		this->sz = sz;
		dim dm(2, sz.d0, sz.d1);

		// Create a matrix
		create(dm);
	}

	algmat::algmat(const msize& sz, const double val) : rows(this->sz.d0), cols(this->sz.d1)
	{
		// Set dimension information
		this->sz = sz;
		dim dm(2, sz.d0, sz.d1);

		// Create a matrix
		create(dm);

		// Set the matrix
		set(val);
	}

	algmat::algmat(const ndmatrix<2>& obj) : rows(sz.d0), cols(sz.d1)
	{
		// Copy the object
		copyObject(obj);
	}

	algmat::algmat(const algmat& obj) : rows(sz.d0), cols(sz.d1)
	{
		// Clone the object
		*this = obj;
	}

	algmat::~algmat()
	{
		// Do nothing
	}

	algmat& algmat::operator=(const algmat& obj)
	{
		// Copy the object
		copyObject(obj);

		return *this;
	}

	void algmat::copyObject(const ndarray& obj)
	{
		// Do down casting
		ndmatrix<2>::copyObject(obj);

		// Copy the parameters
		sz.set(obj.dm[0], obj.dm[1]);
	}

	algmat algmat::to_algmat(const vector<double>& vdata)
	{
		// Check a vector data
		if (vdata.empty() == true)
		{
			return algmat();
		}

		// Convert a vector data to the matrix
		algmat result(msize((int)vdata.size()));
		for (int i = 0; i < result.tlen; i++)
		{
			result.ddata[i] = vdata[i];
		}

		return result;
	}

	algmat algmat::to_algmat(const vector<double>& vdata, const int nrows)
	{
		// Check a vector data
		if (vdata.empty() == true)
		{
			return algmat();
		}

		// Convert a vector data to the matrix
		algmat temp = to_algmat(vdata);

		// Reshape the matrix
		return temp.reshape(nrows);
	}

	algmat algmat::append(const vector<algmat>& vmat, const int axis)
	{
		// Check an input vector
		assert(vmat.empty() == false);

		// Check an axis information
		if (axis == 0)
		{
			// Append a matrix
			return vappend(vmat);
		}
		else
		{
			// Append a matrix
			return happend(vmat);
		}
	}

	const algmat algmat::vappend(const vector<algmat>& vmat)
	{
		// Check dimension information
		const int nmat = (int)vmat.size();
		int rows = 0;
		int cols = vmat[0].cols;
		for (int i = 0; i < nmat; i++)
		{
			if (i != 0)
			{
				assert(vmat[0].cols == vmat[i].cols);
			}
			rows += vmat[i].rows;
		}

		// Vertically, append a matrix
		algmat result(msize(rows, cols));
		for (int i = 0, j = 0; i < nmat; i++)
		{
			for (int k = 0; k < vmat[i].rows; k++)
			{
				for (int l = 0; l < vmat[i].cols; l++)
				{
					result[j + k][l] = vmat[i][k][l];
				}
			}
			j += vmat[i].rows;
		}

		return result;
	}

	const algmat algmat::happend(const vector<algmat>& vmat)
	{
		// Check dimension information
		const int nmat = (int)vmat.size();
		int rows = vmat[0].rows;
		int cols = 0;
		for (int i = 0; i < nmat; i++)
		{
			if (i != 0)
			{
				assert(vmat[0].rows == vmat[i].rows);
			}
			cols += vmat[i].cols;
		}

		// Horizontally, append a matrix
		algmat result(msize(rows, cols));
		for (int i = 0, j = 0; i < nmat; i++)
		{
			for (int k = 0; k < vmat[i].cols; k++)
			{
				for (int l = 0; l < vmat[i].rows; l++)
				{
					result[l][j + k] = vmat[i][l][k];
				}
			}
			j += vmat[i].cols;
		}

		return result;
	}

	algmat algmat::append(const algmat& mat, const int axis) const
	{
		// Set a matrix vector
		vector<algmat> vmat;
		vmat.push_back(*this);
		vmat.push_back(mat);

		// Append a matrix
		return append(vmat, axis);
	}

	algmat algmat::submat(const int idx, const int axis) const
	{
		// Get a selected sub-matrix
		return submat(idx, idx, -1, axis);
	}

	algmat algmat::submat(const int bi, const int ei, const int ii, const int axis) const
	{
		// Check matrix status
		assert(empty() == false);

		// Check an axis information
		if (axis == 0)
		{
			// Get a selected sub-matrix
			return vsubmat(bi, ei, ii);
		}
		else
		{
			// Get a selected sub-matrix
			return hsubmat(bi, ei, ii);
		}
	}

	const algmat algmat::vsubmat(int bi, int ei, int ii) const
	{
		// Check the indices
		assert(bi >= 0 && bi < rows);
		assert(ei >= 0 && ei < rows);
		assert(bi <= ei);
		if (ii < 0)
		{
			ii = 1;
		}

		// Create a selected row matrix
		algmat result(msize((ei - bi) / ii + 1, cols));
		for (int i = bi, k = 0; i <= ei; i += ii, k++)
		{
			for (int j = 0; j < cols; j++)
			{
				result[k][j] = (*this)[i][j];
			}
		}

		return result;
	}

	const algmat algmat::hsubmat(int bi, int ei, int ii) const
	{
		// Check the indices
		assert(bi >= 0 && bi < cols);
		assert(ei >= 0 && ei < cols);
		assert(bi <= ei);
		if (ii < 0)
		{
			ii = 1;
		}

		// Create a selected column matrix
		algmat result(msize(rows, (ei - bi) / ii + 1));
		for (int i = 0; i < rows; i++)
		{
			for (int j = bi, k = 0; j <= ei; j += ii, k++)
			{
				result[i][k] = (*this)[i][j];
			}
		}

		return result;
	}

	algmat algmat::minor(const int idx, const int axis) const
	{
		// Check matrix status
		assert(empty() == false);

		// Check an axis information
		if (axis == 0)
		{
			// Swap the matrix
			return vminor(idx);
		}
		else
		{
			// Swap the matrix
			return hminor(idx);
		}
	}

	algmat algmat::reshape(const int nrows) const
	{
		// Check matrix status
		assert(empty() == false);

		// Reshape a matrix
		return ndmatrix::reshape<2>(dim(2, nrows, -1));
	}

	algmat algmat::t() const
	{
		// Check matrix status
		assert(empty() == false);

		// Transpose a matrix
		algmat result(msize(cols, rows));
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				result[j][i] = (*this)[i][j];
			}
		}

		return result;
	}

	const algmat algmat::vminor(const int idx) const
	{
		// Check an index
		assert(idx >= 0 && idx < rows);

		// Check the row size
		if (rows < 2)
		{
			return *this;
		}

		// Get a minor matrix except the selected index
		algmat result(msize(rows - 1, cols));
		for (int i = 0, k = 0; i < rows; i++)
		{
			// Check an index
			if (i == idx)
			{
				continue;
			}

			// Copy the data
			for (int j = 0; j < cols; j++)
			{
				result[k][j] = (*this)[i][j];
			}
			k++;
		}

		return result;
	}

	const algmat algmat::hminor(const int idx) const
	{
		// Check an index
		assert(idx >= 0 && idx < cols);

		// Check the column size
		if (cols < 2)
		{
			return *this;
		}

		// Get a minor matrix except the selected index
		algmat result(msize(rows, cols - 1));
		for (int i = 0, k = 0; i < cols; i++)
		{
			// Check an index
			if (i == idx)
			{
				continue;
			}

			// Copy the data
			for (int j = 0; j < rows; j++)
			{
				result[j][k] = (*this)[j][i];
			}
			k++;
		}

		return result;
	}

	void algmat::swap(const int fi, const int ti, const int axis) const
	{
		// Check matrix status
		assert(empty() == false);

		// Check an axis information
		if (axis == 0)
		{
			// Swap the matrix
			vswap(fi, ti);
		}
		else
		{
			// Swap the matrix
			hswap(fi, ti);
		}
	}

	void algmat::vswap(const int fi, const int ti) const
	{
		// Check the indices
		assert(fi >= 0 && fi < rows);
		assert(ti >= 0 && ti < rows);

		// Swap the matrix
		algmat temp = submat(fi, 0);
		for (int i = 0; i < cols; i++)
		{
			(*this)[fi][i] = (*this)[ti][i];
			(*this)[ti][i] = temp[0][i];
		}
	}

	void algmat::hswap(const int fi, const int ti) const
	{
		// Check the indices
		assert(fi >= 0 && fi < cols);
		assert(ti >= 0 && ti < cols);

		// Swap the matrix
		algmat temp = submat(fi, 1);
		for (int i = 0; i < cols; i++)
		{
			(*this)[i][fi] = (*this)[i][ti];
			(*this)[i][ti] = temp[i][0];
		}
	}

	void algmat::sort(const bool desc, const int axis) const
	{
		// Check matrix status
		assert(empty() == false);

		// Check an axis information
		if (axis == 0)
		{
			// Sort the matrix
			vsort(desc);
		}
		else
		{
			// Sort the matrix
			hsort(desc);
		}
	}

	void algmat::vsort(const bool desc) const
	{
		// Transpose the matrix
		algmat trans = transpose(dim(2, 1, 0));

		// Sort the matrix
		trans.hsort(desc);
		for (int i = 0; i < trans.rows; i++)
		{
			for (int j = 0; j < trans.cols; j++)
			{
				(*this)[j][i] = trans[i][j];
			}
		}
	}

	void algmat::hsort(const bool desc) const
	{
		// Check a sorting method
		if (desc == true)
		{
			// Sort the matrix
			for (int i = 0; i < rows; i++)
			{
				std::sort(ddata + i * cols, ddata + (i + 1) * cols, greater<double>());
			}
		}
		else
		{
			// Sort the matrix
			for (int i = 0; i < rows; i++)
			{
				std::sort(ddata + i * cols, ddata + (i + 1) * cols);
			}
		}
	}

	algmat algmat::count(const double val, const int axis) const
	{
		// Check matrix status
		assert(empty() == false);

		// Check an axis information
		if (axis == 0)
		{
			// Find a value in the matrix
			return vcount(val);
		}
		else if (axis == 1)
		{
			// Find a value in the matrix
			return hcount(val);
		}
		else
		{
			return algmat::sum(vcount(val));
		}
	}

	const algmat algmat::vcount(const double val) const
	{
		// Find a value in the matrix
		algmat result(msize(1, cols), 0.0);
		for (int i = 0; i < cols; i++)
		{
			for (int j = 0; j < rows; j++)
			{
				if ((*this)[j][i] == val)
				{
					result[0][i] += 1.0;
				}
			}
		}

		return result;
	}

	const algmat algmat::hcount(const double val) const
	{
		// Find a value in the matrix
		algmat result(msize(rows, 1), 0.0);
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				if ((*this)[i][j] == val)
				{
					result[i][0] += 1.0;
				}
			}
		}

		return result;
	}

	algmat algmat::min(const algmat& mat, const int axis)
	{
		// Check matrix status
		assert(mat.empty() == false);

		// Check an axis information
		if (axis == 0)
		{
			// Find min values
			algmat loc(msize(1, mat.cols));

			return vmin(mat, loc);
		}
		else if (axis == 1)
		{
			// Find min values
			algmat loc(msize(mat.rows, 1));

			return hmin(mat, loc);
		}
		else
		{
			// Find min values
			algmat loc0(msize(1, mat.cols));
			algmat loc1(msize(1, 1));

			return hmin(vmin(mat, loc0), loc1);
		}
	}

	const algmat algmat::vmin(const algmat& mat, algmat& loc)
	{
		// Find min values
		algmat result(msize(1, mat.cols));
		for (int i = 0; i < mat.cols; i++)
		{
			int arg = 0;
			double val = mat[0][i];
			for (int j = 1; j < mat.rows; j++)
			{
				if (val > mat[j][i])
				{
					arg = j;
					val = mat[j][i];
				}
			}
			loc[0][i] = arg;
			result[0][i] = val;
		}

		return result;
	}

	const algmat algmat::hmin(const algmat& mat, algmat& loc)
	{
		// Find min values
		algmat result(msize(mat.rows, 1));
		for (int i = 0; i < mat.rows; i++)
		{
			int arg = 0;
			double val = mat[i][0];
			for (int j = 1; j < mat.cols; j++)
			{
				if (val > mat[i][j])
				{
					arg = j;
					val = mat[i][j];
				}
			}
			loc[i][0] = arg;
			result[i][0] = val;
		}

		return result;
	}

	algmat algmat::argmin(const algmat& mat, const int axis)
	{
		// Check matrix status
		assert(mat.empty() == false);

		// Check an axis information
		if (axis == 0)
		{
			// Find min arguments
			algmat loc(msize(1, mat.cols));
			vmin(mat, loc);

			return loc;
		}
		else if (axis == 1)
		{
			// Find min arguments
			algmat loc(msize(mat.rows, 1));
			hmin(mat, loc);

			return loc;
		}
		else
		{
			// Find min arguments
			algmat loc0(msize(1, mat.cols));
			algmat loc1(msize(1, 1));
			hmin(vmin(mat, loc0), loc1);

			return loc1;
		}
	}

	algmat algmat::max(const algmat& mat, const int axis)
	{
		// Check matrix status
		assert(mat.empty() == false);

		// Check an axis information
		if (axis == 0)
		{
			// Find max values
			algmat loc(msize(1, mat.cols));

			return vmax(mat, loc);
		}
		else if (axis == 1)
		{
			// Find max values
			algmat loc(msize(mat.rows, 1));

			return hmax(mat, loc);
		}
		else
		{
			// Find max values
			algmat loc0(msize(1, mat.cols));
			algmat loc1(msize(1, 1));

			return hmax(vmax(mat, loc0), loc1);
		}
	}

	const algmat algmat::vmax(const algmat& mat, algmat& loc)
	{
		// Find max values
		algmat result(msize(1, mat.cols));
		for (int i = 0; i < mat.cols; i++)
		{
			int arg = 0;
			double val = mat[0][i];
			for (int j = 1; j < mat.rows; j++)
			{
				if (val < mat[j][i])
				{
					arg = j;
					val = mat[j][i];
				}
			}
			loc[0][i] = arg;
			result[0][i] = val;
		}

		return result;
	}

	const algmat algmat::hmax(const algmat& mat, algmat& loc)
	{
		// Find max values
		algmat result(msize(mat.rows, 1));
		for (int i = 0; i < mat.rows; i++)
		{
			int arg = 0;
			double val = mat[i][0];
			for (int j = 1; j < mat.cols; j++)
			{
				if (val < mat[i][j])
				{
					arg = j;
					val = mat[i][j];
				}
			}
			loc[i][0] = arg;
			result[i][0] = val;
		}

		return result;
	}

	algmat algmat::argmax(const algmat& mat, const int axis)
	{
		// Check matrix status
		assert(mat.empty() == false);

		// Check an axis information
		if (axis == 0)
		{
			// Find max arguments
			algmat loc(msize(1, mat.cols));
			vmax(mat, loc);

			return loc;
		}
		else if (axis == 1)
		{
			// Find max arguments
			algmat loc(msize(mat.rows, 1));
			hmax(mat, loc);

			return loc;
		}
		else
		{
			// Find max arguments
			algmat loc0(msize(1, mat.cols));
			algmat loc1(msize(1, 1));
			hmax(vmax(mat, loc0), loc1);

			return loc1;
		}
	}

	algmat algmat::sum(const algmat& mat, const int axis)
	{
		// Check matrix status
		assert(mat.empty() == false);

		// Check an axis information
		if (axis == 0)
		{
			// Calculate sum values
			return vsum(mat);
		}
		else if (axis == 1)
		{
			// Calculate sum values
			return hsum(mat);
		}
		else
		{
			// Calculate sum values
			return hsum(vsum(mat));
		}
	}

	const algmat algmat::vsum(const algmat& mat)
	{
		// Calculate sum values
		algmat result(msize(1, mat.cols), 0.0);
		for (int i = 0; i < mat.cols; i++)
		{
			for (int j = 0; j < mat.rows; j++)
			{
				result[0][i] += mat[j][i];
			}
		}

		return result;
	}

	const algmat algmat::hsum(const algmat& mat)
	{
		// Calculate sum values
		algmat result(msize(mat.rows, 1), 0.0);
		for (int i = 0; i < mat.rows; i++)
		{
			for (int j = 0; j < mat.cols; j++)
			{
				result[i][0] += mat[i][j];
			}
		}

		return result;
	}

	algmat algmat::mean(const algmat& mat, const int axis)
	{
		// Check matrix status
		assert(mat.empty() == false);

		// Check an axis information
		if (axis == 0)
		{
			// Calculate mean values
			return vsum(mat) / (double)mat.rows;
		}
		else if (axis == 1)
		{
			// Calculate mean values
			return hsum(mat) / (double)mat.cols;
		}
		else
		{
			// Calculate mean values
			return hsum(vsum(mat)) / (double)mat.dm.length();
		}
	}

	algmat algmat::var(const algmat& mat, const int axis)
	{
		// Check matrix status
		assert(mat.empty() == false);

		// Calculate var values
		algmat mmat = mean(mat, axis);
		algmat dmat = mmat - mat;
		algmat vmat = mean(dmat * dmat, axis);

		return vmat;
	}

	algmat algmat::std(const algmat& mat, const int axis)
	{
		// Check matrix status
		assert(mat.empty() == false);

		// Calculate std values
		algmat vmat = var(mat, axis);
		algmat smat = sqrt(vmat);

		return smat;
	}

	const bool algmat::almost(const double val, const double target, const int round)
	{
		// Compare the values
		double power = std::pow(10.0, round);
		long long _val = (long long)(val * power + 0.5);
		long long _target = (long long)(target * power + 0.5);
		if (_val == _target)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	algmat algmat::diag() const
	{
		// Check matrix status
		assert(empty() == false);
		assert(isSquared() == true);

		// Get a diagonal matrix
		algmat result(msize(rows, cols), 0.0);
		for (int i = 0; i < rows; i++)
		{
			result[i][i] = (*this)[i][i];
		}

		return result;
	}

	const bool algmat::isSquared() const
	{
		// Check a size
		if (rows == cols)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	algmat algmat::values(const msize& sz, const double val)
	{
		// Get a values matrix
		algmat result(sz, val);
		for (int i = 0; i < result.tlen; i++)
		{
			result.ddata[i] = val;
		}

		return result;
	}

	algmat algmat::zeros(const msize& sz)
	{
		// Get a zeros matrix
		return values(sz, 0.0);
	}

	algmat algmat::ones(const msize& sz)
	{
		// Get an ones matrix
		return values(sz, 1.0);
	}

	algmat algmat::eyes(const int len)
	{
		// Get an eyes matrix
		algmat result(msize(len, len), 0.0);
		for (int i = 0; i < len; i++)
		{
			result[i][i] = 1.0;
		}

		return result;
	}

	algmat algmat::permut(const int len, const int i0, const int i1)
	{
		// Check the indices
		assert(i0 >= 0 && i0 < len);
		assert(i1 >= 0 && i1 < len);

		// Get a permutation matrix
		algmat result = eyes(len);
		result[i0][i0] = 0.0;
		result[i1][i1] = 0.0;
		result[i0][i1] = 1.0;
		result[i1][i0] = 1.0;

		return result;
	}

	algmat algmat::uniq() const
	{
		// Check matrix status
		assert(empty() == false);

		// Fine unique values on the matrix
		vector<double> uv;
		for (int i = 0; i < tlen; i++)
		{
			// Check a redundant value or not
			bool redundant = false;
			for (int j = 0; j < (int)uv.size(); j++)
			{
				if (ddata[i] == uv[j])
				{
					redundant = true;
					break;
				}
			}
			if (redundant == false)
			{
				uv.push_back(ddata[i]);
			}
		}

		// Get an unique matrix
		return to_algmat(uv);
	}

	algmat algmat::dot(const algmat& mat) const
	{
		// Check matrix status
		assert(empty() == false);
		assert(mat.empty() == false);
		assert(cols == mat.rows);

		// Get a dot product matrix
		algmat result(msize(rows, mat.cols));
		for (int i = 0; i < result.rows; i++)
		{
			for (int j = 0; j < result.cols; j++)
			{
				double val = 0.0;
				for (int k = 0; k < cols; k++)
				{
					val += (*this)[i][k] * mat[k][j];
				}
				result[i][j] = val;
			}
		}

		return result;
	}

	algmat algmat::cof() const
	{
		// Check matrix status
		assert(empty() == false);
		assert(isSquared() == true);

		// Get a cofactor expansion matrix
		algmat result(msize(rows, cols));
		for (int i = 0; i < result.rows; i++)
		{
			for (int j = 0; j < result.cols; j++)
			{
				result[i][j] = std::pow(-1.0, i + j) * expansion(minor(i, 0).minor(j, 1));
			}
		}

		return result;
	}

	const algmat algmat::expansion(const algmat& mat) const
	{
		// Check the size
		algmat result(msize(1), 0.0);
		if (mat.rows == 2 && mat.cols == 2)
		{
			result[0][0] = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
		}
		else
		{
			// Check zero data
			int val = 0, idx = 0;
			for (int i = 0; i < mat.rows; i++)
			{
				int zc = 0;
				for (int j = 0; j < mat.cols; j++)
				{
					if (almost(mat[i][j], 0.0) == true)
					{
						zc++;
					}
				}
				if (val < zc)
				{
					val = zc;
					idx = i;
				}
			}

			// Calculate a determinant
			for (int i = 0; i < mat.cols; i++)
			{
				if (almost(mat[idx][i], 0.0) == false)
				{
					result += std::pow(-1.0, idx + i) * mat[idx][i] * expansion(mat.minor(idx, 0).minor(i, 1));
				}
			}
		}

		return result;
	}

	algmat algmat::inv() const
	{
		// Check matrix status
		assert(empty() == false);
		assert(isSquared() == true);

		// Get an inverse matrix
		algmat result(msize(rows, cols));
		if (result.rows == 2 && result.cols == 2)
		{
			// Check a determinant
			double val = (*this)[0][0] * (*this)[1][1] - (*this)[0][1] * (*this)[1][0];
			assert(almost(val, 0.0) == false);

			// Calculate an inverse matrix
			result[0][0] = (*this)[1][1] / val;
			result[1][1] = (*this)[0][0] / val;
			result[1][0] = -(*this)[1][0] / val;
			result[0][1] = -(*this)[0][1] / val;
		}
		else
		{
			// Get a cofactor matrix
			algmat c = cof();

			// Check a determinant
			double val = 0.0;
			for (int i = 0; i < cols; i++)
			{
				val += c[0][i] * (*this)[0][i];
			}
			assert(almost(val, 0.0) == false);

			// Calculate an inverse matrix
			result = c.t() * (1.0 / val);
		}

		return result;
	}

	algmat algmat::det() const
	{
		// Check matrix status
		assert(empty() == false);
		assert(isSquared() == true);

		// Get an inverse matrix
		algmat result(msize(1), 0.0);
		if (rows == 2 && cols == 2)
		{
			result = expansion(*this);
		}
		else
		{
			// Get a cofactor matrix
			algmat c(msize(1, cols));
			for (int i = 0; i < cols; i++)
			{
				c[0][i] = std::pow(-1.0, i) * expansion(minor(0, 0).minor(i, 1));
			}

			// Calculate a determinant
			for (int i = 0; i < cols; i++)
			{
				result[0][0] += c[0][i] * (*this)[0][i];
			}
		}

		return result;
	}
}