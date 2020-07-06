#include "stdafx.h"
#include "Decomposition.h"

namespace nml
{
	const algmat& decomposition::operator[](const int idx) const
	{
		// Check the index
		assert(idx >= 0 && idx < 2);

		return D[idx];
	}

	void decomposition::checkEmpty(const algmat& mat) const
	{
		// Check the matrix is empty or not
		assert(mat.empty() == false);
	}

	void decomposition::checkSquareness(const algmat& mat) const
	{
		// Check the matrix is square or not
		assert(mat.rows == mat.cols);
	}

	const vector<algmat> decomposition::result() const
	{
		// Get the decomposition result
		return D;
	}

	LUD::LUD()
	{
		// Set an object
		setObject();
	}

	LUD::LUD(const algmat& mat)
	{
		// Set an object
		setObject();

		// Set a partial pivot as a default mode
		pivot(LU_PARTIAL_PIVOT);

		// Decompose the matrix
		decomp(mat);
	}

	LUD::LUD(const algmat& mat, const pivotMode mode)
	{
		// Set an object
		setObject();

		// Set a pivot mode
		pivot(mode);

		// Decompose the matrix
		decomp(mat);
	}

	LUD::LUD(const LUD& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	LUD::~LUD()
	{
		// Clear the object
		clearObject();
	}

	LUD& LUD::operator=(const LUD& obj)
	{
		// Copy the object
		copyObject(obj);

		return *this;
	}

	inline void LUD::setObject()
	{
		// Set the parameters
		setType(*this);
		mode = LU_PARTIAL_PIVOT;

		// Set the memories
		P.clear();
		M.clear();
		D.clear();
	}

	inline void LUD::copyObject(const object& obj)
	{
		// Do down casting
		const LUD* _obj = static_cast<const LUD*>(&obj);

		// Copy the parameters
		mode = _obj->mode;

		// Copy the memories
		P = _obj->P;
		M = _obj->M;
		D = _obj->D;
	}

	inline void LUD::clearObject()
	{
		// Clear the memories
		P.clear();
		M.clear();
		D.clear();
	}

	void LUD::pivot(const pivotMode mode)
	{
		// Set a pivot mode
		this->mode = mode;
	}

	const vector<algmat> LUD::decomp(const algmat& mat)
	{
		// Check the matrix
		checkEmpty(mat);

		// Check the matrix
		checkSquareness(mat);

		// Clear the old memories
		P.clear();
		M.clear();
		D.clear();

		// Copy the input matrix
		algmat A = mat;

		// Factorize the matrix
		factorize(mat);

		// Create a triangle matrix
		vector<algmat> triangle;
		for (int i = mat.rows - 1; i >= 0; i--)
		{
			triangle.push_back(M[i].dot(P[i]));
		}

		// Create a decomposition matrix vector
		algmat temp = triangle[0];
		for (int i = 1; i < (int)triangle.size(); i++)
		{
			temp = temp.dot(triangle[i]);
		}
		D.push_back(temp.inv());
		D.push_back(temp.dot(A));

		return D;
	}

	void LUD::factorize(const algmat& mat)
	{
		// Factorize the matrix with pivoting
		algmat temp = mat;
		for (int i = 0; i < temp.cols; i++)
		{
			// Find the row and column maximum index
			int rowMaxIdx = i;
			int colMaxIdx = i;
			double maxVal = temp(rowMaxIdx, colMaxIdx);
			for (int j = i; j < temp.rows; j++)
			{
				// Check the pivot mode
				if (mode == LU_PARTIAL_PIVOT)
				{
					if (temp(j, i) > maxVal)
					{
						maxVal = temp(j, i);
						rowMaxIdx = j;
						colMaxIdx = i;
					}
				}
				else
				{
					for (int k = i; k < temp.cols; k++)
					{
						if (temp(j, k) > maxVal)
						{
							maxVal = temp(j, k);
							rowMaxIdx = j;
							colMaxIdx = k;
						}
					}
				}
			}

			// Swap the matrix
			temp.swap(i, rowMaxIdx, 0);
			temp.swap(i, colMaxIdx, 1);

			// Save the result
			P.push_back(algmat::permut(temp.rows, i, rowMaxIdx));
			M.push_back(algmat::eyes(temp.rows));
			for (int j = i + 1; j < temp.rows; j++)
			{
				M[M.size() - 1](j, i) = -temp(j, i) / temp(i, i);
				for (int k = 0; k < temp.cols; k++)
				{
					temp(j, k) += temp(i, k) * M[M.size() - 1](j, i);
				}
			}
		}
	}

	QRD::QRD()
	{
		// Set an object
		setObject();
	}

	QRD::QRD(const algmat& mat)
	{
		// Set an object
		setObject();

		// Decompose the matrix
		decomp(mat);
	}

	QRD::QRD(const QRD& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	QRD::~QRD()
	{
		// Clear the object
		clearObject();
	}

	QRD& QRD::operator=(const QRD& obj)
	{
		// Copy the object
		copyObject(obj);

		return *this;
	}

	inline void QRD::setObject()
	{
		// Set the parameters
		setType(*this);

		// Set the memories
		T.clear();
		D.clear();
	}

	inline void QRD::copyObject(const object& obj)
	{
		// Do down casting
		const QRD* _obj = static_cast<const QRD*>(&obj);

		// Copy the memories
		T = _obj->T;
		D = _obj->D;
	}

	inline void QRD::clearObject()
	{
		// Clear the memories
		T.clear();
		D.clear();
	}

	const vector<algmat> QRD::decomp(const algmat& mat)
	{
		// Check Matrix
		checkEmpty(mat);

		// Clear the old memories
		T.clear();
		D.clear();

		// Check the row information of the matrix
		if (mat.rows == 1)
		{
			D.push_back(algmat::ones(msize(1)));
			D.push_back(mat);

			return D;
		}

		// Copy the input matrix
		algmat temp = mat;
		algmat A = mat;

		// Do House Holder Transform on the matrix
		for (int i = 0; i < A.rows - 1; i++)
		{
			int j = (i != 0) ? 1 : 0;
			if (j - temp.cols == 0)
			{
				break;
			}
			algmat S = reduction(temp, j);
			algmat Q = houseHolder(S);
			temp = Q.dot(S);
			T.push_back(extension(Q, i, A.rows));
		}

		// Create a decomposition matrix vector
		const int count = (int)T.size();
		D.push_back(T[0].t());
		for (int i = 1; i < count; i++)
		{
			D[0] = D[0].dot(T[i].t());
		}
		D.push_back(T[count - 1]);
		for (int i = count - 2; i >= 0; i--)
		{
			D[1] = D[1].dot(T[i]);
		}
		D[1] = D[1].dot(A);

		return D;
	}

	const algmat QRD::reduction(const algmat& mat, int flag) const
	{
		// Check the reduction flag
		if (flag == 0)
		{
			// Get the input matrix
			return mat;
		}
		else
		{
			// Get a reduction matrix
			return mat.minor(0, 0).minor(0, 1);
		}
	}

	const algmat QRD::houseHolder(const algmat& mat) const
	{
		// Do House Holder Transform on the matrix
		algmat a = mat.submat(0, 1);
		algmat e = reflection(a, mat.rows);
		algmat u = a - e;
		algmat v = u / 2.0;
		algmat vvt = v.dot(v.t());

		return algmat::eyes(mat.rows) - (vvt * (2.0 / vecsize(vvt)));
	}

	const algmat QRD::reflection(const algmat& mat, int length) const
	{
		// Create a reflection matrix
		algmat result(msize(length, 1), 0.0);
		result(0) = vecsize(mat);

		return result;
	}

	const double QRD::vecsize(const algmat& mat) const
	{
		// Calculate a size of the input vector
		return std::sqrt(algmat::sum(mat.mul(mat))(0));
	}

	const algmat QRD::extension(const algmat& mat, const int offset, const int length) const
	{
		// Get an extension matrix
		algmat result(msize(length, length), 0.0);
		for (int i = 0; i < offset; i++)
		{
			result(i, i) = 1.0;
		}
		for (int i = 0; i < mat.rows; i++)
		{
			for (int j = 0; j < mat.cols; j++)
			{
				result(offset + i, offset + j) = mat(i, j);
			}
		}

		return result;
	}

	EVD::EVD()
	{
		// Set an object
		setObject();
	}

	EVD::EVD(const algmat& mat)
	{
		// Set an object
		setObject();

		// Decompose the matrix
		decomp(mat);
	}

	EVD::EVD(const EVD& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	EVD::~EVD()
	{
		// Clear the object
		clearObject();
	}

	EVD& EVD::operator=(const EVD& obj)
	{
		// Copy the object
		copyObject(obj);

		return *this;
	}

	inline void EVD::setObject()
	{
		// Set the parameters
		setType(*this);

		// Set the memories
		D.clear();
	}

	inline void EVD::copyObject(const object& obj)
	{
		// Do down casting
		const EVD* _obj = static_cast<const EVD*>(&obj);

		// Copy the memories
		D = _obj->D;
	}

	inline void EVD::clearObject()
	{
		// Clear the memories
		D.clear();
	}

	const vector<algmat> EVD::decomp(const algmat& mat)
	{
		// Check the matrix
		checkEmpty(mat);

		// Check the matrix
		checkSquareness(mat);

		// Clear the old memories
		D.clear();

		// Initialize the parameters
		const int maxIter = 100 * max(mat.rows, mat.cols);
		int loopIter = 0;
		bool convFlag = false;

		// Do eigen value decomposition using QRD
		QRD qr(mat);
		algmat eigVec = qr[0];
		algmat eigVal = qr[1].dot(qr[0]);
		algmat tempVal = eigVal;
		while (true)
		{
			// Do QR decomposition
			qr.decomp(eigVal);

			// Calculate an eigen value and an eigen vector
			eigVal = qr[1].dot(qr[0]);
			convFlag = isConvergence(eigVal, tempVal);
			eigVec = eigVec.dot(qr[0]);

			// Check the convergence flag
			if (convFlag == true)
			{
				break;
			}
			else
			{
				// Check the stop condition
				if (loopIter == maxIter)
				{
					break;
				}
				else
				{
					tempVal = eigVal;
					loopIter++;
				}
			}
		}

		// Create a decomposition matrix vector
		D.push_back(eigVec);
		D.push_back(eigVal);

		return D;
	}

	const bool EVD::isConvergence(const algmat& curr, const algmat& prev) const
	{
		// Calculate a delta between matrices
		algmat delta = algmat::sqrt(curr.mul(curr) - prev.mul(prev));

		// Check the matrix is convergence or not
		for (int i = 0; i < delta.length(); i++)
		{
			if (algmat::almost(delta(i), 0.0, 5) == false)
			{
				return false;
			}
		}

		return true;
	}

	SVD::SVD()
	{
		// Set an object
		setObject();
	}

	SVD::SVD(const algmat& mat)
	{
		// Set an object
		setObject();

		// Decompose the matrix
		decomp(mat);
	}

	SVD::SVD(const SVD& obj)
	{
		// Set an object
		setObject();

		// Clone the object
		*this = obj;
	}

	SVD::~SVD()
	{
		// Clear the object
		clearObject();
	}

	SVD& SVD::operator=(const SVD& obj)
	{
		// Copy the object
		copyObject(obj);

		return *this;
	}

	const algmat& SVD::operator[](const int idx) const
	{
		// Check the index
		assert(idx >= 0 && idx < 3);

		return D[idx];
	}

	inline void SVD::setObject()
	{
		// Set the parameters
		setType(*this);
	}

	const vector<algmat> SVD::decomp(const algmat& mat)
	{
		// Check the matrix
		checkEmpty(mat);

		// Clear the old memories
		D.clear();

		// Initialize the parameters
		const int maxIter = 100 * max(mat.rows, mat.cols);
		int loopIter = 0;

		// Do eigen value decomposition using QRD
		QRD qru(mat);
		QRD qrv(qru[1].t());
		algmat T = qrv[1];
		algmat U = qru[0];
		algmat S = qrv[1];
		algmat V = qrv[0];
		while (true)
		{
			// Do QR decomposition
			qru.decomp(S.t());
			qrv.decomp(qru[1].t());
			U = U.dot(qru[0]);
			V = V.dot(qrv[0]);
			S = qrv[1];

			// Check the convergence flag
			if (isConvergence(S, T) == true)
			{
				break;
			}
			else
			{
				// Check the stop condition
				if (loopIter == maxIter)
				{
					break;
				}
				else
				{
					T = S;
					loopIter++;
				}
			}
		}

		// Create a decomposition matrix vector
		S = S.t();
		int k = (mat.cols > mat.rows) ? mat.rows : mat.cols;
		int l = mat.rows;
		for (int i = 0; i < k; i++)
		{
			if (S(i, i) < 0.0)
			{
				S(i, i) = -S(i, i);
				for (int j = 0; j < l; j++)
				{
					U(j, i) = -U(j, i);
				}
			}
		}
		D.push_back(U);
		D.push_back(S);
		D.push_back(V);

		return D;
	}

	algmat SVD::pinv(const algmat& mat)
	{
		// Factorize the matrix
		SVD sv(mat);

		// Inverse the elements on the eigen matrix
		algmat E = sv[1];
		for (int i = 0; i < E.rows; i++)
		{
			for (int j = 0; j < E.cols; j++)
			{
				if (algmat::almost(E(i, j), 0.0, 5) == false)
				{
					E(i, j) = 1.0 / E(i, j);
				}
			}
		}

		return sv[2].dot(E.t()).dot(sv[0].t());
	}
}