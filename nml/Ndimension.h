#pragma once

namespace nml
{
	// The dimension for array
	class dim
	{
		// Variables
	public:
		const std::vector<int>& d;			// dimension information
		const int& N;			// number of dimensions
		const int& len;			// length of dimensions

		// Functions
	public:
		template <typename T> void set(T arg);
		template <typename T, typename... Types> void set(T arg, Types... args);
		void clear();			// clear dimension information

		// Operators
	public:
		dim& operator=(const dim& obj);
		const int& operator[](const int idx) const;

		// Constructor & destructor
	public:
		dim();
		template <typename T, typename... Types> dim(T arg, Types... args);
		dim(const dim& obj);
		virtual ~dim();

		// Variables
	private:
		bool init;			// initial status flag
		std::vector<int> _d;			// dimension information
		int _N;			// number of dimensions
		int _len;			// length of dimensions

		// Functions
	private:
		void initialize();			// initialize dimension information

	};
}

#include "Ndimension.hpp"