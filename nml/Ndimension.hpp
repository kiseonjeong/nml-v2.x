namespace nml
{
	dim::dim() : d(_d), N(_N), len(_len)
	{
		// Set default values
		_N = 0;
		_len = 0;
		init = true;
	}

	template <typename T, typename... Types> dim::dim(T arg, Types... args) : d(_d), N(_N), len(_len)
	{
		// Set default values
		_N = 0;
		_len = 0;
		init = true;

		// Set dimension information
		set(arg, args...);
	}

	dim::dim(const dim& obj) : d(_d), N(_N), len(_len)
	{
		// Set default values
		_N = 0;
		_len = 0;
		init = true;

		// Copy dimension information
		_d = obj._d;
		_N = obj._N;
		_len = obj._len;
	}

	dim::~dim()
	{
		// Clear dimension information
		clear();
	}

	dim& dim::operator=(const dim& obj)
	{
		// Copy from the input object
		_d = obj._d;
		_N = obj._N;
		_len = obj._len;

		return *this;
	}

	const int& dim::operator[](const int idx) const
	{
		// Get a dimension length
		return _d[idx];
	}

	template <typename T> void dim::set(T arg)
	{
		// Initialize and set dimension information
		initialize();
		_d.push_back(arg);
		_N += 1;
		_len += arg;

		// Reset the flag
		init = true;
	}

	template <typename T, typename... Types> void dim::set(T arg, Types... args)
	{
		// Initialize and set dimension information
		initialize();
		_d.push_back(arg);
		_N += 1;
		_len += arg;
		set(args...);
	}

	void dim::initialize()
	{
		// Check the flag
		if (init == true)
		{
			_d.clear();
			_N = 0;
			_len = 0;
			init = false;
		}
	}

	void dim::clear()
	{
		// Clear dimension information
		_d.clear();
		_N = 0;
		_len = 0;
		init = true;
	}
}
