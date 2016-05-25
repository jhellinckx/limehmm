#ifndef __UTILS_HPP
#define __UTILS_HPP

#include <string>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include <limits>
#include <typeinfo>
#include <type_traits>
#include <mach/mach.h>
#include "constants.hpp"

namespace utils {
	extern const double kInf;
	extern const double kNegInf;

	double round_double(double, int = global_config::kDoublePrecision);

	double sum_log_prob(double log_x, double log_y);

	double log_normalize(double log_x, double log_sum);

	std::pair<std::string, std::string> split_first(const std::string& s, char c);

	void mem_info();
}

#endif