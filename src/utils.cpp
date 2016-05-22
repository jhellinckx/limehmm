#include <string>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include <limits>
#include <typeinfo>
#include <type_traits>
#include "utils.hpp"


namespace utils {
	const double kInf =  std::numeric_limits<double>::infinity();
	const double kNegInf = -std::numeric_limits<double>::infinity();

	double round_double(double d, int precision){
		return round(d * pow(10, precision)) / pow(10, precision);
	}

	double sum_log_prob(double log_x, double log_y){
		// prob(x) == inf, prob(y) == inf
		if(log_x == kInf || log_y == kInf) return kInf;
		// prob(x) == 0
		if(log_x == kNegInf) return log_y;
		if(log_y == kNegInf) return log_x;
		return (log_x > log_y) ? log_x + log(1 + exp(log_y - log_x)) : log_y + log(1 + exp(log_x - log_y));
	}

	double log_normalize(double log_x, double log_sum){
		if(log_x == kInf || log_sum == kInf) return kInf;
		if(log_x == kNegInf || log_sum == kNegInf) return kNegInf;
		return log_x - log_sum;
	}

	std::pair<std::string, std::string> split_first(const std::string& s, char c){
		std::size_t split_i = 0;
		bool found = false;
		while(!found && split_i < s.length()){
			if(s[split_i] == c){ found = true; }
			else{ ++split_i; }
		}
		if(split_i == s.length()) { return std::make_pair(s, ""); }
		return std::make_pair(s.substr(0, split_i), s.substr(split_i + 1, std::string::npos));
	}
}
