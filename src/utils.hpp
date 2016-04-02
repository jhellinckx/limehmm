#ifndef __UTILS_HPP
#define __UTILS_HPP

#include <string>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include <limits>
#include <typeinfo>
#include <type_traits>

namespace utils {
	extern const double kInf =  std::numeric_limits<double>::infinity();
	extern const double kNegInf = -std::numeric_limits<double>::infinity();
	extern const int kPrecision = 8;

	double round_double(const double& d){
		return round(d * pow(10, kPrecision)) / pow(10, kPrecision);
	}

	extern double sum_log_prob(double log_x, double log_y){
		// prob(x) == inf, prob(y) == inf
		if(log_x == kInf || log_y == kInf) return kInf;
		// prob(x) == 0
		if(log_x == kNegInf) return log_y;
		if(log_y == kNegInf) return log_x;
		return (log_x > log_y) ? log_x + log(1 + exp(log_y - log_x)) : log_y + log(1 + exp(log_x - log_y));
	}

	template<typename Iterator>
	extern double sum_log_prob(Iterator begin, Iterator end, double init_sum = kNegInf){
		Iterator it = begin;
		double log_sum = init_sum;
		while(it != end){
			log_sum = sum_log_prob(log_sum, *it);
			++it;
		}
		return log_sum;
	}

	extern double log_normalize(double log_x, double log_sum){
		if(log_x == kInf || log_sum == kInf) return kInf;
		if(log_x == kNegInf || log_sum == kNegInf) return kNegInf;
		return log_x - log_sum;
	}

	template<typename Iterator>
	extern void for_each_log_normalize(Iterator begin, Iterator end, double log_sum){
		Iterator it = begin;
		while(it != end){
			*it = log_normalize(*it, log_sum);
			++it;
		}
	}

	template <typename T>
	class to_string {
		template <typename C> static char test(char[sizeof(&C::to_string)]);
		template <typename C> static long test(...);    
	public:
		static bool exists() {
			return sizeof(test<T>(0)) == sizeof(char);	
		}
	};
}

#endif