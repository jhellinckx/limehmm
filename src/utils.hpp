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

	extern double sum_log_prob(double log_x, double log_y){
		// prob(x) == inf, prob(y) == inf
		if(log_x == kInf or log_y == kInf) return kInf;
		// prob(x) == 0
		if(log_x == kNegInf) return log_y;
		if(log_y == kNegInf) return log_x;
		// log(x) > log(y) -> log(x + y) = log(x) + log(1 + exp(log(y) - log(x)))
		return (log_x > log_y) ? log_x + log(1 + exp(log_y - log_x)) : log_y + log(1 + exp(log_x - log_y));
	}

	template<typename CollectionType>
	extern double sum_log_prob(typename CollectionType::iterator begin, typename CollectionType::iterator end, double init_sum = kNegInf){
		typename CollectionType::iterator it = begin;
		double log_sum = init_sum;
		while(it != end){
			log_sum = sum_log_prob(log_sum, *it);
			++it;
		}
		return log_sum;
	}

	template<typename CollectionType>
	extern void log_normalize(typename CollectionType::iterator begin, typename CollectionType::iterator end){
		typename CollectionType::iterator it = begin;
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