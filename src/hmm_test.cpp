#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <vector>
#include <string>
#include <stdlib.h>

#include "hmm.hpp"

#define ASSERT(expr) assert(expr, #expr, __FILE__, __LINE__)
#define ASSERT_VERBOSE(expr, msg) assert(expr, #expr, __FILE__, __LINE__, msg)
#define ASSERT_ABORT(expr, msg) assert(expr, #expr, __FILE__, __LINE__, msg, true)
#define ASSERT_EXCEPT(instruction, except_type) assert_except<except_type>([&](){instruction;}, #instruction, #except_type, __FILE__, __LINE__)

#define TEST_UNIT(name, instructions) run_unit_test(name, [](){instructions});

#define VERBOSE 1
#define BIG_SEPARATOR "======================================================"
#define THIN_SEPARATOR "------------------------------------------------------"
#define RED "\033[31m"
#define GREEN "\033[32m"
#define MAGENTA "\033[35m"
#define RESET "\033[0m"
#define BOLDRED "\033[1;31m"
#define BOLDGREEN "\033[1;32m"
#define BOLDCYAN "\033[1;36m"
#define BOLDMAGENTA "\033[1;35m"

void print_exception(const std::exception& e){ std::cerr << RED << e.what() << RESET << std::endl; }

int assertions = 0;
int units = 0;
int failed = 0;
int successful = 0;

template<typename Functor>
void run_unit_test(const std::string& name, Functor functor){
	++units;
	if(VERBOSE) std::cout << THIN_SEPARATOR << std::endl << "Testing " << MAGENTA << name << RESET << "..."  << std::endl;
	functor();
}

void assert(bool assertion, const char* assertion_c_str, const char* filename, long int line, std::string fail_message = "", bool fail_abort = false){
	++assertions;
	if(VERBOSE) std::cout << assertion_c_str << " ? ";
	if(assertion == true){
		++successful;
		if(VERBOSE) std::cout << BOLDGREEN << "OK" << RESET << std::endl;
	}
	else{
		++failed;
		if(VERBOSE) { 
			fail_message = (fail_message.length() > 0) ? ": " + fail_message : "";
			std::cout << BOLDRED << "FAIL" << RESET
			<< " -> " << filename  << ": " << + "line " << line 
			<< fail_message << std::endl;
		}
		if(fail_abort) abort();
	}
}

template<typename ExceptionType, typename Functor>
void assert_except(Functor functor, const char* instruction_c_str, const char* exception_c_str, const char* filename, long int line){
	try{
		functor();
		assert(false, instruction_c_str, filename, line, std::string(exception_c_str) + " expected but no exception was thrown");
	}
	catch(const ExceptionType& e){
		assert(true, instruction_c_str, filename, line);
	}
	catch(const std::exception& e){
		assert(false, instruction_c_str, filename, line, std::string(exception_c_str) + " expected but another exception was thrown: " + e.what());
	}
}

void tests_init(){
	if(VERBOSE){
		std::cout << BIG_SEPARATOR << std::endl;
		std::cout << BOLDMAGENTA << "Running " << "tests" << "..." << RESET << std::endl;
	}
}

void tests_results(){
	if(VERBOSE){
		std::cout << BIG_SEPARATOR << std::endl;
		std::cout << BOLDMAGENTA << "Ran " << assertions << " assertion(s) for "<< units << " test(s) : " << RESET;
		std::cout << BOLDGREEN << successful << " succeeded " << RESET;
		std::cout << BOLDRED << failed << " failed." << RESET << std::endl;
	}
	if(failed > 0) throw std::runtime_error("Tests failed.");
}

int main(){
	try{

		tests_init();
		
		/* IEEE 754 floating points are required in order to use std::infinity numeric limit. */
		TEST_UNIT(
			"platform type",
			ASSERT_ABORT(std::numeric_limits<double>::is_iec559, "IEEE 754 required");
		)

		TEST_UNIT(
			"state creation/distribution",
			State s1("state");
			State s2("state");
			ASSERT(s1 == s2);
			ASSERT(s1.is_silent());
			ASSERT_EXCEPT(s1.distribution(), StateDistributionException);
			State s3("state", DiscreteDistribution({{"yes", 0.9}, {"no", 0.1}}));
			ASSERT(s3.distribution().is_discrete());
			ASSERT(!s3.distribution().is_continuous());
			ASSERT(s3 == s1);
			State s4(s3);
			ASSERT(s4.distribution() == s3.distribution());
			DiscreteDistribution dist1({{"A", 0.2}, {"G", 0.4}, {"C", 0.1}, {"T", 0.3}});
			DiscreteDistribution dist2 = dist1;
			ASSERT(dist2["A"] == 0.2);
			ASSERT(dist1 == dist2);
			ASSERT(dist1["A"] == 0.2);
			ASSERT((dist1["A"] = 0.5) == 0.5);
			ASSERT((dist1 != dist2));
		)
		
		TEST_UNIT(
			"begin/end state",
			State begin("begin");
			State end("end");
			HiddenMarkovModel hmm = HiddenMarkovModel(begin, end);
			ASSERT(hmm.begin() == begin);
			ASSERT(hmm.end() == end);
			ASSERT(hmm.contains(begin));
			ASSERT(hmm.contains(end));
		)

		TEST_UNIT(
			"add/remove state",
			HiddenMarkovModel hmm = HiddenMarkovModel();
			State s("s");
			ASSERT(!hmm.contains("s"));
			ASSERT(!hmm.contains(s));
			hmm.add_state(s);
			ASSERT(hmm.contains("s"));
			ASSERT(hmm.contains(s));
			ASSERT_EXCEPT(hmm.add_state(s), VertexExistsException);
			hmm.remove_state("s");
			ASSERT(!hmm.contains("s"));
			ASSERT(!hmm.contains(s));
		)

		TEST_UNIT(
			"add/remove transition",
			HiddenMarkovModel hmm = HiddenMarkovModel();
			State s1("s1");
			State s2("s2");

			hmm.add_state(s1);
		)

		/* Emission */

		/* Initial probability aka pi */

		/* Backward */

		/* Forward */

		/* Log Probability */

		/* Decode Viterbi */

		/* Sample */

		/* Train Viterbi */

		/* Train B-W */

		/* Train stochastic EM */


		tests_results();

	} catch(const std::exception& e){
		print_exception(e);
		return 1;
	}
	return 0;
}
