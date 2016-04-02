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
#include <math.h>
#include "utils.hpp"
#include "hmm.hpp" // tested hmm library

#define ASSERT(expr) assert(expr, #expr, __FILE__, __LINE__)
#define ASSERT_VERBOSE(expr, msg) assert(expr, #expr, __FILE__, __LINE__, msg)
#define ASSERT_ABORT(expr, msg) assert(expr, #expr, __FILE__, __LINE__, msg, true)
#define ASSERT_EXCEPT(instruction, except_type) assert_except<except_type>([&](){instruction;}, #instruction, #except_type, __FILE__, __LINE__)

#define TEST_UNIT(name, instructions) run_unit_test(name, [](){instructions});

#define VERBOSE 1
#define BIG_SEPARATOR std::string(60, '=')
#define THIN_SEPARATOR std::string(60, '-')
#define RED "\033[31m"
#define GREEN "\033[32m"
#define MAGENTA "\033[35m"
#define CYAN "\033[36m"
#define BOLDRED "\033[1;31m"
#define BOLDGREEN "\033[1;32m"
#define BOLDMAGENTA "\033[1;35m"
#define BOLDCYAN "\033[1;36m"
#define RESET "\033[0m"


void print_exception(const std::exception& e){ std::cerr << RED << e.what() << RESET << std::endl; }

int assertions = 0;
int units = 0;
int failed = 0;
int successful = 0;

template<typename Runnable>
void run_unit_test(const std::string& name, Runnable runnable){
	++units;
	if(VERBOSE) std::cout << THIN_SEPARATOR << std::endl << "Testing " << MAGENTA << name << RESET << "..."  << std::endl;
	runnable();
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

template<typename ExceptionType, typename Runnable>
void assert_except(Runnable runnable, const char* instruction_c_str, const char* exception_c_str, const char* filename, long int line){
	try{
		runnable();
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
			/* Construct silent states with same name. */
			State s1("state");
			State s2("state");
			ASSERT(s1 == s2);
			/* If no distribution is given at construction, state is silent. */
			ASSERT(s1.is_silent());
			/* Throw exception if access distribution of silent state. */
			ASSERT_EXCEPT(s1.distribution(), StateDistributionException);
			/* An empty distribution makes the state silent. */
			DiscreteDistribution dist1;
			State s3("state", dist1);
			ASSERT(s3.is_silent());
			/* State has own copy of distribution. */
			ASSERT(s3.distribution() == dist1);
			dist1["A"] = 0;
			ASSERT(s3.distribution() != dist1);
			/* A discrete distribution with probabilites summing to 0 is considered empty. */
			dist1["B"] = 0;
			State s4("state", dist1);
			ASSERT(s4.is_silent());
			/* Probabilites summing to i > 0 makes a state not silent. */
			dist1["C"] = 0.4;
			State s5("state", dist1);
			ASSERT(!s5.is_silent());
			/* Check distribution type. */
			ASSERT(s3.distribution().is_discrete());
			ASSERT(!s3.distribution().is_continuous());
			/* Despite having a distribution, s3 is equal to s1 because they have the same name. 
			States equality is currently equivalent to their name equality. */
			ASSERT(s3 == s1);
			/* State copy constructor copies the distribution. */
			State s6(s3);
			ASSERT(s6.distribution() == s3.distribution());
			/* DiscreteDistribution is constructible from an initializer_list of (std::string, double) pairs. */
			DiscreteDistribution dist2({{"A", 0.2}, {"G", 0.4}, {"C", 0.1}, {"T", 0.3}});
			ASSERT(dist2["A"] == 0.2);
			ASSERT((dist2["A"] = 0.5) == 0.5);
			/* Copy constructor/assignment operator for DiscreteDistribution. */
			DiscreteDistribution dist3 = dist2;
			ASSERT(dist3["A"] == 0.5);
			ASSERT(dist2 == dist3);
			/* Accessing a symbol not contained by the distribution will call the default
			constructor of the symbol and add it to the distribution
			(this has ultimately the same behavior as std::map). */
			double default_value = dist2["NotKey"];
			ASSERT(double() == default_value);
			/* Accessing "NotKey" created an entry in the distribution : (NotKey : 0). */
			ASSERT(dist2["NotKey"] == default_value);
			/* Thus, distributions will now differ. */
			ASSERT((dist2 != dist3));
		)
		
		TEST_UNIT(
			"begin/end state",
			State begin("begin");
			State end("end");
			/* Construct hmm by specifying begin and end state. */
			HiddenMarkovModel hmm = HiddenMarkovModel(begin, end);
			ASSERT(hmm.begin() == begin);
			ASSERT(hmm.end() == end);
			/* Check if the begin/end states were added to the hmm. */
			ASSERT(hmm.has_state(begin));
			ASSERT(hmm.has_state(end));
			/* Currently removing the begin state is allowed. */ //TODO
			hmm.remove_state(begin);
			ASSERT(!hmm.has_state(begin));
			/* Accessing begin when it has been removed throws an exception. */
			ASSERT_EXCEPT(hmm.begin(), StateNotFoundException);
		)

		TEST_UNIT(
			"add/remove state",
			HiddenMarkovModel hmm = HiddenMarkovModel();
			State s("s");
			/* Not yet added */
			ASSERT(!hmm.has_state("s"));
			ASSERT(!hmm.has_state(s));
			hmm.add_state(s);
			/* Check if the state was added to the hmm. */
			ASSERT(hmm.has_state("s"));
			ASSERT(hmm.has_state(s));
			/* Try to add an existing state. */
			ASSERT_EXCEPT(hmm.add_state(s), StateExistsException);
			hmm.remove_state("s");
			/* Remove a state not contained by the hmm. */
			ASSERT_EXCEPT(hmm.remove_state("s"), StateNotFoundException);
			ASSERT(!hmm.has_state("s"));
			ASSERT(!hmm.has_state(s));
		)

		TEST_UNIT(
			"add/remove transition",
			HiddenMarkovModel hmm = HiddenMarkovModel();
			State s1("s1");
			State s2("s2");
			hmm.add_state(s1);
			/* Throw an exception when a transition is added with a state not contained by the hmm. */
			ASSERT_EXCEPT(hmm.add_transition(s1, s2, 0.3), StateNotFoundException);
			hmm.add_state(s2);
			/* State is added, transition OK. */
			hmm.add_transition(s1, s2, 0.3);
			/* Check if transition exists. */
			ASSERT(hmm.has_transition(s1, s2));
			ASSERT(!hmm.has_transition(s2, s1));
			/* Removing it twice throws an exception. */
			hmm.remove_transition(s1, s2);
			ASSERT_EXCEPT(hmm.remove_transition(s1, s2), TransitionNotFoundException);
			/* Check successful removal. */
			ASSERT(!hmm.has_transition(s1, s2));
			/* Re-add it. */
			hmm.add_transition(s1, s2, 0.3);
			ASSERT(hmm.has_transition(s1, s2));
			/* Removing a state removes its transitions from and to other states. */
			hmm.remove_state(s1);
			ASSERT(!hmm.has_transition(s1, s2));
		)

		TEST_UNIT(
			"initial probability aka pi",
			State s1("s1");
			HiddenMarkovModel hmm;
			hmm.add_state(s1);
			/* Set initial probability by adding a transition to the being state of the hmm. */
			hmm.add_transition(hmm.begin(), s1, 0.4);
			ASSERT(hmm.has_state(s1));
			ASSERT(hmm.has_transition(hmm.begin(), s1));
			/* Or by directly calling begin_transition. */
			State s2("s2");
			hmm.add_state(s2);
			hmm.begin_transition(s2, 0.5);
			ASSERT(hmm.has_transition(hmm.begin(), s2));
		)

		TEST_UNIT(
			"brew",
			HiddenMarkovModel hmm;
			DiscreteDistribution dist1 = DiscreteDistribution({{"A",0.3}, {"T", 0.2}, {"G", 0.5}});
			hmm.add_state(State("s1", dist1));
			dist1["C"] = 0.2;
			hmm.add_state(State("s2", dist1));
			double s2_t = 0.5;
			hmm.add_transition("s1", "s2", s2_t);
			/* No begin transition throws exception. */
			ASSERT_EXCEPT(hmm.brew(), std::logic_error);
			hmm.add_transition(hmm.begin(), "s1", 1);
			/* Each state needs to have at least one out transition with prob > 0. */
			ASSERT_EXCEPT(hmm.brew(), std::logic_error);
			hmm.add_transition("s2","s1",1);
			/* Now OK ! */
			hmm.brew();
			std::size_t s1_index = hmm.states_indices()["s1"];
			std::size_t s2_index =  hmm.states_indices()["s2"];
			double brewed_transition = hmm.raw_transitions()[s1_index][s2_index];
			ASSERT(brewed_transition == log(1.0));
			hmm.add_state(State("s3", dist1));
			hmm.add_state(State("s4", dist1));
			double s3_t = 0.2;
			double s4_t = 0.3;
			hmm.add_transition("s1","s3", s3_t);
			hmm.add_transition("s1","s4", s4_t);
			/* Set sum of out transitions > 0 for each state. */
			hmm.add_transition("s3","s1", s3_t);
			hmm.add_transition("s4","s1", s4_t);
			hmm.brew();
			ASSERT(hmm.raw_transitions().size() == 4);
			ASSERT(hmm.raw_pi_begin().size() == 4);
			ASSERT(hmm.raw_pi_end().size() == 4);
			ASSERT(hmm.raw_pdfs().size() == 4);
			std::size_t s1_i = hmm.states_indices()["s1"];
			std::size_t s2_i =  hmm.states_indices()["s2"];
			std::size_t s3_i =  hmm.states_indices()["s3"];
			std::size_t s4_i =  hmm.states_indices()["s4"];
			double t_2 = hmm.raw_transitions()[s1_i][s2_i];
			double t_3 = hmm.raw_transitions()[s1_i][s3_i];
			double t_4 = hmm.raw_transitions()[s1_i][s4_i];
			/* Sum to 1, no normalization was needed, same probabilities. */
			ASSERT(t_2 == log(s2_t));
			ASSERT(t_3 == log(s3_t));
			ASSERT(t_4 == log(s4_t));
			ASSERT(utils::round_double(exp((*(hmm.raw_pdfs()[s2_index]))["A"])) == 0.25);
			hmm.add_state(State("s5", dist1));
			hmm.add_state(State("s6", dist1));
			double s5_t = 0.2; 
			double s6_t = 0.6;
			hmm.remove_transition("s2","s1");
			hmm.add_transition("s2","s5",s5_t);
			hmm.add_transition("s2","s6",s6_t);
			hmm.add_transition("s5","s2",s5_t);
			hmm.add_transition("s6","s2",s6_t);
			hmm.brew();
			std::size_t s2_n = hmm.states_indices()["s2"];
			std::size_t s5_n = hmm.states_indices()["s5"];
			std::size_t s6_n = hmm.states_indices()["s6"];
			double t_5 = hmm.raw_transitions()[s2_n][s5_n];
			double t_6 = hmm.raw_transitions()[s2_n][s6_n];
			/* sum is 0.2 + 0.6 thus with normalization 0.2 becomes 0.25, 0.6 become 0.75. */
			ASSERT(utils::round_double(exp(t_5)) == 0.25);
			ASSERT(utils::round_double(exp(t_6)) == 0.75);
		)

		/* Forward */
		
		/* Backward */

		/* Likelihood */

		/* Sample */

		/* Decode Viterbi */

		/* Train B-W */

		/* Train Viterbi */

		/* Train stochastic EM */


		tests_results();

	} catch(const std::exception& e){
		print_exception(e);
		return 1;
	}
	return 0;
}
