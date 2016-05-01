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
#include <utility>
#include "utils.hpp"
#include "hmm.hpp" // tested hmm library

#define ASSERT(expr) assert(expr, #expr, __FILE__, __LINE__)
#define ASSERT_VERBOSE(expr, msg) assert(expr, #expr, __FILE__, __LINE__, msg)
#define ASSERT_ABORT(expr, msg) assert(expr, #expr, __FILE__, __LINE__, msg, true)
#define ASSERT_EXCEPT(instruction, except_type) assert_except<except_type>([&](){instruction;}, #instruction, #except_type, __FILE__, __LINE__)

#define TEST_UNIT(name, instructions) run_unit_test(name, [&](){instructions});

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

unsigned int assertions = 0;
unsigned int units = 0;
unsigned int failed = 0;
unsigned int successful = 0;

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
		/* Create hmm examples. */

		/* Simple fair/biased model. */
		HiddenMarkovModel casino_hmm("casino");
		DiscreteDistribution fair_dist({{"H", 0.5}, {"T", 0.5}});
		DiscreteDistribution biased_dist({{"H", 0.75}, {"T", 0.25}});
		State fair = State("fair", fair_dist);
		State biased = State("biased", biased_dist);
		casino_hmm.add_state(fair);
		casino_hmm.add_transition(casino_hmm.begin(), fair, 0.5);
		casino_hmm.add_state(biased);
		casino_hmm.add_transition(casino_hmm.begin(), biased, 0.5);
		casino_hmm.add_transition(fair, fair, 0.9);
		casino_hmm.add_transition(fair, biased, 0.1);
		casino_hmm.add_transition(biased, biased, 0.9);
		casino_hmm.add_transition(biased, fair, 0.1);
		casino_hmm.brew();
		/* Precomputed casino values. */
		std::vector<std::string> casino_symbols({"T","H","H","T","T","T","H","H"});
		double casino_precomputed_init_fwd_fair = 0.25;
		double casino_precomputed_init_fwd_biased = 0.125;
		double casino_precomputed_mid_fwd_fair = 0.0303;
		double casino_precomputed_mid_fwd_biased = 0.0191;
		double casino_precomputed_end_fwd_fair = 0.0015;
		double casino_precomputed_end_fwd_biased = 0.0013;
		double casino_precomputed_init_bwd_fair = 1;
		double casino_precomputed_init_bwd_biased = 1;
		double casino_precomputed_mid_bwd_fair = 0.0679;
		double casino_precomputed_mid_bwd_biased = 0.0366;
		double casino_precomputed_end_bwd_fair = 0.0075;
		double casino_precomputed_end_bwd_biased = 0.0071;
		double casino_precomputed_likelihood = 0.0028;
		std::vector<std::string> casino_precomputed_viterbi_path_2_states({"fair", "fair", "fair", "fair", "fair", "fair", "fair", "fair"});
		/* Training casino. */
		std::vector<std::vector<std::string>> casino_training_sequences_1 = 
		{{"H", "T", "H"}, {"H", "H", "H", "H"}, {"T", "H", "T", "H", "T"}, 
		{"T", "H", "T", "H", "T", "H", "T", "H", "T", "H"}, 
		{"T", "T", "H", "H", "H", "H", "H", "H", "H", "H", "H"}, 
		{"T", "H", "T", "H", "T", "H", "T", "H", "T", "H", "T", "H", "T", "H", "T", "H"}, 
		{"T", "T", "T"}};	
		double precomputed_casino_viterbi_improvement_no_pseudocount = utils::round_double(0.639841864861836, 4);
		double precomputed_casino_viterbi_improvement_with_pseudocount = utils::round_double(1.7764802457455673, 4);
		double precomputed_casino_bw_improvement_no_pseudocount = utils::round_double(2.6887270225223574);

		double precomputed_casino_bw_improvement_no_pseudocount_2 = utils::round_double(5.05069902785, 4);

		std::vector<std::vector<std::string>> casino_training_sequences_2 = 
		{	{"T", "H", "H", "T"}, {"T", "H", "H", "T"}, {"T", "H", "H", "T"}, 
			{"T", "H", "T", "H"}, {"T", "T", "T", "T"}, {"T", "T", "T", "T"}, 
			{"T", "H", "T", "H"}, {"H", "T", "H", "H"}, {"H", "T", "H", "H"}};

		std::vector<std::vector<std::string>> casino_training_sequences_3 = 
		{{"T", "H", "T", "T", "T", "H"}};

		/* Simple hmm with 3 states emitting nucleobases. */
		HiddenMarkovModel nucleobase_3_states_hmm("nucleobase 3 states");
		DiscreteDistribution dist1({{"A", 0.35}, {"C", 0.20}, {"G", 0.05}, {"T", 0.40}});
		DiscreteDistribution dist2({{"A", 0.25}, {"C", 0.25}, {"G", 0.25}, {"T", 0.25}});
		DiscreteDistribution dist3({{"A", 0.10}, {"C", 0.40}, {"G", 0.40}, {"T", 0.10}});
		State s1("s1", dist1);
		State s2("s2", dist2);
		State s3("s3", dist3);
		nucleobase_3_states_hmm.add_state(s1);
		nucleobase_3_states_hmm.add_state(s2);
		nucleobase_3_states_hmm.add_state(s3);
		nucleobase_3_states_hmm.begin_transition(s1, 0.90);
		nucleobase_3_states_hmm.begin_transition(s2, 0.10);
		nucleobase_3_states_hmm.add_transition(s1, s1, 0.80);
		nucleobase_3_states_hmm.add_transition(s1, s2, 0.20);
		nucleobase_3_states_hmm.add_transition(s2, s2, 0.30);
		nucleobase_3_states_hmm.add_transition(s2, s3, 0.10);
		nucleobase_3_states_hmm.add_transition(s3, s3, 0.70);
		nucleobase_3_states_hmm.end_transition(s3, 0.30);
		nucleobase_3_states_hmm.end_transition(s2, 0.60);
		nucleobase_3_states_hmm.brew(); 
		/* Precomputed values. */
		std::vector<std::string> nucleobase_symbols({"A", "C", "G", "A", "C", "T", "A", "T", "T", "C", "G", "A", "T"});
		double nucleobase_precomputed_viterbi_log_likelihood = utils::round_double(-23.834436455461574, 4);
		std::vector<std::string> nucleobase_precomputed_viterbi_path_3_states({"s1", "s1", "s1", "s1", "s1", "s1", "s1", "s1", "s1", "s1", "s1", "s1", "s2"});
		
		std::vector<std::vector<std::string>> nucleobase_training_sequences;
		nucleobase_training_sequences.push_back(nucleobase_symbols);

		/* Profile hmm with 10 states. */
		HiddenMarkovModel profile_10_states_hmm("profile 10 states");
		DiscreteDistribution i_d({{"A", 0.25}, {"C", 0.25}, {"G", 0.25}, {"T", 0.25}});
		/* Create insert states. */
		State i0 = State("I0", i_d);
		State i1 = State("I1", i_d);
		State i2 = State("I2", i_d);
		State i3 = State("I3", i_d);
		/* Create match states. */
		State m1 = State("M1", DiscreteDistribution({{"A", 0.95},  {"C", 0.01}, {"G", 0.01},  {"T", 0.02 }}));
		State m2 = State("M2", DiscreteDistribution({{"A", 0.003}, {"C", 0.99}, {"G", 0.003}, {"T", 0.004}}));
		State m3 = State("M3", DiscreteDistribution({{"A", 0.01},  {"C", 0.01}, {"G", 0.01},  {"T", 0.97 }}));
		/* Create delete states. */
		State d1 = State("D1");
		State d2 = State("D2");
		State d3 = State("D3");
		/* Add all the states. */
		profile_10_states_hmm.add_state(i0);
		profile_10_states_hmm.add_state(i1);
		profile_10_states_hmm.add_state(i2);
		profile_10_states_hmm.add_state(i3);
		profile_10_states_hmm.add_state(m1);
		profile_10_states_hmm.add_state(m2);
		profile_10_states_hmm.add_state(m3);
		profile_10_states_hmm.add_state(d1);
		profile_10_states_hmm.add_state(d2);
		profile_10_states_hmm.add_state(d3);
		/* Transitions from match states. */
		profile_10_states_hmm.add_transition(profile_10_states_hmm.begin(), m1, 0.9);
		profile_10_states_hmm.add_transition(profile_10_states_hmm.begin(), i0, 0.1);
		profile_10_states_hmm.add_transition(m1, m2, 0.9);
		profile_10_states_hmm.add_transition(m1, i1, 0.05);
		profile_10_states_hmm.add_transition(m1, d2, 0.05);
		profile_10_states_hmm.add_transition(m2, m3, 0.9);
		profile_10_states_hmm.add_transition(m2, i2, 0.05);
		profile_10_states_hmm.add_transition(m2, d3, 0.05);
		profile_10_states_hmm.add_transition(m3, profile_10_states_hmm.end(), 0.9);
		profile_10_states_hmm.add_transition(m3, i3, 0.1);
		/* Transitions from insert states. */
		profile_10_states_hmm.add_transition(i0, i0, 0.70);
		profile_10_states_hmm.add_transition(i0, d1, 0.15);
		profile_10_states_hmm.add_transition(i0, m1, 0.15);
		profile_10_states_hmm.add_transition(i1, i1, 0.70);
		profile_10_states_hmm.add_transition(i1, d2, 0.15);
		profile_10_states_hmm.add_transition(i1, m2, 0.15);
		profile_10_states_hmm.add_transition(i2, i2, 0.70);
		profile_10_states_hmm.add_transition(i2, d3, 0.15);
		profile_10_states_hmm.add_transition(i2, m3, 0.15);
		profile_10_states_hmm.add_transition(i3, i3, 0.85);
		profile_10_states_hmm.add_transition(i3, profile_10_states_hmm.end(), 0.15);
		/* Transitions from delete states. */
		profile_10_states_hmm.add_transition(d1, d2, 0.15);
		profile_10_states_hmm.add_transition(d1, i1, 0.15);
		profile_10_states_hmm.add_transition(d1, m2, 0.70);
		profile_10_states_hmm.add_transition(d2, d3, 0.15);
		profile_10_states_hmm.add_transition(d2, i2, 0.15);
		profile_10_states_hmm.add_transition(d2, m3, 0.70);
		profile_10_states_hmm.add_transition(d3, i3, 0.30);
		profile_10_states_hmm.add_transition(d3, profile_10_states_hmm.end(), 0.70);
		profile_10_states_hmm.brew(false);
		/* Precomputed values. */
		std::vector<std::vector<std::string>> profile_sequences = {{"A","C","T"}, {"G","G","C"},{"G","A","T"},{"A","C","C"}};
		double profile_seq1_viterbi_log_likelihood_precomputed = utils::round_double(-0.513244900357, 6);
		double profile_seq2_viterbi_log_likelihood_precomputed = utils::round_double(-11.0481012413, 6);
		double profile_seq3_viterbi_log_likelihood_precomputed = utils::round_double(-9.12551967402, 6);
		double profile_seq4_viterbi_log_likelihood_precomputed = utils::round_double(-5.08795587886, 6);
		std::vector<std::string> profile_seq1_viterbi_path_precomputed({"M1", "M2", "M3"});
		std::vector<std::string> profile_seq2_viterbi_path_precomputed({"I0", "I0", "D1", "M2", "D3"});
		std::vector<std::string> profile_seq3_viterbi_path_precomputed({"I0", "M1", "D2", "M3"});
		std::vector<std::string> profile_seq4_viterbi_path_precomputed({"M1", "M2", "M3"});
		std::vector<double> precomputed_viterbi_log_likelihoods = {-5.406181012423981, -10.88681993576597, 
		-3.6244718790494277, -3.644880750680635, -10.674332964640293, -10.393824835172445, -8.67126440174503, 
		-16.903451796110275, -16.451699654050792, -0.5132449003570658, -11.048101241343396, -9.125519674022627, 
		-5.0879558788604475};
		for(double& likelihood : precomputed_viterbi_log_likelihoods){ likelihood = utils::round_double(likelihood); }
		std::vector<std::vector<std::string>> profile_viterbi_decode_sequences = 
		{{"A"}, {"G", "A"}, {"A", "C"}, {"A", "T"}, {"A", "T", "C", "C"}, 
		{"A", "C", "G", "T", "G"}, {"A", "T", "T", "T"}, {"T", "A", "C", "C", "C", "T", "C"}, 
		{"T", "G", "T", "C", "A", "A", "C", "A", "C", "T"}, {"A", "C", "T"}, {"G", "G", "C"},
		{"G", "A", "T"}, {"A", "C", "C"}};
		std::vector<std::vector<std::string>> profile_training_sequences = 
		{{"A", "C", "T"}, {"A", "C", "T"}, {"A", "C", "C"}, {"A", "C", "T", "C"}, 
		{"A", "C", "T"}, {"A", "C", "T"}, {"C", "C", "T"}, {"C", "C", "C"}, 
		{"A", "A", "T"}, {"C", "T"}, {"A", "T"}, {"C", "T"}, {"C", "T"}, {"C", "T"}, 
		{"C", "T"}, {"C", "T"}, {"C", "T"}, {"A", "C", "T"}, {"A", "C", "T"}, 
		{"C", "T"}, {"A", "C", "T"}, {"C", "T"}, {"C", "T"}, {"C", "T"}, {"C", "T"}};
		double precomputed_profile_improvement_no_pseudocount = 84.9318;
		double precomputed_profile_improvement_with_pseudocount = 78.9441;
		std::vector<double> precomputed_observation_likelihood = 
		{-0.505202786679,-0.505202786679,-4.85332533917,-5.86447364745,-0.505202786679,-0.505202786679,
		-4.77279340553,-8.94719166767,-5.71623676018,-7.35302419896,-3.61187131339,-7.35302419896,
		-7.35302419896,-7.35302419896,-7.35302419896,-7.35302419896,-7.35302419896,-0.505202786679,
		-0.505202786679,-7.35302419896,-0.505202786679,-7.35302419896,-7.35302419896,-7.35302419896,
		-7.35302419896};
		for(double& likelihood : precomputed_observation_likelihood){ likelihood = utils::round_double(likelihood); }
		std::vector<std::vector<std::string>>& profile_observation_likelihood_sequences = profile_training_sequences;


		tests_init();
		
		/* IEEE 754 floating points are required in order to use std::infinity numeric limit. */
		TEST_UNIT(
			"platform type",
			ASSERT_ABORT(std::numeric_limits<double>::is_iec559, "IEEE 754 required");
		)

		TEST_UNIT(
			"graph",
			/* Topological sort. */
			Graph<std::string> g;
			g.add_vertex("B");
			g.add_vertex("E");
			g.add_vertex("A");
			g.add_vertex("D");
			g.add_vertex("C");
			g.add_edge("A", "B");
			g.add_edge("A", "D");
			g.add_edge("B", "C");
			g.add_edge("C", "D");
			g.add_edge("D", "E");
			g.add_edge("C", "E");
			std::vector<std::string> precomputed_toposort({"A", "B", "C", "D", "E"}); 
			g.topological_sort();
			std::vector<std::string*> vertices = g.get_vertices();
			std::vector<std::string> toposort;
			for(std::string* str_ptr : vertices) {
				toposort.push_back(*str_ptr);
			}
			ASSERT(toposort == precomputed_toposort);
			/* Subgraph. */
			std::vector<std::string> sub_vertices({"C","E"});
			Graph<std::string> subgraph = g.sub_graph(sub_vertices);
			std::vector<std::string> sub_graph_vertices;
			std::vector<std::string*> sub_ptrs = subgraph.get_vertices();
			for(std::string* sub_ptr : sub_ptrs){
				sub_graph_vertices.push_back(*sub_ptr);
			}
			ASSERT(sub_graph_vertices == sub_vertices);
			ASSERT(subgraph.has_edge("C", "E"));
			ASSERT(!subgraph.has_vertex("A"));
			ASSERT(!subgraph.has_vertex("B"));
			ASSERT(!subgraph.has_vertex("D"));
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
			/* Add transition where from == to. */
			hmm.add_transition(s1, s1, 0.9);
			ASSERT(hmm.has_transition(s1, s1));
			/* Remove it. */
			hmm.remove_transition(s1, s1);
			ASSERT(!hmm.has_transition(s1, s1));
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

		TEST_UNIT(
			"forward",
			HiddenMarkovModel hmm = casino_hmm;
			std::vector<std::string> symbols = casino_symbols;
			/* Test init forward aka t = 1 aka first forward column. */
			/* Second parameter gives t. */
			std::vector<double> init_fwd = hmm.forward(symbols, 1);
			ASSERT(init_fwd.size() == 2);
			double init_fwd_fair = utils::round_double(exp(init_fwd[hmm.states_indices()["fair"]]), 2);
			double init_fwd_biased = utils::round_double(exp(init_fwd[hmm.states_indices()["biased"]]), 3);
			ASSERT(init_fwd_fair == casino_precomputed_init_fwd_fair);
			ASSERT(init_fwd_biased == casino_precomputed_init_fwd_biased);
			/* Test middle column. */
			std::vector<double> mid_fwd = hmm.forward(symbols, 4);
			double mid_fwd_fair = utils::round_double(exp(mid_fwd[hmm.states_indices()["fair"]]), 4);
			double mid_fwd_biased = utils::round_double(exp(mid_fwd[hmm.states_indices()["biased"]]), 4);
			ASSERT(mid_fwd_fair == casino_precomputed_mid_fwd_fair);
			ASSERT(mid_fwd_biased == casino_precomputed_mid_fwd_biased);
			/* Test last fwd column. */
			/* T not given, iterate on all the symbols. */
			std::vector<double> fwd_end = hmm.forward(symbols);
			double end_fwd_fair = utils::round_double(exp(fwd_end[hmm.states_indices()["fair"]]), 4);
			double end_fwd_biased = utils::round_double(exp(fwd_end[hmm.states_indices()["biased"]]), 4);
			ASSERT(end_fwd_fair == casino_precomputed_end_fwd_fair);
			ASSERT(end_fwd_biased == casino_precomputed_end_fwd_biased);
		)

		TEST_UNIT(
			"backward",
			/* Same model as forward test. */
			HiddenMarkovModel hmm = casino_hmm;
			std::vector<std::string> symbols = casino_symbols;
			/* Test init backward. */
			std::vector<double> init_bwd = hmm.backward(symbols, symbols.size());
			ASSERT(init_bwd.size() == 2);
			double init_bwd_fair = utils::round_double(exp(init_bwd[hmm.states_indices()["fair"]]), 3);
			double init_bwd_biased = utils::round_double(exp(init_bwd[hmm.states_indices()["biased"]]), 3);
			ASSERT(init_bwd_fair == casino_precomputed_init_bwd_fair);
			ASSERT(init_bwd_biased == casino_precomputed_init_bwd_biased);
			/* Test middle column. */
			std::vector<double> mid_bwd = hmm.backward(symbols, 4);
			double mid_bwd_fair = utils::round_double(exp(mid_bwd[hmm.states_indices()["fair"]]), 4);
			double mid_bwd_biased = utils::round_double(exp(mid_bwd[hmm.states_indices()["biased"]]), 4);
			ASSERT(mid_bwd_fair == casino_precomputed_mid_bwd_fair);
			ASSERT(mid_bwd_biased == casino_precomputed_mid_bwd_biased);
			/* Test last backward column. */
			std::vector<double> bwd_end = hmm.backward(symbols);
			double end_bwd_fair = utils::round_double(exp(bwd_end[hmm.states_indices()["fair"]]), 4);
			double end_bwd_biased = utils::round_double(exp(bwd_end[hmm.states_indices()["biased"]]), 4);
			ASSERT(end_bwd_fair == casino_precomputed_end_bwd_fair);
			ASSERT(end_bwd_biased == casino_precomputed_end_bwd_biased);
		)

		TEST_UNIT(
			"observation likelihood (casino)",
			/* Same model as fwd/bwd. */
			HiddenMarkovModel hmm = casino_hmm;
			/* Test likelihood with forward algorithm. */
			double forward_observation_likelihood = utils::round_double(hmm.likelihood(casino_symbols), 4);
			ASSERT(forward_observation_likelihood == casino_precomputed_likelihood);
			/* Test likelihood with backward algorithm. */
			double backward_observation_likelihood = utils::round_double(hmm.likelihood(casino_symbols, false), 4);
			ASSERT(backward_observation_likelihood == casino_precomputed_likelihood);
		)

		TEST_UNIT(
			"observation likelihood (profile)",
				HiddenMarkovModel hmm = profile_10_states_hmm;
				std::size_t n_tests = 2;
				/* Randomly choose n_tests sequences. Useless to test all 25 sequences. */
				for(std::size_t i = 0; i < n_tests; ++i){
					std::size_t random_sequence = (std::size_t)rand() % profile_observation_likelihood_sequences.size();	
					const std::vector<std::string>& seq = profile_observation_likelihood_sequences[random_sequence];
					double forward_observation_likelihood = utils::round_double(hmm.log_likelihood(seq));
					ASSERT(forward_observation_likelihood == precomputed_observation_likelihood[random_sequence]);
					double backward_observation_likelihood = utils::round_double(hmm.log_likelihood(seq,false));
					ASSERT(backward_observation_likelihood == precomputed_observation_likelihood[random_sequence]);
				}
		)

		TEST_UNIT(
			"viterbi decode (casino)",
			HiddenMarkovModel hmm = casino_hmm;
			std::vector<std::string> symbols = casino_symbols;
			std::vector<std::string> viterbi_path_2_states = hmm.decode(symbols).first;
			ASSERT(viterbi_path_2_states == casino_precomputed_viterbi_path_2_states);
		)

		TEST_UNIT(
			"viterbi decode/likelihood (nucleobase)",
			HiddenMarkovModel hmm = nucleobase_3_states_hmm;
			auto viterbi_decode = hmm.decode(nucleobase_symbols);
			std::vector<std::string> viterbi_path_3_states = viterbi_decode.first;
			double viterbi_log_likelihood = utils::round_double(viterbi_decode.second, 4);
			std::cout << viterbi_log_likelihood << std::endl;
			ASSERT(viterbi_log_likelihood == nucleobase_precomputed_viterbi_log_likelihood);
			ASSERT(viterbi_path_3_states == nucleobase_precomputed_viterbi_path_3_states);
		)

		TEST_UNIT(
			"viterbi decode (profile)",
			HiddenMarkovModel hmm = profile_10_states_hmm;
			auto viterbi_seq1 = hmm.decode(profile_sequences[0]);
			auto viterbi_seq2 = hmm.decode(profile_sequences[1]);
			auto viterbi_seq3 = hmm.decode(profile_sequences[2]);
			auto viterbi_seq4 = hmm.decode(profile_sequences[3]);
			double seq1_viterbi_log_likelihood = utils::round_double(viterbi_seq1.second, 6);
			double seq2_viterbi_log_likelihood = utils::round_double(viterbi_seq2.second, 6);
			double seq3_viterbi_log_likelihood = utils::round_double(viterbi_seq3.second, 6);
			double seq4_viterbi_log_likelihood = utils::round_double(viterbi_seq4.second, 6);
			ASSERT(seq1_viterbi_log_likelihood == profile_seq1_viterbi_log_likelihood_precomputed);
			ASSERT(seq2_viterbi_log_likelihood == profile_seq2_viterbi_log_likelihood_precomputed);
			ASSERT(seq3_viterbi_log_likelihood == profile_seq3_viterbi_log_likelihood_precomputed);
			ASSERT(seq4_viterbi_log_likelihood == profile_seq4_viterbi_log_likelihood_precomputed);
			std::vector<std::string> seq1_viterbi_path = viterbi_seq1.first;
			std::vector<std::string> seq2_viterbi_path = viterbi_seq2.first;
			std::vector<std::string> seq3_viterbi_path = viterbi_seq3.first;
			std::vector<std::string> seq4_viterbi_path = viterbi_seq4.first;
			ASSERT(seq1_viterbi_path == profile_seq1_viterbi_path_precomputed);
			ASSERT(seq2_viterbi_path == profile_seq2_viterbi_path_precomputed);
			ASSERT(seq3_viterbi_path == profile_seq3_viterbi_path_precomputed);
			ASSERT(seq4_viterbi_path == profile_seq4_viterbi_path_precomputed);
			std::size_t n_tests = 2;
			for(std::size_t i = 0; i < n_tests ; ++i){
				std::size_t random_sequence = (std::size_t)rand() % profile_viterbi_decode_sequences.size();
				double viterbi_log_likelihood = utils::round_double(hmm.decode(profile_viterbi_decode_sequences[random_sequence]).second);
				ASSERT(viterbi_log_likelihood == precomputed_viterbi_log_likelihoods[random_sequence]);
			}
		)

		TEST_UNIT(
			"viterbi training (casino)",
			HiddenMarkovModel hmm = casino_hmm;
			double viterbi_improvement = utils::round_double(hmm.train_viterbi(casino_training_sequences_1), 4);
			ASSERT(viterbi_improvement == precomputed_casino_viterbi_improvement_no_pseudocount);
		)

		TEST_UNIT(
			"viterbi training with pseudocounts (casino)",
			HiddenMarkovModel hmm = casino_hmm;
			double viterbi_improvement = utils::round_double(hmm.train_viterbi(casino_training_sequences_1, 1.0), 4);
			ASSERT(viterbi_improvement == precomputed_casino_viterbi_improvement_with_pseudocount);
		)

		TEST_UNIT(
			"viterbi training with silent states (profile)",
			HiddenMarkovModel hmm = profile_10_states_hmm;
			double viterbi_improvement = utils::round_double(hmm.train_viterbi(profile_training_sequences), 4);
			ASSERT(viterbi_improvement == precomputed_profile_improvement_no_pseudocount);
		)

		TEST_UNIT(
			"viterbi training with pseudocounts and with silent states (profile)",
			HiddenMarkovModel hmm = profile_10_states_hmm;
			double viterbi_improvement = utils::round_double(hmm.train_viterbi(profile_training_sequences, 1.0), 4);
			ASSERT(viterbi_improvement == precomputed_profile_improvement_with_pseudocount);
		)

		TEST_UNIT(
			"baum-welch training 1 sequence (casino)",
			HiddenMarkovModel hmm = casino_hmm;
			// print_transitions(hmm.raw_transitions(), hmm.states_indices(), false);
			// print_pi_begin(hmm.raw_pi_begin(), hmm.states_names(), false);
			// print_pi_end(hmm.raw_pi_end(), hmm.states_names(), false);
			// print_distributions(hmm.raw_pdfs(), hmm.states_names(), false);
			// hmm.train_baum_welch(casino_training_sequences_3);
			// print_transitions(hmm.raw_transitions(), hmm.states_indices(), true);
			// print_pi_begin(hmm.raw_pi_begin(), hmm.states_names(), true);
			// print_pi_end(hmm.raw_pi_end(), hmm.states_names(), true);
			// print_distributions(hmm.raw_pdfs(), hmm.states_names(), false);
		)

		TEST_UNIT(
			"baum-welch training batch of sequences (casino)",
			HiddenMarkovModel hmm = casino_hmm;
			print_transitions(hmm.raw_transitions(), hmm.states_indices(), true);
			print_pi_begin(hmm.raw_pi_begin(), hmm.states_names(), true);
			print_pi_end(hmm.raw_pi_end(), hmm.states_names(), true);
			print_distributions(hmm.raw_pdfs(), hmm.states_names(), false);
			hmm.train_baum_welch(casino_training_sequences_3, 0.0, 2);
			print_transitions(hmm.raw_transitions(), hmm.states_indices(), false);
			print_pi_begin(hmm.raw_pi_begin(), hmm.states_names(), false);
			print_pi_end(hmm.raw_pi_end(), hmm.states_names(), false);
			print_distributions(hmm.raw_pdfs(), hmm.states_names(), false);
		)

		TEST_UNIT(
			"baum-welch training with end state (nucleobase)",
			// HiddenMarkovModel hmm = nucleobase_3_states_hmm;
			// print_transitions(hmm.raw_transitions(), hmm.states_indices(), false);
			// print_pi_begin(hmm.raw_pi_begin(), hmm.states_names(), false);
			// print_pi_end(hmm.raw_pi_end(), hmm.states_names(), false);
			// print_distributions(hmm.raw_pdfs(), hmm.states_names(), false);
			// hmm.train_baum_welch(nucleobase_training_sequences);
			// print_transitions(hmm.raw_transitions(), hmm.states_indices(), true);
			// print_pi_begin(hmm.raw_pi_begin(), hmm.states_names(), true);
			// print_pi_end(hmm.raw_pi_end(), hmm.states_names(), true);
			// print_distributions(hmm.raw_pdfs(), hmm.states_names(), false);
		)

		TEST_UNIT(
			"baum-welch training with silent states (profile)",
			// HiddenMarkovModel hmm = profile_10_states_hmm;
			// print_transitions(hmm.raw_transitions(), hmm.states_indices(), false);
			// print_pi_begin(hmm.raw_pi_begin(), hmm.states_names(), false);
			// print_pi_end(hmm.raw_pi_end(), hmm.states_names(), false);
			// print_distributions(hmm.raw_pdfs(), hmm.states_names(), false);
			// hmm.train_baum_welch(profile_training_sequences);
			// print_transitions(hmm.raw_transitions(), hmm.states_indices(), true);
			// print_pi_begin(hmm.raw_pi_begin(), hmm.states_names(), true);
			// print_pi_end(hmm.raw_pi_end(), hmm.states_names(), true);
			// print_distributions(hmm.raw_pdfs(), hmm.states_names(), false);
		)



		/* Test fix / free parameters */

		/* Test update from raw model */

		/* Train stochastic EM */

		/* Profile HMM */

		/* MLE */

		/* Randomized params */

		tests_results();

	} catch(const std::exception& e){
		print_exception(e);
		return 1;
	}
	return 0;
}
