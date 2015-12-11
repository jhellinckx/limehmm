#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <vector>
#include <string>

#include "LMDHMM.hpp"

static const std::string RED = "\033[31m";
static const std::string GREEN = "\033[32m";
static const std::string RESET = "\033[0m";
static const std::string BOLDRED = "\033[1;31m";
static const std::string BOLDGREEN = "\033[1;32m";
static const std::string BOLDCYAN = "\033[1;36m";
static const std::string BOLDMAGENTA = "\033[1;35m";

void printException(const std::exception& e){ std::cerr << RED << "*** " << e.what() << " ***" << RESET << std::endl; }

static const std::string separator = "*********************************************************************";

int assertions = 0;
int failed = 0;
int successful = 0;

void ASSERT(std::string description, bool assertion){
	++assertions;
	std::cout << description << " ? ... ";
	if(assertion == true){
		++successful;
		std::cout << BOLDGREEN << "OK" << RESET << " !" << std::endl;
	}
	else{
		++failed;
		std::cout << BOLDRED << "FAIL" << RESET << " !" << std::endl;
	}
}

void printInit(){
	std::cout << separator << std::endl;
	std::cout << BOLDMAGENTA << "Running " << "tests" << "..." << RESET << std::endl;
}

void printResults(){
	std::cout << separator << std::endl;
	std::cout <<  BOLDMAGENTA << "Ran " << assertions << " test(s) : " << RESET;
	std::cout << BOLDGREEN << successful << " succeeded " << RESET;
	std::cout << BOLDRED << failed << " failed" << RESET << std::endl;
	
}

int main(){
	try{
		std::vector<std::string> observations({"walk", "shop", "clean"});
		std::vector<std::string> states({"rainy", "sunny"});
		std::vector<double> initial_state_probabilities({0.6,0.4});
		
		std::vector<double> transition_probabilities({	0.7, 0.3,
														0.4, 0.6 	});
		
		std::vector<double> emission_probabilities({	0.1, 0.4, 0.5,
														0.6, 0.3, 0.1	});

		/* Create a Discrete Hidden Markov Model with Linear Memory */
		LMDHMM<std::string, std::string> hmm(observations, emission_probabilities, states, 
			transition_probabilities, initial_state_probabilities);

		
		printInit();
		/* Testing constructors */
		ASSERT("Observations number", observations.size() == hmm.observations());
		ASSERT("States number",states.size() == hmm.states());
		ASSERT("Emission probability by object",hmm.emi_p_by_object("rainy","clean") == 0.5);
		ASSERT("Emission probability by object", hmm.emi_p_by_object("sunny","shop") == 0.3);
		ASSERT("Transition probability by object",hmm.trans_p_by_object("rainy","sunny") == 0.3);
		ASSERT("Transition probability by object",hmm.trans_p_by_object("sunny","sunny") == 0.6);
		ASSERT("Initial state probability by object",hmm.init_p_by_object("rainy") == 0.6);
		
		/* Testing digraph */
		ASSERT("Successors initialization",hmm.successors("rainy") == states);
		hmm.addSuccessor("rainy","rainy");
		ASSERT("Adding existing successor",hmm.successors("rainy") == states);


		/* Testing decoding */
		//TODO

		/* Testing training */
		//TODO
		printResults();

		hmm.to_dot_file();
		
	} catch(const std::exception& e){
		printException(e);
		return 1;
	}
	return 0;
}
