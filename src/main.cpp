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
	std::cout << BOLDRED << failed << " failed " << RESET;
	std::cout << BOLDGREEN << successful << " succeeded" << RESET << std::endl;
}

int main(){
	try{
		std::vector<char> observations({'a', 'b'});
		std::vector<int> states({1, 2});
		std::vector<double> transP({0.3,0.7,
									0.7,0.3});
		std::vector<double> initP({0.0,1});
		std::vector<double> emiP({	0.1,0.9,
									0.9,0.1});
		LMDHMM<int, char> hmm(observations, emiP, states, transP, initP);

		
		printInit();
		ASSERT("observations.size() == hmm.observations()", observations.size() == hmm.observations());
		ASSERT("states.size() == hmm.states()",states.size() == hmm.states());
		printResults();
		
	} catch(const std::exception& e){
		printException(e);
		return 1;
	}
	return 0;
}
