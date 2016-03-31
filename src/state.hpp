#ifndef __STATE_HPP
#define __STATE_HPP

#include <vector>
#include <exception>
#include "constants.hpp"
#include "distribution.hpp"

/* <-------- Exceptions --------> */

class StateException : public std::logic_error {
protected:
	StateException(const std::string& message) :
		std::logic_error(message) {}
};

class StateDistributionException : public StateException {
public:
	template<typename T>
	StateDistributionException(const T& t) : 
		StateException(error_message::format("StateDistributionException: " + error_message::kSilentStateHasNoDistribution, t)) {}
};

/* <----------------------------> */

class State{
private:
	/* State name */
	std::string _name;
	/* Emission probabilities */
	Distribution* _distribution;

public:
	State(const std::string& name) : _name(name), _distribution(nullptr) {}

	State(const char* c_str) : _name(), _distribution(nullptr){
		if(c_str != nullptr){
			_name = c_str;
		}
	}

	template<typename DistributionType>
	explicit State(const std::string& name, const DistributionType& distribution) : 
		_name(name), _distribution(nullptr) {
			_distribution = new DistributionType(distribution);
	} 

	State(const State& other) : _name(other._name), _distribution(nullptr) {
		if(other._distribution != nullptr){
			_distribution = other._distribution->clone();
		}
	}

	State(State&& other) : _name(other._name), _distribution(other._distribution) {
		other._distribution = nullptr;
	}

	State& operator=(const State& other){
		if(this != &other){
			_name = other._name;
			if(_distribution != nullptr){
				delete _distribution;
			}
			_distribution = other._distribution->clone();
		}
		return *this;
	}

	State& operator=(State&& other){
		if(this != &other){
			_name = other._name;
			if(_distribution != nullptr){
				delete _distribution;
			}
			_distribution = other._distribution;
			other._distribution = nullptr;
		}
		return *this;
	}

	inline bool operator==(const State& other) const {
		return _name == other.name();
	}

	Distribution& distribution() const {
		if(_distribution != nullptr){
			return *_distribution;	
		}
		else{
			throw StateDistributionException(*this);
		}
	}

	std::string name() const { return _name; }

	bool is_silent() const { return _distribution == nullptr; }

	virtual ~State(){
		delete _distribution;
	}

};

std::ostream& operator<<(std::ostream& out, const State& state){
	out << state.name() << "("
		<< (state.is_silent() ? "silent" : state.distribution().name())
		<< ")";
	return out;
}


#endif