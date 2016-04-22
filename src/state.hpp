#ifndef __STATE_HPP
#define __STATE_HPP

#include <vector>
#include <exception>
#include "constants.hpp"
#include "distributions.hpp"

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
	/* Free emission */
	bool _free_emission;
	/* Free transition */
	bool _free_transition;

public:
	State(const std::string& name) : _name(name), _distribution(nullptr), 
		_free_emission(hmm_config::kDefaultFreeEmission), 
		_free_transition(hmm_config::kDefaultFreeTransition) {}

	State(const char* c_str) : State(std::string(c_str)) {}

	explicit State(const std::string& name, const Distribution& distribution) : 
		_name(name), _distribution(nullptr),
		_free_emission(hmm_config::kDefaultFreeEmission), 
		_free_transition(hmm_config::kDefaultFreeTransition) {
			_distribution = distribution.clone();
	} 

	State(const State& other) : _name(other._name), _distribution(nullptr), 
	_free_emission(other._free_emission), _free_transition(other._free_transition) {
		if(other._distribution != nullptr){
			_distribution = other._distribution->clone();
		}
	}

	State(State&& other) : _name(std::move(other._name)), _distribution(other._distribution),
	_free_emission(other._free_emission), _free_transition(other._free_transition) {
		other._distribution = nullptr;
	}

	State& operator=(const State& other){
		if(this != &other){
			_name = other._name;
			_free_emission = other._free_emission;
			_free_transition = other._free_transition;
			if(_distribution != nullptr){
				delete _distribution;
			}
			_distribution = other._distribution->clone();
		}
		return *this;
	}

	State& operator=(State&& other){
		if(this != &other){
			_name = std::move(other._name);
			_free_emission = other._free_emission;
			_free_transition = other._free_transition;
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

	inline bool operator!=(const State& other) const {
		return ! operator==(other);
	}

	std::string to_string() const {
		std::string repr = _name + "(";
		repr += (is_silent() ? "silent" : distribution().name()) + ")";
		return repr;
	}

	bool has_free_emission() const { return _free_emission; }
	bool has_free_transition() const { return _free_transition; }
	void fix_emission() { _free_emission = false; }
	void fix_transition() { _free_transition = false; }
	void free_emission() { _free_emission = true; }
	void free_transition() { _free_transition = true; }

	Distribution& distribution() const {
		if(_distribution != nullptr){
			return *_distribution;	
		}
		else{
			throw StateDistributionException(*this);
		}
	}

	std::string name() const { return _name; }
	void set_name(const std::string& name) { _name = name; }   

	bool is_silent() const { 
		return _distribution == nullptr || _distribution->empty();
	}

	virtual ~State(){
		delete _distribution;
	}

};

std::ostream& operator<<(std::ostream& out, const State& state){
	out << state.to_string();
	return out;
}


#endif