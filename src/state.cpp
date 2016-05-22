#include <vector>
#include <exception>
#include "constants.hpp"
#include "distributions.hpp"
#include "state.hpp"

/* <-------- Exceptions --------> */

StateException::StateException(const std::string& message) :
		std::logic_error(message) {}


StateDistributionException::StateDistributionException(const std::string& t) : 
		StateException(error_message::format("StateDistributionException: " + error_message::kSilentStateHasNoDistribution, t)) {}

/* <----------------------------> */

State::State(const std::string& name) : _name(name), _distribution(nullptr), 
	_free_emission(hmm_config::kDefaultFreeEmission), 
	_free_transition(hmm_config::kDefaultFreeTransition) {}

State::State(const char* c_str) : State(std::string(c_str)) {}

State::State(const std::string& name, const Distribution& distribution) : 
	_name(name), _distribution(nullptr),
	_free_emission(hmm_config::kDefaultFreeEmission), 
	_free_transition(hmm_config::kDefaultFreeTransition) {
		_distribution = distribution.clone();
} 

State::State(const State& other) : _name(other._name), _distribution(nullptr), 
_free_emission(other._free_emission), _free_transition(other._free_transition) {
	if(other._distribution != nullptr){
		_distribution = other._distribution->clone();
	}
}

State::State(State&& other) : _name(std::move(other._name)), _distribution(other._distribution),
_free_emission(other._free_emission), _free_transition(other._free_transition) {
	other._distribution = nullptr;
}

State& State::operator=(const State& other){
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

State& State::operator=(State&& other){
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

bool State::operator==(const State& other) const {
	return _name == other.name();
}

bool State::operator!=(const State& other) const {
	return ! operator==(other);
}

std::string State::to_string() const {
	return _name;
}

bool State::has_free_emission() const { return _free_emission; }
bool State::has_free_transition() const { return _free_transition; }
void State::fix_emission() { _free_emission = false; }
void State::fix_transition() { _free_transition = false; }
void State::free_emission() { _free_emission = true; }
void State::free_transition() { _free_transition = true; }

Distribution& State::distribution() const {
	if(_distribution != nullptr){
		return *_distribution;	
	}
	else{
		throw StateDistributionException(this->to_string());
	}
}

void State::set_distribution(const Distribution& dist){
	if(_distribution != nullptr){
		delete _distribution;
	}
	_distribution = dist.clone();
}

std::string State::name() const { return _name; }
void State::set_name(const std::string& name) { _name = name; }   

bool State::is_silent() const { 
	return _distribution == nullptr || _distribution->empty();
}

State::~State(){
	if(_distribution != nullptr) { delete _distribution; }
}

std::ostream& operator<<(std::ostream& out, const State& state){
	out << state.to_string();
	return out;
}