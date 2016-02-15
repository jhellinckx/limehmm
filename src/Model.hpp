#ifndef __MODEL_HPP
#define __MODEL_HPP

#include <vector>
#include <iostream>
#include <sstream>

class Model{
private:
	std::vector<State> _states;
	std::vector<Symbol> _alphabet;

public:
	Model() : _states(std::vector<State>()), _alphabet(){}

	Model(const std::vector<State>& states) : 
		_states(states.size()), 
		_alphabet(){
			std::copy(states.begin(), states;end(), _states.begin());
	}

	Model(const std::vector<Symbol>& alphabet) :
		_states(),
		_alphabet(alphabet.size()){
			std::copy(alphabet.begin(), alphabet.end(), _alphabet.begin());
		}

	Model(const std::vector<State>& states, const std::vector<Symbol>& alphabet) :
		_states(states.size()),
		_alphabet(alphabet.size()){
			std::copy(states.begin(), states.end(), _states.begin());
			std::copy(alphabet.begin(), alphabet.end(), _alphabet.begin());
		}

	Model(const std::vector<std::string>& statesIdentifiers, const std::vector<std::string>& alphabetIdentifiers) :
		_states(statesIdentifiers.size()),
		_alphabet(alphabetIdentifiers.size()){
			for(std::size_t i = 0; i < statesIdentifiers.size(); ++i){
				_states[i] = State(statesIdentifiers[i], statesIdentifiers.size(), 
					alphabetIdentifiers.size());
			}
			for(std::size_t i = 0; i < alphabetIdentifiers.size(); ++i){
				_alphabet[i] = Symbol(alphabetIdentifiers[i]);
			}
		}

	std::size_t states() const { return _states.size(); }
	std::size_t alphabet() const { return _alphabet.size(); }

	/* Overloaded state getter methods */
	virtual State* get_state(std::size_t i) {
		if(i >= 0 && i < _states.size())
			return &_states[i];
		return nullptr;
	}

	virtual State* get_state(const std::string& id) {
		for(State& state : _states){
			if(state.identifier() == id)
				return &state;
		}
		return nullptr;
	}
};


#endif