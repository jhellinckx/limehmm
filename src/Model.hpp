#ifndef __MODEL_HPP
#define __MODEL_HPP

class Model{
private:
	std::vector<State> _states;
	std::vector<Symbol> _alphabet;

public:
	Model() : _states(std::vector<State>()), _alphabet(std::vector<Symbol>){

	}

	Model(const std::vector<State>& states) : 
		_states(std::vector<State>(states.size())), 
		_alphabet(std::vector<Symbol>()){
			std::copy(states.begin(), states;end(), _states.begin());
	}

	Model(const std::vector<Symbol>& alphabet) :
		_states(std::vector<State>()),
		_alphabet(std::vector<Symbol>(alphabet.size())){
			std::copy(alphabet.begin(), alphabet.end(), _alphabet.begin());
		}

	Model(const std::vector<State>& states, const std::vector<Symbol>& alphabet) :
		_states(std::vector<State>(states.size())),
		_alphabet(std::vector<Symbol>(alphabet.size())){
			std::copy(states.begin(), states.end(), _states.begin());
			std::copy(alphabet.begin(), alphabet.end(), _alphabet.begin());
		}

	virtual std::size_t states() const { return _states.size(); }
	virtual std::size_t symbols() const { return _alphabet.size(); }

	virtual State* getState(std::size_t i) {
		if(i >= 0 && i < _states.size())
			return &_states[i];
		return nullptr;
	}

	virtual State* getState(const std::string& id) {
		for(State& state : _states){
			if(state.identifier() == id)
				return &state;
		}
		return nullptr;
	}
};


#endif