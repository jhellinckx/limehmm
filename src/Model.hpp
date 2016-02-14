#ifndef __MODEL_HPP
#define __MODEL_HPP

class Model{
private:
	std::vector<State> _states;
	std::vector<std::string> _alphabet;

public:
	Model() : _states(std::vector<State>()), _observations(std::vector<Observation>){

	}

	Model(const std::vector<State>& states) : 
		_states(std::vector<State>(states.size())), 
		_alphabet(std::vector<std::string>()){
			std::copy(states.begin(), states;end(), _states.begin());
	}

	Model(const std::vector<Observation>& observations) :
		_states(std::vector<State>()),
		_observations(std::vector<Observation>(observations.size())){
			std::copy(observations.begin(), observations.end(), _observations.begin());
		}

	Model(const std::vector<State>& states, std::vector<Observation>& observations) :
		_states(std::vector<State>(states.size())),
		_observations(std::vector<Observation>(observations.size())){
			std::copy(states.begin(), states.end(), _states.begin());
			std::copy(observations.begin(), observations.end(), _observations.begin());
		}

	virtual std::size_t statesNumber() const { return _states.size(); }
	virtual std::size_t observationsNumber() const { return _observations.size(); }

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


	virtual double[] trans_p_matrix
};


#endif