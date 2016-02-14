#ifndef __STATE_HPP
#define __STATE_HPP

class State{
private:
	std::vector<double> _emissionProbabilities;
	std::vector<double> _transitionProbabilities;
	double _initProbability;
	std::string _identifier;

public:
	State(std::size_t statesNumber, std::size_t alphabetSize) : 
		_emissionProbabilities(alphabetSize),
		_transitionProbabilities(statesNumber),
		_initProbability() {

		}

	State(const std::vector<double>& emiP, const std::vector<double>& transP)
	std::string identifier() { return _identifier; }

	void newState() {

	}

};


#endif