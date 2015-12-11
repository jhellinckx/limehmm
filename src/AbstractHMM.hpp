#ifndef __HMM_HPP
#define __HMM_HPP 

#include <exception>
#include <stdexcept>
#include <vector>
#include <map>

template<typename S, typename O> 
class AbstractHMM {
private:
	std::vector<S> _states;
	std::vector<double> _transP;
	std::vector<double> _initP;
	std::map<S, std::vector<S>> _outStates;

protected:
	AbstractHMM(
		const std::vector<S>&, 
		const std::vector<double>&,
		const std::vector<double>&);

	virtual bool state_inBounds(const std::size_t i) const{
		return (i>0 && i<states());
	}

	virtual bool observation_inBounds(const std::size_t) const = 0;

public:
	/* Transition probabilities getters */
	virtual double trans_p(std::size_t, std::size_t) const;
	virtual double trans_p_by_object(S, S) const;

	/* Emission probabilities getters, defined in derived classes */
	virtual double emi_p(std::size_t, std::size_t) const = 0;
	virtual double emi_p_by_object(S, O) const = 0;

	/* Initial state probabilites getters */
	virtual double init_p(std::size_t) const;
	virtual double init_p_by_object(S) const; 

	/* Decoding */
	virtual void viterbi(const std::vector<O>&) = 0;

	/* Training procedures */
	virtual void train_viterbi(const std::vector<O>&) = 0;
	virtual void train_baumWelch(const std::vector<O>&) = 0;
	virtual void train_stochasticEM(const std::vector<O>&) = 0;

	virtual ~AbstractHMM(){}

	virtual std::size_t find_state(const S state) const{
		for(std::size_t i=0;i<_states.size();++i){
			if(_states[i] == state){
				return i;
			}
		}
		throw std::out_of_range("State not found");
	}

	virtual std::size_t find_observation(const O obs) const = 0;

	virtual std::size_t states() const { return _states.size(); };
	virtual std::size_t observations() const = 0;


	/* Digraph methods */
	virtual bool hasSuccessor(S from, S to){
		find_state(from); find_state(to);
		for(S outState : _outStates[from]){
			if(outState == to)
				return true;
		}
		return false;
	}

	virtual void addSuccessor(S from, S to){
		find_state(from); find_state(to);
		if(!hasSuccessor(from, to))
			_outStates[from].push_back(to);
	}

	virtual void setSuccessors(S from, const std::vector<S>& successors){
		for(S state : successors)
			find_state(state);
		_outStates[from].resize(successors.size());
		std::copy(successors.begin(), successors.end(), _outStates[from].begin());
	}

	std::vector<S> successors(S from) {
		return _outStates[from];
	}

	virtual void to_dot_file() const{

	}


};

template<typename S, typename O>
AbstractHMM<S, O>::AbstractHMM(	
	const std::vector<S>& states,
	const std::vector<double>& transP,
	const std::vector<double>& initP):
		_states(std::vector<S>(states.size())), 
		_transP(std::vector<double>(states.size()*states.size())), 
		_initP(std::vector<double>(states.size())),
		_outStates(std::map<S, std::vector<S>>()){
			for(std::size_t i=0;i<states.size();++i){
				try{
					find_state(states[i]);
					throw std::logic_error("Two identical states were given in states initializer vector");
				} catch(std::out_of_range e){
					_states[i] = states[i];
				}
				
			}
			if(transP.size() != states.size()*states.size())
				throw std::logic_error("Wrong vector size for transitions probabilities");
			if(initP.size() != states.size())
				throw std::logic_error("Wrond vector size for initial probabilities");
			std::copy(transP.begin(), transP.end(), _transP.begin());
			std::copy(initP.begin(), initP.end(), _initP.begin());
			/* By default, outdegree(state) = nstate */
			for(S state : _states){
				_outStates[state] = std::vector<S>(_states.size());
				std::copy(_states.begin(), _states.end(), _outStates[state].begin());
			}
		}

template<typename S, typename O>
double AbstractHMM<S, O>::trans_p(std::size_t i, std::size_t j) const{
	if(!state_inBounds(i) || !state_inBounds(j))
		throw std::out_of_range("State out of bounds");
	return _transP[i*states() + j];
}

template<typename S, typename O>
double AbstractHMM<S, O>::trans_p_by_object(S firstState, S secondState) const{
	return _transP[find_state(firstState)*states() + find_state(secondState)];
}

template<typename S, typename O>
double AbstractHMM<S, O>::init_p(std::size_t i) const{
	if(!state_inBounds(i))
		throw std::out_of_range("State out of bounds");
	return _initP[i];
}

template<typename S, typename O>
double AbstractHMM<S, O>::init_p_by_object(S state) const{
	return _initP[find_state(state)];
} 

#endif