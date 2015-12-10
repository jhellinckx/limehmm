#ifndef __HMM_HPP
#define __HMM_HPP 

#include <exception>
#include <stdexcept>
#include <vector>

template<typename S, typename O> 
class AbstractHMM {
private:
	S* _states;
	std::size_t _nStates;
	double* _transP;
	double* _initP;

protected:
	AbstractHMM(
		const std::vector<S>&, 
		const std::vector<double>&,
		const std::vector<double>&);

	virtual bool state_inBounds(const std::size_t i) const{
		return (i>0 && i<_nStates);
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

	virtual ~AbstractHMM(){ delete[] _states; delete[] _transP; delete[] _initP; }

	virtual std::size_t find_state(const S state) const{
		for(std::size_t i=0;i<_nStates;++i){
			if(_states[i] == state){
				return i;
			}
		}
		throw std::out_of_range("State not found");
	}

	virtual std::size_t find_observation(const O obs) const = 0;

	virtual std::size_t states() const { return _nStates; };
	virtual std::size_t observations() const = 0;


};

template<typename S, typename O>
AbstractHMM<S, O>::AbstractHMM(	
	const std::vector<S>& states,
	const std::vector<double>& transP,
	const std::vector<double>& initP):
		_states(new S[states.size()]), _nStates(states.size()), 
		_transP(new double[states.size()*states.size()]), 
		_initP(new double[states.size()]){
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
			std::copy(transP.begin(), transP.end(), _transP);
			std::copy(std::begin(initP), std::end(initP), _initP);
		}

template<typename S, typename O>
double AbstractHMM<S, O>::trans_p(std::size_t i, std::size_t j) const{
	if(!state_inBounds(i) || !state_inBounds(j))
		throw std::out_of_range("State out of bounds");
	return _transP[i*_nStates + j];
}

template<typename S, typename O>
double AbstractHMM<S, O>::trans_p_by_object(S firstState, S secondState) const{
	return _transP[find_state(firstState)*_nStates + find_state(secondState)];
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