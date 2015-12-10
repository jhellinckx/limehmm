#ifndef __HMM_HPP
#define __HMM_HPP 

#include <exception>
#include <stdexcept>

template<typename S, typename O> 
class AbstractHMM {
private:
	S* _states;
	std::size_t _nStates;
	float* _transP;
	float* _initP;

	AbstractHMM(const std::initializer_list<S>& states);

protected:
	AbstractHMM(
		const std::initializer_list<S>&, 
		const std::initializer_list<std::initializer_list<float>>&,
		const std::initializer_list<float>&);

	virtual bool state_inBounds(const std::size_t i) const{
		return (i>0 && i<_nStates);
	}

	virtual bool observation_inBounds(const std::size_t) const = 0;

public:
	/* Transition probabilities getters */
	virtual float trans_p(std::size_t, std::size_t) const;
	virtual float trans_p_by_object(S, S) const;

	/* Emission probabilities getters, defined in derived classes */
	virtual float emi_p(std::size_t, std::size_t) const = 0;
	virtual float emi_p_by_object(S, O) const = 0;

	/* Initial state probabilites getters */
	virtual float init_p(std::size_t) const;
	virtual float init_p_by_object(S) const; 

	/* Decoding */
	virtual void viterbi(const std::array<O>&) = 0;

	/* Training procedures */
	virtual void train_viterbi(const std::array<O>&) = 0;
	virtual void train_baumWelch(const std::array<O>&) = 0;
	virtual void train_stochasticEM(const std::array<O>&) = 0;

	virtual ~AbstractHMM(){ delete[] _states; delete[] _transP; delete[] _initP; }

	virtual std::size_t find_state(const S state) const{
		for(std::size_t i=0;i<_nStates;++i){
			if(_states[i] == state){
				return i
			}
		}
		return nullptr;
	}

	virtual std::size_t find_observation(const O obs) const = 0;

	virtual std::size_t states() const { return _nStates };
	virtual std::size_t observations() const = 0;


};

template<typename S, typename O> 
AbstractHMM<S, O>::AbstractHMM(const std::initializer_list<S>& states): 
	_states(new S[states.size()]), _nStates(states.size()){
		for(std::size_t i=0;i<states.size();++i){
			if(find_state(states[i]) != nullptr)
				throw std::logic_error("Two identical states were given in states initializer array");
			_states[i] = states[i];
		}
}

template<typename S, typename O>
AbstractHMM<S, O>::AbstractHMM(	
	const std::initializer_list<S>& states,
	const std::initializer_list<std::initializer_list<float>>& transP,
	const std::initializer_list<float>& initP):
		AbstractHMM(states), _transP(new float[states.size()*states.size()]), 
		_initP(new float[states.size()]){
			if(transP.size() != states.size() || transP[0].size() != states.size())
				throw std::logic_error("Wrong initializer_list size for transitions probabilities");
			if(initP.size() != states.size())
				throw std::logic_error("Wrond initializer_list size for initial probabilities")
			std::copy(std::begin(transP), std::end(transP), _transP);
			std::copy(std::begin(initP), std::end(initP), _initP);
		}

template<typename S, typename O>
float AbstractHMM<S, O>::trans_p(std::size_t i, std::size_t j) const{
	if(!state_inBounds(i) || !state_inBounds(j))
		throw std::out_of_range("State out of bounds");
	return _transp[i*_nStates + j];
}

template<typename S, typename O>
float AbstractHMM<S, O>::trans_p_by_object(S firstState, S secondState) const{
	std::size_t i,j;
	if(i=find_state(firstState) == nullptr || j=find_state(secondState) == nullptr)
		throw std::logic_error("State not found");
	return _transp[i*_nStates + j];
}

template<typename S, typename O>
float AbstractHMM<S, O>::init_p(std::size_t i) const{
	if(!state_inBounds(i))
		throw std::out_of_range("State out of bounds");
	return _initP[i]
}

template<typename S, typename O>
AbstracHMM<S, O>::init_p_by_object(S state) const{
	std::size_t i;
	if(i=find_state(state) == nullptr)
		throw std::logic_error("State not found");
	return _initP[i];
} 

#endif