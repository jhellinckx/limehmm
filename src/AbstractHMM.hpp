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


protected:
	AbstractHMM(const std::initializer_list<S>& states);
	AbstractHMM(
		const std::initializer_list<S>&, 
		const std::initializer_list<std::initializer_list<float>>&);

	bool stateInBounds(const std::size_t& i) const{
		return (i>0 && i<_nStates)
	}

public:

	float transp(std::size_t, std::size_t) const;
	float transp(S, S) const;
	float emip(std::size_t, std::size_t) const;
	float emip(S, O) const;

	/* Decoding */
	virtual void viterbi(const std::array<O>&) = 0;

	/* Training procedures */
	virtual void train_viterbi(const std::array<O>&) = 0;
	virtual void train_baumWelch(const std::array<O>&) = 0;
	virtual void train_stochasticEM(const std::array<O>&) = 0;

	virtual ~AbstractHMM(){ delete[] _states; delete[] _transP; }

	std::size_t find(const S state) const{
		for(std::size_t i=0;i<_nStates;++i){
			if(_states[i] == state){
				return i
			}
		}
		return nullptr;
	}


};

template<typename S, typename O> 
AbstractHMM<S, O>::AbstractHMM(const std::initializer_list<S>& states): 
	_states(new S[states.size()]), _nStates(states.size()){
		for(std::size_t i=0;i<states.size();++i){
			if(find(states[i]) != nullptr)
				throw std::logic_error("Two identical states were given in states initializer array");
			_states[i] = states[i];
		}
}

template<typename S, typename O>
AbstractHMM<S, O>::AbstractHMM(	
	const std::initializer_list<S>& states,
	const std::initializer_list<std::initializer_list<float>>& initTransitions):
		AbstractHMM(states), _transP(new S[states.size()*states.size()]){
			if(initTransitions.size() != states.size() || initTransitions[0].size() != states.size())
				throw std::logic_error("Wrong size of init transitions probabilities");
			std::copy(std::begin(initTransitions), std::end(initTransitions), _transP);
		}

template<typename S, typename O>
float AbstractHMM<S, O>::transp(std::size_t i, std::size_t j){
	if(!stateInBounds(i) || !stateInBounds(j))
		throw std::out_of_range("State out of bounds");
	return _transp[i*j];
}

template<typename S, typename O>
float AbstractHMM<S, O>::transp(S firstState, S secondState){
	std::size_t i,j;
	if(i=find(firstState) == nullptr || j=find(secondState) == nullptr)
		throw logic_error("State not found");
	return _transp[i*j];
}

#endif