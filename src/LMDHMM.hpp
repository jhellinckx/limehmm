#ifndef __LMDHMM_HPP
#define __LMDHMM_HPP

#include "AbstractLMHMM.hpp"
#include "AbstractDHMM.hpp"

template<typename S, typename O>
class LMDHMM : public AbstractLMHMM<S, O>, public AbstractDHMM<S, O>{
public:
	LMDHMM(
		const std::initializer_list<O>& observations,
		const std::initializer_list<std::initializer_list<float>>& initEmissions,
		const std::initializer_list<S>& states,
		const std::initializer_list<std::initializer_list<float>>& initTransitions):
			AbstractHMM(states, initTransitions), AbstractLMHMM(), 
			AbstractDHMM(observations, initEmissions){}
};

#endif