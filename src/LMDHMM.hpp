#ifndef __LMDHMM_HPP
#define __LMDHMM_HPP

#include "AbstractHMM.hpp"
#include "AbstractLMHMM.hpp"
#include "AbstractDHMM.hpp"

#include <vector>

template<typename S, typename O>
class LMDHMM : public AbstractLMHMM<S, O>, public AbstractDHMM<S, O>{
public:
	LMDHMM(
		const std::vector<O>& observations,
		const std::vector<double>& emiP,
		const std::vector<S>& states,
		const std::vector<double>& transP,
		const std::vector<double>& initP):
			AbstractHMM<S, O>(states, transP, initP), AbstractLMHMM<S, O>(), 
			AbstractDHMM<S, O>(observations, emiP, states.size()){}
};

#endif