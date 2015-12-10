#ifndef __ABSTRACTLMHMM_HPP
#define __ABSTRACTLMHMM_HPP

#include <exception>
#include <stdexcept>

#include "AbstractHMM.hpp"

template<typename S, typename O> 
class AbstractLMHMM : public virtual AbstractHMM<S, O>{

protected:
	AbstractLMHMM(): AbstractHMM<S, O>({},{{}}) {}

public:
	virtual void viterbi(std::array<O>&);
	virtual void train_viterbi(const std::array<O>&);
	virtual void train_baumWelch(const std::array<O>&);
	virtual void train_stochasticEM(const std::array<O>&);
};


#endif