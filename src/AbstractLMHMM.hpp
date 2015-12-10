#ifndef __ABSTRACTLMHMM_HPP
#define __ABSTRACTLMHMM_HPP

#include <exception>
#include <stdexcept>
#include <vector>

#include "AbstractHMM.hpp"

template<typename S, typename O> 
class AbstractLMHMM : public virtual AbstractHMM<S, O>{

protected:
	AbstractLMHMM(): AbstractHMM<S, O>({}, {}, {}) {}

public:
	virtual void viterbi(const std::vector<O>&){}
	virtual void train_viterbi(const std::vector<O>&){}
	virtual void train_baumWelch(const std::vector<O>&){}
	virtual void train_stochasticEM(const std::vector<O>&){}
};


#endif