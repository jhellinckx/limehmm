#ifndef __ABSTRACTLMHMM_HPP
#define __ABSTRACTLMHMM_HPP

#include <exception>
#include <stdexcept>

#include "AbstractHMM.hpp"

template<typename S, typename O> 
class AbstractLMHMM : virtual public AbstractHMM<S, O>{

	void viterbi(std::array<O>&);
	void train_viterbi(const std::array<O>&);
	void train_baumWelch(const std::array<O>&);
	void train_stochasticEM(const std::array<O>&);
};


#endif