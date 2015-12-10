#ifndef __LMCHMM_HPP
#define __LMCHMM_HPP

#include "AbstractLMHMM.hpp"
#include "AbstractCHMM.hpp"

template<typename S, typename O>
class LMCHMM : public AbstractLMHMM<S, O>, public AbstractCHMM<S, O>{


};

#endif