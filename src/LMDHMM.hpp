#ifndef __LMDHMM_HPP
#define __LMDHMM_HPP

#include "AbstractLMHMM.hpp"
#include "AbstractDHMM.hpp"

template<typename S, typename O>
class LMDHMM : public AbstractLMHMM<S, O>, public AbstractDHMM<S, O>{

};

#endif