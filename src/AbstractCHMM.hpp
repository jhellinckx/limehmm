#ifndef __ABSTRACTCHMM_HPP
#define __ABSTRACTCHMM_HPP

#include <exception>
#include <stdexcept>

#include "AbstractHMM.hpp"

template<typename S, typename O>
class AbstractCHMM: public virtual AbstractHMM<S, O>{

};

#endif