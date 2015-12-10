#ifndef __ABSTRACTDHMM_HPP
#define __ABSTRACTDHMM_HPP

#include <exception>
#include <stdexcept>

#include "AbstractHMM.hpp"

template<typename S, typename O>
class AbstractDHMM: public virtual AbstractHMM<S, O>{
private:
	O* _observations;
	std::size_t _nObservations;
	float* _emiP;

protected:
	AbstractDHMM(
		const std::initializer_list<O>& observations,
		const std::initializer_list<std::initializer_list<float>>& initEmissions,
		std::size_t nStates) :
			AbstractHMM<S, O>({}, {{}}), _observations(new O[observations.size()]), 
			_nObservations(observations.size()), _emiP(new float[nStates*observations.size()]){
				for(size_t i=0;i<observations.size();++i){
					if(find_observation(observations[i]) != nullptr
						throw std::logic_error("Two identical observations were given in observations initializer array");
					_observations[i] = observations[i];
				}
				if(initEmissions.size() != observations.size() || initEmissions[0].size() != observations.size())
					throw std::logic_error("Wrong size of init observation probabilities");
				std::copy(std::begin(initEmissions), std::end(initEmissions), _emiP);
			}

	bool observation_inBounds(const std::size_t i) const{
		return (i>0 && i<_nObservations);
	}

public:
	std::size_t find_observation(const O obs) const{
		for(std::size_t i=0;i<_nObservations;++i){
			if(_observations[i] == obs){
				return i
			}
		}
		return nullptr;
	}
};

#endif