#ifndef __ABSTRACTDHMM_HPP
#define __ABSTRACTDHMM_HPP

#include <exception>
#include <stdexcept>
#include <vector>

#include "AbstractHMM.hpp"

template<typename S, typename O>
class AbstractDHMM: public virtual AbstractHMM<S, O>{
private:
	std::vector<O> _observations;
	std::vector<double> _emiP;

protected:
	AbstractDHMM(
		const std::vector<O>& observations,
		const std::vector<double>& emiP,
		std::size_t nStates) :
			AbstractHMM<S, O>({}, {}, {}), _observations(std::vector<O>(observations.size())), 
			_emiP(std::vector<double>(nStates*observations.size())){				
				for(size_t i=0;i<observations.size();++i){
					try{
						find_observation(observations[i]);
						throw std::logic_error("Two identical observations were given in observations initializer array");
					} catch(std::out_of_range e){
						_observations[i] = observations[i];
					}
				}
				if(emiP.size() != nStates*observations.size())
					throw std::logic_error("Wrong size of init observation probabilities");
				std::copy(emiP.begin(), emiP.end(), _emiP.begin());
			}

	bool observation_inBounds(const std::size_t i) const{
		return (i>0 && i<observations());
	}

public:

	virtual double emi_p(std::size_t state, std::size_t obs) const{
		if(!this->state_inBounds(state))
			throw std::out_of_range("State index out of bounds");
		if(!observation_inBounds(obs))
			throw std::out_of_range("Observation index out of bounds");
		return _emiP[state*this->states() + obs];
	}
	virtual double emi_p_by_object(S state, O obs) const{
		return _emiP[this->find_state(state)*this->states() + find_observation(obs)];
	}

	std::size_t find_observation(const O obs) const{
		for(std::size_t i=0;i<observations();++i){
			if(_observations[i] == obs){
				return i;
			}
		}
		throw std::out_of_range("Observation not found");
	}

	std::size_t observations() const { return _observations.size(); }

	virtual ~AbstractDHMM(){}
};

#endif