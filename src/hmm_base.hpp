#ifndef __HMM_BASE_HPP
#define __HMM_BASE_HPP

#include "distributions.hpp"
#include <string>
#include <utility>
#include <vector>

typedef std::vector<std::vector<double>> Matrix;

struct RawModel{
	std::map<std::string, std::size_t> states_indices;
	std::vector<std::string> states_names;
	Matrix A;
	std::vector<Distribution*> B;
	std::vector<double> pi_begin;
	std::vector<double> pi_end;
	bool is_finite;
	std::size_t silent_states_index;
	std::vector<std::string> alphabet;
	std::vector<std::size_t> free_pi_begin;
	std::vector<std::size_t> free_pi_end;
	std::vector<std::pair<std::size_t, std::size_t>> free_transitions;
	/* Only discrete ! */
	std::vector<std::pair<std::size_t, std::string>> free_emissions; //TODO : For now, free/fixed parameters PER state, do it for every parameter. 

	RawModel();
	RawModel(const RawModel&);
	RawModel(RawModel&&);
	RawModel& operator=(const RawModel&);
	RawModel& operator=(RawModel&&);
	virtual void clean();
	virtual ~RawModel();
};

#endif