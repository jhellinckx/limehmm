#include <vector>
#include <string>
#include <utility>

#include "distributions.hpp"
#include "hmm_base.hpp"

RawModel::RawModel() {}

RawModel::RawModel(const RawModel& other) : 
	states_indices(other.states_indices), states_names(other.states_names),
	A(other.A), B(other.B.size()), pi_begin(other.pi_begin), pi_end(other.pi_end),
	is_finite(other.is_finite), silent_states_index(other.silent_states_index),
	alphabet(other.alphabet), free_pi_begin(other.free_pi_begin), 
	free_pi_end(other.free_pi_end), free_transitions(other.free_transitions),
	free_emissions(other.free_emissions) {
		for(std::size_t i = 0; i < other.B.size(); ++i){
			B[i] = (other.B[i] == nullptr) ? nullptr : other.B[i]->clone();
		}
	}

RawModel::RawModel(RawModel&& other) :
	states_indices(std::move(other.states_indices)), 
	states_names(std::move(other.states_names)), A(std::move(other.A)), B(std::move(other.B)), 
	pi_begin(std::move(other.pi_begin)), pi_end(std::move(other.pi_end)), is_finite(other.is_finite), 
	silent_states_index(std::move(other.silent_states_index)), 
	alphabet(std::move(other.alphabet)), free_pi_begin(std::move(other.free_pi_begin)), 
	free_pi_end(std::move(other.free_pi_end)), free_transitions(std::move(other.free_transitions)),
	free_emissions(std::move(other.free_emissions)) {}

RawModel& RawModel::operator=(const RawModel& other) {
	if(this != &other){
		clean();
		states_indices = other.states_indices;
		states_names = other.states_names;
		A = other.A;
		B.resize(other.B.size());
		for(std::size_t i = 0; i < other.B.size(); ++i){
			B[i] = (other.B[i] == nullptr) ? nullptr : other.B[i]->clone();
		}
		pi_begin = other.pi_begin;
		pi_end = other.pi_end;
		is_finite = other.is_finite;
		silent_states_index = other.silent_states_index;
		alphabet = other.alphabet;
		free_pi_begin = other.free_pi_begin;
		free_pi_end = other.free_pi_end;
		free_transitions = other.free_transitions;
		free_emissions = other.free_emissions;
	}
	return *this;
}

RawModel& RawModel::operator=(RawModel&& other){
	if(this != &other){
		states_indices = std::move(other.states_indices);
		states_names = std::move(other.states_names);
		A = std::move(other.A);
		B = std::move(other.B);
		pi_begin = std::move(other.pi_begin);
		pi_end = std::move(other.pi_end);
		is_finite = other.is_finite;
		silent_states_index = other.silent_states_index;
		alphabet = std::move(other.alphabet);
		free_pi_begin = std::move(other.free_pi_begin);
		free_pi_end = std::move(other.free_pi_end);
		free_transitions = std::move(other.free_transitions);
		free_emissions = std::move(other.free_emissions);
	}
	return *this;
}

void RawModel::clean() {
	for(Distribution* dist : B){
		if(dist != nullptr) delete dist;
	}
	states_indices.clear();
	states_names.clear();
	A.clear();
	B.clear();
	pi_begin.clear();
	pi_end.clear();
	is_finite = false;
	silent_states_index = std::size_t();
	alphabet.clear();
	free_pi_begin.clear();
	free_pi_end.clear();
	free_transitions.clear();
	free_emissions.clear();	
}

RawModel::~RawModel() {
	for(Distribution* dist : B){
		if(dist != nullptr) delete dist;
	}
}