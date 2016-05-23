#include <sstream>
#include <unordered_map>
#include <exception>
#include <stdexcept>
#include <vector>
#include <string>
#include <algorithm>	
#include <math.h>	// log, exp
#include <utility>	// std::pair
#include <memory> // std::shared_ptr
#include <iomanip> // std::setprecision
#include <tuple>
#include <iostream>
#include "utils.hpp"
#include "constants.hpp"
#include "hmm_algorithms.hpp"
#include "hmm_base.hpp"


/* ===================== BASE CLASSES ===================== */


HMMAlgorithm::HMMAlgorithm(const std::string& name, RawModel* model) : _name(name), _model(model) {}
std::string HMMAlgorithm::name() const { return _name; }
std::string HMMAlgorithm::type() const { return name(); }
void HMMAlgorithm::set_model(RawModel* model) { _model = model; }
HMMAlgorithm::~HMMAlgorithm() {}

ForwardAlgorithm::ForwardAlgorithm(const std::string& name, RawModel* model) : HMMAlgorithm(name, model) {}
ForwardAlgorithm::~ForwardAlgorithm() {}

BackwardAlgorithm::BackwardAlgorithm(const std::string& name, RawModel* model) : HMMAlgorithm(name, model) {}
BackwardAlgorithm::~BackwardAlgorithm() {}

DecodingAlgorithm::DecodingAlgorithm(const std::string& name, RawModel* model) : HMMAlgorithm(name, model) {}
DecodingAlgorithm::~DecodingAlgorithm() {}

TrainingAlgorithm::TrainingAlgorithm(const std::string& name, RawModel* model) : HMMAlgorithm(name, model) {}
TrainingAlgorithm::~TrainingAlgorithm() {}


/* ===================== LINEAR MEMORY FORWARD ===================== */

LinearMemoryForwardAlgorithm::LinearMemoryForwardAlgorithm(RawModel* model) : ForwardAlgorithm(hmm_config::kLinearMemoryForwardAlgorithmName, model) {}
LinearMemoryForwardAlgorithm* LinearMemoryForwardAlgorithm::clone() const { return new LinearMemoryForwardAlgorithm(*this); }
LinearMemoryForwardAlgorithm::~LinearMemoryForwardAlgorithm() {}

std::vector<double> LinearMemoryForwardAlgorithm::forward(const std::vector<std::string>& sequence, std::size_t t_max) {
	if(t_max == 0) t_max = sequence.size();
	if(sequence.size() == 0) throw std::logic_error("forward on empty sequence");
	else{
		std::vector<double> alpha = forward_init(sequence);
		for(std::size_t t = 1; t < std::min(sequence.size(), t_max); ++t) {
			alpha = forward_step(sequence, alpha, t);
		}
			return alpha;
	}
}

std::vector<double> LinearMemoryForwardAlgorithm::forward_init(const std::vector<std::string>& sequence){
	std::vector<double> alpha_0(_model->A.size(), utils::kNegInf);
	/* First iterate over the silent states to compute the probability of
	passing through silent states before emitting the first symbol. */
	for(std::size_t i = _model->silent_states_index; i < _model->A.size(); ++i){
		alpha_0[i] = _model->pi_begin[i];
		for(std::size_t j = _model->silent_states_index; j < i; ++j){
			alpha_0[i] = utils::sum_log_prob(alpha_0[i], _model->A[j][i] + alpha_0[j]);
		}
	}
	/* Fill alpha_0 for non-silent states. To compute alpha_0, we need to 
	sum the probability to directly begin at non-silent state i (pi)
	with the probabilities to transit from all the silent states which
	have a begin probability > 0 to non-silent state i. */
	for(std::size_t i = 0; i < _model->silent_states_index; ++i){
		alpha_0[i] = _model->pi_begin[i];
		for(std::size_t j = _model->silent_states_index; j < _model->A.size(); ++j){
			alpha_0[i] = utils::sum_log_prob(alpha_0[i], _model->A[j][i] + alpha_0[j]);
		}
		
	}
	/* We can now compute alpha_1. */
	std::vector<double> alpha_1(_model->A.size(), utils::kNegInf);
	/* First iterate over non-silent states. */
	for(std::size_t i = 0; i < _model->silent_states_index; ++i){
		alpha_1[i] = alpha_0[i] + (*_model->B[i])[sequence[0]];
	}
	/* Then silent states, in toporder. */
	for(std::size_t i = _model->silent_states_index; i < _model->A.size(); ++i){
		alpha_1[i] = utils::kNegInf;
		for(std::size_t j = 0; j < i; ++j){
			alpha_1[i] = utils::sum_log_prob(alpha_1[i], _model->A[j][i] + alpha_1[j]);
		}
	}
	return alpha_1;
}

std::vector<double> LinearMemoryForwardAlgorithm::forward_step(const std::vector<std::string>& sequence, const std::vector<double>& alpha_prev_t, std::size_t t) {
	std::vector<double> alpha_t(_model->A.size(), utils::kNegInf);
	/* Normal states. */
	for(std::size_t i = 0; i < _model->silent_states_index; ++i){
		alpha_t[i] = utils::kNegInf;
		for(std::size_t j = 0; j < _model->A.size(); ++j){
			alpha_t[i] = utils::sum_log_prob(alpha_t[i], alpha_prev_t[j] + _model->A[j][i]);
		}
		alpha_t[i] = alpha_t[i] + (*_model->B[i])[sequence[t]];
	}
	/* Silent states. */
	for(std::size_t i = _model->silent_states_index; i < _model->A.size(); ++i){
		alpha_t[i] = utils::kNegInf;
		for(std::size_t j = 0; j < i; ++j){
			alpha_t[i] = utils::sum_log_prob(alpha_t[i], alpha_t[j] + _model->A[j][i]);
		}
	}
	return alpha_t;
}

std::pair<std::vector<double>, double> LinearMemoryForwardAlgorithm::forward_terminate(const std::vector<double>& alpha_T){
	double log_prob = utils::kNegInf;
	std::vector<double> alpha_end(_model->A.size(), utils::kNegInf);
	if(_model->is_finite){
		/* Sum all and add end transitions. */
		for(std::size_t i = 0; i < alpha_T.size(); ++i){
			alpha_end[i] = alpha_T[i] + _model->pi_end[i];
			log_prob = utils::sum_log_prob(log_prob, alpha_end[i]);
		}
	}
	else{
		/* Non finite hmm end in non-silent states. */
		for(std::size_t i = 0; i < _model->silent_states_index; ++i){
			alpha_end[i] = alpha_T[i];
			log_prob = utils::sum_log_prob(log_prob, alpha_end[i]);
		}
		for(std::size_t i = _model->silent_states_index; i < _model->A.size(); ++i){
			alpha_end[i] = utils::kNegInf;
		}
	}
	return std::make_pair(alpha_end, log_prob);
}

double LinearMemoryForwardAlgorithm::log_likelihood(const std::vector<std::string>& sequence){
	return forward_terminate(forward(sequence, sequence.size())).second;	
}

double LinearMemoryForwardAlgorithm::log_likelihood(const std::vector<std::vector<std::string>>& sequences){
	double likelihood = 0;
	for(const std::vector<std::string>& sequence : sequences){
		likelihood += log_likelihood(sequence);	
	}
	return likelihood;
}


/* ===================== LINEAR MEMORY BACKWARD ===================== */


LinearMemoryBackwardAlgorithm::LinearMemoryBackwardAlgorithm(RawModel* model) : BackwardAlgorithm(hmm_config::kLinearMemoryBackwardAlgorithmName, model) {}
LinearMemoryBackwardAlgorithm* LinearMemoryBackwardAlgorithm::clone() const { return new LinearMemoryBackwardAlgorithm(*this); }
LinearMemoryBackwardAlgorithm::~LinearMemoryBackwardAlgorithm() {}


std::vector<double> LinearMemoryBackwardAlgorithm::backward(const std::vector<std::string>& sequence, std::size_t t_min) {
	if(t_min > 0) --t_min;
	if(sequence.size() == 0) throw std::runtime_error("backward on empty sequence");
	else{
		std::vector<double> beta = backward_init();
		for(std::size_t t = sequence.size() - 2; t >= t_min && t < sequence.size(); --t){
			beta = backward_step(beta, sequence, t);
		}
		return beta;
	}
}

std::vector<double> LinearMemoryBackwardAlgorithm::backward_init() {
	std::vector<double> beta_T(_model->A.size());
	if(_model->is_finite){
		for(std::size_t i = _model->A.size() - 1; i >= _model->silent_states_index; --i){
			beta_T[i] = _model->pi_end[i];
			for(std::size_t j = _model->A.size() - 1; j > i; --j){
				beta_T[i] = utils::sum_log_prob(beta_T[i], _model->A[i][j] + beta_T[j]);
			}
		}
		for(std::size_t i = 0; i < _model->silent_states_index; ++i){
			beta_T[i] = _model->pi_end[i];
			for(std::size_t j = _model->silent_states_index; j < _model->A.size(); ++j){
				beta_T[i] = utils::sum_log_prob(beta_T[i], _model->A[i][j] + beta_T[j]); 
			}
		}
	}
	else{
		for(std::size_t i = 0; i < _model->silent_states_index; ++i){
			beta_T[i] = 0.0;
		}
		for(std::size_t i = _model->silent_states_index; i < _model->A.size(); ++i){
			beta_T[i] = utils::kNegInf;
		}
	}
	return beta_T;
};

std::vector<double> LinearMemoryBackwardAlgorithm::backward_step(const std::vector<double>& beta_previous_t, const std::vector<std::string>& sequence, std::size_t t) {
	std::vector<double> beta_t(_model->A.size());
	for(std::size_t i = _model->A.size(); i-- > 0;){
		beta_t[i] = utils::kNegInf;
		/* Consider previous step non-silent states. */
		for(std::size_t j = 0; j < _model->silent_states_index; j++){
			beta_t[i] = utils::sum_log_prob(beta_t[i], beta_previous_t[j] + _model->A[i][j] + (*_model->B[j])[sequence[t + 1]]);
		}
		/* Consider current step silent states. 
		If i is a silent state (i.e. i > _silent_state_index), only iterate for each j > i (topological order !). 
		Else if i is a non-silent state, iterate over all the silent states. */
		for(std::size_t j = std::max(i + 1, _model->silent_states_index); j < _model->A.size(); j++){
			beta_t[i] = utils::sum_log_prob(beta_t[i], beta_t[j] + _model->A[i][j]);
		}
	}
	return beta_t;
};

std::tuple<std::vector<double>, std::vector<double>, double> LinearMemoryBackwardAlgorithm::backward_terminate(const std::vector<double>& beta_1, const std::vector<std::string>& sequence){
	std::vector<double> beta_0(_model->A.size());
	for(std::size_t i = _model->A.size() - 1; i >= _model->silent_states_index; --i){
		beta_0[i] = utils::kNegInf;
		/* Consider previous step non-silent states. */
		for(std::size_t j = 0; j < _model->silent_states_index; j++){
			beta_0[i] = utils::sum_log_prob(beta_0[i], beta_1[j] + _model->A[i][j] + (*_model->B[j])[sequence[0]]);
		}
		/* Consider current step silent states. */
		for(std::size_t j = i + 1; j < _model->A.size(); j++){
			beta_0[i] = utils::sum_log_prob(beta_0[i], beta_0[j] + _model->A[i][j]);
		}
	}
	std::vector<double> beta_end(_model->A.size());
	double log_prob = utils::kNegInf;
	for(std::size_t i = 0; i < _model->silent_states_index; ++i){
		beta_end[i] = _model->pi_begin[i] + (*_model->B[i])[sequence[0]] + beta_1[i];
		log_prob = utils::sum_log_prob(log_prob, beta_end[i]);
	}
	for(std::size_t i = _model->silent_states_index; i < _model->A.size(); ++i){
		beta_end[i] = _model->pi_begin[i] + beta_0[i];
		log_prob = utils::sum_log_prob(log_prob, beta_end[i]);
	}
	return std::make_tuple(beta_0, beta_end, log_prob);
}

double LinearMemoryBackwardAlgorithm::log_likelihood(const std::vector<std::string>& sequence){
	return std::get<2>(backward_terminate(backward(sequence, 0), sequence));
}

double LinearMemoryBackwardAlgorithm::log_likelihood(const std::vector<std::vector<std::string>>& sequences){
	double likelihood = 0;
	for(const std::vector<std::string>& sequence : sequences){
		likelihood += log_likelihood(sequence);	
	}
	return likelihood;
}


/* ===================== LINEAR MEMORY VITERBI DECODE ===================== */

/* ------------- TRACEBACK -------------  */

void LinearMemoryViterbiDecodingAlgorithm::Traceback::_init_previous() {
	for(std::size_t i = 0; i < _nodes; ++i){
		_previous_nodes[i] = NodePtr(new Node(i));
	}
}
void LinearMemoryViterbiDecodingAlgorithm::Traceback::_init_current() {
	for(std::size_t i = 0; i < _nodes; ++i){
			_current_nodes[i] = NodePtr(new Node(i));
	}
}

LinearMemoryViterbiDecodingAlgorithm::Traceback::Traceback(std::size_t num_nodes) : 
	_nodes(num_nodes), _previous_nodes(_nodes), _current_nodes(_nodes) {
		_init_previous();
		_init_current();
	}

void LinearMemoryViterbiDecodingAlgorithm::Traceback::add_link(std::size_t previous, std::size_t current, bool link_to_current) {
	_current_nodes[current]->set_previous((link_to_current) ? _current_nodes[previous] : _previous_nodes[previous]);
}

void LinearMemoryViterbiDecodingAlgorithm::Traceback::next_column() {
	_previous_nodes = _current_nodes;
	_init_current();
}

void LinearMemoryViterbiDecodingAlgorithm::Traceback::reset() {
	_init_previous();
	_init_current();
}

std::vector<std::size_t> LinearMemoryViterbiDecodingAlgorithm::Traceback::from(std::size_t k){
	std::vector<std::size_t> traceback;
	NodePtr node_ptr = _previous_nodes[k];
	traceback.push_back(node_ptr->value);
	while(node_ptr->previous){
		node_ptr = node_ptr->previous; 
		traceback.push_back(node_ptr->value);
	}
	std::reverse(traceback.begin(), traceback.end());
	return traceback;
}

std::string LinearMemoryViterbiDecodingAlgorithm::Traceback::to_string() const {
	std::ostringstream oss;
	for(NodePtr p_node : _current_nodes){
		oss << p_node->value << " -> ";
		if(p_node->previous){
			oss << p_node->previous->value;
		}
		else{
			oss << "END";
		}
		oss << " / ";
	}
	return oss.str();
}

/* ------------- DECODE -------------  */

LinearMemoryViterbiDecodingAlgorithm::LinearMemoryViterbiDecodingAlgorithm(RawModel* model) : DecodingAlgorithm(hmm_config::kLinearMemoryViterbiDecodeAlgorithmName, model) {}
LinearMemoryViterbiDecodingAlgorithm* LinearMemoryViterbiDecodingAlgorithm::clone() const { return new LinearMemoryViterbiDecodingAlgorithm(*this); }
LinearMemoryViterbiDecodingAlgorithm::~LinearMemoryViterbiDecodingAlgorithm() {}

std::vector<double> LinearMemoryViterbiDecodingAlgorithm::viterbi_init(Traceback& psi, const std::vector<std::string>& sequence) {
	std::vector<double> phi_0(_model->A.size(), utils::kNegInf);
	/* First iterate over the silent states to compute the max probability of
	passing through silent states before emitting the first symbol. */
	double max_phi;
	double current_phi;
	std::size_t max_psi;
	for(std::size_t i = _model->silent_states_index; i < _model->A.size(); ++i){
		max_phi = _model->pi_begin[i];
		max_psi = _model->A.size();
		for(std::size_t j = _model->silent_states_index; j < i; ++j){
			current_phi = _model->A[j][i] + phi_0[j];
			if(current_phi > max_phi){
				max_phi = current_phi;
				max_psi = j;
			}
		}
		if(max_phi != utils::kNegInf){
			phi_0[i] = max_phi;
		}
		if(max_psi < _model->A.size()){
			psi.add_link(max_psi, i, true);	
		}
	}
	psi.next_column();
	std::vector<double> phi_1(_model->A.size(), utils::kNegInf);
	/* Fill phi_1 for non-silent states. */
	for(std::size_t i = 0; i < _model->silent_states_index; ++i){
		max_phi = _model->pi_begin[i];
		max_psi = _model->A.size();
		for(std::size_t j = _model->silent_states_index; j < _model->A.size(); ++j){
			current_phi = _model->A[j][i] + phi_0[j];
			if(current_phi > max_phi){
				max_phi = current_phi;
				max_psi = j;
			}
		}
		if(max_phi != utils::kNegInf){
			phi_1[i] = max_phi + (*_model->B[i])[sequence[0]];
		}
		if(max_psi < _model->A.size()){
			psi.add_link(max_psi, i);
		}
	}
	/* Then silent states, in toporder. */
	for(std::size_t i = _model->silent_states_index; i < _model->A.size(); ++i){
		max_phi = utils::kNegInf;
		max_psi = _model->A.size();
		for(std::size_t j = 0; j < i; ++j){
			current_phi = _model->A[j][i] + phi_1[j];
			if(current_phi > max_phi){
				max_phi = current_phi;
				max_psi = j;
			}
		}
		if(max_phi != utils::kNegInf && max_psi < _model->A.size()){
			phi_1[i] = max_phi;
			psi.add_link(max_psi, i, true);
		}
	}
	psi.next_column();
	return phi_1;
}

std::vector<double> LinearMemoryViterbiDecodingAlgorithm::viterbi_step(const std::vector<double>& phi_prev_t, Traceback& psi, std::size_t t, const std::vector<std::string>& sequence) {
	std::vector<double> phi_t(_model->A.size(), utils::kNegInf);
		double max_phi;
		double current_phi;
		std::size_t max_psi;
		/* Normal states. */
		for(std::size_t i = 0; i < _model->silent_states_index; ++i){
			max_phi = utils::kNegInf;
			max_psi = _model->A.size();
			for(std::size_t j = 0; j < _model->A.size(); ++j){
				current_phi = phi_prev_t[j] + _model->A[j][i];
				if(current_phi > max_phi){
					max_phi = current_phi;
					max_psi = j;
				}
			}
			if(max_phi != utils::kNegInf && max_psi != _model->A.size()){
				phi_t[i] = max_phi + (*_model->B[i])[sequence[t]];
				psi.add_link(max_psi, i);
			}
		}
		/* Silent states. */
		for(std::size_t i = _model->silent_states_index; i < _model->A.size(); ++i){
			max_phi = utils::kNegInf;
			max_psi = _model->A.size();
			for(std::size_t j = 0; j < i; ++j){
				current_phi = phi_t[j] + _model->A[j][i];
				if(current_phi > max_phi){
					max_phi = current_phi;
					max_psi = j;
				}
			}
			if(max_phi != utils::kNegInf && max_psi != _model->A.size()){
				phi_t[i] = max_phi;
				psi.add_link(max_psi, i, true);
			}
		}
		psi.next_column();
		return phi_t;
}

std::size_t LinearMemoryViterbiDecodingAlgorithm::viterbi_terminate(std::vector<double>& phi_T){
	double max_phi_T = utils::kNegInf;
	std::size_t max_state_index = _model->A.size();
	if(_model->is_finite){
		/* Add end transitions.*/
		for(std::size_t i = 0; i < _model->A.size(); ++i){
			phi_T[i] = phi_T[i] + _model->pi_end[i];
			if(phi_T[i] > max_phi_T){
				max_phi_T = phi_T[i];
				max_state_index = i;
			}
		}
	}
	else{
		/* Only consider normal states. */
		for(std::size_t i = 0; i < _model->silent_states_index; ++i){
			if(phi_T[i] > max_phi_T){
				max_phi_T = phi_T[i];
				max_state_index = i;
			}
		}
	}
	return max_state_index;
}

std::pair<std::vector<std::string>, double> LinearMemoryViterbiDecodingAlgorithm::decode(const std::vector<std::string>& sequence, std::size_t t_max) {
	if(t_max == 0) t_max = sequence.size();
	if(sequence.size() == 0) throw std::logic_error("viterbi on empty sequence");
	else{
		Traceback psi(_model->A.size());
		std::vector<double> phi = viterbi_init(psi, sequence);
		for(std::size_t t = 1; t < std::min(sequence.size(), t_max); ++t) {
			phi = viterbi_step(phi, psi, t, sequence);
		}
		std::size_t max_state_index = viterbi_terminate(phi);
		double max_phi_T = phi[max_state_index];
		if(max_phi_T != utils::kNegInf && max_state_index < _model->A.size()){
			std::vector<std::size_t> path_indices = psi.from(max_state_index);
			std::vector<std::string> path;
			path.reserve(path_indices.size());
			for(std::size_t path_index : path_indices){
				path.push_back(_model->states_names[path_index]);
			}
			return std::make_pair(path, max_phi_T);
		}
		else{
			/* Sequence is impossible. */
			return std::make_pair(std::vector<std::string>(), utils::kNegInf);
		}
	}
}

/* ===================== LINEAR MEMORY TRAINING ===================== */

LinearMemoryTrainingAlgorithm::LinearMemoryTrainingAlgorithm(const std::string& name, RawModel* model) : TrainingAlgorithm(name, model) {}
LinearMemoryTrainingAlgorithm::~LinearMemoryTrainingAlgorithm() {}

unsigned int LinearMemoryTrainingAlgorithm::delta(std::size_t i, std::size_t j){
	return (unsigned int)(i == j);
}

unsigned int LinearMemoryTrainingAlgorithm::delta(std::string i, std::string j){
	return (unsigned int)(i == j);
}

unsigned int LinearMemoryTrainingAlgorithm::any_of_transitions(const std::vector<std::size_t>& traceback, std::size_t i, std::size_t j){
	unsigned int delta_sum = 0;
	for(std::size_t l = 0; l < traceback.size() - 1; ++l){
		delta_sum += delta(traceback[l], i) * delta(traceback[l+1], j);
	}
	return delta_sum;
}

double LinearMemoryTrainingAlgorithm::log_score(std::string first_symbol, std::string second_symbol) {
	return (first_symbol == second_symbol) ? 0 : utils::kNegInf;
}

double LinearMemoryTrainingAlgorithm::log_delta(std::size_t i, std::size_t j) {
	return (i == j) ? 0 : utils::kNegInf;
}

void LinearMemoryTrainingAlgorithm::print_scores(const TransitionScore& score, std::string from_str, bool from_all, std::size_t from, bool log_prob){
	if(from_all){
		for(std::size_t m = 0; m < _model->A.size(); ++m){
			std::cout << score.to_string(m, _model->states_names, from_str, log_prob) << std::endl;
		}
	}
	else{
		std::cout << score.to_string(from, _model->states_names, from_str, log_prob) << std::endl;
	}
}

void LinearMemoryTrainingAlgorithm::print_scores(const EmissionScore& score, std::string from_str, bool from_all, std::size_t from, bool log_prob){
	if(from_all){
		for(std::size_t m = 0; m < _model->A.size(); ++m){
			std::cout << score.to_string(m, _model->states_names, from_str, log_prob) << std::endl;
		}
	}
	else{
		std::cout << score.to_string(from, _model->states_names, from_str, log_prob) << std::endl;
	}
}

void LinearMemoryTrainingAlgorithm::print_total_scores(const TransitionScore& score, bool log_prob){
	print_scores(score, "total", false, 0, log_prob);
}

void LinearMemoryTrainingAlgorithm::print_all_scores(const EmissionScore& score, bool log_prob){
	print_scores(score, "", true, 0, log_prob);
}

/* ------------- TRANSITION SCORE -------------  */

LinearMemoryTrainingAlgorithm::TransitionScore::TransitionScore(const std::vector<std::pair<std::size_t, std::size_t>>& free_transitions,
	const std::vector<std::size_t>& free_pi_begin,
	const std::vector<std::size_t>& free_pi_end,
	std::size_t num_states,
	double default_score) :
		_transitions_scores(num_states, std::vector<double>(free_transitions.size(), default_score)),
		_pi_begin_scores(num_states, std::vector<double>(free_pi_begin.size(), default_score)), 
		_pi_end_scores(num_states, std::vector<double>(free_pi_end.size(), default_score)),
		_free_transitions(&free_transitions),
		_free_pi_begin(&free_pi_begin),
		_free_pi_end(&free_pi_end),
		_default_score(default_score) {}

LinearMemoryTrainingAlgorithm::TransitionScore& LinearMemoryTrainingAlgorithm::TransitionScore::operator=(const TransitionScore& other) {
	if(this != &other){
		for(std::size_t m = 0; m < _transitions_scores.size(); ++m){
			for(std::size_t id = 0; id < _transitions_scores[m].size(); ++id){
				_transitions_scores[m][id] = other._transitions_scores[m][id];
			}
			for(std::size_t id = 0; id < _pi_begin_scores[m].size(); ++id){
				_pi_begin_scores[m][id] = other._pi_begin_scores[m][id];
			}
			for(std::size_t id = 0; id < _pi_end_scores[m].size(); ++id){
				_pi_end_scores[m][id] = other._pi_end_scores[m][id];
			}
		}	
	}
	return *this;
}

void LinearMemoryTrainingAlgorithm::TransitionScore::add(const TransitionScore& other, std::size_t m, std::size_t l){
	for(std::size_t id = 0; id < _transitions_scores[m].size(); ++id){
		_transitions_scores[m][id] += other._transitions_scores[l][id];
	}
	for(std::size_t id = 0; id < _pi_begin_scores[m].size(); ++id){
		_pi_begin_scores[m][id] += other._pi_begin_scores[l][id];
	}
	for(std::size_t id = 0; id < _pi_end_scores[m].size(); ++id){
		_pi_end_scores[m][id] += other._pi_end_scores[l][id];
	}
}

double LinearMemoryTrainingAlgorithm::TransitionScore::score(std::size_t m, std::size_t free_transition_id) const {
	return _transitions_scores[m][free_transition_id];
}

double LinearMemoryTrainingAlgorithm::TransitionScore::score_begin(std::size_t m, std::size_t free_transition_id) const {
	return _pi_begin_scores[m][free_transition_id];
}

double LinearMemoryTrainingAlgorithm::TransitionScore::score_end(std::size_t m, std::size_t end_transition_id) const {
	return _pi_end_scores[m][end_transition_id];
}

std::size_t LinearMemoryTrainingAlgorithm::TransitionScore::num_free_transitions() const { return _free_transitions->size(); }
std::size_t LinearMemoryTrainingAlgorithm::TransitionScore::num_free_begin_transitions() const { return _free_pi_begin->size(); }
std::size_t LinearMemoryTrainingAlgorithm::TransitionScore::num_free_end_transitions() const { return _free_pi_end->size(); }

void LinearMemoryTrainingAlgorithm::TransitionScore::set_begin_score(std::size_t m, std::size_t free_begin_transition_id, double score){
	_pi_begin_scores[m][free_begin_transition_id] = score;
}

void LinearMemoryTrainingAlgorithm::TransitionScore::set_score(std::size_t m, std::size_t free_transition_id, double score){
	_transitions_scores[m][free_transition_id] = score;
}
void LinearMemoryTrainingAlgorithm::TransitionScore::set_end_score(std::size_t m, std::size_t free_end_transition_id, double score){
	_pi_end_scores[m][free_end_transition_id] = score;
}

void LinearMemoryTrainingAlgorithm::TransitionScore::copy_begin(const TransitionScore& other, std::size_t l, std::size_t m){
	for(std::size_t begin_transition_id = 0; begin_transition_id < _pi_begin_scores[m].size(); ++begin_transition_id){
		_pi_begin_scores[m][begin_transition_id] = other.score_begin(l, begin_transition_id); 
	}
}

std::size_t LinearMemoryTrainingAlgorithm::TransitionScore::get_from_state_id(std::size_t free_transition_id) const {
	return (*_free_transitions)[free_transition_id].first;
}

std::size_t LinearMemoryTrainingAlgorithm::TransitionScore::get_to_state_id(std::size_t free_transition_id) const {
	return (*_free_transitions)[free_transition_id].second;	
}

std::size_t LinearMemoryTrainingAlgorithm::TransitionScore::get_state_id_from_begin(std::size_t free_begin_transition_id) const {
	return (*_free_pi_begin)[free_begin_transition_id];
}

std::size_t LinearMemoryTrainingAlgorithm::TransitionScore::get_state_id_to_end(std::size_t free_end_transition_id) const {
	return (*_free_pi_end)[free_end_transition_id];
}

void LinearMemoryTrainingAlgorithm::TransitionScore::reset(double reset_score){
	for(std::size_t m = 0; m < _transitions_scores.size(); ++m){
		for(std::size_t id = 0; id < _transitions_scores[m].size(); ++id){
			_transitions_scores[m][id] = reset_score;
		}
		for(std::size_t id = 0; id < _pi_begin_scores[m].size(); ++id){
			_pi_begin_scores[m][id] = reset_score;
		}
		for(std::size_t id = 0; id < _pi_end_scores[m].size(); ++id){
			_pi_end_scores[m][id] = reset_score;
		}
	}
}

void LinearMemoryTrainingAlgorithm::TransitionScore::reset(){
	reset(_default_score);
}

std::string LinearMemoryTrainingAlgorithm::TransitionScore::to_string(std::size_t m, const std::vector<std::string>& names, const std::string& from, bool log_prob) const {
	std::ostringstream oss; double score;
	std::string name = from.empty() ? names[m] : from;
	oss << "From state " << name << std::endl;
	oss << "Begin scores : " << std::endl;
	for(std::size_t begin_transition_id = 0; begin_transition_id < _pi_begin_scores[m].size(); ++begin_transition_id){
		score = log_prob ? _pi_begin_scores[m][begin_transition_id] : exp(_pi_begin_scores[m][begin_transition_id]);
		oss << "(" << names[(*_free_pi_begin)[begin_transition_id]] << " = " << score << ") "; 
	}
	oss << std::endl << "Mid scores : " << std::endl;
	for(std::size_t transition_id = 0; transition_id < _transitions_scores[m].size(); ++transition_id){
		score = log_prob ? _transitions_scores[m][transition_id] : exp(_transitions_scores[m][transition_id]);
		oss << "(" << names[(*_free_transitions)[transition_id].first] << "->" << 
		names[(*_free_transitions)[transition_id].second] << " = " << score << ") "; 
	}
	oss << std::endl << "End scores : " << std::endl;
	for(std::size_t end_transition_id = 0; end_transition_id < _pi_end_scores[m].size(); ++end_transition_id){
		score = log_prob ?  _pi_end_scores[m][end_transition_id] : exp(_pi_end_scores[m][end_transition_id]);
		oss << "(" << names[(*_free_pi_end)[end_transition_id]] << " = " << score << ") "; 
	}
	oss << std::endl;
	return oss.str();
}

/* ------------- EMISSION SCORE -------------  */

LinearMemoryTrainingAlgorithm::EmissionScore::EmissionScore(
	const std::vector<std::pair<std::size_t, std::string>>& free_emissions, 
	std::size_t num_states, double default_score) :
		_emissions_scores(num_states, std::vector<double>(free_emissions.size(), default_score)),
		_free_emissions(&free_emissions),
		_default_score(default_score) {}

LinearMemoryTrainingAlgorithm::EmissionScore& LinearMemoryTrainingAlgorithm::EmissionScore::operator=(const EmissionScore& other) {
	if(this != &other){
		for(std::size_t m = 0; m < _emissions_scores.size(); ++m){
			for(std::size_t id = 0; id < _emissions_scores[m].size(); ++id){
				_emissions_scores[m][id] = other._emissions_scores[m][id];
			}
		}	
	}
	return *this;
}

std::size_t LinearMemoryTrainingAlgorithm::EmissionScore::get_state_id(std::size_t free_emission_id) const {
	return (*_free_emissions)[free_emission_id].first;
}

std::string LinearMemoryTrainingAlgorithm::EmissionScore::get_symbol(std::size_t free_emission_id) const {
	return (*_free_emissions)[free_emission_id].second;
}

double LinearMemoryTrainingAlgorithm::EmissionScore::score(std::size_t m, std::size_t free_emission_id) const {
	return _emissions_scores[m][free_emission_id];
}

void LinearMemoryTrainingAlgorithm::EmissionScore::set_score(std::size_t m, std::size_t free_emission_id, double score){
	_emissions_scores[m][free_emission_id] = score;
}

std::size_t LinearMemoryTrainingAlgorithm::EmissionScore::num_free_emissions() const {
	return _free_emissions->size();
}

void LinearMemoryTrainingAlgorithm::EmissionScore::add(const EmissionScore& other, std::size_t m, std::size_t l){
	for(std::size_t id = 0; id < _emissions_scores[m].size(); ++id){
		_emissions_scores[m][id] += other._emissions_scores[l][id];
	}	
}

void LinearMemoryTrainingAlgorithm::EmissionScore::reset(double reset_score){
	for(std::size_t m = 0; m < _emissions_scores.size(); ++m){
		for(std::size_t id = 0; id < _emissions_scores[m].size(); ++id){
			_emissions_scores[m][id] = reset_score;
		}
	}
}

void LinearMemoryTrainingAlgorithm::EmissionScore::reset(){
	reset(_default_score);
}

std::string LinearMemoryTrainingAlgorithm::EmissionScore::to_string(std::size_t m, const std::vector<std::string>& names, const std::string& from, bool log_prob) const {
	std::ostringstream oss; double score;
	std::string name = from.empty() ? names[m] : from;
	oss << "From state " << name << std::endl;
	oss << "Emissions scores : " << std::endl;
	for(std::size_t emission_id = 0; emission_id < _emissions_scores[m].size(); ++emission_id){
		score = log_prob ? _emissions_scores[m][emission_id] : exp(_emissions_scores[m][emission_id]);
		oss << "(" << names[(*_free_emissions)[emission_id].first] << "->" << (*_free_emissions)[emission_id].second << " = " << score << ") ";
	}
	oss << std::endl;
	return oss.str();
}

/* ------------- COUNTS HELPERS -------------  */

std::size_t LinearMemoryTrainingAlgorithm::last_non_silent_state(const std::vector<std::size_t>& traceback){
	for(std::size_t i = traceback.size(); i-- > 0;){
		if(traceback[i] < _model->silent_states_index){ 
			return traceback[i];
		}
	}
	return _model->A.size(); //Not found sentinel value. Should never happen though.
}

void LinearMemoryTrainingAlgorithm::update_emissions(const EmissionScore& previous_counts, EmissionScore& current_counts, const std::vector<std::size_t>& traceback, std::string symbol){
	if(!traceback.empty()){
		std::size_t l = traceback[0]; std::size_t m = traceback[traceback.size() - 1];
		std::size_t transmitter = last_non_silent_state(traceback);
		if(transmitter == _model->A.size()) { return; } // This should not happen. 
		std::size_t i; std::string gamma;
		for(std::size_t free_emission_id = 0; free_emission_id < current_counts.num_free_emissions(); ++free_emission_id){
			i = current_counts.get_state_id(free_emission_id);
			gamma = current_counts.get_symbol(free_emission_id);
			current_counts.set_score(m, free_emission_id, previous_counts.score(l, free_emission_id) + delta(transmitter, i) * delta(gamma, symbol));
		}	
	}
}

void LinearMemoryTrainingAlgorithm::update(const TransitionScore& previous_counts, TransitionScore& current_counts, const std::vector<std::size_t>& traceback){
	if(traceback.size() >= 2) {
		std::size_t l = traceback[0]; std::size_t m = traceback[traceback.size() - 1];
		std::size_t i,j;
		current_counts.copy_begin(previous_counts, l, m);
		for(std::size_t free_transition_id = 0; free_transition_id < current_counts.num_free_transitions(); ++free_transition_id){
			i = current_counts.get_from_state_id(free_transition_id);
			j = current_counts.get_to_state_id(free_transition_id);
			current_counts.set_score(m, free_transition_id, previous_counts.score(l, free_transition_id) + any_of_transitions(traceback, i, j));
		}	
	}
}

	

void LinearMemoryTrainingAlgorithm::update_begin(TransitionScore& counts, const std::vector<std::size_t>& traceback){
	if(!traceback.empty()){
		std::size_t l = traceback[0]; std::size_t m = traceback[traceback.size() - 1];
		std::size_t i,j;
		for(std::size_t begin_transition_id = 0; begin_transition_id < counts.num_free_begin_transitions(); ++begin_transition_id){
			j = counts.get_state_id_from_begin(begin_transition_id);
			counts.set_begin_score(m, begin_transition_id, delta(l, j));
		}
		if(traceback.size() >= 2){
			for(std::size_t free_transition_id = 0; free_transition_id < counts.num_free_transitions(); ++free_transition_id){
				i = counts.get_from_state_id(free_transition_id);
				j = counts.get_to_state_id(free_transition_id);
				counts.set_score(m, free_transition_id, any_of_transitions(traceback, i, j));
			}
		}
	}
}

void LinearMemoryTrainingAlgorithm::update_end(TransitionScore& counts, std::size_t m){
	for(std::size_t end_transition_id = 0; end_transition_id < counts.num_free_end_transitions(); ++end_transition_id){
		if(counts.get_state_id_to_end(end_transition_id) == m){
			counts.set_end_score(m, end_transition_id, counts.score_end(m, end_transition_id) + 1.0);
		}
	}
}

/* ===================== LINEAR MEMORY VITERBI TRAINING ===================== */

LinearMemoryViterbiTraining::LinearMemoryViterbiTraining(RawModel* model) : 
	LinearMemoryTrainingAlgorithm(hmm_config::kLinearMemoryViterbiTrainingAlgorithmName, model), 
	_decoding_algorithm(model), _forward_algorithm(model) {}
LinearMemoryViterbiTraining* LinearMemoryViterbiTraining::clone() const { return new LinearMemoryViterbiTraining(*this); }
void LinearMemoryViterbiTraining::set_model(RawModel* model) { 
	_model = model; 
	_decoding_algorithm.set_model(model); 
	_forward_algorithm.set_model(model);
}

LinearMemoryViterbiTraining::~LinearMemoryViterbiTraining() {}

double LinearMemoryViterbiTraining::train(const std::vector<std::vector<std::string>>& sequences, 
	double transition_pseudocount, double convergence_threshold, unsigned int min_iterations, unsigned int max_iterations){

	/* This holds all the counts for the batch of sequences. */
	TransitionScore total_transition_count(_model->free_transitions, _model->free_pi_begin, _model->free_pi_end, 1);
	EmissionScore total_emission_count(_model->free_emissions, 1);
	/* This hold the counts for each sequence. */
	TransitionScore previous_transition_count(_model->free_transitions, _model->free_pi_begin, _model->free_pi_end, _model->A.size());
	TransitionScore current_transition_count(_model->free_transitions, _model->free_pi_begin, _model->free_pi_end, _model->A.size());
	EmissionScore previous_emission_count(_model->free_emissions, _model->A.size());
	EmissionScore current_emission_count(_model->free_emissions, _model->A.size());
	unsigned int iteration = 0;
	/* Use likelihood to determine convergence. */
	double delta = utils::kInf;
	double initial_likelihood = _forward_algorithm.log_likelihood(sequences);
	double previous_likelihood = initial_likelihood;
	double current_likelihood;
	while((iteration < min_iterations || delta > convergence_threshold) 
		&& iteration < max_iterations) {
		/* Iterate over each sequence and compute the counts. */
		for(const std::vector<std::string>& sequence : sequences){
			/* If sequence is empty, go to next sequence. */
			if(sequence.size() == 0) { continue; }
			LinearMemoryViterbiDecodingAlgorithm::Traceback psi(_model->A.size());
			/* The initial step is a special case, since we use initial transition probabilities which
			are not stored in the raw A matrix. */
			std::vector<double> phi = _decoding_algorithm.viterbi_init(psi, sequence);
			/* First iterate only on normal states since emission count are only needed for such states. */
			for(std::size_t m = 0; m < _model->A.size(); ++m){
				std::vector<std::size_t> traceback_m = psi.from(m);
				update_begin(current_transition_count, traceback_m);
				update_emissions(previous_emission_count, current_emission_count, traceback_m, sequence[0]);
			}
			previous_transition_count = current_transition_count;
			previous_emission_count = current_emission_count;
			/* Resetting the traceback since we only need the traceback of current viterbi step. */
			psi.reset();
			/* Main loop for current sequence. */
			for(std::size_t k = 1; k < sequence.size(); ++k){
				phi = _decoding_algorithm.viterbi_step(phi, psi, k, sequence);
				for(std::size_t m = 0; m < _model->A.size(); ++m){
					std::vector<std::size_t> traceback_m = psi.from(m);
					update(previous_transition_count, current_transition_count, traceback_m);
					update_emissions(previous_emission_count, current_emission_count, traceback_m, sequence[k]);
				}
				psi.reset();
				previous_transition_count = current_transition_count;
				previous_emission_count = current_emission_count;
			}
			std::size_t max_state_index = _decoding_algorithm.viterbi_terminate(phi);
			/* Test wether the sequence is possible. */
			if(max_state_index < _model->A.size()){
				/* Add 1 to the end transition count of the max state index if model has end state. */
				if(_model->is_finite){
					update_end(current_transition_count, max_state_index);
				}
				/* Update the total counts. */
				total_transition_count.add(current_transition_count, 0, max_state_index);
				total_emission_count.add(current_emission_count, 0, max_state_index);
			}
			/* Reset counts. */
			current_transition_count.reset();
			previous_transition_count.reset();
			current_emission_count.reset();
			previous_emission_count.reset();
		}
		update_model_from_scores(total_transition_count, total_emission_count, transition_pseudocount);
		total_transition_count.reset();
		total_emission_count.reset();
		current_likelihood = _forward_algorithm.log_likelihood(sequences);
		delta = current_likelihood - previous_likelihood;
		previous_likelihood = current_likelihood;
		++iteration;
	}
	/* Return total improvement. */
	return current_likelihood - initial_likelihood;
}

void LinearMemoryViterbiTraining::update_model_from_scores(const TransitionScore& transitions_scores, 
	const EmissionScore& emissions_scores, double transition_pseudocount){
		update_model_transitions_from_scores(transitions_scores, transition_pseudocount);
		update_model_emissions_from_scores(emissions_scores);
}

void LinearMemoryViterbiTraining::update_model_transitions_from_scores(const TransitionScore& transitions_counts, double transition_pseudocount){
	/* Update begin transitions. */
	/* First, sum all the begin transitions counts. */
	double begin_transitions_count = 0;
	for(std::size_t begin_transition_id = 0; begin_transition_id < _model->free_pi_begin.size(); ++begin_transition_id){
		begin_transitions_count += transitions_counts.score_begin(0, begin_transition_id) + transition_pseudocount;
	}
	/* Then, normalize the count of each begin transition by using the total count. */
	std::size_t state_id;
	for(std::size_t begin_transition_id = 0; begin_transition_id < _model->free_pi_begin.size(); ++begin_transition_id){
		state_id = _model->free_pi_begin[begin_transition_id];
		if(begin_transitions_count > 0){
			_model->pi_begin[state_id] = log((transitions_counts.score_begin(0, begin_transition_id) + transition_pseudocount) / begin_transitions_count);
		}
	}

	/* Update other transitions (don't forget to take end transitions into account). */
	/* Similarly to begin transitions, sum, for each state, all its out transitions counts. */
	std::unordered_map<std::size_t, double> out_transitions_counts;
	for(std::size_t transition_id = 0; transition_id < _model->free_transitions.size(); ++transition_id){
		state_id = _model->free_transitions[transition_id].first;
		/* If state i not yet in map, init count for state i at 0. */
		if(out_transitions_counts.find(state_id) == out_transitions_counts.end()){
			out_transitions_counts[state_id] = 0;
		}
		out_transitions_counts[state_id] += transitions_counts.score(0, transition_id) + transition_pseudocount;
	}
	/* Also add end transitions counts to the the sum. */
	for(std::size_t end_transition_id = 0; end_transition_id < _model->free_pi_end.size(); ++end_transition_id){
		state_id = _model->free_pi_end[end_transition_id];
		if(out_transitions_counts.find(state_id) == out_transitions_counts.end()){
			out_transitions_counts[state_id] = 0;
		}
		out_transitions_counts[state_id] += transitions_counts.score_end(0, end_transition_id) + transition_pseudocount;
	}
	/* Normalize each transition by using the sum. */
	std::size_t from_state, to_state;
	for(std::size_t transition_id = 0; transition_id < _model->free_transitions.size(); ++transition_id){
		from_state = _model->free_transitions[transition_id].first; to_state = _model->free_transitions[transition_id].second;
		if(out_transitions_counts[from_state] > 0){
			_model->A[from_state][to_state] =  log((transitions_counts.score(0, transition_id) + transition_pseudocount) / out_transitions_counts[from_state]);
		}
	}
	/* Don't forget to update the end transitions ! */
	for(std::size_t end_transition_id = 0; end_transition_id < _model->free_pi_end.size(); ++end_transition_id){
		state_id = _model->free_pi_end[end_transition_id];
		if(out_transitions_counts[state_id] > 0){
			_model->pi_end[state_id] = log((transitions_counts.score_end(0, end_transition_id) + transition_pseudocount) / out_transitions_counts[state_id]);
		}
	}
}

void LinearMemoryViterbiTraining::update_model_emissions_from_scores(const EmissionScore& emissions_counts){
	std::unordered_map<std::size_t, double> all_emissions_counts;
	std::size_t state_id;
	for(std::size_t emission_id = 0; emission_id < _model->free_emissions.size(); ++emission_id){
		state_id = _model->free_emissions[emission_id].first;
		if(all_emissions_counts.find(state_id) == all_emissions_counts.end()){
			all_emissions_counts[state_id] = 0;
		}
		all_emissions_counts[state_id] += emissions_counts.score(0, emission_id);
	}
	std::string symbol;
	for(std::size_t emission_id = 0; emission_id < _model->free_emissions.size(); ++emission_id){
		state_id = _model->free_emissions[emission_id].first;
		symbol = _model->free_emissions[emission_id].second;
		if(all_emissions_counts[state_id] > 0) {
			(*(_model->B[state_id]))[symbol] = log(emissions_counts.score(0, emission_id) / all_emissions_counts[state_id]);
		}
	}
}

/* ===================== LINEAR MEMORY BAUM WELCH TRAINING ===================== */

LinearMemoryBaumWelchTraining::LinearMemoryBaumWelchTraining(RawModel* model) : 
	LinearMemoryTrainingAlgorithm(hmm_config::kLinearMemoryBaumWelchTrainingAlgorithmName, model),
	_backward_algorithm(LinearMemoryBackwardAlgorithm(model)) {}

LinearMemoryBaumWelchTraining* LinearMemoryBaumWelchTraining::clone() const { return new LinearMemoryBaumWelchTraining(*this); }
void LinearMemoryBaumWelchTraining::set_model(RawModel* model) { _model = model; _backward_algorithm.set_model(model); }
LinearMemoryBaumWelchTraining::~LinearMemoryBaumWelchTraining() {}

void LinearMemoryBaumWelchTraining::update_model_from_log_scores(const TransitionScore& transitions_scores, 
	const EmissionScore& emissions_scores){
		update_model_transitions_from_log_scores(transitions_scores);
		update_model_emissions_from_log_scores(emissions_scores);
}

void LinearMemoryBaumWelchTraining::update_model_transitions_from_log_scores(const TransitionScore& transitions_scores){
	/* Update begin transitions. */
	/* First, sum all the begin transitions scores. */
	double begin_transitions_score = utils::kNegInf;
	for(std::size_t begin_transition_id = 0; begin_transition_id < _model->free_pi_begin.size(); ++begin_transition_id){
		begin_transitions_score = utils::sum_log_prob(begin_transitions_score, transitions_scores.score_begin(0, begin_transition_id));
	}
	/* Then, normalize the score of each begin transition by using the total score. */
	std::size_t state_id;
	for(std::size_t begin_transition_id = 0; begin_transition_id < _model->free_pi_begin.size(); ++begin_transition_id){
		state_id = _model->free_pi_begin[begin_transition_id];
		if(begin_transitions_score != utils::kNegInf){
			_model->pi_begin[state_id] = transitions_scores.score_begin(0, begin_transition_id) - begin_transitions_score;
		}
	}

	/* Update other transitions (don't forget to take end transitions into account). */
	/* Similarly to begin transitions, sum, for each state, all its out transitions scores. */
	std::unordered_map<std::size_t, double> out_transitions_scores;
	for(std::size_t transition_id = 0; transition_id < _model->free_transitions.size(); ++transition_id){
		state_id = _model->free_transitions[transition_id].first;
		/* If state i not yet in map, init count for state i at 0. */
		if(out_transitions_scores.find(state_id) == out_transitions_scores.end()){
			out_transitions_scores[state_id] = utils::kNegInf;
		}
		out_transitions_scores[state_id] = utils::sum_log_prob(out_transitions_scores[state_id], transitions_scores.score(0, transition_id));
	}
	/* Also add end transitions scores to the the sum. */
	for(std::size_t end_transition_id = 0; end_transition_id < _model->free_pi_end.size(); ++end_transition_id){
		state_id = _model->free_pi_end[end_transition_id];
		if(out_transitions_scores.find(state_id) == out_transitions_scores.end()){
			out_transitions_scores[state_id] = utils::kNegInf;
		}
		out_transitions_scores[state_id] = utils::sum_log_prob(out_transitions_scores[state_id], transitions_scores.score_end(0, end_transition_id));
	}
	/* Normalize each transition by using the sum. */
	std::size_t from_state, to_state;
	for(std::size_t transition_id = 0; transition_id < _model->free_transitions.size(); ++transition_id){
		from_state = _model->free_transitions[transition_id].first; to_state = _model->free_transitions[transition_id].second;
		if(out_transitions_scores[from_state] != utils::kNegInf){
			_model->A[from_state][to_state] =  transitions_scores.score(0, transition_id) - out_transitions_scores[from_state];
		}
	}
	/* Don't forget to update the end transitions ! */
	for(std::size_t end_transition_id = 0; end_transition_id < _model->free_pi_end.size(); ++end_transition_id){
		state_id = _model->free_pi_end[end_transition_id];
		if(out_transitions_scores[state_id] != utils::kNegInf){
			_model->pi_end[state_id] = transitions_scores.score_end(0, end_transition_id) - out_transitions_scores[state_id];
		}
	}
}

void LinearMemoryBaumWelchTraining::update_model_emissions_from_log_scores(const EmissionScore& emissions_scores){
	std::unordered_map<std::size_t, double> all_emissions_scores;
	std::size_t state_id;
	for(std::size_t emission_id = 0; emission_id < _model->free_emissions.size(); ++emission_id){
		state_id = _model->free_emissions[emission_id].first;
		if(all_emissions_scores.find(state_id) == all_emissions_scores.end()){
			all_emissions_scores[state_id] = utils::kNegInf;
		}
		all_emissions_scores[state_id] = utils::sum_log_prob(all_emissions_scores[state_id], emissions_scores.score(0, emission_id));
	}
	std::string symbol;
	for(std::size_t emission_id = 0; emission_id < _model->free_emissions.size(); ++emission_id){
		state_id = _model->free_emissions[emission_id].first;
		symbol = _model->free_emissions[emission_id].second;
		if(all_emissions_scores[state_id] != utils::kNegInf) {
			(*(_model->B[state_id]))[symbol] = emissions_scores.score(0, emission_id) - all_emissions_scores[state_id];
		}
	}
}

double LinearMemoryBaumWelchTraining::train(const std::vector<std::vector<std::string>>& sequences, 
	double transition_pseudocount, double convergence_threshold, unsigned int min_iterations, unsigned int max_iterations){

	if(transition_pseudocount > 0) { std::cout << "Warning : baum-welch algorithm does not add pseudocounts ! "; }

	TransitionScore total_transition_score(_model->free_transitions, _model->free_pi_begin, _model->free_pi_end, 1, utils::kNegInf);
	EmissionScore total_emission_score(_model->free_emissions, 1, utils::kNegInf);

	TransitionScore previous_transition_score(_model->free_transitions, _model->free_pi_begin, _model->free_pi_end, _model->A.size(), utils::kNegInf);
	TransitionScore current_transition_score(_model->free_transitions, _model->free_pi_begin, _model->free_pi_end, _model->A.size(), utils::kNegInf);
	EmissionScore previous_emission_score(_model->free_emissions, _model->A.size(), utils::kNegInf);
	EmissionScore current_emission_score(_model->free_emissions, _model->A.size(), utils::kNegInf);
	
	unsigned int iteration = 0;
	double delta = utils::kInf;
	double initial_likelihood = _backward_algorithm.log_likelihood(sequences);
	double previous_likelihood = initial_likelihood;
	double current_likelihood;
	std::vector<double> previous_beta, beta, beta_end;
	std::size_t i, j, state_id;
	double score;
	std::string gamma;
	while((iteration < min_iterations || delta > convergence_threshold) 
		&& iteration < max_iterations) {
			/* Iterate over each sequence and compute the counts. */
		for(const std::vector<std::string>& sequence : sequences){
			/* If sequence is empty, go to current sequence. */
			if(sequence.size() == 0) { continue; }

			/* Initialization. */
			beta = _backward_algorithm.backward_init();
			for(std::size_t m = _model->A.size(); m-- > 0;){
				for(std::size_t free_emission_id = 0; free_emission_id < current_emission_score.num_free_emissions(); ++free_emission_id){
					state_id = current_emission_score.get_state_id(free_emission_id);
					gamma = current_emission_score.get_symbol(free_emission_id);
					score = beta[state_id] + log_score(sequence[sequence.size() - 1], gamma) + log_delta(state_id, m);
					current_emission_score.set_score(m, free_emission_id, score);	
				}

				/* Compute the transitions scores for silent states paths to the end state. Same behavior as in backward_init. */
				for(std::size_t free_end_transition_id = 0; free_end_transition_id < current_transition_score.num_free_end_transitions(); ++free_end_transition_id){
					state_id = current_transition_score.get_state_id_to_end(free_end_transition_id);
					score = _model->pi_end[state_id] + log_delta(m, state_id);
					for(std::size_t n = std::max(m + 1, _model->silent_states_index); n < _model->A.size(); ++n){
						score = utils::sum_log_prob(score, current_transition_score.score_end(n, free_end_transition_id) + _model->A[m][n]);
					}
					current_transition_score.set_end_score(m, free_end_transition_id, score);
				}

				for(std::size_t free_transition_id = 0; free_transition_id < current_transition_score.num_free_transitions(); ++free_transition_id){
					i = current_transition_score.get_from_state_id(free_transition_id);
					j = current_transition_score.get_to_state_id(free_transition_id);
					if(j >= _model->silent_states_index){
						score = beta[j] + log_delta(i, m) + _model->A[m][j];
						for(std::size_t n = std::max(m + 1, _model->silent_states_index); n < _model->A.size(); ++n){
							score = utils::sum_log_prob(score, current_transition_score.score(n, free_transition_id) + _model->A[m][n]);
						}
						current_transition_score.set_score(m, free_transition_id, score);
					}
				}
			}
			previous_beta = beta;
			previous_transition_score = current_transition_score; 
			previous_emission_score = current_emission_score;
			current_transition_score.reset();
			current_emission_score.reset();
			/* Recurrence. */
			for(std::size_t t = sequence.size() - 1; t-- > 0;){
				beta = _backward_algorithm.backward_step(previous_beta, sequence, t);
				for(std::size_t m = _model->A.size(); m-- > 0;){
					/* Compute transitions scores for current step. */
					for(std::size_t free_transition_id = 0; free_transition_id < current_transition_score.num_free_transitions(); ++free_transition_id){
						i = current_transition_score.get_from_state_id(free_transition_id);
						j = current_transition_score.get_to_state_id(free_transition_id);
						score = (j < _model->silent_states_index) ? previous_beta[j] + _model->A[m][j] + (*_model->B[j])[sequence[t + 1]] + log_delta(i, m) : beta[j] + _model->A[m][j] + log_delta(i, m);
						/* Consider previous step non-silent states. */
						for(std::size_t n = 0; n < _model->silent_states_index; ++n){
							score = utils::sum_log_prob(score, previous_transition_score.score(n, free_transition_id) + _model->A[m][n] + (*_model->B[n])[sequence[t + 1]]);
						}
						/* Consider current step silent states. */
						for(std::size_t n = std::max(m + 1, _model->silent_states_index); n < _model->A.size(); ++n){
							score = utils::sum_log_prob(score,  current_transition_score.score(n, free_transition_id) + _model->A[m][n]);
						}
						current_transition_score.set_score(m, free_transition_id, score);
					}
					/* Compute end transitions scores. */
					for(std::size_t free_end_transition_id = 0; free_end_transition_id < current_transition_score.num_free_end_transitions(); ++free_end_transition_id){
						state_id = current_transition_score.get_state_id_to_end(free_end_transition_id);
						score = utils::kNegInf;
						/* Consider previous step non-silent states. */
						for(std::size_t n = 0; n < _model->silent_states_index; ++n){
							score = utils::sum_log_prob(score, previous_transition_score.score_end(n, free_end_transition_id) + _model->A[m][n] + (*_model->B[n])[sequence[t + 1]]);
						}
						/* Consider current step silent states. */
						for(std::size_t n = std::max(m + 1, _model->silent_states_index); n < _model->A.size(); ++n){
							score = utils::sum_log_prob(score,  current_transition_score.score_end(n, free_end_transition_id) + _model->A[m][n]);
						}
						current_transition_score.set_end_score(m, free_end_transition_id, score);
					}
					/* Compute emissions score for current step. */
					for(std::size_t free_emission_id = 0; free_emission_id < current_emission_score.num_free_emissions(); ++free_emission_id){
						state_id = current_emission_score.get_state_id(free_emission_id);
						gamma = current_emission_score.get_symbol(free_emission_id);
						score = beta[m] + log_score(sequence[t], gamma) + log_delta(m, state_id);
						/* Consider previous step non-silent states. */
						for(std::size_t n = 0; n < _model->silent_states_index; ++n){
							score = utils::sum_log_prob(score, previous_emission_score.score(n, free_emission_id) + _model->A[m][n] + (*_model->B[n])[sequence[t + 1]]);
						}
						/* Consider current step silent states. */
						for(std::size_t n = std::max(m + 1, _model->silent_states_index); n < _model->A.size(); ++n){
							score = utils::sum_log_prob(score,  current_emission_score.score(n, free_emission_id) + _model->A[m][n]);
						}
						current_emission_score.set_score(m, free_emission_id, score);
					}
				}
				previous_beta = beta;
				previous_transition_score = current_transition_score;
				previous_emission_score = current_emission_score;
				current_transition_score.reset();
				current_emission_score.reset();
			}
			/* Termination. */
			std::tie(beta, beta_end, std::ignore) = _backward_algorithm.backward_terminate(beta, sequence);

			/* Compute the transitions scores for silent states paths to the begin state. 
			This essentially uses the same loop as the first loop in backward_terminate. */
			for(std::size_t m = _model->A.size(); m-- > _model->silent_states_index;){
				for(std::size_t free_transition_id = 0; free_transition_id < current_transition_score.num_free_transitions(); ++free_transition_id){
					i = current_transition_score.get_from_state_id(free_transition_id);
					j = current_transition_score.get_to_state_id(free_transition_id);
					score = (j < _model->silent_states_index) ? previous_beta[j] + _model->A[m][j] + (*_model->B[j])[sequence[0]] + log_delta(i, m) : beta[j] + _model->A[m][j] + log_delta(i, m);
					/* Consider previous step non-silent states. */
					for(std::size_t n = 0; n < _model->silent_states_index; ++n){
						score = utils::sum_log_prob(score, previous_transition_score.score(n, free_transition_id) + _model->A[m][n] + (*_model->B[n])[sequence[0]]);
					}
					/* Consider current step silent states. */
					for(std::size_t n = std::max(m + 1, _model->silent_states_index); n < _model->A.size(); ++n){
						score = utils::sum_log_prob(score,  current_transition_score.score(n, free_transition_id) + _model->A[m][n]);
					}
					current_transition_score.set_score(m, free_transition_id, score);
				}
				for(std::size_t free_end_transition_id = 0; free_end_transition_id < current_transition_score.num_free_end_transitions(); ++free_end_transition_id){
					state_id = current_transition_score.get_state_id_to_end(free_end_transition_id);
					score = utils::kNegInf;
					/* Consider previous step non-silent states. */
					for(std::size_t n = 0; n < _model->silent_states_index; ++n){
						score = utils::sum_log_prob(score, previous_transition_score.score_end(n, free_end_transition_id) + _model->A[m][n] + (*_model->B[n])[sequence[0]]);
					}
					/* Consider current step silent states. */
					for(std::size_t n = std::max(m + 1, _model->silent_states_index); n < _model->A.size(); ++n){
						score = utils::sum_log_prob(score,  current_transition_score.score_end(n, free_end_transition_id) + _model->A[m][n]);
					}
					current_transition_score.set_end_score(m, free_end_transition_id, score);
				}
				for(std::size_t free_emission_id = 0; free_emission_id < current_emission_score.num_free_emissions(); ++free_emission_id){
					state_id = current_emission_score.get_state_id(free_emission_id);
					gamma = current_emission_score.get_symbol(free_emission_id);
					score = utils::kNegInf;
					/* Consider previous step non-silent states. */
					for(std::size_t n = 0; n < _model->silent_states_index; ++n){
						score = utils::sum_log_prob(score, previous_emission_score.score(n, free_emission_id) + _model->A[m][n] + (*_model->B[n])[sequence[0]]);
					}
					/* Consider current step silent states. */
					for(std::size_t n = std::max(m + 1, _model->silent_states_index); n < _model->A.size(); ++n){
						score = utils::sum_log_prob(score,  current_emission_score.score(n, free_emission_id) + _model->A[m][n]);
					}
					current_emission_score.set_score(m, free_emission_id, score);
				}

			}
			for(std::size_t m = 0; m < _model->silent_states_index; ++m){
				for(std::size_t free_transition_id = 0; free_transition_id < current_transition_score.num_free_transitions(); ++free_transition_id){
					current_transition_score.set_score(m, free_transition_id, previous_transition_score.score(m, free_transition_id));
				}
				for(std::size_t free_end_transition_id = 0; free_end_transition_id < current_transition_score.num_free_end_transitions(); ++free_end_transition_id){
					current_transition_score.set_end_score(m, free_end_transition_id, previous_transition_score.score_end(m, free_end_transition_id));
				}
				for(std::size_t free_emission_id = 0; free_emission_id < previous_emission_score.num_free_emissions(); ++free_emission_id){
					current_emission_score.set_score(m, free_emission_id, previous_emission_score.score(m, free_emission_id));
				}
			}

			/* Begin transitions. */
			for(std::size_t free_begin_transition_id = 0; free_begin_transition_id < previous_transition_score.num_free_begin_transitions(); ++free_begin_transition_id){
				state_id = previous_transition_score.get_state_id_from_begin(free_begin_transition_id);
				current_transition_score.set_begin_score(0, free_begin_transition_id, beta_end[state_id]);
			}

			for(std::size_t m = 0; m < _model->A.size(); ++m){
				score = (m < _model->silent_states_index) ? _model->pi_begin[m] + (*_model->B[m])[sequence[0]] :  _model->pi_begin[m];
				for(std::size_t free_transition_id = 0; free_transition_id < current_transition_score.num_free_transitions(); ++free_transition_id){
					current_transition_score.set_score(m, free_transition_id, current_transition_score.score(m, free_transition_id) + score);
				}
				for(std::size_t free_end_transition_id = 0; free_end_transition_id < previous_transition_score.num_free_end_transitions(); ++free_end_transition_id){
					current_transition_score.set_end_score(m, free_end_transition_id, current_transition_score.score_end(m, free_end_transition_id) + score);
				}
				for(std::size_t free_emission_id = 0; free_emission_id < previous_emission_score.num_free_emissions(); ++free_emission_id){
					current_emission_score.set_score(m, free_emission_id, current_emission_score.score(m, free_emission_id) + score);
				}
			}

			/* Update total scores. */
			double seq_log_likelihood = _backward_algorithm.log_likelihood(sequence);
			/* Transitions. */
			log_update_transition_score(current_transition_score, total_transition_score, seq_log_likelihood);

			/* Emissions. */
			log_update_emission_score(current_emission_score, total_emission_score, seq_log_likelihood);

			current_transition_score.reset();
			previous_transition_score.reset();
			current_emission_score.reset();
			previous_emission_score.reset();
		}

		/* No pseudocount for b-w training ! */
		update_model_from_log_scores(total_transition_score, total_emission_score);
		total_transition_score.reset();
		total_emission_score.reset();
		current_likelihood = _backward_algorithm.log_likelihood(sequences);
		delta = current_likelihood - previous_likelihood;
		previous_likelihood = current_likelihood;
		++iteration;
	}
	return current_likelihood - initial_likelihood;
}

void LinearMemoryBaumWelchTraining::log_update_transition_score(const TransitionScore& current_transition_score, TransitionScore& total_transition_score, double seq_log_likelihood){
	double score;
	/* Begin transitions. */
	for(std::size_t free_begin_transition_id = 0; free_begin_transition_id < current_transition_score.num_free_begin_transitions(); ++free_begin_transition_id){
		score = utils::sum_log_prob(total_transition_score.score_begin(0, free_begin_transition_id), current_transition_score.score_begin(0, free_begin_transition_id) - seq_log_likelihood);
		total_transition_score.set_begin_score(0, free_begin_transition_id, score);
	}

	/* Mid transitions. */
	for(std::size_t free_transition_id = 0; free_transition_id < current_transition_score.num_free_transitions(); ++free_transition_id){
		score = utils::kNegInf;
		for(std::size_t m = 0; m < _model->A.size(); ++m){
			score = utils::sum_log_prob(score, current_transition_score.score(m, free_transition_id));
		}

		score = utils::sum_log_prob(score - seq_log_likelihood, total_transition_score.score(0, free_transition_id));
		total_transition_score.set_score(0, free_transition_id, score);
	}
	
	/* End transitions. */
	for(std::size_t free_end_transition_id = 0; free_end_transition_id < current_transition_score.num_free_end_transitions(); ++free_end_transition_id){
		score = utils::kNegInf;
		for(std::size_t m = 0; m < _model->A.size(); ++m){
			score = utils::sum_log_prob(score, current_transition_score.score_end(m, free_end_transition_id));
		}
		score = utils::sum_log_prob(score - seq_log_likelihood, total_transition_score.score_end(0, free_end_transition_id));
		total_transition_score.set_end_score(0, free_end_transition_id, score);
	}
}

void LinearMemoryBaumWelchTraining::log_update_emission_score(const EmissionScore& current_emission_score, EmissionScore& total_emission_score, double seq_log_likelihood){
	double score;
	for(std::size_t free_emission_id = 0; free_emission_id < current_emission_score.num_free_emissions(); ++free_emission_id){
		score = utils::kNegInf;
		for(std::size_t m = 0; m < _model->A.size(); ++m){
			score = utils::sum_log_prob(score, current_emission_score.score(m, free_emission_id));
		}
		score = utils::sum_log_prob(score - seq_log_likelihood, total_emission_score.score(0, free_emission_id));
		total_emission_score.set_score(0, free_emission_id, score);
	}
}


