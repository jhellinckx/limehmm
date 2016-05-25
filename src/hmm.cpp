#include <iostream>
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
#include <fstream>
#include "constants.hpp"
#include "state.hpp"
#include "graph.hpp"
#include "utils.hpp"
#include "hmm_algorithms.hpp"
#include "hmm_base.hpp"
#include "hmm.hpp"

#define CYAN "\033[36m"
#define RESET "\033[0m"

HiddenMarkovModel::HiddenMarkovModel() : HiddenMarkovModel(hmm_config::kDefaultHMMName) {}

HiddenMarkovModel::HiddenMarkovModel(const std::string& name) : 
	HiddenMarkovModel(name, 
		LinearMemoryForwardAlgorithm(nullptr),
		LinearMemoryBackwardAlgorithm(nullptr), 
		LinearMemoryViterbiDecodingAlgorithm(nullptr),
		LinearMemoryViterbiTraining(nullptr)) {}

HiddenMarkovModel::HiddenMarkovModel(
	const std::string& name, 
	const ForwardAlgorithm& forward, const BackwardAlgorithm& backward, 
	const DecodingAlgorithm& decode, const TrainingAlgorithm& train) :
		_name(name), _begin(nullptr), _end(nullptr), _graph(), _model(new RawModel()),
		_forward_algorithm(forward.clone()), _backward_algorithm(backward.clone()),
		_decoding_algorithm(decode.clone()), _training_algorithm(train.clone())
			{	
				_graph.add_vertex(hmm_config::kDefaultStartStateLabel);
				_begin = _graph.get_vertex(hmm_config::kDefaultStartStateLabel);
				_graph.add_vertex(hmm_config::kDefaultEndStateLabel);
				_end = _graph.get_vertex(hmm_config::kDefaultEndStateLabel);
				_forward_algorithm->set_model(_model); _backward_algorithm->set_model(_model);
				_decoding_algorithm->set_model(_model); _training_algorithm->set_model(_model);
			}

HiddenMarkovModel::HiddenMarkovModel(const HiddenMarkovModel& other) :
	_name(other._name), _begin(), _end(), _graph(other._graph), _model(new RawModel(*(other._model))),
	_forward_algorithm(other._forward_algorithm->clone()), _backward_algorithm(other._backward_algorithm->clone()),
	_decoding_algorithm(other._decoding_algorithm->clone()), _training_algorithm(other._training_algorithm->clone())
		{
			_begin = _graph.get_vertex(*other._begin); _end = _graph.get_vertex(*other._end);
			_forward_algorithm->set_model(_model); _backward_algorithm->set_model(_model);
			_decoding_algorithm->set_model(_model); _training_algorithm->set_model(_model);
		}

HiddenMarkovModel::HiddenMarkovModel(HiddenMarkovModel&& other) : 
	_name(std::move(other._name)), _begin(std::move(other._begin)), _end(std::move(other._end)), 
	_graph(std::move(other._graph)),  _model(std::move(other._model)), _forward_algorithm(std::move(other._forward_algorithm)), 
	_backward_algorithm(std::move(other._backward_algorithm)), _decoding_algorithm(std::move(other._decoding_algorithm)), 
	_training_algorithm(std::move(other._training_algorithm)) 
		{
			_forward_algorithm->set_model(_model); _backward_algorithm->set_model(_model);
			_decoding_algorithm->set_model(_model); _training_algorithm->set_model(_model);
		}	

HiddenMarkovModel& HiddenMarkovModel::operator=(const HiddenMarkovModel& other){
	if(this != &other){
		_name = other._name;
		_graph = other._graph;
		_begin = _graph.get_vertex(*other._begin);
		_end = _graph.get_vertex(*other._end);
		delete _forward_algorithm; delete _backward_algorithm; delete _decoding_algorithm; delete _training_algorithm;
		_forward_algorithm = other._forward_algorithm->clone();
		_backward_algorithm = other._backward_algorithm->clone();
		_decoding_algorithm = other._decoding_algorithm->clone();
		_training_algorithm = other._training_algorithm->clone();
		*_model = *other._model;
		_forward_algorithm->set_model(_model); _backward_algorithm->set_model(_model);
		_decoding_algorithm->set_model(_model); _training_algorithm->set_model(_model);
	}
	return *this;
}

HiddenMarkovModel& HiddenMarkovModel::operator=(HiddenMarkovModel&& other){
	if(this != &other){
		_model->clean();
		_name = std::move(other._name);
		_begin = std::move(other._begin);
		_end = std::move(other._end);
		_graph = std::move(other._graph);
		_forward_algorithm = std::move(other._forward_algorithm);
		_backward_algorithm = std::move(other._backward_algorithm);
		_decoding_algorithm = std::move(other._decoding_algorithm);
		_training_algorithm = std::move(other._training_algorithm);
		_model = std::move(other._model);
		_forward_algorithm->set_model(_model); _backward_algorithm->set_model(_model);
		_decoding_algorithm->set_model(_model); _training_algorithm->set_model(_model);
	}
	return *this;
}

void HiddenMarkovModel::print_transitions(bool log_prob){
	__print_transitions(_model->A, _model->states_indices, log_prob);
}

void HiddenMarkovModel::print_distributions(bool log_prob){
	__print_distributions(_model->B, _model->states_names, log_prob);
}

void HiddenMarkovModel::print_pi_begin(bool log_prob){
	__print_pi_begin(_model->pi_begin, _model->states_names, log_prob);
}

void HiddenMarkovModel::print_pi_end(bool log_prob){
	__print_pi_end(_model->pi_end, _model->states_names, log_prob);
}

std::string HiddenMarkovModel::name() const { return _name; }
void HiddenMarkovModel::set_name(const std::string& name) { _name = name; } 
std::size_t HiddenMarkovModel::num_states() const { return _graph.num_vertices(); }
std::size_t HiddenMarkovModel::num_transitions() const { return _graph.num_edges(); }	

Graph<State> HiddenMarkovModel::get_graph() { return _graph; }

bool HiddenMarkovModel::has_state(const State& state) const {
	return _graph.has_vertex(state);
}

bool HiddenMarkovModel::has_transition(const State& from_state, const State& to_state) const {
	return _graph.has_edge(from_state, to_state);
}

State& HiddenMarkovModel::begin() { 
	if(_begin != nullptr){
		return *_begin;
	}
	else{
		throw StateNotFoundException(error_message::kHMMHasNoBeginState);
	}
}

State& HiddenMarkovModel::end() {
	if(_end != nullptr){
		return *_end;
	}
	else{
		throw StateNotFoundException(error_message::kHMMHasNoEndState);
	}
}

State& HiddenMarkovModel::get_state(const State& state) {
	try{
		return *_graph.get_vertex(state);
	}
	catch(const VertexNotFoundException<State>& e){
		throw StateNotFoundException(e.trigger(), error_message::kHMMGetStateNotFound);
	}
}

void HiddenMarkovModel::add_state(const State& state){
	try{
		_graph.add_vertex(state);	
	}
	catch(const VertexExistsException<State>& e){
		throw StateExistsException(e.trigger(), error_message::kHMMAddStateExists);
	}
}

void HiddenMarkovModel::remove_state(const State& state){
	if(state == *_begin) _begin = nullptr;
	else if(state == *_end) _end = nullptr;
	try{
		_graph.remove_vertex(state);	
	}
	catch(const VertexNotFoundException<State>& e){
		throw StateNotFoundException(e.trigger(), error_message::kHMMRemoveStateNotFound);
	}
}

std::string HiddenMarkovModel::transition_string(const State& from, const State& to) const {
	return from.to_string() + " -> " + to.to_string();
}

void HiddenMarkovModel::add_transition(const State& from, const State& to, double probability){
	if(from == end()) throw TransitionLogicException(transition_string(from, to), error_message::kHMMAddedTransitionFromEndState);
	if(to == begin()) throw TransitionLogicException(transition_string(from, to), error_message::kHMMAddedTransitionToBeginState);
	if(probability < 0) throw TransitionLogicException(transition_string(from, to), error_message::kHMMTransitionNegativeProbability);
	try{
		_graph.add_edge(from, to, probability);	
	} 
	catch(const EdgeExistsException<Edge<State>>& e){
		throw TransitionExistsException(transition_string(*(e.trigger().from()), *(e.trigger().to())), error_message::kHMMAddTransitionExists);
	}
	catch(const IncidentVertexNotFoundException<State>& e){
		throw StateNotFoundException(e.trigger(), error_message::kHMMAddTransitionStateNotFound);
	}
}

double HiddenMarkovModel::get_transition(const State& from, const State& to){
	try{
		Edge<State>* edge = _graph.get_edge(from, to);
		double* weight = edge->weight();
		if(weight == nullptr) { 
			throw TransitionLogicException(transition_string(*(edge->from()), *(edge->to())), error_message::kHMMGetTransitionNullWeight);
		}
		else{ return *weight; }
	}
	catch(const EdgeNotFoundException<Edge<State>>& e){
		throw TransitionNotFoundException(transition_string(*(e.trigger().from()), *(e.trigger().to())), error_message::kHMMGetTransitionNotFound);
	}
}

void HiddenMarkovModel::set_transition(const State& from, const State& to, double probability){
	if(probability < 0) throw TransitionLogicException(transition_string(from, to), error_message::kHMMTransitionNegativeProbability);
	try{
		_graph.set_weight(from, to, probability);	
	}
	catch(const IncidentVertexNotFoundException<State>& e){
		throw StateNotFoundException(e.trigger(), error_message::kHMMGetTransitionNotFound);
	}
}

void HiddenMarkovModel::begin_transition(const State& state, double probability) {
	add_transition(begin(), state, probability);
}

void HiddenMarkovModel::end_transition(const State& state, double probability) {
	add_transition(state, end(), probability);
}

void HiddenMarkovModel::remove_transition(const State& from, const State& to){
	try{
		_graph.remove_edge(from, to);	
	}
	catch(const EdgeNotFoundException<Edge<State>>& e){
		throw TransitionNotFoundException(transition_string(*(e.trigger().from()), *(e.trigger().to())), error_message::kHMMRemoveTransitionNotFound);
	}
	
}


void HiddenMarkovModel::brew(bool normalize) {
	/* Get rid of previous data. */
	_model->clean();

	/* Get the states from graph. */
	std::vector<State*> states = _graph.get_vertices();

	/* Remove begin and end states. */
	states.erase(std::remove_if(states.begin(), states.end(), [this](State* p_state){ return (*p_state) == begin() || (*p_state) == end(); }), states.end());
	std::vector<State*> silent_states;
	std::size_t num_states = states.size();

	/* Keep track of the matrix index of each state. */
	std::map<std::string, std::size_t> states_indices;
	std::vector<std::string> states_names(num_states);

	/* Init size of raw transition matrix. */
	Matrix A(num_states);
	std::vector<double> pi_begin(num_states, utils::kNegInf);
	std::vector<double> pi_end(num_states, utils::kNegInf);

	/* Transitions to end state exist and are not empty. */
	bool finite = false;

	/* Check if silent/end states are indeed silent. */
	if(!begin().is_silent()) { throw std::logic_error("begin state has to be silent."); }
	if(!end().is_silent()) { throw std::logic_error("end state has to be silent."); }

	std::size_t normal_states_index = 0;
	for(State* p_state : states){
		if(p_state->is_silent()) {
			silent_states.push_back(p_state);
		}
		else{
			A[normal_states_index] = std::vector<double>(num_states, utils::kNegInf);
			/* Map state name to row index. */
			states_indices[p_state->name()] = normal_states_index;
			states_names[normal_states_index] = p_state->name();
			++normal_states_index;
		}
	}
	std::size_t num_silent_states = silent_states.size();
	/* This points to the beginning of the silent states array in A and B. */
	std::size_t silent_states_index = normal_states_index;
	/* We use a topological sort on the silent states sub graph in order to adapt 
	the HMM to silent states. */
	/* See p. 71 at http://www.upch.edu.pe/facien/fc/dbmbqf/zimic/ubioinfo/bks/Bioinformatics/Biological%20Sequence%20Analysis%20Hmm%20Bioinformatics%20(Durbin).pdf */
	/* Start by dereferencing silent states pointers to pass them to subgraph. */
	std::vector<State> silent_states_values;
	silent_states_values.reserve(num_silent_states);
	for(State* p_state : silent_states) { silent_states_values.push_back(*p_state); }
	Graph<State> subgraph = _graph.sub_graph(silent_states_values);
	subgraph.topological_sort();
	/* Get toposorted silent states. */
	silent_states = subgraph.get_vertices();
	/* Init the toposorted silent states rows in the matrix. */
	for(State* p_silent_state : silent_states){
		A[silent_states_index] = std::vector<double>(num_states, utils::kNegInf);
		states_indices[p_silent_state->name()] = silent_states_index;
		states_names[silent_states_index] = p_silent_state->name();
		++silent_states_index;
	}

	/* Fill transitions with log probabilities and check whether a normalization is needed. */
	auto fill_normalize = [this, &pi_end, &states_indices] (const std::vector<Edge<State>*>& edges, std::vector<double>& prob_vec_to_fill, bool normalize) {
		std::vector<double> vec_to_normalize;
		vec_to_normalize.reserve(edges.size());
		double prob_sum = 0;
		double prob;
		for(Edge<State>* edge : edges){ 
			prob = (edge->weight() == nullptr) ? 0 : *(edge->weight());
			prob_sum += prob;
			vec_to_normalize.push_back(log(prob));
		}
		if(prob_sum != 1.0 && normalize){
			std::vector<double>::iterator it = vec_to_normalize.begin();
			while(it != vec_to_normalize.end()){
				*it = utils::log_normalize(*it, log(prob_sum));
				++it;
			}
		}
		std::size_t i = 0;
		for(double& log_prob : vec_to_normalize){
			if(*(edges[i]->to()) == end()){
				pi_end[states_indices[edges[i]->from()->name()]] = log_prob;
			}
			else{
				prob_vec_to_fill[states_indices[edges[i]->to()->name()]] = log_prob;
			}
			i++;
		}
		return prob_sum;
	};

	/* Add the begin state transitions. */
	State& begin_state = begin();
	if(_graph.get_in_edges(begin_state).size() > 0) { throw std::logic_error("begin state cannot have predecessors"); }
	std::vector<Edge<State>*> out_edges = _graph.get_out_edges(begin_state);
	double prob_sum = fill_normalize(out_edges, pi_begin, normalize);
	if(prob_sum == 0.0) { throw std::logic_error("hmm has no begin transition"); }
	
	/* Check if end state has out edges. */
	State& end_state = end();
	if(_graph.get_out_edges(end_state).size() > 0) { throw std::logic_error("end state cannot have successors"); }

	/* Iterate through all the other states and add them to the matrix. */
	for(State* p_state : states){
		/* Fill normal transitions aka matrix A. */
		std::vector<Edge<State>*> out_edges = _graph.get_out_edges(*p_state);
		double prob_sum = fill_normalize(out_edges, A[states_indices[p_state->name()]], normalize);
		if(prob_sum == 0.0) { throw std::logic_error("hmm has no transition from " + p_state->to_string()); }
	}

	/* Determine whether the hmm is finite by summing the end state in transitions probabilities. */
	double prob_sum_to_end = 0;
	for(std::size_t i = 0; i < A.size(); ++i) { prob_sum_to_end += exp(pi_end[i]); }
	if(prob_sum_to_end > 0.0) { finite = true; }

	/* Fill emission matrix with the states PDFs. */
	std::vector<Distribution*> B(num_states);
	for(State* p_state : states){
		Distribution* distribution = nullptr;
		if(! p_state->is_silent()) {
			distribution = p_state->distribution().clone();
			distribution->log_probabilities(true);
			if(normalize){
				distribution->log_normalize();	
			}
		}
		B[states_indices[p_state->name()]] = distribution;
	}

	/* Get alphabet. Only discrete ! */
	std::vector<std::string> alphabet;
	for(const State* p_state : states){
		if(! p_state->is_silent()){
			std::vector<std::string> dist_symbols = static_cast<DiscreteDistribution*>(&p_state->distribution())->symbols();
			for(const std::string& symbol : dist_symbols){
				if(std::find(alphabet.begin(), alphabet.end(), symbol) == alphabet.end()){
					alphabet.push_back(symbol);
				}
			}
		}
	}

	/* Set free emissions. Only discrete ! */
	std::vector<std::pair<std::size_t, std::string>> free_emissions;
	for(const State* p_state : states){
		if((! p_state->is_silent()) && p_state->has_free_emission()){
			std::size_t state_id = states_indices[p_state->name()];
			for(const std::string& symbol : alphabet){
				free_emissions.push_back(std::make_pair(state_id, symbol));	
			}
		}
	}
	
	/* Set free transitions. */
	std::vector<std::pair<std::size_t, std::size_t>> free_transitions;
	std::vector<std::size_t> free_pi_begin;
	std::vector<std::size_t> free_pi_end;
	for(const State* p_state : states){
		if(p_state->has_free_transition()){
			std::size_t state_id = states_indices[p_state->name()];
			auto out_edges = _graph.get_out_edges(*p_state);
			for(auto& edge: out_edges){
				if(*(edge->to()) == end()){ free_pi_end.push_back(state_id); }
				else { free_transitions.push_back(std::make_pair(state_id, states_indices[edge->to()->name()])); }
			}
		}
	}
	if(begin().has_free_transition()){
		auto begin_out_edges = _graph.get_out_edges(begin());
		for(auto edge: begin_out_edges){
			free_pi_begin.push_back(states_indices[edge->to()->name()]);
		}	
	}
	
	/* Set fields for the hmm raw values. */
	_model->A = std::move(A);
	_model->B = std::move(B);
	_model->pi_begin = std::move(pi_begin);
	_model->pi_end = std::move(pi_end);
	_model->states_indices = std::move(states_indices);
	_model->states_names = std::move(states_names);
	_model->is_finite = finite;
	_model->silent_states_index = normal_states_index;
	_model->alphabet = std::move(alphabet); // DISCRETE ONLY !!
	_model->free_emissions = std::move(free_emissions);
	_model->free_transitions = std::move(free_transitions);
	_model->free_pi_begin = std::move(free_pi_begin);
	_model->free_pi_end = std::move(free_pi_end);
}

Matrix HiddenMarkovModel::raw_transitions() { return _model->A; }
std::vector<double> HiddenMarkovModel::raw_pi_begin() { return _model->pi_begin; }
std::vector<double> HiddenMarkovModel::raw_pi_end() { return _model->pi_end; }
std::vector<Distribution*> HiddenMarkovModel::raw_pdfs() { return _model->B; }
std::map<std::string, std::size_t> HiddenMarkovModel::states_indices() { return _model->states_indices; }
std::vector<std::string> HiddenMarkovModel::states_names() { return _model->states_names; }

void HiddenMarkovModel::set_forward(const ForwardAlgorithm& forward) {
	delete _forward_algorithm; 
	_forward_algorithm = forward.clone();
	_forward_algorithm->set_model(_model);
}

void HiddenMarkovModel::set_backward(const BackwardAlgorithm& backward) {
	delete _backward_algorithm; 
	_backward_algorithm = backward.clone();
	_backward_algorithm->set_model(_model);
}

void HiddenMarkovModel::set_decoding(const DecodingAlgorithm& decode) {
	delete _decoding_algorithm; 
	_decoding_algorithm = decode.clone();
	_decoding_algorithm->set_model(_model);
}

void HiddenMarkovModel::set_training(const TrainingAlgorithm& training) {
	delete _training_algorithm; 
	_training_algorithm = training.clone();
	_training_algorithm->set_model(_model);
}

std::string HiddenMarkovModel::forward_type() const 	{ return _forward_algorithm->type(); }
std::string HiddenMarkovModel::backward_type() const 	{ return _backward_algorithm->type(); }
std::string HiddenMarkovModel::decoding_type() const 	{ return _decoding_algorithm->type(); }
std::string HiddenMarkovModel::training_type() const 	{ return _training_algorithm->type(); }

std::vector<double> HiddenMarkovModel::forward(const std::vector<std::string>& sequence, std::size_t t_max){
	return _forward_algorithm->forward(sequence, t_max);
}

std::vector<double> HiddenMarkovModel::backward(const std::vector<std::string>& sequence, std::size_t t_min){
	return _backward_algorithm->backward(sequence, t_min);
}

double HiddenMarkovModel::log_likelihood(const std::vector<std::string>& sequence, bool do_fwd){
	if(do_fwd){
		return _forward_algorithm->log_likelihood(sequence);
	}
	else{
		return _backward_algorithm->log_likelihood(sequence);
	}
}

double HiddenMarkovModel::log_likelihood(const std::vector<std::vector<std::string>>& sequences, bool do_fwd){
	if(do_fwd){
		return _forward_algorithm->log_likelihood(sequences);
	}
	else{
		return _backward_algorithm->log_likelihood(sequences);
	}
}

double HiddenMarkovModel::likelihood(const std::vector<std::string>& sequence, bool do_fwd){
	return exp(log_likelihood(sequence, do_fwd));
}

double HiddenMarkovModel::likelihood(const std::vector<std::vector<std::string>>& sequences, bool do_fwd){
	return exp(log_likelihood(sequences, do_fwd));
}

std::pair<std::vector<std::string>, double> HiddenMarkovModel::decode(const std::vector<std::string>& sequence, std::size_t t_max){
	return _decoding_algorithm->decode(sequence, t_max);
}

double HiddenMarkovModel::train(const std::vector<std::vector<std::string>>& sequences,
	double transition_pseudocount, double convergence_threshold,
	unsigned int min_iterations, unsigned int max_iterations){

		double improvement = _training_algorithm->train(sequences, transition_pseudocount, convergence_threshold, min_iterations, max_iterations);
		_update_from_raw();
		return improvement;
}

void HiddenMarkovModel::_update_from_raw(){
	/* Update transitions. Since we use log probabilities in the raw data, don't forget to exp() the log prob. */
	std::string from_state_name, to_state_name;
	std::size_t from_state_id, to_state_id;
	double log_probability;
	/* Update begin transitions. */
	for(std::size_t begin_transition_id = 0; begin_transition_id < _model->free_pi_begin.size(); ++begin_transition_id){
		to_state_id = _model->free_pi_begin[begin_transition_id];
		log_probability = _model->pi_begin[to_state_id];
		to_state_name = _model->states_names[to_state_id];
		_graph.set_weight(begin(), State(to_state_name), exp(log_probability));
	}
	/* Update mid transitions. */
	for(std::size_t transition_id = 0; transition_id < _model->free_transitions.size(); ++transition_id){
		from_state_id = _model->free_transitions[transition_id].first;
		to_state_id = _model->free_transitions[transition_id].second;
		log_probability = _model->A[from_state_id][to_state_id];
		from_state_name = _model->states_names[from_state_id];
		to_state_name = _model->states_names[to_state_id];
		_graph.set_weight(State(from_state_name), State(to_state_name), exp(log_probability));
	}
	/* Update end transitions. */
	for(std::size_t end_transition_id = 0; end_transition_id < _model->free_pi_end.size(); ++end_transition_id){
		from_state_id = _model->free_pi_end[end_transition_id];
		log_probability = _model->pi_end[from_state_id];
		from_state_name = _model->states_names[from_state_id];
		_graph.set_weight(State(from_state_name), end(), exp(log_probability));
	}

	/* Update emissions. */
	//DISCRETE ONLY !!
	std::string symbol;
	std::size_t state_id;
	std::string state_name;
	for(std::size_t emission_id = 0; emission_id < _model->free_emissions.size(); ++emission_id){
		state_id = _model->free_emissions[emission_id].first;
		state_name = _model->states_names[state_id];
		symbol = _model->free_emissions[emission_id].second;
		_graph.get_vertex(State(state_name))->distribution()[symbol] = exp((*_model->B[state_id])[symbol]);
	}
}


void HiddenMarkovModel::save() {
	save(_name);
}

void HiddenMarkovModel::save(const std::string& filename, const std::string& extension){
	std::ofstream savefile(filename + "." + extension);
	if(savefile.is_open()){
		std::vector<State*> states = _graph.get_vertices();
		states.erase(std::remove_if(states.begin(), states.end(), [this](State* p_state){ return (*p_state) == begin() || (*p_state) == end(); }), states.end());
		std::vector<Edge<State>*> edges;
		/* Name */
		savefile << _name << std::endl;
		/* Algorithms */
		savefile << _forward_algorithm->type() << std::endl;
		savefile << _backward_algorithm->type() << std::endl;
		savefile << _decoding_algorithm->type() << std::endl;
		savefile << _training_algorithm->type() << std::endl;
		/* Begin / end states names */
		savefile << begin().name() << std::endl;
		savefile << end().name() << std::endl;
		/* Number of non begin/end states */
		savefile << states.size() << std::endl;
		for(std::size_t i = 0; i < states.size(); ++i){
			savefile << states[i]->name() << std::endl;
		}
		/* Begin transitions */
		edges = _graph.get_out_edges(begin());
		savefile << begin().name() << std::endl;
		savefile << edges.size() << std::endl;
		for(std::size_t j = 0; j < edges.size(); ++j){
			Edge<State>& edge = *(edges[j]);
			savefile << edge.to()->name() << global_config::kProbabilitySeparator << (edge.weight() == nullptr ? global_config::kNullValue : std::to_string(*(edge.weight()))) << std::endl;
		}
		/* All transitions */
		for(std::size_t i = 0; i < states.size(); ++i){
			State& s = *(states[i]);
			savefile << s.name() << std::endl;
			edges = _graph.get_out_edges(s);
			savefile << edges.size() << std::endl;
			for(std::size_t j = 0; j < edges.size(); ++j){
				Edge<State>& edge = *(edges[j]);
				savefile << edge.to()->name() << global_config::kProbabilitySeparator << (edge.weight() == nullptr ? global_config::kNullValue : std::to_string(*(edge.weight()))) << std::endl;
			}
			if(! s.is_silent()){
				Distribution& dist = s.distribution();
				savefile << dist.name() << std::endl;
				dist.save(savefile);
			}
			else{
				savefile << global_config::kNullValue << std::endl;
			}
		}
		savefile.close();
	}
	else{
		throw std::runtime_error("Could not open save file.");
	}
}

void HiddenMarkovModel::load(const std::string& filename, const std::string& extension){
	std::ifstream loadfile(filename + "." + extension);
	if(loadfile.is_open()){
		std::string line;
		/* Get name */
		std::getline(loadfile, line);
		set_name(line);
		/* Algorithms */
		std::string algo_type;
		std::getline(loadfile, algo_type);
		if(algo_type == hmm_config::kLinearMemoryForwardAlgorithmName){
			set_forward(LinearMemoryForwardAlgorithm(_model));
		}
		else{
			std::cout << "Warning : unknown forward algorithm type. Defaults to linear memory forward." << std::endl;
			set_forward(LinearMemoryForwardAlgorithm(_model));
		}
		std::getline(loadfile, algo_type);
		if(algo_type == hmm_config::kLinearMemoryBackwardAlgorithmName){
			set_backward(LinearMemoryBackwardAlgorithm(_model));
		}
		else{
			std::cout << "Warning : unknown backward algorithm type. Defaults to linear memory backward." << std::endl;
			set_backward(LinearMemoryBackwardAlgorithm(_model));
		}
		std::getline(loadfile, algo_type);
		if(algo_type == hmm_config::kLinearMemoryViterbiDecodeAlgorithmName){
			set_decoding(LinearMemoryViterbiDecodingAlgorithm(_model));
		}
		else{
			std::cout << "Warning : unknown decoding algorithm type. Defaults to linear memory viterbi." << std::endl;
			set_decoding(LinearMemoryViterbiDecodingAlgorithm(_model));
		}
		std::getline(loadfile, algo_type);
		if(algo_type == hmm_config::kLinearMemoryViterbiTrainingAlgorithmName){
			set_training(LinearMemoryViterbiTraining(_model));
		}
		else if(algo_type == hmm_config::kLinearMemoryBaumWelchTrainingAlgorithmName){
			set_training(LinearMemoryBaumWelchTraining(_model));
		}
		else{
			std::cout << "Warning : unknown decoding algorithm type. Defaults to linear memory viterbi." << std::endl;
			set_training(LinearMemoryViterbiTraining(_model));
		}
		/* Begin / end states names */
		std::getline(loadfile, line);
		begin().set_name(line);
		std::getline(loadfile, line);
		end().set_name(line);
		/* Number of non begin/end states */
		std::size_t num_states;
		std::getline(loadfile, line);
		num_states = (std::size_t) std::stoi(line);
		for(std::size_t i = 0; i < num_states; ++i){
			std::getline(loadfile, line);
			add_state(State(line));
		}
		/* Transitions / Emissions */
		std::size_t num_transitions;
		std::string from_state;
		std::string state_name;
		std::string prob_str;
		/* Begin name */
		std::getline(loadfile, line);
		/* Num begin transitions */
		std::getline(loadfile, line);
		num_transitions = (std::size_t) std::stoi(line);
		/* Begin transitions */
		for(std::size_t j = 0; j < num_transitions; ++j){
			std::getline(loadfile, line);
			std::tie(state_name, prob_str) = utils::split_first(line, global_config::kProbabilitySeparator);
			if(prob_str != global_config::kNullValue && ! prob_str.empty()){
				add_transition(begin(), state_name, std::stod(prob_str));
			}
		}
		
		for(std::size_t i = 0; i < num_states; ++i){
			std::getline(loadfile, line);
			from_state = line;
			std::getline(loadfile, line);
			num_transitions = (std::size_t) std::stoi(line);
			/* All transitions */
			for(std::size_t j = 0; j < num_transitions; ++j){
				std::getline(loadfile, line);
				std::tie(state_name, prob_str) = utils::split_first(line, global_config::kProbabilitySeparator);
				if(prob_str != global_config::kNullValue && ! prob_str.empty()){
					add_transition(from_state, state_name, std::stod(prob_str));
				}
			}
			/* All emissions */
			std::getline(loadfile, line);
			std::string dist_name = line;
			if(dist_name == global_config::kNullValue){ continue; }
			else if(dist_name == distribution_config::kDiscreteDistributionName){
				DiscreteDistribution dist = DiscreteDistribution();
				dist.load(loadfile);
				get_state(from_state).set_distribution(dist);
			}
		}
		loadfile.close();
	}
	else{
		throw std::runtime_error("Could not open load file.");
	}
}

HiddenMarkovModel::~HiddenMarkovModel(){
	delete _forward_algorithm; delete _backward_algorithm; 
	delete _decoding_algorithm; delete _training_algorithm;
	delete _model;
}


/* ========================== UTILS ========================== */

void __print_transitions(const Matrix& matrix, const std::map<std::string, std::size_t>& indices, bool log_prob){
	std::size_t longest_string = 0;
	for(std::size_t i = 0; i < matrix.size(); ++i){
		for(std::size_t j = 0; j < matrix[i].size(); ++j){
			std::string double_string = (log_prob) ? std::to_string(matrix[i][j]) : std::to_string(exp(matrix[i][j]));
			if(double_string.length() > longest_string) longest_string = double_string.length();
		}
	}
	std::vector<std::string> sorted_names(indices.size());
	std::ostringstream out;
	for(auto& pair : indices) {
		sorted_names[pair.second] = pair.first;
	}
	out << std::string(longest_string + 1, ' ');
	for(std::string& name : sorted_names) {
		out << std::string(longest_string - name.length(), ' ');
		out << CYAN << name << RESET;
		out << ' ';
	}
	out << std::endl;
	for(std::size_t i = 0; i < matrix.size(); ++i){
		out << std::string(longest_string - sorted_names[i].length(), ' ');
		out << CYAN << sorted_names[i] << RESET;
		out << ' ';
		for(std::size_t j = 0; j < matrix[i].size(); ++j){
			std::string double_string = (log_prob) ? std::to_string(matrix[i][j]) : std::to_string(exp(matrix[i][j]));
			out << std::string(longest_string - double_string.length(), ' ');
			out << double_string;
			out << ' ';
		}
		out << std::endl;
	}
	std::cout << out.str() << std::endl;
}

void __print_distributions(const std::vector<Distribution*>& dists, const std::vector<std::string>& names, bool log_prob){
	std::ostringstream oss;
	for(std::size_t state_id = 0; state_id < dists.size(); ++state_id){
		Distribution* dist = dists[state_id];
		if(dist == nullptr) oss << "Silent";
		else{
			oss << names[state_id] << " : ";
			bool used_log_prob = dist->uses_log_probabilities();
			dist->log_probabilities(log_prob);
			oss << *dist;
			dist->log_probabilities(used_log_prob);
		}
		
		oss << std::endl;
	}
	std::cout << oss.str() << std::endl;
}

void __print_distributions(std::vector<DiscreteDistribution>& dists, const std::vector<std::string>& names, bool log_prob){
	std::vector<Distribution*> tmp_dists;
	for(auto& dist : dists){ tmp_dists.push_back(&dist); }
	__print_distributions(tmp_dists, names, log_prob);
}

void __print_names(const std::vector<std::size_t>& ids, const std::vector<std::string>& names){
	std::ostringstream oss;
	for(std::size_t id : ids){
		if(id < names.size()){
			oss << names[id] << " ";
		}
	}
	std::cout << oss.str() << std::endl;
}

void __print_pi_begin(const std::vector<double>& pi, const std::vector<std::string>& names, bool log_prob){
	std::ostringstream oss;
	oss << "Pi begin : ";
	for(std::size_t state_id = 0; state_id < pi.size(); ++state_id){
		oss << names[state_id] << "(";
		if(log_prob) oss << pi[state_id];
		else oss << exp(pi[state_id]);
		oss << ") ";
	}
	std::cout << oss.str() << std::endl;
}

void __print_pi_end(const std::vector<double>& pi, const std::vector<std::string>& names, bool log_prob){
	std::ostringstream oss;
	oss << "Pi end : ";
	for(std::size_t state_id = 0; state_id < pi.size(); ++state_id){
		oss << names[state_id] << "(";
		if(log_prob) oss << pi[state_id];
		else oss << exp(pi[state_id]);
		oss << ") ";
	}
	oss << std::endl;
	std::cout << oss.str() << std::endl;
}

void print_prob(const std::vector<double>& probs, bool log_prob){
	std::ostringstream oss;
	for(double d : probs){
		if(log_prob){oss << d << " ";}
		else oss << exp(d) << " ";
	}
	std::cout << oss.str() << std::endl;
}

std::ostream& operator<<(std::ostream& out, const std::vector<double>& vec){
	for(double d : vec){
		out << exp(d) << " ";
	}
	out << std::endl;
	return out;
}

std::ostream& operator<<(std::ostream& out, const std::vector<std::size_t>& vec){
	for(std::size_t ul : vec){
		out << ul << " ";
	}
	out << std::endl;
	return out;
}

std::ostream& operator<<(std::ostream& out, const std::vector<std::string>& vec){
	for(const std::string& s : vec){
		out << s << " ";
	}
	out << std::endl;
	return out;
}

std::ostream& operator<<(std::ostream& out, const std::vector<State>& vec){
	out << "States :" << std::endl;
	for(const State& s : vec){
		out << s << std::endl;
	}
	out << "----------------------------------------------" << std::endl;
	return out;
}

std::ostream& operator<<(std::ostream& out, const std::vector<State*>& vec){
	out << "States :" << std::endl;
	for(const State* s : vec){
		out << *s << std::endl;
	}
	out << "----------------------------------------------" << std::endl;
	return out;
}
