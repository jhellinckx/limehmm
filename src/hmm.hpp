#ifndef __HIDDENMARKOVMODEL_HPP
#define __HIDDENMARKOVMODEL_HPP

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

#define CYAN "\033[36m"
#define RESET "\033[0m"

/* <-------- HMM Exceptions --------> */

class HMMException : public std::logic_error {
protected:
	HMMException(const std::string& message) :
		std::logic_error(message) {}
};

class StateNotFoundException : public HMMException {
public:
	template<typename T>
	StateNotFoundException(const T& t, const std::string& msg) : 
		HMMException(error_message::format("StateNotFoundException: " + msg, t)) {}

	StateNotFoundException(const std::string& msg) : 
		HMMException("StateNotFoundException: " + msg) {}
};

class StateExistsException : public HMMException {
public:
	template<typename T>
	StateExistsException(const T& t, const std::string& msg) : 
		HMMException(error_message::format("StateExistsException: " + msg, t)) {}
};

class TransitionNotFoundException : public HMMException {
public:
	template<typename T>
	TransitionNotFoundException(const T& t, const std::string& msg) : 
		HMMException(error_message::format("TransitionNotFoundException: " + msg, t)) {}

	TransitionNotFoundException(const std::string& msg) : 
		HMMException("TransitionNotFoundException: " + msg) {}
};

class TransitionExistsException : public HMMException {
public:
	template<typename T>
	TransitionExistsException(const T& t, const std::string& msg) : 
		HMMException(error_message::format("TransitionExistsException: " + msg, t)) {}
};

class TransitionLogicException : public HMMException {
public:
	template<typename T>
	TransitionLogicException(const T& t, const std::string& msg) : 
		HMMException(error_message::format("TransitionLogicException: " + msg, t)) {}
};

/* <----------------------------> */

class HiddenMarkovModel {
	/* Name of this hmm. */
	std::string _name;

	/* Begin and end states. */
	State* _begin;
	State* _end;

	/* Holds the states and the transitions when building the hmm. */
	Graph<State> _graph;	

	/* Generated via brew() */
	RawModel* _model;

	/* Algorithms */
	ForwardAlgorithm* _forward_algorithm;
	BackwardAlgorithm* _backward_algorithm;
	DecodingAlgorithm* _decoding_algorithm;
	TrainingAlgorithm* _training_algorithm;

	/* Helper method. Used by train() to update the HMM values (i.e. its graph and PDFs) from the RawModel. */
	void _update_from_raw();

public:
	/* Constructors */
	HiddenMarkovModel();
	HiddenMarkovModel(const std::string& name);
	HiddenMarkovModel(
		const std::string& name, 
		const ForwardAlgorithm& forward, const BackwardAlgorithm& backward, 
		const DecodingAlgorithm& decode, const TrainingAlgorithm& train);
	HiddenMarkovModel(const HiddenMarkovModel& other);
	HiddenMarkovModel(HiddenMarkovModel&& other);
	HiddenMarkovModel& operator=(const HiddenMarkovModel& other);
	HiddenMarkovModel& operator=(HiddenMarkovModel&& other);

	/* Name getter / setter */
	std::string name() const;
	void set_name(const std::string& name);

	/* Returns the hmm graph */
	Graph<State> get_graph();

	/* Get the number of states of the hmm */
	std::size_t num_states() const;
	/* Get the number of transitions of the hmm */
	std::size_t num_transitions() const;

	/* True if the hmm contains the given state */
	bool has_state(const State& state) const;
	/* True if the hmm contains a transition from from_state to to_state*/
	bool has_transition(const State& from_state, const State& to_state) const;

	/* Returns the begin state of the hmm. Throws an exception if it is null (when removed). */
	State& begin();
	/* Returns the end state of the hmm. Throws an exception if it is null (when removed). */
	State& end();
	/* Returns the state with the given state name. Note that this worke by simply giving a string. */
	State& get_state(const State& state);
	/* Add the given state to the hmm if it is not yet contained by the hmm. Throws an exception if the 
	hmm already has the given state (by name (string) comparison between states) */
	void add_state(const State& state);
	/* Removes the given state of the hmm and all its transitions. Throws an exception if state not found. */
	void remove_state(const State& state);
	
	/* Adds a transition with given probability into this hmm. Throws an exception if one of the given state 
	is not contained by the hmm or if the given probability is < 0. */
	void add_transition(const State& from, const State& to, double probability);

	/* Returns the probability of the transition from the state from to the state to. Throws an exception 
	if transition not found or the probability was set to null (which should never happen) */
	double get_transition(const State& from, const State& to);

	/* Update the transition probability between from and to states. Throws if the transition is not found. */
	void set_transition(const State& from, const State& to, double probability);

	/* Sets a transition from the begin silent state and the given state with the given probability. Throws if state not
	contained by the hmm. */
	void begin_transition(const State& state, double probability);

	/* Sets a transition from the given state to the end silent state with the given probability. Throws if the state
	is not contained by the hmm. */
	void end_transition(const State& state, double probability);

	/* Remove a transition if it is contained by the graph else throws an exception. */
	void remove_transition(const State& from, const State& to);

	/* Prepares the hmm before calling algorithms on it. This method MUST always be called before any
	algorithm is called. This is because all the hmm algorithms use the RawModel attribute that is initialized
	by this method. If normalize is true then the probabilies will be normalized, ex if State s has 2 out transitions
	witch each a probability of 2, each is set to have a probability of 0.5. */
	void brew(bool normalize = true);

	/* RawModel getters */
	Matrix raw_transitions();
	std::vector<double> raw_pi_begin();
	std::vector<double> raw_pi_end();
	std::vector<Distribution*> raw_pdfs();
	std::map<std::string, std::size_t> states_indices();
	std::vector<std::string> states_names();

	/* Setters for the algorithms. */
	void set_forward(const ForwardAlgorithm& forward);
	void set_backward(const BackwardAlgorithm& backward);
	void set_decoding(const DecodingAlgorithm& decode);
	void set_training(const TrainingAlgorithm& training);

	/* Get the type for each algorithm. */
	std::string forward_type() const;
	std::string backward_type() const;
	std::string decoding_type() const;
	std::string training_type() const;

	/* Interface calling the algorithms */
	
	/* Calls the forward algorithm on given sequence. t_max is the t at which the forward 
	recursion will stop. If t is set to 0 or is greater than the sequence length, the recursion 
	will go all the way to t=T, the sequence length. */
	std::vector<double> forward(const std::vector<std::string>& sequence, std::size_t t_max = 0);

	/* Calls the backward algorithm on given sequence. t_min is the t at which the backward 
	recursion will stop. If t is set to 0, the recursion will go to t=1 (complete). If t
	is greater than the sequence length only the backward initialisation will take place. */
	std::vector<double> backward(const std::vector<std::string>& sequence, std::size_t t_min = 0);

	/* Returns the log likelihood by using the forward algorithm if do_fwd is true, else 
	the backward algorithm is used. */
	double log_likelihood(const std::vector<std::string>& sequence, bool do_fwd = true);
	double log_likelihood(const std::vector<std::vector<std::string>>& sequences, bool do_fwd = true);
	double likelihood(const std::vector<std::string>& sequence, bool do_fwd = true);
	double likelihood(const std::vector<std::vector<std::string>>& sequences, bool do_fwd = true);

	/* Calls the decoding algorithm on the given sequence by executing the recursion until t_max. 
	If t_max is set to 0 or is greater than the sequence length, the recursion will be complete.
	Retunrs the optimal state path and its likelihood. */
	std::pair<std::vector<std::string>, double> decode(const std::vector<std::string>& sequence, std::size_t t_max = 0);

	/* Calls the training algorithm on the given set of training sequences. Return the obtained improvement. */
	double train(const std::vector<std::vector<std::string>>& sequences,
		double transition_pseudocount = hmm_config::kDefaultTransitionPseudocount,
		double convergence_threshold = hmm_config::kDefaultConvergenceThreshold,
		unsigned int min_iterations = hmm_config::kDefaultMinIterations, 
		unsigned int max_iterations = hmm_config::kDefaultMaxIterations);

	/* IO operations */
	/* Save the hmm. The file name is the HMM name with the default hmm extension. */
	void save();
	/* Saves the hmm in the given filename with given extension. */
	void save(const std::string& filename, const std::string& extension = global_config::kDefaultFileExtension);
	/* Sets this HMM to the values of the HMM contained in the given filegit . */
	void load(const std::string& filename, const std::string& extension = global_config::kDefaultFileExtension);

	/* Printing */
	std::string transition_string(const State& from, const State& to) const;
	void print_transitions(bool log_prob = true);
	void print_distributions(bool log_prob = true);
	void print_pi_begin(bool log_prob = true);
	void print_pi_end(bool log_prob = true);

	virtual ~HiddenMarkovModel();
};

/* <-------- PRINT UTILS --------> */

void __print_transitions(const Matrix& matrix, const std::map<std::string, std::size_t>& indices, bool log_prob = true);
void __print_distributions(const std::vector<Distribution*>& dists, const std::vector<std::string>& names, bool log_prob=true);
void __print_distributions(std::vector<DiscreteDistribution>& dists, const std::vector<std::string>& names, bool log_prob = true);
void __print_names(const std::vector<std::size_t>& ids, const std::vector<std::string>& names);
void __print_pi_begin(const std::vector<double>& pi, const std::vector<std::string>& names, bool log_prob=true);
void __print_pi_end(const std::vector<double>& pi, const std::vector<std::string>& names, bool log_prob=true);
void print_prob(const std::vector<double>& probs, bool log_prob=true);
std::ostream& operator<<(std::ostream& out, const std::vector<double>& vec);
std::ostream& operator<<(std::ostream& out, const std::vector<std::size_t>& vec);
std::ostream& operator<<(std::ostream& out, const std::vector<std::string>& vec);
std::ostream& operator<<(std::ostream& out, const std::vector<State>& vec);
std::ostream& operator<<(std::ostream& out, const std::vector<State*>& vec);

/* <----------------------------> */

#endif