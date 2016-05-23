#ifndef __HMMALGORITHMS_HPP
#define __HMMALGORITHMS_HPP

#include <vector>
#include <string>
#include <utility>
#include "state.hpp"
#include "distributions.hpp"
#include "hmm_base.hpp"

/* ===================== BASE CLASSES ===================== */

class HMMAlgorithm {
	std::string _name;
protected:
	RawModel* _model;
	HMMAlgorithm(const std::string&, RawModel*);
public:
	std::string name() const;
	std::string type() const;
	virtual void set_model(RawModel*);
	virtual HMMAlgorithm* clone() const = 0;
	virtual ~HMMAlgorithm();
};

class ForwardAlgorithm : public HMMAlgorithm {
protected:
	ForwardAlgorithm(const std::string&, RawModel*);
public:
	virtual ForwardAlgorithm* clone() const = 0;
	virtual std::vector<double> forward(const std::vector<std::string>&, std::size_t) = 0;
	virtual double log_likelihood(const std::vector<std::string>&) = 0;
	virtual double log_likelihood(const std::vector<std::vector<std::string>>&) = 0;

	virtual ~ForwardAlgorithm();
};

class BackwardAlgorithm : public HMMAlgorithm {
protected:
	BackwardAlgorithm(const std::string&, RawModel*);
public:
	virtual BackwardAlgorithm* clone() const = 0;
	virtual std::vector<double> backward(const std::vector<std::string>&, std::size_t) = 0;
	virtual double log_likelihood(const std::vector<std::string>&) = 0;
	virtual double log_likelihood(const std::vector<std::vector<std::string>>&) = 0;
	virtual ~BackwardAlgorithm();
};

class DecodingAlgorithm : public HMMAlgorithm {
protected:
	DecodingAlgorithm(const std::string&, RawModel*);
public:
	virtual DecodingAlgorithm* clone() const = 0;
	virtual std::pair<std::vector<std::string>, double> decode(const std::vector<std::string>&, std::size_t) = 0;
	virtual ~DecodingAlgorithm();
};

class TrainingAlgorithm : public HMMAlgorithm {
protected:
	TrainingAlgorithm(const std::string&, RawModel*);
public:
	virtual TrainingAlgorithm* clone() const = 0;
	virtual double train(const std::vector<std::vector<std::string>>&, double, double, unsigned int, unsigned int) = 0;
	virtual ~TrainingAlgorithm();
};

/* ===================== LINEAR MEMORY FORWARD ===================== */

class LinearMemoryForwardAlgorithm : public ForwardAlgorithm {
public:
	LinearMemoryForwardAlgorithm(RawModel*);
	LinearMemoryForwardAlgorithm* clone() const;

	std::vector<double> forward(const std::vector<std::string>&, std::size_t);
	double log_likelihood(const std::vector<std::string>&);
	double log_likelihood(const std::vector<std::vector<std::string>>&);

	std::vector<double> forward_init(const std::vector<std::string>&);
	std::vector<double> forward_step(const std::vector<std::string>&, const std::vector<double>&, std::size_t t);
	std::pair<std::vector<double>, double> forward_terminate(const std::vector<double>&);

	virtual ~LinearMemoryForwardAlgorithm();
};

/* ===================== LINEAR MEMORY BACKWARD ===================== */

class LinearMemoryBackwardAlgorithm : public BackwardAlgorithm {
public:
	LinearMemoryBackwardAlgorithm(RawModel*);
	LinearMemoryBackwardAlgorithm* clone() const;

	std::vector<double> backward(const std::vector<std::string>& sequence, std::size_t);
	double log_likelihood(const std::vector<std::string>&);
	double log_likelihood(const std::vector<std::vector<std::string>>&);

	std::vector<double> backward_init();
	std::vector<double> backward_step(const std::vector<double>&, const std::vector<std::string>&, std::size_t);
	std::tuple<std::vector<double>, std::vector<double>, double> backward_terminate(const std::vector<double>&, const std::vector<std::string>&);

	virtual ~LinearMemoryBackwardAlgorithm();
};

/* ===================== LINEAR MEMORY VITERBI DECODE ===================== */

class LinearMemoryViterbiDecodingAlgorithm : public DecodingAlgorithm{
public:
	class Traceback {
		struct Node; // forward declaration.
		typedef std::shared_ptr<Node> NodePtr;
		struct Node{
			NodePtr previous;
			std::size_t value;
			Node(std::size_t v) : previous(), value(v) {}
			Node(std::size_t v, NodePtr p) : previous(p), value(v) {}
			void set_previous(const NodePtr& p){ previous = p; }
		};
		std::size_t _nodes;
		std::vector<NodePtr> _previous_nodes;
		std::vector<NodePtr> _current_nodes;

		void _init_previous();
		void _init_current();

	public:
		Traceback(std::size_t);
		void add_link(std::size_t, std::size_t, bool = false);
		void next_column();
		void reset();
		std::vector<std::size_t> from(std::size_t);
		std::string to_string() const;

	};

	LinearMemoryViterbiDecodingAlgorithm(RawModel*);
	LinearMemoryViterbiDecodingAlgorithm* clone() const;

	std::pair<std::vector<std::string>, double> decode(const std::vector<std::string>&, std::size_t);

	std::vector<double> viterbi_init(Traceback&, const std::vector<std::string>&);
	std::vector<double> viterbi_step(const std::vector<double>&, Traceback&, std::size_t, const std::vector<std::string>&);
	std::size_t viterbi_terminate(std::vector<double>&);

	virtual ~LinearMemoryViterbiDecodingAlgorithm();
};

/* ===================== LINEAR MEMORY TRAINING ===================== */

class LinearMemoryTrainingAlgorithm : public TrainingAlgorithm {
protected:
	LinearMemoryTrainingAlgorithm(const std::string&, RawModel*);

	/* Test wether a transition from i to j occurs in the given traceback.
		This method should only return 0 or 1. */
	static unsigned int any_of_transitions(const std::vector<std::size_t>&, std::size_t, std::size_t);

	class TransitionScore{
		std::vector<std::vector<double>> _transitions_scores;
		std::vector<std::vector<double>> _pi_begin_scores;
		std::vector<std::vector<double>> _pi_end_scores;
		const std::vector<std::pair<std::size_t, std::size_t>>* _free_transitions;
		const std::vector<std::size_t>* _free_pi_begin;
		const std::vector<std::size_t>* _free_pi_end;
		double _default_score;
	public:

		/* ===================== TRANSITION SCORE ===================== */
		TransitionScore(const std::vector<std::pair<std::size_t, std::size_t>>&,
			const std::vector<std::size_t>&, const std::vector<std::size_t>&, std::size_t, double = 0.0);
		
		TransitionScore& operator=(const TransitionScore&);

		void add(const TransitionScore&, std::size_t, std::size_t);
		/* Returns the transitions score of given transition for a path finishing at state m. */
		double score(std::size_t, std::size_t) const;
		double score_begin(std::size_t, std::size_t) const;
		double score_end(std::size_t, std::size_t) const;
		std::size_t num_free_transitions() const;
		std::size_t num_free_begin_transitions() const;
		std::size_t num_free_end_transitions() const;

		void set_begin_score(std::size_t, std::size_t, double);
		void set_score(std::size_t, std::size_t, double);
		void set_end_score(std::size_t, std::size_t, double);

		void copy_begin(const TransitionScore&, std::size_t, std::size_t);
		std::size_t get_from_state_id(std::size_t) const;
		std::size_t get_to_state_id(std::size_t) const;
		std::size_t get_state_id_from_begin(std::size_t) const;
		std::size_t get_state_id_to_end(std::size_t) const;
		void reset(double);
		void reset();
		std::string to_string(std::size_t, const std::vector<std::string>&, const std::string& = "", bool = true) const;
	};

	/* ===================== EMISSION SCORE ===================== */
	class EmissionScore{
		std::vector<std::vector<double>> _emissions_scores;
		const std::vector<std::pair<std::size_t, std::string>>* _free_emissions;
		double _default_score;
	public:
		EmissionScore(const std::vector<std::pair<std::size_t, std::string>>&, std::size_t, double = 0.0);
		EmissionScore& operator=(const EmissionScore&);
		std::size_t get_state_id(std::size_t) const;
		std::string get_symbol(std::size_t) const;
		double score(std::size_t, std::size_t) const;
		void set_score(std::size_t, std::size_t, double);
		std::size_t num_free_emissions() const;

		/* Adds the scores for arriving at state m of other EmissionScore to the scores of arriving 
		at state 0 of this EmissionScore. Both scores should have the same sizes. */
		void add(const EmissionScore&, std::size_t, std::size_t);
		void reset(double reset_score);
		void reset();

		std::string to_string(std::size_t, const std::vector<std::string>&, const std::string& = "", bool = true) const;
	};

	static unsigned int delta(std::size_t, std::size_t);
	static unsigned int delta(std::string, std::string);
	static double log_score(std::string, std::string);
	static double log_delta(std::size_t i, std::size_t j);

	void print_scores(const TransitionScore& score, std::string from_str = "", bool from_all = true, std::size_t from = 0, bool log_prob = true);
	void print_scores(const EmissionScore& score, std::string from_str = "", bool from_all = true, std::size_t from = 0, bool log_prob = true);
	void print_total_scores(const TransitionScore& score, bool log_prob = true);
	void print_total_scores(const EmissionScore& score, bool log_prob = true);
	void print_all_scores(const TransitionScore& score, bool log_prob = true);
	void print_all_scores(const EmissionScore& score, bool log_prob = true);

	std::size_t last_non_silent_state(const std::vector<std::size_t>&);
	void update_emissions(const EmissionScore&, EmissionScore&, const std::vector<std::size_t>&, std::string);
	/* Updates all the Tij counts for paths finishing at m by using the previous Tij 
		counts for paths finishing at l and increments it if i == l and j == m. */
	void update(const TransitionScore&, TransitionScore&, const std::vector<std::size_t>&);
	void update_begin(TransitionScore&, const std::vector<std::size_t>&);
	/* Adds 1 to the end transition count of m for path arriving at m. */
	void update_end(TransitionScore&, std::size_t);


public:
	virtual LinearMemoryTrainingAlgorithm* clone() const = 0;
	virtual ~LinearMemoryTrainingAlgorithm();
};

/* ===================== LINEAR MEMORY VITERBI TRAINING ===================== */

class LinearMemoryViterbiTraining : public LinearMemoryTrainingAlgorithm{
	LinearMemoryViterbiDecodingAlgorithm _decoding_algorithm;
	LinearMemoryForwardAlgorithm _forward_algorithm;
public:
	LinearMemoryViterbiTraining(RawModel*);
	LinearMemoryViterbiTraining* clone() const;
	virtual void set_model(RawModel*);
	double train(const std::vector<std::vector<std::string>>& sequences, double transition_pseudocount,
		double convergence_threshold, unsigned int min_iterations, unsigned int max_iterations);

	void update_model_from_scores(const TransitionScore&, const EmissionScore&, double);
	void update_model_transitions_from_scores(const TransitionScore&, double);
	void update_model_emissions_from_scores(const EmissionScore&); 

	virtual ~LinearMemoryViterbiTraining();
};

/* ===================== LINEAR MEMORY BAUM-WELCH TRAINING ===================== */

class LinearMemoryBaumWelchTraining : public LinearMemoryTrainingAlgorithm{
	LinearMemoryBackwardAlgorithm _backward_algorithm;
public:
	LinearMemoryBaumWelchTraining(RawModel*);
	LinearMemoryBaumWelchTraining* clone() const;
	virtual void set_model(RawModel*);

	void update_model_from_log_scores(const TransitionScore&, const EmissionScore&);
	void update_model_transitions_from_log_scores(const TransitionScore&);
	void update_model_emissions_from_log_scores(const EmissionScore&);

	double train(const std::vector<std::vector<std::string>>& sequences, double transition_pseudocount, 
		double convergence_threshold, unsigned int min_iterations, unsigned int max_iterations);

	void log_update_transition_score(const TransitionScore&, TransitionScore&, double);
	void log_update_emission_score(const EmissionScore&, EmissionScore&, double);

	virtual ~LinearMemoryBaumWelchTraining();
};

#endif