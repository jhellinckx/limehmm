#ifndef __HMMFACTORY_HPP
#define __HMMFACTORY_HPP

enum Algorithm {
	DECODING_LINEAR_VITERBI,
	TRAINING_LINEAR_VITERBI,
	TRAINING_LINEAR_BAUM_WELCH,
	TRAINING_LINEAR_EM
};

class HMMFactory{
private:
	HMMFactory(){}
public:
	static std::unique_ptr<Model> make(	const std::vector<S>& states,
									const std::vector<O>& observations,
									const std::vector<double>& transP,
									const std::vector<double>& emiP,
									const std::vector<double>& initP,
									int trainingAlgorithm){

		
	}
};

#endif