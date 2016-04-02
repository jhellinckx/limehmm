#ifndef __HMMALGORITHMS_HPP
#define __HMMALGORITHMS_HPP

#include <vector>
#include "distributions.hpp"

class HMMAlgorithm{};

class ForwardAlgorithm : public HMMAlgorithm {

	template<typename Symbol>
	std::vector<double> forward(const std::vector<Symbol>& symbols){

	} 

};

class BackwardAlgorithm : public HMMAlgorithm {
	template<typename Symbol>
	std::vector<double> backward(const std::vector<Symbol>& symbols){
		
	} 
};

class DecodingAlgorithm : public HMMAlgorithm {
public:
	virtual std::vector<std::string> decode(const std::vector<Symbol>&, Model&) = 0;
};

class LinearViterbiDecoding : public DecodingAlgorithm{
	virtual std::vector<std::string> decode(const std::vector<Symbol>& decoding_set, Model& model){

	}

};

class TrainingAlgorithm : public HMMAlgorithm {
public:
	virtual void train(const std::vector<Symbol>&, Model&) = 0;
};

class LinearViterbiTraining : public TrainingAlgorithm{
	virtual void train(const std::vector<Symbol>& training_set, Model& trained_model){

	}

};

class LinearBaumWelchTraining : public TrainingAlgorithm{
	virtual void train(const std::vector<Symbol>& training_set, Model& trained_model){

	}
};
	
class LinearEMTraining : public TrainingAlgorithm{
	virtual void train(const std::vector<Symbol>& training_set, Model& trained_model){

	}
};

#endif