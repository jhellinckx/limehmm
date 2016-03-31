#ifndef __HMMALGORITHMS_HPP
#define __HMMALGORITHMS_HPP

#include <vector>
#include "model.hpp"
#include "state.hpp"
#include "symbol.hpp"

class DecodingAlgorithm{
public:
	virtual std::vector<std::string> decode(const std::vector<Symbol>&, Model&) = 0;
};

class LinearViterbiDecoding : public DecodingAlgorithm{
	virtual std::vector<std::string> decode(const std::vector<Symbol>& decoding_set, Model& model){

	}

};

class TrainingAlgorithm{
public:
	virtual void trainModel(const std::vector<Symbol>&, Model&) = 0;
};

class LinearViterbiTraining : public TrainingAlgorithm{
	virtual void trainModel(const std::vector<Symbol>& training_set, Model& trained_model){

	}

};

class LinearBaumWelchTraining : public TrainingAlgorithm{
	virtual void trainModel(const std::vector<Symbol>& training_set, Model& trained_model){

	}
};
	
class LinearEMTraining : public TrainingAlgorithm{
	virtual void trainModel(const std::vector<Symbol>& training_set, Model& trained_model){

	}
};

#endif