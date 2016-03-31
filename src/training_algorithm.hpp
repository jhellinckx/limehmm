#ifndef __TRAININGALGORITHM_HPP
#define __TRAININGALGORITHM_HPP

#include <vector>
#include "model.hpp"
#include "state.hpp"
#include "symbol.hpp"

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