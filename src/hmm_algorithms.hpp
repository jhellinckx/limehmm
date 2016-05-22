#ifndef __HMMALGORITHMS_HPP
#define __HMMALGORITHMS_HPP

#include <vector>
#include <string>
#include "state.hpp"
#include "distributions.hpp"


class HMMAlgorithm {
	std::string _name;
protected:
	HMMAlgorithm(const std::string& name) : _name(name) {}
public:
	std::string name() const { return _name; }
	std::string type() const { return name(); }

	virtual HMMAlgorithm* clone() const = 0;

	virtual ~HMMAlgorithm() {}
};

class ForwardAlgorithm : public HMMAlgorithm {
protected:
	ForwardAlgorithm(const std::string& name) : HMMAlgorithm(name) {}
public:
	virtual ForwardAlgorithm* clone() const = 0;

	virtual std::vector<double> forward(const std::vector<std::string>&) = 0;

	virtual ~ForwardAlgorithm() {}
};

class LinearMemoryForwardAlgorithm : public ForwardAlgorithm {
public:
	LinearMemoryForwardAlgorithm() : ForwardAlgorithm(hmm_config::kLinearMemoryForwardAlgorithmName) {}

	LinearMemoryForwardAlgorithm* clone() const { return new LinearMemoryForwardAlgorithm(*this); }

	std::vector<double> forward(const std::vector<std::string>& sequence);

	virtual ~LinearMemoryForwardAlgorithm() {}
};

class BackwardAlgorithm : public HMMAlgorithm {
protected:
	BackwardAlgorithm(const std::string& name) : HMMAlgorithm(name) {}
public:
	virtual BackwardAlgorithm* clone() const = 0;

	virtual std::vector<double> backward(const std::vector<std::string>&) = 0;

	virtual ~BackwardAlgorithm() {}
};

class LinearMemoryBackwardAlgorithm : public BackwardAlgorithm {
public:
	LinearMemoryBackwardAlgorithm() : BackwardAlgorithm(hmm_config::kLinearMemoryBackwardAlgorithmName) {}

	LinearMemoryBackwardAlgorithm* clone() const { return new LinearMemoryBackwardAlgorithm(*this); }

	std::vector<double> backward(const std::vector<std::string>& sequence){}

	virtual ~LinearMemoryBackwardAlgorithm() {}
};

class DecodingAlgorithm : public HMMAlgorithm {
protected:
	DecodingAlgorithm(const std::string& name) : HMMAlgorithm(name) {}
public:
	virtual DecodingAlgorithm* clone() const = 0;

	virtual std::vector<std::string> decode(const std::vector<std::string>&) = 0;

	virtual ~DecodingAlgorithm() {}
};

class LinearMemoryViterbiDecodingAlgorithm : public DecodingAlgorithm{
public:
	LinearMemoryViterbiDecodingAlgorithm() : DecodingAlgorithm(hmm_config::kLinearMemoryViterbiDecodeAlgorithmName) {}
	
	LinearMemoryViterbiDecodingAlgorithm* clone() const { return new LinearMemoryViterbiDecodingAlgorithm(*this); }

	std::vector<std::string> decode(const std::vector<std::string>& sequence){

	}

	virtual ~LinearMemoryViterbiDecodingAlgorithm() {}

};

class TrainingAlgorithm : public HMMAlgorithm {
protected:
	TrainingAlgorithm(const std::string& name) : HMMAlgorithm(name) {}
public:
	virtual TrainingAlgorithm* clone() const = 0;

	virtual void train(const std::vector<std::vector<std::string>>&) = 0;

	virtual ~TrainingAlgorithm() {}
};

class LinearMemoryTrainingAlgorithm : public TrainingAlgorithm {
protected:
	LinearMemoryTrainingAlgorithm(const std::string& name) : TrainingAlgorithm(name) {}
public:
	virtual LinearMemoryTrainingAlgorithm* clone() const = 0;

	virtual ~LinearMemoryTrainingAlgorithm() {}
};

class LinearMemoryViterbiTraining : public LinearMemoryTrainingAlgorithm{
public:
	LinearMemoryViterbiTraining() : LinearMemoryTrainingAlgorithm(hmm_config::kLinearMemoryViterbiTrainingAlgorithmName){}

	LinearMemoryViterbiTraining* clone() const { return new LinearMemoryViterbiTraining(*this); }

	void train(const std::vector<std::vector<std::string>>&){

	}

	virtual ~LinearMemoryViterbiTraining() {}
};

class LinearMemoryBaumWelchTraining : public LinearMemoryTrainingAlgorithm{
	LinearMemoryBaumWelchTraining() : LinearMemoryTrainingAlgorithm(hmm_config::kLinearMemoryBaumWelchTrainingAlgorithmName){}

	LinearMemoryBaumWelchTraining* clone() const { return new LinearMemoryBaumWelchTraining(*this); }
	
	void train(const std::vector<std::vector<std::string>>&){

	}

	virtual ~LinearMemoryBaumWelchTraining() {}
};

#endif