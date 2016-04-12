#ifndef __DISTRIBUTION_HPP
#define __DISTRIBUTION_HPP

#include <string>
#include <algorithm>
#include <map>
#include <numeric>
#include <functional>
#include "constants.hpp"
#include "utils.hpp"


/* <-------- Exceptions --------> */

class DistributionException : public std::logic_error {
protected:
	DistributionException(const std::string& message) :
		std::logic_error(message) {}
};

class DistributionSymbolNotFoundException : public DistributionException {
public:
	template<typename T>
	DistributionSymbolNotFoundException(const T& t) : 
		DistributionException(error_message::format("DistributionSymbolNotFoundException: " + error_message::kDistributionSymbolNotFound, t)) {}
};

/* <----------------------------> */

class Distribution {
private:
	std::string _name;
	bool _log;
public:
	Distribution() : 
		Distribution(distribution_config::kDistributionName) {}
	Distribution(const std::string& name) : _name(name), _log(distribution_config::kDefaultLogUse) {}

	virtual std::string name() const { return _name; }
	virtual bool is_discrete() const { return false; }
	virtual bool is_continuous() const { return false; }
	virtual bool uses_log_probabilities() const { return _log; }
	virtual void log_probabilities(bool use_log) { _log = use_log; }
	/* Pure virtual methods */
	/* Get probabilities with operator[] */
	virtual bool empty() const = 0;
	virtual std::string to_string() const {
		return _name;
	}
	virtual void log_normalize() = 0;
	virtual double& operator[] (const std::string&) = 0;
	virtual double& operator[] (double) = 0;
	virtual bool operator==(const Distribution& other) const = 0;
	virtual bool operator!=(const Distribution& other) const = 0;
	/* For polymorphic use */
	virtual Distribution* clone() const = 0;
	virtual ~Distribution() {}
};

std::ostream& operator<<(std::ostream& out, const Distribution& dist){
	out << dist.to_string();
	return out;
}


/* This class is basically a wrapper around std::map... */
class DiscreteDistribution : public Distribution{
private:
	std::map<std::string, double> _distribution;
public:
	DiscreteDistribution() : 
		DiscreteDistribution(distribution_config::kDiscreteDistributionName) {}

	DiscreteDistribution(const std::string& name) :
		Distribution(name), _distribution() {}

	DiscreteDistribution(std::initializer_list<std::pair<const std::string, double>> distribution) : 
		Distribution(distribution_config::kDiscreteDistributionName), _distribution(distribution) {}

	DiscreteDistribution(const DiscreteDistribution& other) :
		Distribution(other.name()), _distribution(other._distribution) {}

	/* Covariant return type */
	virtual DiscreteDistribution* clone() const {
		return new DiscreteDistribution(*this);
	}

	double prob_sum() const {
		double init_sum;
		std::function<double(double, const std::pair<std::string, double>&)> prob_adder;
		if(this->uses_log_probabilities()){
			init_sum = utils::kNegInf;
			prob_adder = [](const double previous, const std::pair<const std::string, double>& entry) { 
							return utils::sum_log_prob(previous, entry.second);
						};
		}
		else{
			init_sum = double(0);
			prob_adder = [](const double previous, const std::pair<const std::string, double>& entry) { 
							return previous + entry.second;
						};
		}
		return std::accumulate(_distribution.begin(), _distribution.end(), init_sum, prob_adder);
	}

	std::string to_string() const {
		std::string str = Distribution::to_string() + ": ";
		if(uses_log_probabilities()){
			std::for_each(_distribution.begin(), _distribution.end(), 
				[&str](const std::pair<const std::string,double>& entry){
					str += entry.first + "(" + std::to_string(exp(entry.second)) + ") ";
				});
			str += "-> sum " + std::to_string(exp(prob_sum()));	
		}
		else{
			std::for_each(_distribution.begin(), _distribution.end(), 
			[&str](const std::pair<const std::string,double>& entry){
				str += entry.first + "(" + std::to_string(entry.second) + ") ";
			});
			str += "-> sum " + std::to_string(prob_sum());
		}
		
		return str;
	}

	virtual void log_probabilities(bool use_log) {
		/* Only make changes if distriubtion not yet using what is asked. */
		if(use_log != this->uses_log_probabilities()){
			Distribution::log_probabilities(use_log);
			if(use_log){
				std::for_each(_distribution.begin(), _distribution.end(), 
					[](std::pair<const std::string, double>& entry) {
						entry.second = log(entry.second);
					});
			}
			else{
				std::for_each(_distribution.begin(), _distribution.end(), 
					[](std::pair<const std::string, double>& entry) {
						entry.second = exp(entry.second);
					});
			}
		}
	}

	virtual void log_normalize() {
		if(! this->uses_log_probabilities()) {
			log_probabilities(true);
		}
		double probabilities_sum = prob_sum();
		if(exp(probabilities_sum) != 1.0){
			std::for_each(_distribution.begin(), _distribution.end(),
				[&probabilities_sum](std::pair<const std::string, double>& entry) {
					entry.second = utils::log_normalize(entry.second, probabilities_sum);
				});
		}
	}

	bool empty() const {
		return _distribution.size() == 0 || prob_sum() == 0;
	}

	bool is_discrete() const { return true; }

	bool contains(const std::string& symbol) const {
		return _distribution.find(symbol) != _distribution.end();
	}
	bool contains(const std::string& symbol) {
		return _distribution.find(symbol) != _distribution.end();
	}

	std::vector<std::string> symbols() const {
		std::vector<std::string> keys;
		keys.reserve(_distribution.size());
		for(auto entry : _distribution){
			keys.push_back(entry.first);
		}
		return keys;
	}

	virtual double& operator[] (const std::string& symbol) {
		if(!contains(symbol)){
			_distribution[symbol] = (this->uses_log_probabilities()) ? utils::kNegInf : 0.0;
		}
		return _distribution[symbol];
	}

	virtual double& operator[] (double symbol) {
		std::cout<<"call double"<<std::endl;
		return operator[](std::to_string(symbol));
	}

	virtual bool operator==(const DiscreteDistribution& other) const {
		return (other.name() == this->name()) && (other.uses_log_probabilities() == other.uses_log_probabilities()) && (other._distribution == _distribution);
	}

	virtual bool operator==(const Distribution& other) const {
		if(other.is_discrete()){
			return operator==((DiscreteDistribution&) other);
		}
		return false;
	}

	virtual bool operator!=(const Distribution& other) const{
		return ! operator==(other);
	}

	virtual ~DiscreteDistribution() {}
};


class ContinuousDistribution : public Distribution {
public:
	ContinuousDistribution() : 
		ContinuousDistribution(distribution_config::kContinuousDistributionName) {}
	ContinuousDistribution(const std::string& name) : 
		Distribution(name) {}

	bool is_continuous() const { return true; }
	virtual ~ContinuousDistribution() {}
};

class NormalDistribution : public ContinuousDistribution {
public:
	NormalDistribution() : ContinuousDistribution(distribution_config::kNormalDistributionName) {}
	virtual ~NormalDistribution() {}
};

class UniformDistribution : public ContinuousDistribution {
public:
	UniformDistribution() : ContinuousDistribution(distribution_config::kUniformDistributionName) {}
	virtual ~UniformDistribution() {}
};

#endif