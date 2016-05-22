#include <string>
#include <algorithm>
#include <map>
#include <numeric>
#include <functional>
#include <vector>
#include <fstream>
#include "constants.hpp"
#include "distributions.hpp"
#include "utils.hpp"

/* <-------- Exceptions --------> */
DistributionException::DistributionException(const std::string& msg) : std::logic_error(msg) {}

DistributionSymbolNotFoundException::DistributionSymbolNotFoundException(const std::string& t) : 
		DistributionException(error_message::format("DistributionSymbolNotFoundException: " + error_message::kDistributionSymbolNotFound, t)) {}
/* <----------------------------> */

std::ostream& operator<<(std::ostream& out, const Distribution& dist){
	out << dist.to_string();
	return out;
}


Distribution::Distribution() : 
	Distribution(distribution_config::kDistributionName) {}

Distribution::Distribution(const std::string& name) : _name(name), _log(distribution_config::kDefaultLogUse) {}
Distribution::Distribution(const Distribution& other) : _name(other._name), _log(other._log) {}


std::string Distribution::name() const { return _name; }
bool Distribution::is_discrete() const { return false; }
bool Distribution::is_continuous() const { return false; }
bool Distribution::uses_log_probabilities() const { return _log; }
void Distribution::log_probabilities(bool use_log) { _log = use_log; }

std::string Distribution::to_string() const {
	return _name;
}

Distribution::~Distribution() {}


DiscreteDistribution::DiscreteDistribution(const std::string& name) :
	Distribution(name), _distribution() {}

DiscreteDistribution::DiscreteDistribution() : 
	DiscreteDistribution(distribution_config::kDiscreteDistributionName) {}

DiscreteDistribution::DiscreteDistribution(std::initializer_list<std::pair<const std::string, double>> distribution) : 
	Distribution(distribution_config::kDiscreteDistributionName), _distribution(distribution) {}

DiscreteDistribution::DiscreteDistribution(const DiscreteDistribution& other) :
	Distribution(other), _distribution(other._distribution) {}

DiscreteDistribution* DiscreteDistribution::clone() const {
	return new DiscreteDistribution(*this);
}

void DiscreteDistribution::round(int precision){
	for(auto& entry : _distribution){
		entry.second = utils::round_double(entry.second, precision);
	}
}

double DiscreteDistribution::prob_sum() const {
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

std::string DiscreteDistribution::to_string() const {
	std::string str = Distribution::to_string() + ": ";
	std::for_each(_distribution.begin(), _distribution.end(), 
	[&str](const std::pair<const std::string,double>& entry){
		str += entry.first + "(" + std::to_string(entry.second) + ") ";
	});
	str += "-> sum " + std::to_string(prob_sum());
	return str;
}

void DiscreteDistribution::save(std::ofstream& out) {
	bool used_log = uses_log_probabilities();
	log_probabilities(false);
	out << _distribution.size() << std::endl;
	for(auto& entry : _distribution){
		out << entry.first << '>' << entry.second << std::endl;
	}

	log_probabilities(used_log);
}

void DiscreteDistribution::load(std::ifstream& in) {
	log_probabilities(false);
	std::string line;
	std::string symbol;
	std::string prob_str;
	std::size_t num_symbols;
	std::getline(in, line);
	num_symbols = (std::size_t) stoi(line);
	for(std::size_t i = 0; i < num_symbols; ++i){
		std::getline(in, line);
		std::tie(symbol, prob_str) = utils::split_first(line, global_config::kProbabilitySeparator);
		_distribution[symbol] = std::stod(prob_str);
	}
}

void DiscreteDistribution::log_probabilities(bool use_log) {
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

void DiscreteDistribution::log_normalize() {
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

bool DiscreteDistribution::empty() const {
	return _distribution.size() == 0 || prob_sum() == 0;
}

bool DiscreteDistribution::is_discrete() const { return true; }

bool DiscreteDistribution::contains(const std::string& symbol) const {
	return _distribution.find(symbol) != _distribution.end();
}
bool DiscreteDistribution::contains(const std::string& symbol) {
	return _distribution.find(symbol) != _distribution.end();
}

std::vector<std::string> DiscreteDistribution::symbols() const {
	std::vector<std::string> keys;
	keys.reserve(_distribution.size());
	for(auto entry : _distribution){
		keys.push_back(entry.first);
	}
	return keys;
}

double& DiscreteDistribution::operator[] (const std::string& symbol) {
	if(!contains(symbol)){
		_distribution[symbol] = (this->uses_log_probabilities()) ? utils::kNegInf : 0.0;
	}
	return _distribution[symbol];
}

double& DiscreteDistribution::operator[] (double symbol) {
	return operator[](std::to_string(symbol));
}

bool DiscreteDistribution::operator==(const DiscreteDistribution& other) const {
	return (uses_log_probabilities() == other.uses_log_probabilities()) && (other._distribution == _distribution);
}

bool DiscreteDistribution::operator==(const Distribution& other) const {
	if(other.is_discrete()){
		return operator==((DiscreteDistribution&) other);
	}
	return false;
}

bool DiscreteDistribution::operator!=(const Distribution& other) const{
	return ! operator==(other);
}

DiscreteDistribution::~DiscreteDistribution() {}



ContinuousDistribution::ContinuousDistribution(const std::string& name) : 
	Distribution(name) {}

ContinuousDistribution::ContinuousDistribution() : 
	ContinuousDistribution(distribution_config::kContinuousDistributionName) {}

bool ContinuousDistribution::is_continuous() const { return true; }
bool ContinuousDistribution::is_normal() const { return false; }
bool ContinuousDistribution::is_uniform() const { return false; }
ContinuousDistribution::~ContinuousDistribution() {}

NormalDistribution::NormalDistribution() : ContinuousDistribution(distribution_config::kNormalDistributionName) {}
bool NormalDistribution::is_normal() const { return true; }
NormalDistribution::~NormalDistribution() {}

UniformDistribution::UniformDistribution() : ContinuousDistribution(distribution_config::kUniformDistributionName) {}
bool UniformDistribution::is_uniform() const { return true; }
UniformDistribution::~UniformDistribution() {}
