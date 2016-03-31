#ifndef __DISTRIBUTION_HPP
#define __DISTRIBUTION_HPP

#include <string>
#include <algorithm>
#include <map>
#include "constants.hpp"

class Distribution {
private:
	std::string _name;

public:
	Distribution() : 
		Distribution(distribution_config::kDistributionName) {}
	Distribution(const std::string& name) : _name(name) {}

	virtual std::string name() const { return this->_name; }
	virtual bool is_discrete() const { return false; }
	virtual bool is_continuous() const { return false; }
	
	/* Pure virtual methods */
	/* Get probabilities with operator[] */
	virtual double& operator[] (const std::string&) = 0;
	virtual double& operator[] (const double&) = 0;
	virtual bool operator==(const Distribution& other) const = 0;
	virtual void prepare() = 0;
	/* For polymorphic use */
	virtual Distribution* clone() const = 0;
	virtual ~Distribution() {}
};

class DiscreteDistribution : public Distribution, public std::map<std::string, double> {
private:

public:
	DiscreteDistribution() : 
		DiscreteDistribution(distribution_config::kDiscreteDistributionName) {}

	DiscreteDistribution(const std::string& name) :
		Distribution(name), std::map<std::string, double>() {}

	DiscreteDistribution(std::initializer_list<std::pair<const std::string, double>> distribution) : 
		Distribution(distribution_config::kDiscreteDistributionName), std::map<std::string, double>(distribution) {}

	DiscreteDistribution(const DiscreteDistribution& other) :
		Distribution(other.name()), std::map<std::string, double>(other) {}

	/* Covariant return type */
	virtual DiscreteDistribution* clone() const {
		return new DiscreteDistribution(*this);
	}

	bool is_discrete() const { return true; }
	bool contains(const std::string& symbol) const {
		return this->find(symbol) != this->end();
	}
	bool contains(const std::string& symbol) {
		return this->find(symbol) != this->end();
	}

	virtual double& operator[] (const std::string& symbol){
		return std::map<std::string, double>::operator[](symbol);
	}

	virtual double& operator[] (const double& symbol){
		return std::map<std::string, double>::operator[](std::to_string(symbol));
	}

	virtual bool operator==(const DiscreteDistribution& other) const {
		return (other.name() == this->name()) && ((std::map<std::string, double>&)other == ((std::map<std::string, double>&)*this));
	}

	virtual bool operator==(const Distribution& other) const {
		if(other.is_discrete()){
			return operator==((DiscreteDistribution&) other);
		}
		return false;
	}

	void prepare(){

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