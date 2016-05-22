#ifndef __DISTRIBUTION_HPP
#define __DISTRIBUTION_HPP

#include <string>
#include <algorithm>
#include <map>
#include <numeric>
#include <functional>
#include <fstream>
#include "constants.hpp"
#include "utils.hpp"

#include <string>
#include <algorithm>
#include <map>
#include <numeric>
#include <functional>
#include <fstream>
#include "constants.hpp"
#include "distributions.hpp"
#include "utils.hpp"


/* <-------- Exceptions --------> */

class DistributionException : public std::logic_error {
protected:
	DistributionException(const std::string&);
};

class DistributionSymbolNotFoundException : public DistributionException {
public:
	DistributionSymbolNotFoundException(const std::string&);
};

/* <----------------------------> */

class Distribution {
private:
	std::string _name;
	bool _log;

protected:
	Distribution();
	Distribution(const std::string&);
	Distribution(const Distribution&);

public:
	virtual std::string name() const;
	virtual bool is_discrete() const;
	virtual bool is_continuous() const;
	virtual bool uses_log_probabilities() const;
	virtual void log_probabilities(bool use_log);
	virtual std::string to_string() const;

	/* Pure virtual methods */
	/* Get probabilities with operator[] */
	virtual bool empty() const = 0;
	virtual void log_normalize() = 0;
	virtual double& operator[] (const std::string&) = 0;
	// virtual double& operator[] (double) = 0;
	virtual bool operator==(const Distribution& other) const = 0;
	virtual bool operator!=(const Distribution& other) const = 0;
	virtual void save(std::ofstream&) = 0;
	virtual void load(std::ifstream&) = 0;
	/* For polymorphic use */
	virtual Distribution* clone() const = 0;
	virtual ~Distribution();
};

std::ostream& operator<<(std::ostream&, const Distribution&);

class DiscreteDistribution : public Distribution{
private:
	std::map<std::string, double> _distribution;
protected:
	DiscreteDistribution(const std::string&);
public:
	DiscreteDistribution();
	DiscreteDistribution(std::initializer_list<std::pair<const std::string, double>>);
	DiscreteDistribution(const DiscreteDistribution&);

	/* Covariant return type */
	virtual DiscreteDistribution* clone() const;

	void round(int = global_config::kDoublePrecision);
	double prob_sum() const;
	std::string to_string() const;
	void save(std::ofstream& out);
	void load(std::ifstream& in);
	void log_probabilities(bool use_log);
	void log_normalize();
	bool empty() const;
	bool is_discrete() const;
	bool contains(const std::string&) const;
	bool contains(const std::string&);
	std::vector<std::string> symbols() const;

	double& operator[] (const std::string& symbol);
	double& operator[] (double symbol);
	bool operator==(const DiscreteDistribution& other) const;
	bool operator==(const Distribution& other) const;
	bool operator!=(const Distribution& other) const;

	virtual ~DiscreteDistribution();
};


class ContinuousDistribution : public Distribution {
protected:
	ContinuousDistribution(const std::string&);

public:
	ContinuousDistribution();

	bool is_continuous() const;
	virtual bool is_normal() const;
	virtual bool is_uniform() const;
	virtual ~ContinuousDistribution();
};

class NormalDistribution : public ContinuousDistribution {
public:
	NormalDistribution();
	bool is_normal() const;
	virtual ~NormalDistribution();
};

class UniformDistribution : public ContinuousDistribution {
public:
	UniformDistribution();
	bool is_uniform() const;
	virtual ~UniformDistribution();
};

#endif