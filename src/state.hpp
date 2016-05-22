#ifndef __STATE_HPP
#define __STATE_HPP

#include <vector>
#include <exception>
#include "constants.hpp"
#include "distributions.hpp"

/* <-------- Exceptions --------> */

class StateException : public std::logic_error {
protected:
	StateException(const std::string&);
};

class StateDistributionException : public StateException {
public:
	StateDistributionException(const std::string&);
};

/* <----------------------------> */

class State{
private:
	/* State name */
	std::string _name;
	/* Emission probabilities */
	Distribution* _distribution;
	/* Free emission */
	bool _free_emission;
	/* Free transition */
	bool _free_transition;

public:
	State(const std::string&);
	State(const char*);
	State(const std::string&, const Distribution&);
	State(const State&);
	State(State&&);
	State& operator=(const State&);
	State& operator=(State&&);
	bool operator==(const State&) const;
	bool operator!=(const State&) const;
	std::string to_string() const;
	bool has_free_emission() const;
	bool has_free_transition() const;
	void fix_emission();
	void fix_transition();
	void free_emission();
	void free_transition();
	Distribution& distribution() const;
	void set_distribution(const Distribution&);
	std::string name() const;
	void set_name(const std::string&);
	bool is_silent() const;
	virtual ~State();
};

std::ostream& operator<<(std::ostream& out, const State& state);

#endif