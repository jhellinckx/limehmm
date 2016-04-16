#ifndef __CONSTANTS_HPP
#define __CONSTANTS_HPP

#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <typeinfo>
#include <type_traits>

namespace error_message {
	/* Graph */
	extern const std::string kGetVertexNotFound = "tried to get a vertex but it was not found in the graph";
	extern const std::string kGetEdgeNotFound = "tried to get an edge but it was not found in the graph";
	extern const std::string kGetOutEdgesVertexNotFound = "tried to get out edges for vertex but vertex was not found in the graph";
	extern const std::string kGetInEdgesVertexNotFound = "tried to get in edges for vertex but vertex was not found in the graph";
	extern const std::string kRemoveVertexNotFound = "tried to remove a vertex but it was not found in the graph";
	extern const std::string kVertexNotFound = "vertex was not found in graph";
	extern const std::string kRemoveEdgeNotFound = "tried to remove an edge but it was not found in the graph";
	extern const std::string kEdgeNotFound = "edge was not found in graph";
	extern const std::string kIncidentVertexNotFound = "tried to add an edge but one of its incident vertex was not found in the graph";
	extern const std::string kAddedVertexExists = "tried to add a vertex but an equal vertex was found in the graph";
	extern const std::string kAddedEdgeExists = "tried to add an edge but an equal edge was found in the graph";

	/* State */
	extern const std::string kSilentStateHasNoDistribution = "tried to get the emission probability of a silent state; but a silent state has no distribution";

	/* Distribution */
	extern const std::string kDistributionSymbolNotFound = "symbol not found in distribution";

	/* HMM */
	extern const std::string kHMMHasNoBeginState = "no begin state was found, maybe it has been removed ?";
	extern const std::string kHMMHasNoEndState = "no end state was found, maybe it has been removed ?";
	extern const std::string kHMMRemoveStateNotFound = "tried to remove a state not contained by the hmm";
	extern const std::string kHMMRemoveTransitionNotFound = "tried to remove a transition not contained by the hmm";
	extern const std::string kHMMAddStateExists = "tried to add a state already contained by the hmm";
	extern const std::string kHMMAddTransitionExists = "tried to add a transition already contained by the hmm";
	extern const std::string kAddTransitionStateNotFound = "tried to add a transition with a state not contained by the hmm";
	extern const std::string kAddedTransitionFromEndState = "tried to add a transition from an end state";
	extern const std::string kAddedTransitionToBeginState = "tried to add a transition to a begin state";
	extern const std::string kAddedTransitionNegativeProbability = "tried to add a transition with a negative probability";

	template<typename T>
	extern std::string format(const std::string& error, const T& t) {
		std::ostringstream oss;
		oss << t << " : " << error << std::endl;
		return oss.str();
	}
}

namespace global_config{
	extern const int kDoublePrecision = 8;	
}

namespace hmm_config {
	extern const bool kFreeEmission = true;
	extern const bool kFreeTransition = true;
	extern const int kPiBegin = 0;
	extern const int kPiEnd = 0;
	extern const int kTransitionProbability = 0;
	extern const int kEmissionProbability = 0;
	extern const std::string kStartStateLabel = "begin_state";
	extern const std::string kEndStateLabel = "end_state";

	extern const double kPseudocount = 1;
	extern const double kConvergenceThreshold = 1e-9;
	extern const unsigned int kMaxIterationsViterbi = 1e8;
	extern const unsigned int kMinIterationsViterbi = 0;

	extern const int kAutoStateLabelCountStart = 1;
	extern const std::string kAutoStateLabelString = "state_";
}

namespace distribution_config {
	extern const std::string kDistributionName = "Distribution";
	extern const std::string kDiscreteDistributionName = "Discrete distribution";
	extern const std::string kContinuousDistributionName = "Continuous distribution";
	extern const std::string kNormalDistributionName = "Normal distribution";
	extern const std::string kUniformDistributionName = "Uniform distribution";

	extern const bool kDefaultLogUse = false;
}

#endif