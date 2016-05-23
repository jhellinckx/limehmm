#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <typeinfo>
#include <type_traits>
#include "constants.hpp"

namespace error_message {
	/* Graph */
	const std::string kGetVertexNotFound = "tried to get a vertex but it was not found in the graph";
	const std::string kGetEdgeNotFound = "tried to get an edge but it was not found in the graph";
	const std::string kGetOutEdgesVertexNotFound = "tried to get out edges for vertex but vertex was not found in the graph";
	const std::string kGetInEdgesVertexNotFound = "tried to get in edges for vertex but vertex was not found in the graph";
	const std::string kRemoveVertexNotFound = "tried to remove a vertex but it was not found in the graph";
	const std::string kVertexNotFound = "vertex was not found in graph";
	const std::string kRemoveEdgeNotFound = "tried to remove an edge but it was not found in the graph";
	const std::string kEdgeNotFound = "edge was not found in graph";
	const std::string kIncidentVertexNotFound = "tried to add an edge but one of its incident vertex was not found in the graph";
	const std::string kAddedVertexExists = "tried to add a vertex but an equal vertex was found in the graph";
	const std::string kAddedEdgeExists = "tried to add an edge but an equal edge was found in the graph";

	/* State */
	const std::string kSilentStateHasNoDistribution = "tried to get the emission probability of a silent state; but a silent state has no distribution";

	/* Distribution */
	const std::string kDistributionSymbolNotFound = "symbol not found in distribution";

	/* HMM */
	const std::string kHMMGetStateNotFound = "tried to get a state not contained by the hmm";
	const std::string kHMMGetTransitionNotFound = "tried to get a transition not contained by the hmm";
	const std::string kHMMGetTransitionNullWeight = "tried to get a transition which exists but has a probability set to null";
	const std::string kHMMHasNoBeginState = "no begin state was found, maybe it has been removed ?";
	const std::string kHMMHasNoEndState = "no end state was found, maybe it has been removed ?";
	const std::string kHMMRemoveStateNotFound = "tried to remove a state not contained by the hmm";
	const std::string kHMMRemoveTransitionNotFound = "tried to remove a transition not contained by the hmm";
	const std::string kHMMAddStateExists = "tried to add a state already contained by the hmm";
	const std::string kHMMAddTransitionExists = "tried to add a transition already contained by the hmm";
	const std::string kHMMAddTransitionStateNotFound = "tried to add a transition with a state not contained by the hmm";
	const std::string kHMMAddedTransitionFromEndState = "tried to add a transition from an end state";
	const std::string kHMMAddedTransitionToBeginState = "tried to add a transition to a begin state";
	const std::string kHMMTransitionNegativeProbability = "tried to set a transition with a negative probability";

}

namespace global_config{
	const int kDoublePrecision = 8;	

	const char kProbabilitySeparator = '>';
	const std::string kDefaultFileExtension = "hmm";
	const std::string kNullValue = "null";
}

namespace hmm_config {
	const bool kDefaultFreeEmission = true;
	const bool kDefaultFreeTransition = true;
	const int kDefaultPiBegin = 0;
	const int kDefaultPiEnd = 0;
	const int kDefaultTransitionProbability = 0;
	const int kDefaultEmissionProbability = 0;

	const double kDefaultTransitionPseudocount = 0;
	const double kDefaultConvergenceThreshold = 1e-9;
	const unsigned int kDefaultMaxIterations = 1e8;
	const unsigned int kDefaultMinIterations = 0;

	const std::string kDefaultHMMName = "HiddenMarkovModel";
	const std::string kDefaultStartStateLabel = "begin_state";
	const std::string kDefaultEndStateLabel = "end_state";

	const int kDefaultStateLabelCountStart = 1;
	const std::string kDefaultStateLabelString = "state_";

	const std::string kLinearMemoryForwardAlgorithmName = "Linear Memory Forward";
	const std::string kLinearMemoryBackwardAlgorithmName = "Linear Memory Backward";
	const std::string kLinearMemoryViterbiDecodeAlgorithmName = "Linear Memory Viterbi Decode";
	const std::string kLinearMemoryViterbiTrainingAlgorithmName = "Linear Memory Viterbi Training";
	const std::string kLinearMemoryBaumWelchTrainingAlgorithmName = "Linear Memory Baum-Welch Training";
}

namespace distribution_config {
	const std::string kDistributionName = "Distribution";
	const std::string kDiscreteDistributionName = "Discrete distribution";
	const std::string kContinuousDistributionName = "Continuous distribution";
	const std::string kNormalDistributionName = "Normal distribution";
	const std::string kUniformDistributionName = "Uniform distribution";

	const bool kDefaultLogUse = false;
}