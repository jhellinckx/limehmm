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
	extern const std::string kGetVertexNotFound;
	extern const std::string kGetEdgeNotFound;
	extern const std::string kGetOutEdgesVertexNotFound;
	extern const std::string kGetInEdgesVertexNotFound;
	extern const std::string kRemoveVertexNotFound;
	extern const std::string kVertexNotFound;
	extern const std::string kRemoveEdgeNotFound;
	extern const std::string kEdgeNotFound;
	extern const std::string kIncidentVertexNotFound;
	extern const std::string kAddedVertexExists;
	extern const std::string kAddedEdgeExists;

	/* State */
	extern const std::string kSilentStateHasNoDistribution;

	/* Distribution */
	extern const std::string kDistributionSymbolNotFound;

	/* HMM */
	extern const std::string kHMMGetStateNotFound;
	extern const std::string kHMMGetTransitionNotFound;
	extern const std::string kHMMGetTransitionNullWeight;
	extern const std::string kHMMHasNoBeginState;
	extern const std::string kHMMHasNoEndState;
	extern const std::string kHMMRemoveStateNotFound;
	extern const std::string kHMMRemoveTransitionNotFound;
	extern const std::string kHMMAddStateExists;
	extern const std::string kHMMAddTransitionExists;
	extern const std::string kHMMAddTransitionStateNotFound;
	extern const std::string kHMMAddedTransitionFromEndState;
	extern const std::string kHMMAddedTransitionToBeginState;
	extern const std::string kHMMTransitionNegativeProbability;

	template<typename T>
	static std::string format(const std::string& error, const T& t) {
		std::ostringstream oss;
		oss << t << " : " << error << std::endl;
		return oss.str();
	}
}

namespace global_config{
	extern const int kDoublePrecision;

	extern const char kProbabilitySeparator;
	extern const std::string kDefaultFileExtension;
	extern const std::string kNullValue;
}

namespace hmm_config {
	extern const bool kDefaultFreeEmission;
	extern const bool kDefaultFreeTransition;
	extern const int kDefaultPiBegin;
	extern const int kDefaultPiEnd;
	extern const int kDefaultTransitionProbability;
	extern const int kDefaultEmissionProbability;

	extern const double kDefaultTransitionPseudocount;
	extern const double kDefaultConvergenceThreshold;
	extern const unsigned int kDefaultMaxIterationsViterbi;
	extern const unsigned int kDefaultMinIterationsViterbi;

	extern const std::string kDefaultHMMName;
	extern const std::string kDefaultStartStateLabel;
	extern const std::string kDefaultEndStateLabel;

	extern const int kDefaultStateLabelCountStart;
	extern const std::string kDefaultStateLabelString;

	extern const std::string kLinearMemoryForwardAlgorithmName;
	extern const std::string kLinearMemoryBackwardAlgorithmName;
	extern const std::string kLinearMemoryViterbiDecodeAlgorithmName;
	extern const std::string kLinearMemoryViterbiTrainingAlgorithmName;
	extern const std::string kLinearMemoryBaumWelchTrainingAlgorithmName;
}

namespace distribution_config {
	extern const std::string kDistributionName;
	extern const std::string kDiscreteDistributionName;
	extern const std::string kContinuousDistributionName;
	extern const std::string kNormalDistributionName;
	extern const std::string kUniformDistributionName;

	extern const bool kDefaultLogUse;
}

#endif