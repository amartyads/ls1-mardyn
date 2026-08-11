/*
 * GuidedDomainDecomposition.h
 *
 *  Created on: 11 August 2026
 *      Author: amartyads
 */

#include <array>
#include <sstream>

#include "utils/mardyn_assert.h"
#include "utils/xmlfileUnits.h"

struct StaticDDAtTime {
	unsigned int timestep = 0;
	std::array<std::vector<unsigned int>, 3> subdomainWeights{{{}, {}, {}}};

	void readXML(XMLfileUnits& xmlconfig) {
		xmlconfig.getNodeValue("@t", timestep);
		const std::array<std::string, 3> axes = {"x", "y", "z"};
		for (int i = 0; i < axes.size(); i++) {
			subdomainWeights[i].clear();
			const std::string weights = xmlconfig.getNodeValue_string(axes.at(i));
			if (!weights.empty()) {
				std::stringstream ss(weights);
				// Parse the weights, until the stringstream has chars and extraction
				// doesn't fail, and no EOF or linebreaks etc
				while (ss.good()) {
					int temp;
					// Extraction from stream into int type fails if token is not an int
					// We check for this failure, and additionally check for positive
					// integer
					if (!(ss >> temp) || temp <= 0) {
						std::ostringstream error_message;
						error_message << "Weights in " << axes.at(i)
									  << " axis have a non-natural number! Only integer weights > 0 allowed, please "
										 "check XML file!"
									  << std::endl;
						MARDYN_EXIT(error_message.str());
					}
					subdomainWeights[i].push_back(temp);
					if (ss.peek() == ',' || ss.peek() == ' ')  // skip commas and spaces
						ss.ignore();
				}
			} else {
				// No decomposition requested -> subdomain spans whole domain length
				subdomainWeights[i].push_back(1);
			}
		}
	}
};
