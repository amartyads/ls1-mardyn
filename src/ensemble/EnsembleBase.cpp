#include "EnsembleBase.h"

#include "molecules/mixingrules/MixingRuleBase.h"
#include "molecules/mixingrules/LorentzBerthelot.h"
#include "utils/xmlfileUnits.h"
#include "utils/Logger.h"
#include "Simulation.h"
#include "Domain.h"

#include <vector>
#include <array>

using namespace std;
using Log::global_log;

Ensemble::~Ensemble() {
	delete _domain;
	_domain = nullptr;
	for(auto& m : _mixingrules)
		delete (m);
}

void Ensemble::readXML(XMLfileUnits& xmlconfig) {
	long numComponents = 0;
	XMLfile::Query query = xmlconfig.query("components/moleculetype");
	numComponents = query.card();
	global_log->info() << "Number of components: " << numComponents << endl;
	_components.resize(numComponents);
	XMLfile::Query::const_iterator componentIter;
	string oldpath = xmlconfig.getcurrentnodepath();
	for(componentIter = query.begin(); componentIter; componentIter++) {
		xmlconfig.changecurrentnode(componentIter);
		unsigned int cid = 0;
		xmlconfig.getNodeValue("@id", cid);
		_components[cid - 1].readXML(xmlconfig);
		_componentnamesToIds[_components[cid - 1].getName()] = cid - 1;
		global_log->debug() << _components[cid - 1].getName() << " --> " << cid - 1 << endl;
	}
	xmlconfig.changecurrentnode(oldpath);

	/* mixing rules */
	query = xmlconfig.query("components/mixing/rule");
	XMLfile::Query::const_iterator mixingruletIter;
	uint32_t numMixingrules = 0;
	numMixingrules = query.card();
	global_log->info() << "Found " << numMixingrules << " mixing rules." << endl;
	if(numMixingrules > 0) {
		global_log->info() << "Mixing rules are applied symmetrically: specifying a rule for AB also sets BA." << endl;
	}

	for(mixingruletIter = query.begin(); mixingruletIter; mixingruletIter++) {
		xmlconfig.changecurrentnode(mixingruletIter);
		MixingRuleBase* mixingrule = nullptr;
		string mixingruletype;

		xmlconfig.getNodeValue("@type", mixingruletype);
		global_log->info() << "Mixing rule type: " << mixingruletype << endl;
		if("LB" == mixingruletype) {
			mixingrule = new LorentzBerthelotMixingRule();

		} else {
			global_log->error() << "Unknown mixing rule " << mixingruletype << endl;
			global_log->error() << "Only rule type LB supported." << endl;
			Simulation::exit(1);
		}
		mixingrule->readXML(xmlconfig);
		_mixingrules.push_back(mixingrule);
	}

	checkMixingRules();
	setVectorOfMixingCoefficientsForComp2Param();

	xmlconfig.changecurrentnode(oldpath);
	setComponentLookUpIDs();
}


void Ensemble::checkMixingRuleCID(unsigned cid) const {
	int internalCID = static_cast<int>(cid) - 1;
	if (internalCID < 0) {
		global_log->error() << "Found mixing rule for component ID " << cid << ", which is not allowed." << endl;
		global_log->error() << "Component IDs start from 1." << endl;
		global_log->error() << "Aborting." << endl;
		Simulation::exit(2019);
	}

	int numComponents = _components.size();
	if (internalCID >= numComponents) {
		global_log->error() << "Found mixing rule for component ID " << cid << ", which is not allowed." << endl;
		global_log->error() << "Maximal component ID: " << numComponents << "." << endl;
		global_log->error() << "Aborting." << endl;
		Simulation::exit(2019);
	}
}

void Ensemble::checkMixingRules() const {
	using std::vector;

	int numMixingRules = _mixingrules.size();
	int numComponents = _components.size();

	//1 Comps: 0 coeffs; 2 Comps: 1 coeffs; 3 Comps: 3 coeffs; 4 Comps 6 coeffs : C * (C-1) / 2
	int exactNumMixingRules = numComponents * (numComponents - 1) / 2;

	if (numMixingRules < exactNumMixingRules) {
		global_log->info() << "Found " << numMixingRules << " mixing rules, which is smaller than exact number of mixing rules: " << exactNumMixingRules << "." << endl;
		global_log->info() << "The remaining mixing rules are populated with default values of (xi, eta) = (1.0, 1.0)." << endl;
	} else if (numMixingRules > exactNumMixingRules) {
		global_log->error() << "Found " << numMixingRules << " mixing rules, which is larger than exact number of mixing rules: " << exactNumMixingRules << "." << endl;
		global_log->error() << "Aborting." << endl;
		Simulation::exit(2019);
	}


	// check that:
	// II doesn't appear,
	// if IJ appears, then JI doesn't appear.
	vector<vector<int>> flags(numComponents, vector<int>(numComponents, 0));

	for(auto m = _mixingrules.begin(); m != _mixingrules.end(); ++m) {
		unsigned outerCID1 = (*m)->getCid1();
		unsigned outerCID2 = (*m)->getCid2();
		checkMixingRuleCID(outerCID1);
		checkMixingRuleCID(outerCID2);

		unsigned compID1 = outerCID1-1;
		unsigned compID2 = outerCID2-1;

		if (compID1 == compID2) {
			global_log->error() << "Specified rule for same component ids: " << compID1+1 << " and " << compID2+1 <<", which is not allowed." << endl;
			global_log->error() << "Aborting." << endl;
			Simulation::exit(2019);
		}

		flags[compID1][compID2] ++;
	}

	for (int cid1 = 0; cid1 < numComponents; ++cid1) {
		for (int cid2 = cid1 + 1; cid2 < numComponents; ++cid2) {
			int flag = flags[cid1][cid2];
			if (flag > 1) {
				global_log->error() << "Specified rule for component ids: " << cid1+1 << " and " << cid2+1 <<" more than once, which is not allowed." << endl;
				global_log->error() << "If you specify " << cid1+1 << " and " << cid2+1 <<", then you must not specify " << cid2+1 << " and " << cid1+1 << "." << endl;
				global_log->error() << "Aborting." << endl;
				Simulation::exit(2019);
			}
		}
	}
}

void Ensemble::setVectorOfMixingCoefficientsForComp2Param() const {
	using std::vector;
	using std::array;

	// data structure for mixing coefficients of domain class (still in use!!!)

	/*
	 * Mixing coefficients
	 *
	 * TODO: information of mixing rules (eta, xi) is stored in Domain class and its actually in use
	 * --> we need to decide where this information should be stored in future, in ensemble class,
	 * in the way it is done above?
	 *
	 */

	int numComponents = _components.size();

	std::vector<double>& dmixcoeff = global_simulation->getDomain()->getmixcoeff();
	dmixcoeff.clear();

	// two dimensional vector of Xi and Eta values
	// NOTE: initialise non-specified ones with default values of 1.0, 1.0
	vector<vector<array<double, 2>>> values(numComponents, vector<array<double, 2>>(numComponents,{1.0, 1.0}));

	for(auto m = _mixingrules.begin(); m != _mixingrules.end(); ++m) {
		// cast to LB mixing rule
		LorentzBerthelotMixingRule & rule = dynamic_cast<LorentzBerthelotMixingRule &>(**m);
		unsigned compID1 = rule.getCid1()-1;
		unsigned compID2 = rule.getCid2()-1;

		values[compID1][compID2][0] = rule.getEta();
		values[compID1][compID2][1] = rule.getXi();
	}

	// follow the precise way of initialising shit in Comp2Param::initialize
	// i.e. sort-of symmetrically initialised
	for (int cid1 = 0; cid1 < numComponents; ++cid1) {
		for (int cid2 = cid1 + 1; cid2 < numComponents; ++cid2) {
			double eta = values[cid1][cid2][0];
			double xi = values[cid1][cid2][1];
			dmixcoeff.push_back(eta);
			dmixcoeff.push_back(xi);
		}
	}
}

void Ensemble::setComponentLookUpIDs() {
	// Get the maximum Component ID.
	unsigned maxID = 0;
	for(auto c = _components.begin(); c != _components.end(); ++c)
		maxID = std::max(maxID, c->ID());

	// we need a look-up table for the Lennard-Jones centers in the Vectorized*CellProcessors
	// (for no other centers)
	unsigned centers = 0;
	for(auto c = _components.begin(); c != _components.end(); ++c) {
		c->setLookUpId(centers);
		centers += c->numLJcenters();
	}
}
