/*
 * GuidedDomainDecomposition.h
 *
 *  Created on: 10 August 2026
 *      Author: amartyads
 */

#pragma once
#include <algorithm>

#include "DomainDecompMutable.h"
#include "StaticDDAtTime.h"

class GuidedDDList {
public:
	bool hasMoreRebalancings() { return _currentPositionInList < static_cast<int>(_weightList.size()) - 1; }
	void goToNextRebalancing() { _currentPositionInList++; }
	int nextRebalancingTime() {
		if (hasMoreRebalancings())
			return _weightList[_currentPositionInList + 1].timestep;
		return -1;
	}
	StaticDDAtTime getCurrentDD() { return _weightList[_currentPositionInList]; }
	bool isValid(unsigned int numProcs) {
		if (_weightList.size() <= 0)
			return false;
		bool check = true;
		int prevTime = _weightList[0].timestep;
		for (size_t i = 0; i < _weightList.size(); i++) {
			check &= (_weightList[i].numRanksInDim(0) * _weightList[i].numRanksInDim(1) *
					  _weightList[i].numRanksInDim(2)) == numProcs;
			if (i != 0) {
				check &= _weightList[i].timestep > prevTime;
				prevTime = _weightList[i].timestep;
			}
		}
		return check;
	}
	void sortList() {
		std::sort(_weightList.begin(), _weightList.end(),
				  [](StaticDDAtTime a, StaticDDAtTime b) { return a.timestep < b.timestep; });
	}
	void insertWeights(StaticDDAtTime dd) { _weightList.push_back(dd); }
	void reset() { _currentPositionInList = 0; }
	void reserve(unsigned int s) { _weightList.reserve(s); }
	void print() {
		std::stringstream ss;
		for (size_t i = 0; i < _weightList.size(); i++) {
			ss << "Timestep " << _weightList[i].timestep << std::endl;
			for (int j = 0; j < 3; j++) {
				ss << "Weights for axis " << j << ": ";
				for (auto w : _weightList[i].subdomainWeights[j]) {
					ss << w << " ";
				}
				ss << std::endl;
			}
		}
		Log::global_log->info() << ss.str() << std::endl;
	}

private:
	std::vector<StaticDDAtTime> _weightList{};
	int _currentPositionInList = 0;
};

class GuidedDomainDecomposition : public DomainDecompMutable {
public:
	/**
	 * @brief Construct a new Guided Domain Decomposition object
	 *
	 */
	GuidedDomainDecomposition(Domain* domain);
	GuidedDomainDecomposition(Domain* domain, MPI_Comm comm);

	void readXML(XMLfileUnits& xmlconfig) override;

	static std::tuple<std::array<double, 3>, std::array<double, 3>> getBoxBounds(
		const StaticDDAtTime& staticDD, const std::array<double, 3>& domainLength,
		const std::array<size_t, 3>& gridCoords);

	void balanceAndExchange(double lastTraversalTime, bool forceRebalancing, ParticleContainer* moleculeContainer,
							Domain* domain) override;
	void rebalance(ParticleContainer* moleculeContainer, Domain* domain, std::array<double, 3> newBoxMin,
				   std::array<double, 3> newBoxMax);

	// returns a vector of the neighbour ranks in x y and z direction (only neighbours connected by an area to local
	// area)
	std::vector<int> getNeighbourRanks() override {
		throw std::runtime_error("GeneralDomainDecomposition::getNeighbourRanks() not yet implemented");
	}

	// returns a vector of all 26 neighbour ranks in x y and z direction
	std::vector<int> getNeighbourRanksFullShell() override {
		throw std::runtime_error("GeneralDomainDecomposition::getNeighbourRanksFullShell() not yet implemented");
	}

	// documentation in base class
	void prepareNonBlockingStage(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
								 unsigned int stageNumber) override {
		throw std::runtime_error("GeneralDomainDecomposition::prepareNonBlockingStage() not yet implemented");
	}

	// documentation in base class
	void finishNonBlockingStage(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
								unsigned int stageNumber) override {
		throw std::runtime_error("GeneralDomainDecomposition::prepareNonBlockingStage() not yet implemented");
	}

	// documentation in base class
	bool queryBalanceAndExchangeNonBlocking(bool forceRebalancing, ParticleContainer* moleculeContainer, Domain* domain,
											double etime) override {
		throw std::runtime_error(
			"GeneralDomainDecomposition::queryBalanceAndExchangeNonBlocking() not yet implemented");
	}

	std::vector<CommunicationPartner> getNeighboursFromHaloRegion(Domain* domain, const HaloRegion& haloRegion,
																  double cutoff) override {
		throw std::runtime_error("GeneralDomainDecomposition::getNeighboursFromHaloRegion() not yet implemented");
	}

private:
	GuidedDDList _guidedDDList;
	std::array<size_t, DIMgeom> _gridSize;	//!< Number of processes in each dimension of the MPI process grid
	std::array<size_t, DIMgeom> _coords;	//!< Coordinate of the process in the MPI process grid
};
