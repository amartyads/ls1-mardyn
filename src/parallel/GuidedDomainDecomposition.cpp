/*
 * GuidedDomainDecomposition.cpp
 *
 *  Created on: 10 August 2026
 *      Author: amartyads
 */

#include "GuidedDomainDecomposition.h"

#include "Domain.h"

GuidedDomainDecomposition::GuidedDomainDecomposition(Domain* domain)
	: GuidedDomainDecomposition(domain, MPI_COMM_WORLD) {}

GuidedDomainDecomposition::GuidedDomainDecomposition(Domain* domain, MPI_Comm comm)
	: DomainDecompMutable(comm, {0., 0., 0.}, {0., 0., 0.},
						  {domain->getGlobalLength(0), domain->getGlobalLength(1), domain->getGlobalLength(2)}) {}

void GuidedDomainDecomposition::readXML(XMLfileUnits& xmlconfig) {
	_forceDirectPP = true;
	DomainDecompMPIBase::readXML(xmlconfig);
#ifdef MARDYN_AUTOPAS
	Log::global_log->info() << "AutoPas only supports FS, so setting it." << std::endl;
	setCommunicationScheme("direct-pp", "fs");
#endif
	std::string oldpath = xmlconfig.getcurrentnodepath();
	XMLfile::Query query = xmlconfig.query("subdomainWeights");
	auto numWeights = query.card();
	bool initWeightsGiven = false;
	_guidedDDList.reserve(numWeights);
	for (auto subWeightsXML = query.begin(); subWeightsXML; ++subWeightsXML) {
		xmlconfig.changecurrentnode(subWeightsXML);
		StaticDDAtTime temp;
		temp.readXML(xmlconfig);
		if (temp.timestep < 0) {
			MARDYN_EXIT("NEG");
		}
		if (temp.timestep == 0)
			initWeightsGiven = true;
		_guidedDDList.insertWeights(temp);
	}
	if (!_guidedDDList.isValid(_numProcs)) {
		MARDYN_EXIT("INV");
	}
	_guidedDDList.reset();
	if (numWeights < 1 || !initWeightsGiven) {
		// default behaviour, create default config
		_gridSize = getOptimalGrid(_domainLength, this->getNumProcs());
		_coords = getCoordsFromRank(_gridSize, _rank);
		std::tie(_boxMin, _boxMax) = initializeRegularGrid(_domainLength, _gridSize, _coords);
	} else {
		auto initDD = _guidedDDList.getCurrentDD();
		mardyn_assert(initDD.timestep == 0);
		for (int i = 0; i < 3; i++) {
			_gridSize[i] = static_cast<int>(initDD.numRanksInDim(i));
		}
		_coords = getCoordsFromRank(_gridSize, _rank);
		std::tie(_boxMin, _boxMax) = getBoxBounds(initDD, _domainLength, _coords);
	}
	xmlconfig.changecurrentnode(oldpath);
}

std::tuple<std::array<double, 3>, std::array<double, 3>> GuidedDomainDecomposition::getBoxBounds(
	const StaticDDAtTime& staticDD, const std::array<double, 3>& domainLength,
	const std::array<size_t, 3>& gridCoords) {
	std::array<double, 3> boxMin, boxMax;
	for (int i = 0; i < 3; i++) {
		const auto backWeight =
			std::reduce(staticDD.subdomainWeights[i].begin(), staticDD.subdomainWeights[i].begin() + gridCoords[i], 0u);
		const auto totalWeight = std::reduce(staticDD.subdomainWeights[i].begin() + gridCoords[i],
											 staticDD.subdomainWeights[i].end(), backWeight);

		// calculate box bounds from cumulative weights of previous ranks, and the
		// weight of the current rank
		boxMin[i] = static_cast<double>(backWeight) * domainLength[i] / totalWeight;
		boxMax[i] = boxMin[i] +
					(static_cast<double>(staticDD.subdomainWeights[i][gridCoords[i]]) * domainLength[i] / totalWeight);
	}
	return std::make_tuple(boxMin, boxMax);
}

void GuidedDomainDecomposition::balanceAndExchange(double lastTraversalTime, bool forceRebalancing,
												   ParticleContainer* moleculeContainer, Domain* domain) {
	if (_steps == 0) {
		// ensure that there are no outer particles
		moleculeContainer->deleteOuterParticles();
		// init communication partners
		initCommPartners(moleculeContainer, domain);
		DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);
	} else {
		if (_steps == _guidedDDList.nextRebalancingTime()) {
			// rebalance
			_guidedDDList.goToNextRebalancing();
			auto curDD = _guidedDDList.getCurrentDD();
			for (int i = 0; i < 3; i++) {
				_gridSize[i] = static_cast<int>(curDD.numRanksInDim(i));
			}
			_coords = getCoordsFromRank(_gridSize, _rank);
			auto [newBoxMin, newBoxMax] = getBoxBounds(curDD, _domainLength, _coords);
			rebalance(moleculeContainer, domain, newBoxMin, newBoxMax);	 // sets _boxMin and _boxMax
		} else {
			if (sendLeavingWithCopies()) {
				DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_AND_HALO_COPIES);
			} else {
				DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_ONLY);
#ifndef MARDYN_AUTOPAS
				moleculeContainer->deleteOuterParticles();
#endif
				DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);
			}
		}
	}
	++_steps;
}

void GuidedDomainDecomposition::rebalance(ParticleContainer* moleculeContainer, Domain* domain,
										  std::array<double, 3> newBoxMin, std::array<double, 3> newBoxMax) {
	// first transfer leaving particles
	DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_ONLY);

	// ensure that there are no outer particles
	moleculeContainer->deleteOuterParticles();

	// rebalance
	Log::global_log->info() << "rebalancing..." << std::endl;
	// migrate the particles, this will rebuild the moleculeContainer!
	Log::global_log->info() << "migrating particles" << std::endl;
	migrateParticles(domain, moleculeContainer, newBoxMin, newBoxMax);

#ifndef MARDYN_AUTOPAS
	// The linked cells container needs this (I think just to set the cells to valid...)
	moleculeContainer->update();
#endif

	// set new boxMin and boxMax
	_boxMin = newBoxMin;
	_boxMax = newBoxMax;

	// init communication partners
	Log::global_log->info() << "updating communication partners" << std::endl;
	initCommPartners(moleculeContainer, domain);
	Log::global_log->info() << "rebalancing finished" << std::endl;
	DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);

	_boundaryHandler.setLocalRegion(_boxMin.data(), _boxMax.data());
	_boundaryHandler.updateGlobalWallLookupTable();
}
