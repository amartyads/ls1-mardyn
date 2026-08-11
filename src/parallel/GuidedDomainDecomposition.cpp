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
	_guidedDDList.print();
	if (numWeights < 1 || !initWeightsGiven) {
		// default behaviour, create default config
		_gridSize = getOptimalGrid(_domainLength, this->getNumProcs());
		_coords = getCoordsFromRank(_gridSize, _rank);
		std::tie(_boxMin, _boxMax) = initializeRegularGrid(_domainLength, _gridSize, _coords);
	}
	else {
		auto initWeights = _guidedDDList.getCurrentDD();
		mardyn_assert(initWeights.timestep == 0);
		for (int i = 0; i < 3; i++) {
			_gridSize[i] = static_cast<int>(initWeights.numRanksInDim(i));
		}
		_coords = getCoordsFromRank(_gridSize, _rank);
		std::tie(_boxMin, _boxMax) = getBoxBounds(initWeights, _domainLength, _coords);
	}
	for (int i = 0; i < 3; i++) {
		Log::global_log->set_mpi_output_all();
		Log::global_log->info() << "_boxMin[" << i << "] :" << _boxMin[i] << std::endl;
		Log::global_log->info() << "_boxMax[" << i << "] :" << _boxMax[i] << std::endl;
		Log::global_log->set_mpi_output_root(0);
	}
	xmlconfig.changecurrentnode(oldpath);
}

std::tuple<std::array<double, 3>, std::array<double, 3>> GuidedDomainDecomposition::getBoxBounds(const StaticDDAtTime& staticDD,
	const std::array<double, 3>& domainLength, 
	const std::array<size_t, 3>& gridCoords) {
	std::array<double, 3> boxMin, boxMax;
	for (int i = 0; i < 3; i++) {
		const auto backWeight =
			std::reduce(staticDD.subdomainWeights[i].begin(),
						staticDD.subdomainWeights[i].begin() + gridCoords[i], 0u);
		const auto totalWeight =
			std::reduce(staticDD.subdomainWeights[i].begin() + gridCoords[i],
						staticDD.subdomainWeights[i].end(), backWeight);

		// calculate box bounds from cumulative weights of previous ranks, and the
		// weight of the current rank
		boxMin[i] =
			static_cast<double>(backWeight) * domainLength[i] / totalWeight;
		boxMax[i] =
			boxMin[i] + (static_cast<double>(staticDD.subdomainWeights[i][gridCoords[i]]) *
						domainLength[i] / totalWeight);
	}
	return std::make_tuple(boxMin, boxMax);
}

void GuidedDomainDecomposition::balanceAndExchange(double lastTraversalTime, bool forceRebalancing,
												   ParticleContainer* moleculeContainer, Domain* domain) {
	MARDYN_EXIT("NOP");
}

void GuidedDomainDecomposition::rebalance(std::array<double, 3> newBoxMin, std::array<double, 3> newboxmax) {
	MARDYN_EXIT("NOP");
}
