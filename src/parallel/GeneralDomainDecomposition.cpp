/**
 * @file GeneralDomainDecomposition.cpp
 * @author seckler
 * @date 11.04.19
 */

#include "GeneralDomainDecomposition.h"
#ifdef ENABLE_ALLLBL
#include "ALLLoadBalancer.h"
#endif
#include "Domain.h"
#include "NeighborAcquirer.h"
#include "NeighbourCommunicationScheme.h"

#include "utils/String_utils.h"
#include "utils/mardyn_assert.h"

#include <numeric>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <tuple>
#include <sstream>

GeneralDomainDecomposition::GeneralDomainDecomposition(double interactionLength, Domain* domain, bool forceGrid)
	: DomainDecompMutable( {0.}, {0.}, {domain->getGlobalLength(0), domain->getGlobalLength(1), domain->getGlobalLength(2)}),
	  _interactionLength{interactionLength},
	  _forceLatchingToLinkedCellsGrid{forceGrid} {}

void GeneralDomainDecomposition::initializeALL() {
	Log::global_log->info() << "initializing ALL load balancer..." << std::endl;
	auto gridSize = getOptimalGrid(_domainLength, this->getNumProcs());
	auto gridCoords = getCoordsFromRank(gridSize, _rank);
	Log::global_log->info() << "gridSize:" << gridSize[0] << ", " << gridSize[1] << ", " << gridSize[2] << std::endl;
	Log::global_log->info() << "gridCoords:" << gridCoords[0] << ", " << gridCoords[1] << ", " << gridCoords[2] << std::endl;
	std::tie(_boxMin, _boxMax) = initializeRegularGrid(_domainLength, gridSize, gridCoords);
	if (_forceLatchingToLinkedCellsGrid and not _gridSize.has_value()) {
		std::array<double, 3> forcedGridSize{};
		for(size_t dim = 0; dim < 3; ++dim){
			size_t numCells = _domainLength[dim] / _interactionLength;
			forcedGridSize[dim] = _domainLength[dim] / numCells;
		}
		_gridSize = forcedGridSize;
	}
	if (_gridSize.has_value()) {
		std::tie(_boxMin, _boxMax) = latchToGridSize(_boxMin, _boxMax);
	}
#ifdef ENABLE_ALLLBL
	// Increased slightly to prevent rounding errors.
	const double safetyFactor = 1. + 1.e-10;
	const std::array<double, 3> minimalDomainSize =
		_gridSize.has_value() ? *_gridSize
							  : std::array{_interactionLength * safetyFactor, _interactionLength * safetyFactor,
										   _interactionLength * safetyFactor};

	_loadBalancer = std::make_unique<ALLLoadBalancer>(_boxMin, _boxMax, 4 /*gamma*/, this->getCommunicator(), gridSize,
													  gridCoords, minimalDomainSize);
#else
	std::ostringstream error_message;
	error_message << "ALL load balancing library not enabled. Aborting." << std::endl;
	MARDYN_EXIT(error_message.str());
#endif
	Log::global_log->info() << "GeneralDomainDecomposition initial box: [" << _boxMin[0] << ", " << _boxMax[0] << "] x ["
					   << _boxMin[1] << ", " << _boxMax[1] << "] x [" << _boxMin[2] << ", " << _boxMax[2] << "]"
					   << std::endl;
}

GeneralDomainDecomposition::~GeneralDomainDecomposition() = default;

bool GeneralDomainDecomposition::queryRebalancing(size_t step, size_t updateFrequency, size_t initPhase,
												  size_t initUpdateFrequency, double /*lastTraversalTime*/) {
	return step <= initPhase ? step % initUpdateFrequency == 0 : step % updateFrequency == 0;
}

void GeneralDomainDecomposition::balanceAndExchange(double lastTraversalTime, bool forceRebalancing,
													ParticleContainer* moleculeContainer, Domain* domain) {
	const bool rebalance =
		queryRebalancing(_steps, _rebuildFrequency, _initPhase, _initFrequency, lastTraversalTime) or forceRebalancing;
	if (_steps == 0) {
		// ensure that there are no outer particles
		moleculeContainer->deleteOuterParticles();
		// init communication partners
		initCommPartners(moleculeContainer, domain);
		DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, HALO_COPIES);
	} else {
		if (rebalance) {
			// first transfer leaving particles
			DomainDecompMPIBase::exchangeMoleculesMPI(moleculeContainer, domain, LEAVING_ONLY);

			// ensure that there are no outer particles
			moleculeContainer->deleteOuterParticles();

			// rebalance
			Log::global_log->info() << "rebalancing..." << std::endl;

			Log::global_log->set_mpi_output_all();
			Log::global_log->debug() << "work:" << lastTraversalTime << std::endl;
			Log::global_log->set_mpi_output_root(0);
			auto [newBoxMin, newBoxMax] = _loadBalancer->rebalance(lastTraversalTime);
			if (_gridSize.has_value()) {
				std::tie(newBoxMin, newBoxMax) = latchToGridSize(newBoxMin, newBoxMax);
			}
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
		_boundaryHandler.setLocalRegion(_boxMin.data(),_boxMax.data());
		_boundaryHandler.updateGlobalWallLookupTable();
	}
	++_steps;
}

void GeneralDomainDecomposition::initCommPartners(ParticleContainer* moleculeContainer,
												  Domain* domain) {  // init communication partners
	auto coversWholeDomain = _loadBalancer->getCoversWholeDomain();
	for (int d = 0; d < DIMgeom; ++d) {
		// this needs to be updated for proper initialization of the neighbours
		_neighbourCommunicationScheme->setCoverWholeDomain(d, coversWholeDomain[d]);
	}
	_neighbourCommunicationScheme->initCommunicationPartners(moleculeContainer->getCutoff(), domain, this,
															 moleculeContainer);
}

void GeneralDomainDecomposition::readXML(XMLfileUnits& xmlconfig) {
	// Ensures that the readXML() call to DomainDecompMPIBase forces the direct-pp communication scheme.
	_forceDirectPP = true;

	DomainDecompMPIBase::readXML(xmlconfig);

#ifdef MARDYN_AUTOPAS
	Log::global_log->info() << "AutoPas only supports FS, so setting it." << std::endl;
	setCommunicationScheme("direct-pp", "fs");
#endif

	xmlconfig.getNodeValue("updateFrequency", _rebuildFrequency);
	Log::global_log->info() << "GeneralDomainDecomposition update frequency: " << _rebuildFrequency << std::endl;

	xmlconfig.getNodeValue("initialPhaseTime", _initPhase);
	Log::global_log->info() << "GeneralDomainDecomposition time for initial rebalancing phase: " << _initPhase << std::endl;

	xmlconfig.getNodeValue("initialPhaseFrequency", _initFrequency);
	Log::global_log->info() << "GeneralDomainDecomposition frequency for initial rebalancing phase: " << _initFrequency
					   << std::endl;

	std::string gridSizeString;
	if (xmlconfig.getNodeValue("gridSize", gridSizeString)) {
		Log::global_log->info() << "GeneralDomainDecomposition grid size: " << gridSizeString << std::endl;

		if (gridSizeString.find(',') != std::string::npos) {
			auto strings = string_utils::split(gridSizeString, ',');
			if (strings.size() != 3) {
				std::ostringstream error_message;
				error_message
					<< "GeneralDomainDecomposition's gridSize should have three entries if a list is given, but has "
					<< strings.size() << "!" << std::endl;
				MARDYN_EXIT(error_message.str());
			}
			_gridSize = {std::stod(strings[0]), std::stod(strings[1]), std::stod(strings[2])};
		} else {
			double gridSize = std::stod(gridSizeString);
			_gridSize = {gridSize, gridSize, gridSize};
		}
		for (auto gridSize : *_gridSize) {
			if (gridSize < _interactionLength) {
				std::ostringstream error_message;
				error_message << "GeneralDomainDecomposition's gridSize (" << gridSize
									<< ") is smaller than the interactionLength (" << _interactionLength
									<< "). This is forbidden, as it leads to errors! " << std::endl;
				MARDYN_EXIT(error_message.str());
			}
		}
	}

	if (xmlconfig.changecurrentnode("loadBalancer")) {
		std::string loadBalancerString = "None";
		xmlconfig.getNodeValue("@type", loadBalancerString);
		Log::global_log->info() << "Chosen Load Balancer: " << loadBalancerString << std::endl;

		std::transform(loadBalancerString.begin(), loadBalancerString.end(), loadBalancerString.begin(), ::tolower);

		if (loadBalancerString.find("all") != std::string::npos) {
			initializeALL();
		} else {
			std::ostringstream error_message;
			error_message << "GeneralDomainDecomposition: Unknown load balancer " << loadBalancerString
								<< ". Aborting! Please select a valid option! Valid options: ALL";
			MARDYN_EXIT(error_message.str());
		}
		_loadBalancer->readXML(xmlconfig);
	} else {
		std::ostringstream error_message;
		error_message << "loadBalancer section missing! Aborting!" << std::endl;
		MARDYN_EXIT(error_message.str());
	}
	xmlconfig.changecurrentnode("..");
}
