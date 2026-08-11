/*
 * DomainDecompMutable.cpp
 *
 *  Created on: 10 August 2026
 *      Author: amartyads, seckler
 */

#include "DomainDecompMutable.h"

#include "Domain.h"
#include "NeighborAcquirer.h"
#include "NeighbourCommunicationScheme.h"
#include "utils/Math.h"	 // isNearRel

DomainDecompMutable::DomainDecompMutable()
	: DomainDecompMutable(MPI_COMM_WORLD, {0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}) {}

DomainDecompMutable::DomainDecompMutable(MPI_Comm comm)
	: DomainDecompMutable(comm, {0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}) {}

DomainDecompMutable::DomainDecompMutable(std::array<double, 3> boxMin, std::array<double, 3> boxMax,
										 std::array<double, 3> domainLength)
	: DomainDecompMutable(MPI_COMM_WORLD, boxMin, boxMax, domainLength) {}

DomainDecompMutable::DomainDecompMutable(MPI_Comm comm, std::array<double, 3> boxMin, std::array<double, 3> boxMax,
										 std::array<double, 3> domainLength)
	: DomainDecompMPIBase(comm), _boxMin(boxMin), _boxMax(boxMax), _domainLength(domainLength) {}

DomainDecompMutable::~DomainDecompMutable() = default;

double DomainDecompMutable::getBoundingBoxMin(int dimension, Domain* /*domain*/) { return _boxMin[dimension]; }

double DomainDecompMutable::getBoundingBoxMax(int dimension, Domain* /*domain*/) { return _boxMax[dimension]; }

/**
 * Get the ordering of the input data.
 * The ordering will contain indices of elements of data, starting with the smallest going to the biggest.
 * e.g.:
 * returns {2, 0, 3, 1} for data={5, 16, 4, 7}
 * @tparam ArrayType should be some array type, e.g., vector of int
 * @param data the data to define the order
 * @return The ordering.
 */
template <typename ArrayType>
std::vector<size_t> getOrdering(const ArrayType& data) {
	std::vector<size_t> index(data.size());
	std::iota(index.begin(), index.end(), 0);
	std::sort(index.begin(), index.end(), [&](const size_t& a, const size_t& b) { return (data[a] < data[b]); });
	return index;
}

std::array<size_t, 3> DomainDecompMutable::getOptimalGrid(const std::array<double, 3>& domainLength, int numProcs) {
	// generate default grid
	std::array<int, 3> gridSize{0};
	MPI_CHECK(MPI_Dims_create(numProcs, 3, gridSize.data()));
	std::sort(gridSize.begin(), gridSize.end());

	// we want the grid to actually resemble the domain a bit (long sides -> many processes!)
	std::array<size_t, 3> grid{0};
	auto ordering = getOrdering(domainLength);
	for (size_t i = 0; i < 3; ++i) {
		grid[ordering[i]] = gridSize[i];
	}
	return grid;
}

std::array<size_t, 3> DomainDecompMutable::getCoordsFromRank(const std::array<size_t, 3>& gridSize, int rank) {
	auto yzSize = gridSize[1] * gridSize[2];
	auto zSize = gridSize[2];
	auto x = rank / yzSize;
	auto y = (rank - x * yzSize) / zSize;
	auto z = (rank - x * yzSize - y * zSize);
	return {x, y, z};
}

std::tuple<std::array<double, 3>, std::array<double, 3>> DomainDecompMutable::initializeRegularGrid(
	const std::array<double, 3>& domainLength, const std::array<size_t, 3>& gridSize,
	const std::array<size_t, 3>& gridCoords) {
	std::array<double, 3> boxMin{0.};
	std::array<double, 3> boxMax{0.};
	// initialize it as regular grid!
	for (size_t dim = 0; dim < 3; ++dim) {
		boxMin[dim] = gridCoords[dim] * domainLength[dim] / gridSize[dim];
		boxMax[dim] = (gridCoords[dim] + 1) * domainLength[dim] / gridSize[dim];
		if (gridCoords[dim] == gridSize[dim] - 1) {
			// ensure that the upper domain boundaries match.
			// lower domain boundaries always match, because they are 0.
			boxMax[dim] = domainLength[dim];
		}
	}
	return std::make_tuple(boxMin, boxMax);
}

void DomainDecompMutable::initCommPartners(ParticleContainer* moleculeContainer,
										   Domain* domain) {  // init communication partners
	bool coversWholeDomain;
	for (int d = 0; d < DIMgeom; ++d) {
		coversWholeDomain = isNearRel(_boxMin[d], 0.0) && isNearRel(_boxMax[d], _domainLength[d]);
		// this needs to be updated for proper initialization of the neighbours
		_neighbourCommunicationScheme->setCoverWholeDomain(d, coversWholeDomain);
	}
	_neighbourCommunicationScheme->initCommunicationPartners(moleculeContainer->getCutoff(), domain, this,
															 moleculeContainer);
}

void DomainDecompMutable::migrateParticles(Domain* domain, ParticleContainer* particleContainer,
										   std::array<double, 3> newMin, std::array<double, 3> newMax) {
	std::array<double, 3> oldBoxMin{particleContainer->getBoundingBoxMin(0), particleContainer->getBoundingBoxMin(1),
									particleContainer->getBoundingBoxMin(2)};
	std::array<double, 3> oldBoxMax{particleContainer->getBoundingBoxMax(0), particleContainer->getBoundingBoxMax(1),
									particleContainer->getBoundingBoxMax(2)};

	HaloRegion ownDomain{}, newDomain{};
	for (size_t i = 0; i < 3; ++i) {
		ownDomain.rmin[i] = oldBoxMin[i];
		newDomain.rmin[i] = newMin[i];
		ownDomain.rmax[i] = oldBoxMax[i];
		newDomain.rmax[i] = newMax[i];
		ownDomain.offset[i] = 0;
		newDomain.offset[i] = 0;
	}
	Log::global_log->set_mpi_output_all();
	Log::global_log->debug() << "migrating from"
							 << " [" << oldBoxMin[0] << ", " << oldBoxMax[0] << "] x"
							 << " [" << oldBoxMin[1] << ", " << oldBoxMax[1] << "] x"
							 << " [" << oldBoxMin[2] << ", " << oldBoxMax[2] << "] " << std::endl;
	Log::global_log->debug() << "to"
							 << " [" << newMin[0] << ", " << newMax[0] << "] x"
							 << " [" << newMin[1] << ", " << newMax[1] << "] x"
							 << " [" << newMin[2] << ", " << newMax[2] << "]." << std::endl;
	Log::global_log->set_mpi_output_root(0);
	std::vector<HaloRegion> desiredDomain{newDomain};
	std::vector<CommunicationPartner> sendNeighbors{}, recvNeighbors{};

	std::array<double, 3> globalDomainLength{domain->getGlobalLength(0), domain->getGlobalLength(1),
											 domain->getGlobalLength(2)};
	// 0. skin, as it is not needed for the migration of particles!
	std::tie(recvNeighbors, sendNeighbors) =
		NeighborAcquirer::acquireNeighbors(globalDomainLength, &ownDomain, desiredDomain, _comm);

	std::vector<Molecule> dummy;
	for (auto& sender : sendNeighbors) {
		sender.initSend(particleContainer, _comm, _mpiParticleType, LEAVING_ONLY, dummy,
						false /*don't use invalid particles*/, true /*do halo position change*/,
						true /*removeFromContainer*/);
	}
	// TODO: copying own molecules out and reinserting them can be done within autopas more efficiently
	std::vector<Molecule> ownMolecules{};
	ownMolecules.reserve(particleContainer->getNumberOfParticles());
	for (auto iter = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); iter.isValid(); ++iter) {
		ownMolecules.push_back(*iter);
		// TODO: This check should be in debug mode only
		if (not iter->inBox(newMin.data(), newMax.data())) {
			std::ostringstream error_message;
			error_message << "Particle still in domain that should have been migrated."
						  << "BoxMin: " << particleContainer->getBoundingBoxMin(0) << ", "
						  << particleContainer->getBoundingBoxMin(1) << ", " << particleContainer->getBoundingBoxMin(2)
						  << "\n"
						  << "BoxMax: " << particleContainer->getBoundingBoxMax(0) << ", "
						  << particleContainer->getBoundingBoxMax(1) << ", " << particleContainer->getBoundingBoxMax(2)
						  << "\n"
						  << "Particle: \n"
						  << *iter << std::endl;
			MARDYN_EXIT(error_message.str());
		}
	}
	particleContainer->clear();
	particleContainer->rebuild(newMin.data(), newMax.data());
	particleContainer->addParticles(ownMolecules);
	bool allDone = false;
	double waitCounter = 30.0;
	double deadlockTimeOut = 360.0;
	double startTime = MPI_Wtime();
	while (not allDone) {
		allDone = true;

		// "kickstart" processing of all Isend requests
		for (auto& sender : sendNeighbors) {
			allDone &= sender.testSend();
		}

		// unpack molecules
		for (auto& recv : recvNeighbors) {
			allDone &= recv.iprobeCount(this->getCommunicator(), this->getMPIParticleType());
			allDone &= recv.testRecv(particleContainer, false);
		}

		// catch deadlocks
		double waitingTime = MPI_Wtime() - startTime;
		if (waitingTime > waitCounter) {
			Log::global_log->warning() << "KDDecomposition::migrateParticles: Deadlock warning: Rank " << _rank
									   << " is waiting for more than " << waitCounter << " seconds" << std::endl;
			waitCounter += 1.0;
			for (auto& sender : sendNeighbors) {
				sender.deadlockDiagnosticSend();
			}
			for (auto& recv : recvNeighbors) {
				recv.deadlockDiagnosticRecv();
			}
		}

		if (waitingTime > deadlockTimeOut) {
			Log::global_log->error() << "KDDecomposition::migrateParticles: Deadlock error: Rank " << _rank
									 << " is waiting for more than " << deadlockTimeOut << " seconds" << std::endl;
			for (auto& sender : sendNeighbors) {
				sender.deadlockDiagnosticSend();
			}
			for (auto& recv : recvNeighbors) {
				recv.deadlockDiagnosticRecv();
			}
			break;
		}
	}
}