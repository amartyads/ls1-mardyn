/*
 * DomainDecompMutable.cpp
 *
 *  Created on: 10 August 2026
 *      Author: amartyads, seckler
 */

#include "DomainDecompMutable.h"

DomainDecompMutable::DomainDecompMutable() :  DomainDecompMutable(MPI_COMM_WORLD, {0.}, {0.}, {0.}){}

DomainDecompMutable::DomainDecompMutable(MPI_Comm comm) : DomainDecompMutable(comm, {0.}, {0.}, {0.}) {}

DomainDecompMutable::DomainDecompMutable(std::array<double, 3> boxMin, std::array<double, 3> boxMax,
										 std::array<double, 3> domainLength) : DomainDecompMutable(MPI_COMM_WORLD, boxMin, boxMax, domainLength){}

DomainDecompMutable::DomainDecompMutable(MPI_Comm comm, std::array<double, 3> boxMin, std::array<double, 3> boxMax,
										 std::array<double, 3> domainLength) : DomainDecompMPIBase(comm), _boxMin(boxMin), _boxMax(boxMax), _domainLength(domainLength) {}

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

std::array<size_t, 3> DomainDecompMutable::getOptimalGrid(const std::array<double, 3>& domainLength,
																 int numProcs) {
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


