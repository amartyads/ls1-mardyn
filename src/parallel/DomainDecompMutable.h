/*
 * DomainDecompMutable.h
 *
 *  Created on: 10 August 2026
 *      Author: amartyads, seckler
 */

#include "DomainDecompMPIBase.h"

#pragma once

/**
 * @brief Class containing methods that are shared between all domain decompositions that change during runtime.
 * 
 */
class DomainDecompMutable : public DomainDecompMPIBase {
public:
	DomainDecompMutable();
	DomainDecompMutable(MPI_Comm comm);
	DomainDecompMutable(std::array<double, 3> boxMin, std::array<double, 3> boxMax, std::array<double, 3> domainLength);
	DomainDecompMutable(MPI_Comm comm, std::array<double, 3> boxMin, std::array<double, 3> boxMax, std::array<double, 3> domainLength);
	virtual ~DomainDecompMutable();

	// functions shifted here from GeneralDomainDecomposition.h
	// documentation see father class (DomainDecompBase.h)
	double getBoundingBoxMin(int dimension, Domain* domain) override;

	// documentation see father class (DomainDecompBase.h)
	double getBoundingBoxMax(int dimension, Domain* domain) override;
	
	/**
	 * Get the optimal grid for the given dimensions of the box and the number of processes.
	 * The grid is produced, s.t., the number of grid[0] * grid[1] * grid[2] == numProcs
	 * The edge lengths of the grid will resemble the lengths of the domain, i.e., the longest edge of the domain will
	 * also have the largest amount of grid points.
	 * @param domainLength
	 * @param numProcs
	 * @return
	 */
	static std::array<size_t, 3> getOptimalGrid(const std::array<double, 3>& domainLength, int numProcs);

	/**
	 * Get the coordinates from the rank.
	 * S.t. rank = x * gridSize[1]*gridSize[2] + y * gridSize[2] + z
	 * @param gridSize
	 * @param rank
	 * @return
	 */
	static std::array<size_t, 3> getCoordsFromRank(const std::array<size_t, 3>& gridSize, int rank);

	/**
	 * Returns boxMin and boxMax according to regular grid.
	 * @param domainLength
	 * @param gridSize
	 * @param gridCoords
	 * @return boxMin and boxMax
	 */
	static std::tuple<std::array<double, 3>, std::array<double, 3>> initializeRegularGrid(
		const std::array<double, 3>& domainLength, const std::array<size_t, 3>& gridSize,
		const std::array<size_t, 3>& gridCoords);


protected:
	std::array<double, 3> _boxMin;
	std::array<double, 3> _boxMax;
	std::array<double, 3> _domainLength;
};
