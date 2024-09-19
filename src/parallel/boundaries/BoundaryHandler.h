/*
 * BoundaryHandler.h
 *
 *  Created on: 24 March 2023
 *      Author: amartyads
 */

#pragma once

#include "BoundaryUtils.h"
#include "DimensionUtils.h"
#include "RegionUtils.h"
#include "particleContainer/ParticleContainer.h"

#include <array>
#include <map>
#include <string>
#include <vector>

/**
 * Class to handle boundary conditions, namely leaving and halo particles.
 *
 * The objects of this class store the local and global bounds of the subdomain
 * in every process, and provide functions to deal with leaving particles, and
 * delete halo particles.
 *
 * The internal walls of the subdomain, touching other subdomains are 'local'
 * walls while the walls that are also the limits of the global domain are
 * 'global' walls.
 *
 * Since the behaviour of 'local' walls are unchanged (and is identical to the
 * behaviour of global periodic walls), they are assigned to
 * 'PERIODIC_OR_LOCAL'.
 */

class BoundaryHandler {
  public:
	BoundaryHandler() = default;

	/* Find the boundary type of a global wall for a particular dimension. */
	BoundaryUtils::BoundaryType getGlobalWallType(DimensionUtils::DimensionType dimension) const;

	/* Set the boundary type of a global wall for a particular dimension. */
	void setGlobalWallType(DimensionUtils::DimensionType dimension, BoundaryUtils::BoundaryType value);

	BoundaryUtils::BoundaryType getGlobalWallType(int dimension) const;

	/* Check if any of the global boundaries have invalid types. */
	bool hasInvalidBoundary() const;

	/**
	 *  Check if any of the global boundaries are non-periodic.
	 *
	 *  This check helps bypass all boundary-related code if default behaviour
	 * (all periodic boundaries) is expected.
	 */
	bool hasNonPeriodicBoundary() const;

	/* Set bounds for global subdomain. */
	void setGlobalRegion(const double *start, const double *end);

	/* Set bounds for local subdomain. */
	void setLocalRegion(const double *start, const double *end);
	/**
	 * Determine which walls in the local region are actually global walls.
	 *
	 * Should be called after changing global and local regions (typically after a
	 * rebalance).
	 */
	void updateGlobalWallLookupTable();

	/* Check if the local wall in a particular dimension is actually a global
	 * wall. */
	bool isGlobalWall(DimensionUtils::DimensionType dimension) const;

	/* Check if the local wall in a particular dimension is actually a global
	 * wall. */
	bool isGlobalWall(int dimension) const;

	/**
	 * Processes all particles that would leave the global domain.
	 *
	 * If a subdomain has no global walls, this function does nothing.
	 * For every global wall, the function iterates through all particles that are
	 * within one cutoff distance away from the wall. If these particles would
	 * leave the global box in the next simulation, the following is done:
	 *
	 * PERIODIC_OR_LOCAL - No actions taken (default behaviour).
	 * REFLECTING - The particle's velocity is reversed normal to the wall it's
	 * leaving.
	 * OUTFLOW - The particle is deleted.
	 */
	void processGlobalWallLeavingParticles(ParticleContainer *moleculeContainer, double timestepLength) const;

	/**
	 * Processes all halo particles outside the global domain.
	 *
	 * If a subdomain has no global walls, this function does nothing.
	 * For every global wall, the function iterates through all halo particles
	 * that are within one cutoff distance away from the wall. The following is
	 * done for each particle:
	 *
	 * PERIODIC_OR_LOCAL - No actions taken (default behaviour).
	 * REFLECTING / OUTFLOW - The halo particle is deleted, so that particles
	 * approaching the boundary do not decelerate due to influence from the halo
	 * particles, and preserve their velocities before being bounced/deleted
	 */
	void removeNonPeriodicHalos(ParticleContainer *moleculeContainer) const;

  private:
	/* List of global boundary type per dimension. */
	std::map<DimensionUtils::DimensionType, BoundaryUtils::BoundaryType> _boundaries{
		{DimensionUtils::DimensionType::POSX, BoundaryUtils::BoundaryType::PERIODIC_OR_LOCAL},
		{DimensionUtils::DimensionType::POSY, BoundaryUtils::BoundaryType::PERIODIC_OR_LOCAL},
		{DimensionUtils::DimensionType::POSZ, BoundaryUtils::BoundaryType::PERIODIC_OR_LOCAL},
		{DimensionUtils::DimensionType::NEGX, BoundaryUtils::BoundaryType::PERIODIC_OR_LOCAL},
		{DimensionUtils::DimensionType::NEGY, BoundaryUtils::BoundaryType::PERIODIC_OR_LOCAL},
		{DimensionUtils::DimensionType::NEGZ, BoundaryUtils::BoundaryType::PERIODIC_OR_LOCAL}};

	/* Lookup table to check if a local wall is also global. */
	std::map<DimensionUtils::DimensionType, bool> _isGlobalWall;

	/* Global region start/end. */
	std::array<double, 3> _globalRegionStart, _globalRegionEnd;

	/* Local region start/end. */
	std::array<double, 3> _localRegionStart, _localRegionEnd;
};
