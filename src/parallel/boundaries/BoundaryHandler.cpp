/*
 * BoundaryHandler.cpp
 *
 *  Created on: 24 March 2023
 *      Author: amartyads
 */

#include <algorithm>

#include "BoundaryHandler.h"
#include "BoundaryUtils.h"

#include "Simulation.h"
#include "integrators/Integrator.h"

#include "utils/Logger.h" 

BoundaryHandler::BoundaryHandler() : boundaries{
	{DimensionType::POSX, BoundaryType::PERIODIC}, {DimensionType::POSY, BoundaryType::PERIODIC}, 
	{DimensionType::POSZ, BoundaryType::PERIODIC}, {DimensionType::NEGX, BoundaryType::PERIODIC}, 
	{DimensionType::NEGY, BoundaryType::PERIODIC}, {DimensionType::NEGZ, BoundaryType::PERIODIC},
	{DimensionType::ERROR, BoundaryType::ERROR}
	} {}


BoundaryType BoundaryHandler::getBoundary(DimensionType dimension) const
{
	return boundaries.at(dimension);
}

BoundaryType BoundaryHandler::getBoundary(std::string dimension) const
{
	DimensionType convertedDimension = BoundaryUtils::convertStringToDimension(dimension);
	return getBoundary(convertedDimension);
}

BoundaryType BoundaryHandler::getBoundary(int dimension) const 
{
	return getBoundary(BoundaryUtils::convertLS1DimsToDimensionPos(dimension));
}

void BoundaryHandler::setBoundary(DimensionType dimension, BoundaryType value)
{
	if(dimension != DimensionType::ERROR)
		boundaries[dimension] = value;
}

void BoundaryHandler::setBoundary(std::string dimension, BoundaryType value) 
{
	DimensionType convertedDimension = BoundaryUtils::convertStringToDimension(dimension);
	setBoundary(convertedDimension, value);
}

void BoundaryHandler::setGlobalRegion(double* start, double* end) 
{
	for(short int i = 0; i < 3; i++) {
		_globalRegionStart[i] = start[i];
		_globalRegionEnd[i] = end[i];
	}
}

void BoundaryHandler::setLocalRegion(double* start, double* end) 
{
	for(short int i = 0; i < 3; i++) {
		_localRegionStart[i] = start[i];
		_localRegionEnd[i] = end[i];
	}
}

void BoundaryHandler::setGlobalRegion(std::array<double, 3> start, std::array<double, 3> end) {
	_globalRegionStart = start;
	_globalRegionEnd = end;
}

void BoundaryHandler::setLocalRegion(std::array<double, 3> start, std::array<double, 3> end) {
	_localRegionStart = start;
	_localRegionEnd = end;
}

void BoundaryHandler::findOuterWallsInLocalRegion() 
{
	_isOuterWall[DimensionType::POSX] = (_localRegionEnd[0] == _globalRegionEnd[0]);
	_isOuterWall[DimensionType::NEGX] = (_localRegionStart[0] == _globalRegionStart[0]);
	_isOuterWall[DimensionType::POSY] = (_localRegionEnd[1] == _globalRegionEnd[1]);
	_isOuterWall[DimensionType::NEGY] = (_localRegionStart[1] == _globalRegionStart[1]);
	_isOuterWall[DimensionType::POSZ] = (_localRegionEnd[2] == _globalRegionEnd[2]);
	_isOuterWall[DimensionType::NEGZ] = (_localRegionStart[2] == _globalRegionStart[2]);
}

bool BoundaryHandler::hasInvalidBoundary() const
{
	return boundaries.at(DimensionType::POSX) == BoundaryType::ERROR ||
		boundaries.at(DimensionType::POSY) == BoundaryType::ERROR ||
		boundaries.at(DimensionType::POSZ) == BoundaryType::ERROR ||
		boundaries.at(DimensionType::NEGX) == BoundaryType::ERROR ||
		boundaries.at(DimensionType::NEGY) == BoundaryType::ERROR ||
		boundaries.at(DimensionType::NEGZ) == BoundaryType::ERROR;
}

bool BoundaryHandler::isOuterWall(DimensionType dimension) const
{
	return _isOuterWall.at(dimension);
}

bool BoundaryHandler::isOuterWall(int dimension) const
{
	return isOuterWall(BoundaryUtils::convertLS1DimsToDimensionPos(dimension));
}

bool BoundaryHandler::processOuterWallLeavingParticles()
{
	auto moleculeContainer = global_simulation->getMoleculeContainer(); // :-(
	double timestepLength = (global_simulation->getIntegrator())->getTimestepLength(); 
	double cutoff = moleculeContainer->getCutoff();
	for (auto const& currentWall : _isOuterWall)
	{
		//global_log->info() << "wall number " << BoundaryUtils::convertDimensionToString(currentWall.first) << " : " << currentWall.second << std::endl;
		if(!currentWall.second)
			continue;

		switch(getBoundary(currentWall.first))
		{
			case BoundaryType::PERIODIC:
				//default behaviour
				break;
			
			case BoundaryType::OUTFLOW:
			case BoundaryType::REFLECTING:
			{
				//create region
				std::array<double,3> curWallRegionBegin, curWallRegionEnd;
				std::tie(curWallRegionBegin, curWallRegionEnd) = BoundaryUtils::getInnerBuffer(_localRegionStart, _localRegionEnd, currentWall.first, cutoff);
				//conversion
				const double cstylerbegin[] = {curWallRegionBegin[0], curWallRegionBegin[1], curWallRegionBegin[2]}; 
				const double cstylerend[] = {curWallRegionEnd[0], curWallRegionEnd[1], curWallRegionEnd[2]};
				
				auto particlesInRegion = moleculeContainer->regionIterator(cstylerbegin, cstylerend, ParticleIterator::ONLY_INNER_AND_BOUNDARY);
				for (auto it = particlesInRegion; it.isValid(); ++it)
				{
					Molecule curMolecule = *it;
					global_log->info() << "Boundary particle found " << std::endl;
					if (BoundaryUtils::isMoleculeLeaving(curMolecule, curWallRegionBegin, curWallRegionEnd, currentWall.first, timestepLength))
					{
						global_log->info() << "Boundary particle found leaving" << std::endl;
						int currentDim = BoundaryUtils::convertDimensionToLS1Dims(currentWall.first);
						if(getBoundary(currentWall.first) == BoundaryType::REFLECTING)
						{
							global_log->info() << "Reflection particle found " << std::endl;
							double vel = it->v(currentDim);
							it->setv(currentDim, -vel);
						}
						else
						{
							global_log->info() << "Outflow particle found " << std::endl;
							moleculeContainer->deleteMolecule(it, false);
						}
					}
				}
				break;
			}
			default:
				global_log->error() << "Boundary type error! Received type not allowed!" << std::endl;
				Simulation::exit(1);
		}
	}
	return true;
}

void BoundaryHandler::removeNonPeriodicHalos()
{
	auto moleculeContainer = global_simulation->getMoleculeContainer();
	double cutoff = moleculeContainer->getCutoff();
	for (auto const& currentWall : _isOuterWall)
	{
		//global_log->info() << "wall number " << BoundaryUtils::convertDimensionToString(currentWall.first) << " : " << currentWall.second << std::endl;
		if(!currentWall.second) //not an outer wall
			continue;

		switch(getBoundary(currentWall.first))
		{
			case BoundaryType::PERIODIC:
				//default behaviour
				break;
			
			case BoundaryType::OUTFLOW:
			case BoundaryType::REFLECTING:
			{
				//create region
				std::array<double,3> curWallRegionBegin, curWallRegionEnd;
				std::tie(curWallRegionBegin, curWallRegionEnd) = BoundaryUtils::getOuterBuffer(_localRegionStart, _localRegionEnd, currentWall.first, cutoff);
				//conversion
				const double cstylerbegin[] = {curWallRegionBegin[0], curWallRegionBegin[1], curWallRegionBegin[2]}; 
				const double cstylerend[] = {curWallRegionEnd[0], curWallRegionEnd[1], curWallRegionEnd[2]};

				auto particlesInRegion = moleculeContainer->regionIterator(cstylerbegin, cstylerend, ParticleIterator::ALL_CELLS);
				for (auto it = particlesInRegion; it.isValid(); ++it)
				{
					global_log->info() << "Halo particle found " << std::endl;
					moleculeContainer->deleteMolecule(it, false);
				}
				break;
			}
			default:
				global_log->error() << "Boundary type error! Received type not allowed!" << std::endl;
				Simulation::exit(1);
		}
	}
	#ifndef MARDYN_AUTOPAS
	moleculeContainer->updateBoundaryAndHaloMoleculeCaches();
	#endif
}