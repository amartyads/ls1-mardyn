/*
 * DomainDecompBaseTest.h
 *
 * @Date: 09.05.2012
 * @Author: eckhardw
 */

#ifndef DOMAINDECOMPBASETEST_H_
#define DOMAINDECOMPBASETEST_H_

#include "utils/TestWithSimulationSetup.h"

class DomainDecompBaseTest: public utils::TestWithSimulationSetup {

	TEST_SUITE(DomainDecompBaseTest);
	TEST_METHOD(testNoDuplicatedParticles);
	TEST_METHOD(testNoLostParticles);
	TEST_METHOD(testExchangeMoleculesSimple);
	TEST_METHOD(testExchangeMolecules);
	TEST_SUITE_END();

public:
	DomainDecompBaseTest();

	virtual ~DomainDecompBaseTest();
	void testNoDuplicatedParticlesFilename(const char * filename, double cutoff);
	void testNoDuplicatedParticles();
	void testNoLostParticlesFilename(const char * filename, double cutoff);
	void testNoLostParticles();
	void testExchangeMoleculesSimple();
	void testExchangeMolecules();
};

#endif /* DOMAINDECOMPBASETEST_H_ */
