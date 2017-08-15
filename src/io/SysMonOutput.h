#ifndef SYSMONOUTPUT_H_
#define SYSMONOUTPUT_H_

#include "io/OutputBase.h"


class SysMonOutput : public OutputBase {
public:
	SysMonOutput();
	~SysMonOutput(){}

	void readXML(XMLfileUnits& xmlconfig);

	//! @todo comment
	void initOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);
	//! @todo comment

	void doOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep, std::list<ChemicalPotential>* lmu, std::map<unsigned, CavityEnsemble>* mcav);
	
	//! @todo comment
	void finishOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);

	static std::string getPluginName() {
		return std::string("SysMonOutput");
	}
	static OutputBase* createInstance() { return new SysMonOutput(); }

private:
	//! prefix for the names of all output files
	std::ofstream _resultStream;
	unsigned long _writeFrequency;
};

#endif /* SYSMONOUTPUT_H_ */
