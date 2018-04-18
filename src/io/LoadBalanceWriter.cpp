#include "io/LoadBalanceWriter.h"

#include "Simulation.h"
#include "TimerProfiler.h"
#include "parallel/DomainDecompBase.h"

LoadbalanceWriter::LoadbalanceWriter::LoadbalanceWriter() :
	_writeFrequency(1),
	_outputFilename("mardyn-LB.dat"),
	_defaultTimer(nullptr),
	_timerNames(0),
	_times(0),
	_global_times(0),
	_simsteps(0),
	_warninglevels()
{}

void LoadbalanceWriter::readXML(XMLfileUnits& xmlconfig) {
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << endl;
	xmlconfig.getNodeValue("outputfilename", _outputFilename);
	global_log->info() << "Output filename: " << _outputFilename << endl;

	XMLfile::Query query = xmlconfig.query("timers/timer");
	std::string oldpath = xmlconfig.getcurrentnodepath();
	for(auto timerIter = query.begin(); timerIter; ++timerIter ) {
		xmlconfig.changecurrentnode(timerIter);
		std::string timername;
		xmlconfig.getNodeValue("name", timername);
		if(timername != LB_WRITER_DEFAULT_TIMER_NAME) {
			_timerNames.push_back(timername);
		}
		double warninglevel = 0;
		xmlconfig.getNodeValue("warninglevel", warninglevel);
		_warninglevels[timername] = warninglevel;
		global_log->info() << "Added timer for LB monitoring: " <<  timername << ", warninglevel: " << warninglevel << std::endl;
	}
	xmlconfig.changecurrentnode(oldpath);
}

void LoadbalanceWriter::init(ParticleContainer *particleContainer,
                             DomainDecompBase *domainDecomp, Domain *domain) {
	std::string default_timer_name(LB_WRITER_DEFAULT_TIMER_NAME);
	_defaultTimer = new Timer();
	global_simulation->timers()->registerTimer(default_timer_name, vector<string>{"SIMULATION"}, _defaultTimer);
	_timerNames.push_back(default_timer_name);

	size_t timestep_entry_offset = 2* _timerNames.size();
	_simsteps.reserve(_writeFrequency);
	_times.reserve(timestep_entry_offset * _writeFrequency);
	_global_times.reserve(timestep_entry_offset * _writeFrequency);

	_defaultTimer->reset();
	_defaultTimer->start();

	if(0 == domainDecomp->getRank()) {
		writeOutputFileHeader();
	}
}

void LoadbalanceWriter::endStep(
        ParticleContainer *particleContainer,
        DomainDecompBase *domainDecomp, Domain *domain,
        unsigned long simstep, std::list<ChemicalPotential> *lmu,
        std::map<unsigned, CavityEnsemble> *mcav
)  {
	_defaultTimer->stop();
	recordTimes(simstep);
	_defaultTimer->reset();

	if((simstep % _writeFrequency) == 0) {
		  flush(domainDecomp);
	}
	_defaultTimer->start();
}

void LoadbalanceWriter::writeOutputFileHeader() {
	std::ofstream outputfile(_outputFilename, std::ofstream::out);
	outputfile << "#simstep";
	for(auto timername : _timerNames) {
		outputfile << "\t#" << timername <<"#\t\t";
	}
	outputfile << "\n\t";
	for(auto timername : _timerNames) {
		outputfile << "\tmin\tmax\tf_LB";
	}
	outputfile << std::endl;
	outputfile.close();
}

void LoadbalanceWriter::recordTimes(unsigned long simstep) {
	_simsteps.push_back(simstep);
	for(auto timername : _timerNames) {
		Timer *timer = global_simulation->timers()->getTimer(timername);
		double time = timer->get_etime();
		_times.push_back(time);  // needed for maximum
		_times.push_back(-time); // needed for minimum
	}
}

void LoadbalanceWriter::flush(DomainDecompBase* domainDecomp) {
#ifdef ENABLE_MPI
	//! @todo If this shall become a general LB monitor/manager a MPI_Allreduce will be needed here
	MPI_CHECK(MPI_Reduce( _times.data(), _global_times.data(), _times.size(), MPI_DOUBLE, MPI_MAX, 0, domainDecomp->getCommunicator()));
#else
	//! @todo in case we do not reuse _times later on, we may use here  _global_times = std::move(_times);
	_global_times = _times;
#endif
	if(0 == domainDecomp->getRank()) {
		std::ofstream outputfile(_outputFilename,  std::ofstream::app);
		for(size_t i = 0; i < _simsteps.size(); ++i) {
			writeLBEntry(i, outputfile);
		}
		outputfile.close();
	}
	resetTimes();
}

void LoadbalanceWriter::writeLBEntry(size_t id, std::ofstream &outputfile) {
	outputfile << _simsteps[id];
	size_t timestep_entry_offset = 2* _timerNames.size();
	size_t base_offset = timestep_entry_offset * id;
	for(size_t timer_id = 0; timer_id < _timerNames.size(); ++timer_id) {
		double min, max, f_LB;
		size_t timer_offset = 2 * timer_id;
		size_t offset = base_offset + timer_offset;
		max = _global_times[offset];
		min = -_global_times[offset + 1];
		f_LB = max / min;
		std::string timername = _timerNames[timer_id];
		if(_warninglevels.count(timername) > 0 && _warninglevels[timername] < f_LB) {
			displayWarning(_simsteps[id], timername, f_LB);
		}
		outputfile << "\t" << min << "\t" << max << "\t" << f_LB;
	}
	outputfile << std::endl;
}

void LoadbalanceWriter::displayWarning(unsigned long simstep, std::string timername, double f_LB) {
	global_log->warning() << "Load balance limit exceeded in simstep " << simstep
		<< " for timer " << timername
		<< ", limit: " << _warninglevels[timername]
		<< "  value: " << f_LB << std::endl;
}

void LoadbalanceWriter::resetTimes() {
	_times.clear();
	_global_times.clear();
	_simsteps.clear();
}
