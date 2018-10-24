//
// Created by Kruegener on 10/22/2018.
//

#include "TemperatureProfile.h"
#include "DOFProfile.h"
#include "KineticProfile.h"


void TemperatureProfile::output(string prefix) {
    global_log->info() << "[TemperatureProfile] output" << std::endl;
    _profilePrefix = prefix;
    _profilePrefix += "_kartesian.Temppr";
    ofstream outfile(_profilePrefix.c_str());
    outfile.precision(6);

    // TODO: FIX ACCESS TO ELEMENTS / _yOff ???
    outfile << "//Segment volume: " << _kartProf->segmentVolume << "\n//Accumulated data sets: " << _kartProf->accumulatedDatasets << "\n//Local profile of the temperature. Output file generated by the \"TemperatureProfile\" method, plugins/profiles. \n";
    outfile << "//Temperature expressed by 2Ekin/#DOF\n";
    outfile << "//local temperature profile: Y - Z || X-projection\n";
    //outfile << "//one single matrix of the local number density rho'(phi_i;r_i',h_i') \n//       | r_i'\n// ---------------------\n//   h_i'| rho'(r_i',h_i')\n//       | \n";
    //outfile << "//  T' \t sigma_ii' \t eps_ii' \t yOffset \t DELTA_phi \t DELTA_r2' \t DELTA_h' \t quantities in atomic units are denoted by an apostrophe '\n";
    //outfile << _kartProf->dom->getTargetTemperature(0)<<"\t"<<_kartProf->dom->getSigma(0, 0)<<"\t"<<_kartProf->dom->getepsilonRF()<<"\t"; //changed by Michaela Heier
    outfile << "// \t dX \t dY \t dZ \n";
    outfile << "\t" << 1/_kartProf->universalInvProfileUnit[0] << "\t" << 1/_kartProf->universalInvProfileUnit[1] << "\t" << 1/_kartProf->universalInvProfileUnit[2]<< "\n";
    outfile << "0 \t";

    // Assuming x = 1 -> projection along x axis

    // Z - axis label
    //global_log->info() << "[DensityProfile] output -> units: " << _kartProf->universalProfileUnit[0] << " " << _kartProf->universalProfileUnit[1] << " "<< _kartProf->universalProfileUnit[2] << "\n";
    for(unsigned z = 0; z < _kartProf->universalProfileUnit[2]; z++){
        outfile << (z+0.5) / _kartProf->universalInvProfileUnit[2] <<"  \t"; // Eintragen der z Koordinaten in Header
    }
    outfile << "\n";

    // Y - axis label
    for(unsigned y = 0; y < _kartProf->universalProfileUnit[1]; y++){
        double hval = (y + 0.5) / _kartProf->universalInvProfileUnit[1];
        outfile << hval << "  \t";

        // temperature values
        for(unsigned z = 0; z < _kartProf->universalProfileUnit[2]; z++){
            for(unsigned x = 0; x < _kartProf->universalProfileUnit[0]; x++){
                long unID = (long) (x * _kartProf->universalProfileUnit[0] * _kartProf->universalProfileUnit[2] + y * _kartProf->universalProfileUnit[1] + z);
                int dofs = _dofProfile->getGlobalDOF(unID);
                if(dofs == 0){
                    outfile << 0.0 << "\t";
                }
                else{
                    outfile << (_kineticProfile->getGlobalKineticEnergy(unID) / dofs) << "\t";
                }
            }
        }
        outfile << "\n";
    }
    outfile.close();
}

// TODO: log Ekin / DOF
//_universalDOFProfile
//_universalKineticProfile
