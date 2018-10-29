//
// Created by Kruegener on 8/29/2018.
//

#include "VelocityAbsProfile.h"
#include "plugins/profiles/DensityProfile.h"

void VelocityAbsProfile::output(string prefix) {
    global_log->info() << "[VelocityAbsProfile] output" << std::endl;

    // Setup outfile
    _profilePrefix = prefix + "_kartesian.VAbspr";
    ofstream outfile(_profilePrefix.c_str());
    outfile.precision(6);

    // Write header
    outfile << "//Segment volume: " << _samplInfo.segmentVolume << "\n//Accumulated data sets: " << _samplInfo.accumulatedDatasets << "\n//Local profile of magnitude of velocity. Output file generated by the \"VelocityAbsProfile\" method, plugins/profiles. \n";
    outfile << "//local velocity magnitude profile: Y - Z || X-projection\n";
    // TODO: more info
    outfile << "// \t dX \t dY \t dZ \n";
    outfile << "\t" << 1/_samplInfo.universalInvProfileUnit[0] << "\t" << 1/_samplInfo.universalInvProfileUnit[1] << "\t" << 1/_samplInfo.universalInvProfileUnit[2]<< "\n";
    outfile << "0 \t";

    // Invoke matrix write routine
    writeMatrix(outfile);

    outfile.close();
}

void VelocityAbsProfile::writeDataEntry(unsigned long uID, ofstream &outfile) const {
    // Check for division by 0
    long double vd;
    int numberDensity = _densityProfile->getGlobalNumberDensity(uID);
    if(numberDensity != 0){
        vd = _globalProfile.at(uID) / numberDensity;
    }
    else{
        vd = 0;
    }
    outfile << vd << "\t";
}
