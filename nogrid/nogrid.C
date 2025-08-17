#include <iostream>
#include <fstream>
#include <numeric>

#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Random.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ComponentChargedRing.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/FundamentalConstants.hh"



#define LOG(x) std::cout<<x<<std::endl

using namespace Garfield;

AvalancheMC drift;
AvalancheMicroscopic aval;
ComponentChargedRing rings;

// set time window parameters
double tmin = 0.;
double timestep = 0.0025;

std::vector<double> mean_pos = {0.,0.,0.};
double meanZ,meanX;
int particle_count = 0;

double sft = 0.00001;
double st = 0.00005;

//std::ofstream ion_times;


void userHandleIonisation(double x, double y, double z, double t, int type, int level, Medium * m){
    drift.AddIon(x,y,z,tmin + timestep); // ion added after step.
    //ion_times << (t-tmin)/timestep << "\n";
    //std::cout << "ionisation at t = " << t << " x = " << x << " y = " << y << " z = " << z << "\n";
}

void userHandleAttachment(double x, double y, double z, double t, int type, int level, Medium * m){
    drift.AddNegativeIon(x,y,z,tmin + timestep);
}


std::array<int,3> get_mean(bool & electrons_remaining, bool & ions_remaining, const double & zmax){
    int ne = 0;
    int ni = 0;
    int nNi = 0;
    for (const auto& electron : aval.GetElectrons()) {
        // get electron positions
        double xf = electron.path.back().x;
        double yf = electron.path.back().y;
        double zf = electron.path.back().z;
        int status = electron.status;
        // if still in a drift medium or at the end of timestep:
        if (status==0 || status == -17){
	      ne++;
          electrons_remaining = true;
          if (yf < zmax){
              particle_count++;
              mean_pos[0] += xf;
              mean_pos[1] += yf;
              mean_pos[2] += zf;
            }
        }
    }
    for (const auto& ion : drift.GetIons()) {
        // get electron positions
        double xf = ion.path.back().x;
        double yf = ion.path.back().y;
        double zf = ion.path.back().z;
        int status = ion.status;
        // if still in a drift medium or at the end of timestep:
        if (status==0 || status == -17){
	      ni++;
          ions_remaining = true;
          if (yf < zmax){
              particle_count++;
              mean_pos[0] += xf;
              mean_pos[1] += yf;
              mean_pos[2] += zf;
            }
        }
    }
    for (const auto& ion : drift.GetNegativeIons()) {
        // get electron positions
        double xf = ion.path.back().x;
        double yf = ion.path.back().y;
        double zf = ion.path.back().z;
        int status = ion.status;
        // if still in a drift medium or at the end of timestep:
        if (status==0 || status == -17){
	      nNi++;
          ions_remaining = true;
          if (yf < zmax){
              particle_count++;
              mean_pos[0] += xf;
              mean_pos[1] += yf;
              mean_pos[2] += zf;
            }
        }
    }
    return {ne,ni,nNi};
}

int main(int argc, char * argv[]) {
  
  // Check if the right number of arguments are provided
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <process_index> <voltage>" << std::endl;
    return 1;
  }

  int run_number = std::atoi(argv[1]);
  std::string s_run_number = argv[1];
  int i_dv = std::atoi(argv[2]);
  double dv = std::atof(argv[2]);
  std::string s_dv = std::to_string(i_dv);
  dv = -1*dv;
  s_dv = "_"+s_dv;
  std::ofstream size;
  
  // Convert the process index argument to an integer
  int processIndex = std::atoi(argv[1]);
  
  // Example logic: Print a message including the process index
  LOG("Script: Running job with process index: " << processIndex);
  
  // Setup the gas
  MediumMagboltz gas;
  gas.SetComposition("ar", 93, "co2", 7); // [%]
  gas.SetTemperature(293.15); // [K]
  gas.SetPressure(760.); // [Torr]
  gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");
  gas.Initialise(true);
  gas.EnablePenningTransfer(); // Penning effect. This I forgot to explain... [1]
  
  // Add parallel plate geometry for space-charge
  const double posBottomPlane = 0.; // cm
  const double posTopPlane = 128.e-4; // cm
  
  ComponentAnalyticField mm; // Micromegas field
  mm.SetMedium(&gas);
  mm.AddPlaneY(posBottomPlane, 0.);
  mm.AddPlaneY(posTopPlane, dv);
  
  rings.SetArea(-0.02,posBottomPlane,-0.02,0.02,posTopPlane,0.02);
  rings.SetSpacingTolerance(st);
  rings.SetSelfFieldTolerance(sft);
  rings.SetMedium(&gas);

  
  // Create the sensor
  Sensor sensor;
  sensor.AddComponent(&mm);
  sensor.AddComponent(&rings);
  
  // Create electrons
  aval.SetSensor(&sensor);
  aval.EnableRKNSteps(false);
  
  drift.SetSensor(&sensor);

  const bool enableDebug = false;
  const bool enableSpaceCharge = true;
  const bool plotField = false;
  const bool plotDrift = false;

  //aval.EnableAvalancheSizeLimit(500);
  aval.SetUserHandleIonisation(userHandleIonisation);
  aval.SetUserHandleAttachment(userHandleAttachment);

  std::string path_end = ".txt";
  if (enableSpaceCharge){
      std::string path_begin = "/afs/cern.ch/work/t/tszwarce/nogrid-static/gain_output/on/size_sc_";
    
      size.open(path_begin+s_run_number+s_dv+path_end);
      //ion_times.open("/afs/cern.ch/work/t/tszwarce/nogrid-static/ion_times.txt");
  }
  else{
      std::string path_begin = "/afs/cern.ch/work/t/tszwarce/nogrid-static/gain_output/off/size_no_sc_";
      size.open(path_begin+s_run_number+s_dv+path_end);
  }   

  if (enableDebug) { 
    rings.EnableDebugging();
    aval.EnableDebugging();
    drift.EnableDebugging();
  }
  
    
  const double elCharge = 1.60217663e-19; // C
  
  double midpt = (posTopPlane - posBottomPlane)/2;
  const double x0 = 0.; // cm
  double y0 = 127.9e-4;
  const double z0 = 0.; // cm
  const double t0 = 0.; // ns
  const double e0 = 0.1; // eV

  const unsigned int nFrames = 80; 
  int framecount = 0;
  int num_new_ions;

  aval.AddElectron(x0,y0,z0,0.,0.1); 


  std::array<int,3> particles;
  int iF = -1;

  while (1) {
    ++iF;

    std::cout << "Frame " << iF << "\nTimestep: " << timestep << " ns\n";
    

    meanZ = 0;
    meanX = 0;

    particle_count = 0;
    mean_pos = {0.,0.,0.};

    bool electrons_remaining = false;
    bool ions_remaining = false;

    particles = get_mean(electrons_remaining,ions_remaining,posTopPlane);
    std::cout << "e: " << particles[0] << ", +i: " << particles[1] << ", -i: " << particles[2] << "\n";
    num_new_ions = particles[1] - framecount;
    std::cout << num_new_ions << " new ions\n";
    if (num_new_ions < 10 && particles[0] < 100) {
        timestep = 0.05;
        std::cout << "Too few ionisations: timestep updated to 0.05\n";
    }
    else{timestep = 0.0025;}
    framecount = particles[1];

    meanX = mean_pos[0]/particle_count;
    meanZ = mean_pos[2]/particle_count;
    if(particle_count > 0) rings.UpdateCentre(meanX, meanZ);

    if (enableSpaceCharge){
        rings.ClearActiveRings();

        for (auto & electron:aval.GetElectrons()){
            double xf = electron.path.back().x;
            double yf = electron.path.back().y;
            double zf = electron.path.back().z;
            double tf = electron.path.back().t;
            int status = electron.status;
            if (status==0 || status == -17) rings.AddChargedRing(xf,yf,zf,-1);
        }
        for (auto & ion:drift.GetIons()){
            double xf = ion.path.back().x;
            double yf = ion.path.back().y;
            double zf = ion.path.back().z;
            int status = ion.status;
            if (status==0 || status == -17) rings.AddChargedRing(xf,yf,zf,1);
        }
        for (auto & ion:drift.GetNegativeIons()){
            double xf = ion.path.back().x;
            double yf = ion.path.back().y;
            double zf = ion.path.back().z;
            int status = ion.status;
            if (status==0 || status == -17) rings.AddChargedRing(xf,yf,zf,-1);
        }
      
    }
    int nr;
    rings.GetNumberOfRings(nr);
    std::cout << nr << " rings.\n";

    if (electrons_remaining){
        std::cout << "Microscopic\n";
        aval.SetTimeWindow(tmin, tmin + timestep);
        aval.ResumeAvalanche();

        if (ions_remaining){
            std::cout << "MC\n";
            drift.SetTimeWindow(tmin, tmin + timestep);
            drift.ResumeAvalanche();
            tmin += timestep;
            continue;  
        }
        tmin += timestep;
    }
    else{
        break;
    }
  }
  size << i_dv << "," << particles[1] - particles[2] << "\n"; 
  LOG("Script: Done.\n");
}
