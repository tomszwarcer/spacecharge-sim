#include <array>
#include <iostream>

#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/ComponentChargedRing.hh"
#include "Garfield/Medium.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ComponentGrid2D.hh"

using namespace Garfield;

// These are here so they can be used within the functions defined below
AvalancheMC drift;
AvalancheMicroscopic aval;
ComponentChargedRing rings;
ComponentGrid2D grid;

// set time window parameters
double tmin = 0.;
double timestep = 0.05;

void userHandleIonisation(double x, double y, double z, double t, int type,
                          int level, Medium* m) {
  drift.AddIon(x, y, z,
               tmin + timestep);  // ion added at the start of the next timestep
}

void userHandleAttachment(double x, double y, double z, double t, int type,
                          int level, Medium* m) {
  drift.AddNegativeIon(
      x, y, z, tmin + timestep);  // ion added at the start of the next timestep
}

void get_mean(std::array<double, 3>& mean_pos,
              std::array<int, 3>& particle_counts, int& n_particles,
              bool enableSpaceCharge) {
  // Finds the mean of the particle positions and number
  // of particles in the avalanche. Also add the charged
  // rings if space-charge is enabled.
  //
  // particle counts is {electrons, ions, negative ions}

  mean_pos = {0., 0., 0.};
  particle_counts = {
      0,
      0,
      0,
  };

  for (const auto& electron : aval.GetElectrons()) {
    // get electron positions
    double xf = electron.path.back().x;
    double yf = electron.path.back().y;
    double zf = electron.path.back().z;
    int status = electron.status;
    // if at the end of timestep and within the region of interest:
    // status == 0 allows for the first electron to be added
    if (status == 0 || status == -17) {
      particle_counts[0]++;
      mean_pos[0] += xf;
      mean_pos[1] += yf;
      mean_pos[2] += zf;
      if (enableSpaceCharge) rings.AddChargedRing(xf, yf, zf, -1);
    }
  }
  for (const auto& ion : drift.GetIons()) {
    // get positive ion positions
    double xf = ion.path.back().x;
    double yf = ion.path.back().y;
    double zf = ion.path.back().z;
    int status = ion.status;
    // if at the end of timestep and within the region of interest:
    if (status == 0 || status == -17) {
      particle_counts[1]++;
      mean_pos[0] += xf;
      mean_pos[1] += yf;
      mean_pos[2] += zf;
      if (enableSpaceCharge) rings.AddChargedRing(xf, yf, zf, 1);
    }
  }
  for (const auto& ion : drift.GetNegativeIons()) {
    // get negative ion positions
    double xf = ion.path.back().x;
    double yf = ion.path.back().y;
    double zf = ion.path.back().z;
    int status = ion.status;
    // if at the end of timestep and within the region of interest:
    if (status == 0 || status == -17) {
      particle_counts[2]++;
      mean_pos[0] += xf;
      mean_pos[1] += yf;
      mean_pos[2] += zf;
      if (enableSpaceCharge) rings.AddChargedRing(xf, yf, zf, -1);
    }
  }

  n_particles = particle_counts[0] + particle_counts[1] + particle_counts[2];
  for (int j = 0; j < 3; ++j) {
    mean_pos[j] /= (double)n_particles;
  }
}

int main(int argc, char* argv[]) {

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

  // Setup the gas
  MediumMagboltz gas;
  gas.SetComposition("ar", 93, "co2", 7);  // [%]
  gas.SetTemperature(293.15);              // [K]
  gas.SetPressure(760.);                   // [Torr]
  gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");
  gas.Initialise(true);
  gas.EnablePenningTransfer();

  // Parallel plate configuration
  const double posBottomPlane = 0.;    // cm
  const double posTopPlane = 128.e-4;  // cm

  // Parallel plate field
  ComponentAnalyticField pp;
  pp.SetMedium(&gas);
  pp.AddPlaneY(posBottomPlane, 0.);
  pp.AddPlaneY(posTopPlane, dv);

  rings.SetArea(-0.02, posBottomPlane, -0.02, 0.02, posTopPlane, 0.02);
  rings.SetMedium(&gas);

  grid.SetArea(-0.02,posBottomPlane,-0.02,0.02,posTopPlane,0.02);
  grid.Set2dGrid(posBottomPlane,posTopPlane,200,0.02,80);
  grid.SetMedium(&gas);
  
  // Create and setup the sensor
  Sensor sensor;
  sensor.AddComponent(&pp);
  sensor.AddComponent(&grid);
  aval.SetSensor(&sensor);
  drift.SetSensor(&sensor);

  const bool enableDebug = false;
  const bool enableSpaceCharge = false;
  const bool enableAdaptiveTimestep = true;

  aval.SetUserHandleIonisation(userHandleIonisation);
  aval.SetUserHandleAttachment(userHandleAttachment);

  std::string path_end = ".txt";
  if (enableSpaceCharge){
      std::string path_begin = "/afs/cern.ch/work/t/tszwarce/grid/gain_output/on/size_sc_";
    
      size.open(path_begin+s_run_number+s_dv+path_end);
      //ion_times.open("/afs/cern.ch/work/t/tszwarce/nogrid-static/ion_times.txt");
  }
  else{
      std::string path_begin = "/afs/cern.ch/work/t/tszwarce/grid/gain_output/off/size_no_sc_";
      size.open(path_begin+s_run_number+s_dv+path_end);
  }   

  if (enableDebug) {
    rings.EnableDebugging();
    aval.EnableDebugging();
    drift.EnableDebugging();
    grid.EnableDebugging();
  }

  // Add initial electron
  const double x0 = 0.;        // cm
  const double y0 = 127.9e-4;  // cm
  const double z0 = 0.;        // cm
  const double t0 = 0.;        // ns
  const double e0 = 0.1;       // eV
  const int N_electrons = 200;
  for (int i = 0; i < N_electrons; ++i) aval.AddElectron(x0, y0, z0, 0., 0.1);

  // Variables used in adaptive timestepping
  int current_ion_count = 0;
  int current_electron_count = 0;
  int num_new_ions;

  int frame_number = -1;

  int n_particles;
  std::array<double, 3> mean_pos;
  std::array<int, 3> particle_counts;

  while (1) {
    ++frame_number;

    std::cout << "\nFrame " << frame_number << "\n";

    n_particles = 0;
    if (enableSpaceCharge) rings.ClearActiveRings();

    get_mean(mean_pos, particle_counts, n_particles, enableSpaceCharge);

    std::cout << "Particles in the simulation:\n    e-: " << particle_counts[0]
              << ", i+: " << particle_counts[1]
              << ", i-: " << particle_counts[2] << "\n";

    if (!enableAdaptiveTimestep) timestep = 0.05;          
    else if (current_ion_count > 0) {
      num_new_ions = particle_counts[1] - current_ion_count;
      std::cout << num_new_ions << " new ions created in the last timestep\n";

      double ratio = num_new_ions / (double)current_electron_count;

      // Adaptive time step tolerance
      // This limits the number of ionisations that can happen without the
      // corresponding rings being added to the simulation. Smaller number means
      // shorter dt.
      double tolerance = 0.2;

      if (ratio <= tolerance) {
        timestep = timestep / (1. - ratio);
      } else if (ratio > tolerance) {
        timestep = timestep / (1. + ratio);
      }
    }

    std::cout << "Timestep: " << timestep
              << " ns\n";

    current_ion_count = particle_counts[1];
    current_electron_count = particle_counts[0];

    // Update the symmetry axis of the rings
    if (n_particles > 0) {
      rings.UpdateCentre(mean_pos[0], mean_pos[2]);
      grid.UpdateCentre(mean_pos[0], mean_pos[2]);
    }

    size_t nr;
    nr = rings.GetNumberOfRings();
    std::cout << nr << " rings in the simulation.\n";

    grid.AddElectricField(&rings);

    if (particle_counts[0] > 0) {
      std::cout << "Simulating electrons: AvalancheMicroscopic...\n";
      aval.SetTimeWindow(tmin, tmin + timestep);
      aval.ResumeAvalanche();

      if (particle_counts[1] > 0 || particle_counts[2] > 0) {
        std::cout << "Simulating ions: AvalancheMC...\n";
        drift.SetTimeWindow(tmin, tmin + timestep);
        drift.ResumeAvalanche();
        tmin += timestep;
        continue;
      }
      tmin += timestep;
    } else {
      std::cout << "\nNo electrons remaining: done.\n";
      break;
    }
  }
  get_mean(mean_pos, particle_counts, n_particles, enableSpaceCharge);
  std::cout << "\nFinal numbers of particles:\n           "
            << particle_counts[0] << " electrons\n           "
            << particle_counts[1] << " positive ions\n           "
            << particle_counts[2] << " negative ions\n";
  size << i_dv << "," << particle_counts[1] - particle_counts[2] << "\n"; 
}
