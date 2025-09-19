#include <TApplication.h>
#include <TCanvas.h>
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewField.hh"

#include <array>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentChargedRing.hh"
#include "Garfield/Medium.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ComponentGrid2D.hh"
#include "Garfield/ComponentComsol.hh"

using namespace Garfield;

Sensor sensor;

// set time window parameters
double tmin = 9999.; // this is later updated to the shortest time 
                     // drawn from a distribution (inelegant, I know...)
double timestep = 0.05;

// Parallel plate configuration
const double posBottomPlane = 0.;    // cm
double anode_pos = 1.e-8;
const double posTopPlane = 128.e-4;  // cm
std::vector<std::array<double,4>> start_coords;
int num_initial_electrons = 0;

bool enableResistiveAnode;
bool enableSpaceCharge;
bool enableDebug;


int current_system_index = 0;

struct System {

  // One System object for each ring system

  AvalancheMicroscopic aval;
  AvalancheMC drift;
  ComponentChargedRing rings;
  ComponentChargedRing anode_rings;
  ComponentGrid2D grid;

  std::array<double, 3> mean_pos;
  std::array<int, 3> particle_counts;
  int n_particles;
  int n_anode_electrons;

  System(MediumMagboltz * gas) {

    // setup
    if (enableDebug) {
      //rings.EnableDebugging();
      //anode_rings.EnableDebugging();
      //grid.EnableDebugging();
      aval.EnableDebugging();
      drift.EnableDebugging();
    }

    if (enableSpaceCharge) {
      sensor.AddComponent(&grid);
      grid.SetArea(-0.1,posBottomPlane,-0.1,0.1,posTopPlane,0.1);
      grid.Set2dGrid(posBottomPlane,posTopPlane,200,0.1,80);
      grid.SetMedium(gas);
      rings.SetArea(-0.1, posBottomPlane, -0.1, 0.1, posTopPlane, 0.1);
      rings.SetMedium(gas);
      if (enableResistiveAnode) {
        anode_rings.SetArea(-0.1,posBottomPlane,-0.1,0.1,posTopPlane,0.1);
        anode_rings.SetMedium(gas);
      }
    }

    aval.SetSensor(&sensor);
    drift.SetSensor(&sensor);
  }

  void add_electron_rings() {
    // update mean position and add rings

    for (const auto& electron : aval.GetElectrons()) {
      // get electron positions
      double xf = electron.path.back().x;
      double yf = electron.path.back().y;
      double zf = electron.path.back().z;
      int status = electron.status;
      // if at the end of timestep and within the region of interest:
      // status == 0 allows for the first electron to be added
      if (status == -5 && yf < anode_pos && enableResistiveAnode) {
        n_anode_electrons++;
        particle_counts[0]++;
        mean_pos[0] += xf;
        mean_pos[1] += yf;
        mean_pos[2] += zf;
        if (enableSpaceCharge) anode_rings.AddChargedRing(xf,0.,zf,-1);
      }
      if (status == 0 || status == -17) {
        particle_counts[0]++;
        mean_pos[0] += xf;
        mean_pos[1] += yf;
        mean_pos[2] += zf;
        if (enableSpaceCharge) rings.AddChargedRing(xf, yf, zf, -1);
      }
    }
  }

  void add_ion_rings() {
    // update mean position and add rings

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
  } 

  void get_mean() {
    mean_pos = {0., 0., 0.};
    particle_counts = {0, 0, 0};
    add_electron_rings();
    add_ion_rings();
    n_particles = particle_counts[0] + particle_counts[1] + particle_counts[2];
    for (int j = 0; j < 3; ++j) {
      mean_pos[j] /= (double)n_particles;
    }
  }

};

std::vector<System> systems;

void userHandleIonisation(double x, double y, double z, double t, int type,
                            int level, Medium* m) {
  systems[current_system_index].drift.AddIon(x, y, z, tmin + timestep);  // ion added at the start of the next timestep
}
void userHandleAttachment(double x, double y, double z, double t, int type,
                          int level, Medium* m) {
  systems[current_system_index].drift.AddNegativeIon(x, y, z, tmin + timestep);  // ion added at the start of the next timestep
}

void add_initial_electrons(std::string gasfile_path, double ex, double ey, double ez,std::vector<std::array<double,4>> & start_coords) {
  double vx,vy,vz; // drift velocities
  double dl,dt; // diffusion coeffs
  double interaction_distance,longitudinal_sd,transverse_sd,time_sd,total_distance,mean_time,drift_velocity;

  // setup random
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(0.,0.3);

  // we use a separate object for the gas table as it slows down AvalancheMC
  // when we load it in the original gas object.
  MediumMagboltz gastable;
  gastable.SetComposition("ar", 93, "co2", 7);  // [%]
  gastable.SetTemperature(293.15);              // [K]
  gastable.SetPressure(760.);                   // [Torr]
  gastable.Initialise(true);
  gastable.LoadGasFile(gasfile_path);
  gastable.ElectronVelocity(ex,ey,ez,0.,0.,0.,vx,vy,vz);
  gastable.ElectronDiffusion(ex,ey,ez,0.,0.,0.,dl,dt);

  interaction_distance = dist(gen);
  drift_velocity = std::sqrt(vx*vx+vy*vy+vz*vz);
  transverse_sd = dt * std::sqrt(interaction_distance);
  longitudinal_sd = dl * std::sqrt(interaction_distance);
  time_sd = longitudinal_sd / drift_velocity;
  mean_time = interaction_distance / drift_velocity;
  std::normal_distribution<double> time_norm_dist(mean_time,time_sd);  
  std::normal_distribution<double> position_norm_dist(0.,transverse_sd);
  const double y0 = 150e-4; // cm
  const double e0 = 0.1; // eV
  const int N_electrons = 2; // 212 for Fe-55

  // draw from the distributions
  for (int run_counter = 0; run_counter < N_electrons; ++run_counter) {
    int chosen_value = dist(gen);
    double x0 = position_norm_dist(gen);
    double z0 = position_norm_dist(gen);
    double t0 = time_norm_dist(gen);
    while (t0 < 0.) t0 = time_norm_dist(gen);
    std::cout << "Adding electron at (x,y,z,t) = (" << x0 << "," << y0 << "," << z0 << "," << t0 <<")\n";
    Medium* medium = sensor.GetMedium(x0, y0, z0);
    if (medium != nullptr && medium->IsDriftable()) {
      if (t0 < tmin) tmin = t0;
      start_coords.push_back({x0,y0,z0,t0});
      num_initial_electrons++;
    }
    else {std::cout << "Electron " << run_counter << " is not in a valid medium\n";}
  }
  std::cout << num_initial_electrons << " electrons will be added\nFirst electron arrives at t = " << tmin << "ns\n";
  // make the first time window just before the first electron arrives
  tmin = tmin - 0.01 * timestep;
}

int main(int argc, char* argv[]) {

  // Check if the right number of arguments are provided
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0] << " <process_index>  <voltage>  <space_charge_on>  <resisive_anode_on>" << std::endl;
    return 1;
  }

  TApplication app("app", &argc, argv);

  int run_number = std::atoi(argv[1]);
  std::string s_run_number = argv[1];
  int i_dv = std::atoi(argv[2]);
  double dv = std::atof(argv[2]);
  std::string s_dv = std::to_string(i_dv);
  dv = -1*dv;
  s_dv = "_"+s_dv;
  std::ofstream size;
  enableSpaceCharge = std::atoi(argv[3]);
  enableResistiveAnode = std::atoi(argv[4]);
  if (enableResistiveAnode && !enableSpaceCharge) throw std::runtime_error("WARNING: resistive anode enabled but space charge is disabled.\n");
  
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

  // Defining the field map.
  // Load the field map.
  ComponentComsol fm;
  std::string potential = "Potential";
  std::string path_end = ".txt";
  potential = potential + s_dv + path_end;
  fm.Initialise("mesh.mphtxt","dielectrics.dat",potential, "mm");
  fm.EnablePeriodicityX();
  fm.EnablePeriodicityZ();
  fm.PrintRange();
  const unsigned int nMaterials = fm.GetNumberOfMaterials();
  for (unsigned int i = 0; i < nMaterials; ++i) {
      const double eps = fm.GetPermittivity(i);
      if(eps==1) fm.SetMedium(i, &gas);
  }
  fm.PrintMaterials();
  
  sensor.AddComponent(&fm);

  enableDebug = false;
  bool plotDrift = true;

  // Plot drift lines
  ViewDrift driftView;
  if (plotDrift) {
    driftView.SetCollisionMarkerSize(0.0000000000000000001);
    driftView.SetColourIonisations(7);
    driftView.SetPlane(0, 0, 1, 0, 0, 0);
    driftView.SetArea(-0.03, posBottomPlane, 0.03, posTopPlane);
  }

  // setup the file for exporting the gain
  if (enableSpaceCharge){
      std::string path_begin = "/Users/tomszwarcer/Documents/CERN/MMBasic/Code/gain_output/on/size_sc_";
    
      size.open(path_begin+s_run_number+s_dv+path_end);
  }
  else{
      std::string path_begin = "/Users/tomszwarcer/Documents/CERN/MMBasic/Code/gain_output/off/size_no_sc_";
      size.open(path_begin+s_run_number+s_dv+path_end);
  }   

  // this is used to generate the initial electrons
  std::string gasfile_start = "gas";
  std::string gasfile_path = gasfile_start + s_dv + path_end;
  double ex,ey,ez;
  int stat; // status
  Medium * m = nullptr;
  sensor.ElectricField(0.,0.15 + posTopPlane,0.,ex,ey,ez,m,stat);
  int num_electrons_added = 0;
  add_initial_electrons(gasfile_path,ex,ey,ez,start_coords);

  int frame_number = -1;
  bool electrons_remaining;

  // main loop
  while (1) {
    ++frame_number;

    std::cout << "\nFrame " << frame_number << "\n";

    std::cout << "Time: " << tmin << "\nTimestep: " << timestep
              << " ns\n";

    // add starting electrons
    for (const std::array<double,4> coords : start_coords) {
      if ((tmin <= coords[3]) && (coords[3] < tmin + timestep)) {
        std::cout << "Adding electron " << num_electrons_added << " of " << num_initial_electrons - 1 << " at t = " << coords[3] << "\n";
        // increase the size of the ring system vectors
        systems.emplace_back(&gas);
        // set up user handles
        systems.back().aval.SetUserHandleIonisation(userHandleIonisation);
        systems.back().aval.SetUserHandleAttachment(userHandleAttachment);
        // set up plotting
        if (plotDrift) {
          systems.back().aval.EnablePlotting(&driftView);
          systems.back().drift.EnablePlotting(&driftView);
        }
        // add the electron to the last ring system
        systems.back().aval.AddElectron(coords[0],coords[1],coords[2],coords[3],0.1);
        num_electrons_added++;
      }
    }

    int counter = 0;
    for (System & system : systems) {

      if (enableSpaceCharge) sensor.EnableMagneticField(counter + 1,false); // When I remove this, it crashes... don't know why.
      current_system_index = counter;
      system.rings.ClearActiveRings();
      system.get_mean();

      std::cout << "    Ring system " << counter << ":\n    e-: " << system.particle_counts[0]
                << ", i+: " << system.particle_counts[1]
                << ", i-: " << system.particle_counts[2] << "\n";



      // update symmetry axis          
      if (system.n_particles > 0 && enableSpaceCharge) {
        system.rings.UpdateCentre(system.mean_pos[0], system.mean_pos[2]);
        system.anode_rings.UpdateCentre(system.mean_pos[0], system.mean_pos[2]);
        system.grid.UpdateCentre(system.mean_pos[0], system.mean_pos[2]);

        size_t nr;
        nr = system.rings.GetNumberOfRings();
        std::cout << nr << " rings in system " << counter << ".\nAdding fields of system " << counter << " to grid...\n";

        // add ring fields to grid
        system.grid.AddElectricField(&system.rings);

        if (enableResistiveAnode) {
          std::cout << system.n_anode_electrons << " electrons at anode system " << counter << "\nAdding anode fields of system " << counter << " to grid...\n";
          
          // add anode ring fields to grid
          system.grid.AddElectricField(&system.anode_rings);
        }
      }
      electrons_remaining = false;

      // main simulation
      if (system.particle_counts[0] > 0) {
        std::cout << "    Simulating electrons for ring system " << counter << ": AvalancheMicroscopic...\n";
        system.aval.SetTimeWindow(tmin, tmin + timestep);
        system.aval.ResumeAvalanche();
        electrons_remaining = true;
        std::cout << "    AvalancheMicroscopic done.\n";

        if (system.particle_counts[1] > 0 || system.particle_counts[2] > 0) {
          std::cout << "    Simulating ions for ring system " << counter << ": AvalancheMC...\n";
          system.drift.SetTimeWindow(tmin, tmin + timestep);
          system.drift.ResumeAvalanche();
          std::cout << "    AvalancheMC done.\n";
        }
      }
      counter++;
    }
    // if we've gone through all the systems and there's no electrons
    if (electrons_remaining == false && counter == num_initial_electrons) {
      std::cout << "\nNo electrons remaining: done.\n";
      break;
    }
    tmin += timestep;
  }

  // export gain
  int gain_temp = 0;
  int counter = 0;
  std::cout << "Particles in the simulation:\n";
  for (System & system : systems) {
    system.get_mean();
    std::cout << "    Ring system " << counter << ":\ne-: " << system.particle_counts[0]
              << ", i+: " << system.particle_counts[1]
              << ", i-: " << system.particle_counts[2] << "\n";
    gain_temp += system.particle_counts[1] - system.particle_counts[2];
    counter++;
  }
  size << i_dv << "," << gain_temp << "\n"; 

  // export Magboltz collision infor
  std::ofstream collisions;
  collisions.open("collisions.txt");
  int n_levels = gas.GetNumberOfLevels();
  int ngas,type;
  unsigned int n_coll;
  std::string descr;
  double e;
  for (int i = 0; i < n_levels; ++i) {
    n_coll = gas.GetNumberOfElectronCollisions(i);
    gas.GetLevel(i,ngas,type,descr,e);
    collisions << n_coll << "," << ngas << "," << e << "," << descr << "\n";
  }

  if (plotDrift) {
    TCanvas* cd = new TCanvas();
    driftView.SetCanvas(cd);
    driftView.Plot(true);
  }
  app.Run(true);
}
