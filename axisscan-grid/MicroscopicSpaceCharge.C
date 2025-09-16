#include <iostream>
#include <fstream>
#include <numeric>
#include <TCanvas.h>
#include <TApplication.h>

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ComponentChargedRing.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ComponentGrid2D.hh"


using namespace Garfield;

// These are here so they can be used within the functions defined below
AvalancheMC drift;
AvalancheMicroscopic aval;
ComponentChargedRing rings;
ComponentChargedRing anode_rings;
ComponentGrid2D grid;

// set time window parameters
double tmin = 0.;
double timestep = 0.05;

double anode_pos = 1.e-8;

void userHandleIonisation(double x, double y, double z, double t, int type, int level, Medium * m){
    drift.AddIon(x,y,z,tmin + timestep); // ion added at the start of the next timestep
}

void userHandleAttachment(double x, double y, double z, double t, int type, int level, Medium * m){
    drift.AddNegativeIon(x,y,z,tmin + timestep); // ion added at the start of the next timestep
}

std::ofstream efile_metallic("efile_main.txt");
std::ofstream efile_resistive("efile_combined.txt");
std::ofstream pfile("pfile.txt");
int iF;
int n_anode_rings = 0;


void get_mean(std::array<double,3> & mean_pos, std::array<int,3> & particle_counts, int & n_particles, bool enableSpaceCharge){
    // Finds the mean of the particle positions and number
    // of particles in the avalanche. Also add the charged
    // rings if space-charge is enabled.
    //
    // particle counts is {electrons, ions, negative ions}

    mean_pos = {0., 0. ,0.};
    particle_counts = {0, 0, 0,};

    for (const auto& electron : aval.GetElectrons()) {
        // get electron positions
        double xf = electron.path.back().x;
        double yf = electron.path.back().y;
        double zf = electron.path.back().z;
        int status = electron.status;
        // if at the end of timestep and within the region of interest:
        // status == 0 allows for the first electron to be added
        if (status == -5 && yf < anode_pos) {
            n_anode_rings++;
            particle_counts[0]++;
            mean_pos[0] += xf;
            mean_pos[1] += yf;
            mean_pos[2] += zf;
            if (enableSpaceCharge) anode_rings.AddChargedRing(xf,0.,zf,-1);
        }
        else if (status == 0 || status == -17) {
            particle_counts[0]++;
            mean_pos[0] += xf;
            mean_pos[1] += yf;
            mean_pos[2] += zf;
            if (enableSpaceCharge) rings.AddChargedRing(xf,yf,zf,-1);
            pfile << iF<<",e," << yf << "\n";
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
            if (enableSpaceCharge) rings.AddChargedRing(xf,yf,zf,1);
            pfile << iF <<",i," << yf << "\n";
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
            if (enableSpaceCharge) rings.AddChargedRing(xf,yf,zf,-1);
        }
    }

    n_particles = particle_counts[0] + particle_counts[1] + particle_counts[2];
    for (int j = 0; j < 3; ++j) {
        mean_pos[j] /= (double)n_particles;
    }

}

int main(int argc, char* argv[]) {
  //TApplication app("app", &argc, argv);

  // Setup the gas
  MediumMagboltz gas;
  gas.SetComposition("ar", 93, "co2", 7); // [%]
  gas.SetTemperature(293.15); // [K]
  gas.SetPressure(760.); // [Torr]
  gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");
  gas.Initialise(true);
  gas.EnablePenningTransfer(); 
  
  // Parallel plate configuration
  const double posBottomPlane = 0.; // cm
  const double posTopPlane = 128.e-4; // cm
  
  // Parallel plate field
  ComponentAnalyticField pp; 
  int voltage = -500;
  pp.SetMedium(&gas);
  pp.AddPlaneY(posBottomPlane, 0.);
  pp.AddPlaneY(posTopPlane, voltage);
  
  rings.SetArea(-0.02,posBottomPlane,-0.02,0.02,posTopPlane,0.02);
  rings.SetMedium(&gas);
  rings.SetSpacingTolerance(0.00005);
  rings.SetSelfFieldTolerance(0.00001);
  anode_rings.SetArea(-0.02,posBottomPlane,-0.02,0.02,posTopPlane,0.02);
  anode_rings.SetMedium(&gas);
  anode_rings.SetSpacingTolerance(0.00005);
  anode_rings.SetSelfFieldTolerance(0.00001);

  grid.SetArea(-0.02,posBottomPlane,-0.02,0.02,posTopPlane,0.02);
  grid.Set2dGrid(posBottomPlane,posTopPlane,100,0.02,100);
  grid.SetMedium(&gas);

  // Create and setup the sensor
  Sensor sensor;
  sensor.AddComponent(&pp);
  sensor.AddComponent(&grid);
  aval.SetSensor(&sensor);  
  drift.SetSensor(&sensor);

  const bool enableDebug = false;
  const bool enableSpaceCharge = true;
  const bool plotField = false;
  const bool plotDrift = false;

  aval.SetUserHandleIonisation(userHandleIonisation);
  aval.SetUserHandleAttachment(userHandleAttachment);

  if (enableDebug) { 
    //rings.EnableDebugging();
    aval.EnableDebugging();
    drift.EnableDebugging();
  }
  

  // Add initial electron
  const double x0 = 0.; // cm
  const double y0 = 127.9e-4; // cm
  const double z0 = 0.; // cm
  const double t0 = 0.; // ns
  const double e0 = 0.1; // eV
  for (int i = 0; i < 50; ++i) aval.AddElectron(x0,y0,z0,0.,0.1); 
  

  // Variables used in adaptive timestepping
  int current_ion_count = 0;
  int current_electron_count = 0;
  int num_new_ions;


  int frame_number = -1;

  int n_particles;
  std::array<double,3> mean_pos;
  std::array<int,3> particle_counts;
  
  int axis_samples = 150;

  while (1) {
    ++frame_number;
    iF = frame_number;

    std::cout << "\nFrame " << frame_number << "\nTimestep: " << timestep << " ns\n\n";
    
    n_particles = 0;
    if (enableSpaceCharge) rings.ClearActiveRings();

    get_mean(mean_pos, particle_counts, n_particles, enableSpaceCharge);
    

    std::cout << "Particles in the simulation:\n    e-: "
              << particle_counts[0] << ", i+: " << particle_counts[1]
              << ", i-: " << particle_counts[2] << "\n";

    // Update the symmetry axis of the rings
    if (n_particles > 0) {
        rings.UpdateCentre(mean_pos[0], mean_pos[2]);
        anode_rings.UpdateCentre(mean_pos[0],mean_pos[2]);
        grid.UpdateCentre(mean_pos[0], mean_pos[2]);
    }

    grid.AddElectricField(&anode_rings);
    grid.AddElectricField(&rings);

    size_t nr,nra;
    nr = rings.GetNumberOfRings();
    std::cout << nr << " rings in the simulation.\n";
    std::cout << n_anode_rings << " electrons at anode\n";
    for (int i = 0; i < n_anode_rings; i++) {
        pfile << iF<<",e," << 0.0 << "\n";
    }

    double d = (posTopPlane-posBottomPlane)/axis_samples;
    Medium* m = nullptr;
    double ex,ey,ez,v,y;
    double ex_a,ey_a,ez_a;
    double ex_p,ey_p,ez_p;
    int stat;

    for (int i=0; i<axis_samples;++i){
        y = posBottomPlane + i*d;
        rings.ElectricField(mean_pos[0]+0.0001,y,mean_pos[2],ex,ey,ez,m,stat);
        anode_rings.ElectricField(mean_pos[0]+0.0001,y,mean_pos[2],ex_a,ey_a,ez_a,m,stat);
        pp.ElectricField(mean_pos[0]+0.0001,y,mean_pos[2],ex_p,ey_p,ez_p,m,stat);
        ex += ex_p;
        ey += ey_p;
        ez += ez_p;
        efile_metallic << iF << "," << y << "," << ex << "," << ey << "," << ez << "," << std::sqrt(ex*ex+ey*ey+ez*ez) << "\n";
        ex += ex_a;
        ey += ey_a;
        ez += ez_a;
        efile_resistive << iF << "," << y << "," << ex << "," << ey << "," << ez << "," << std::sqrt(ex*ex+ey*ey+ez*ez) << "\n";
    }

    if (particle_counts[0] > 0){
        std::cout << "Simulating electrons: AvalancheMicroscopic...\n";
        aval.SetTimeWindow(tmin, tmin + timestep);
        aval.ResumeAvalanche();

        if (particle_counts[1] > 0 || particle_counts[2] > 0){
            std::cout << "Simulating ions: AvalancheMC...\n";
            drift.SetTimeWindow(tmin, tmin + timestep);
            drift.ResumeAvalanche();
            tmin += timestep;
            continue;  
        }
        tmin += timestep;
    }
    else{
        std::cout << "\nNo electrons remaining: done.\n";
        break;
    }
  }
  get_mean(mean_pos, particle_counts, n_particles, enableSpaceCharge);
  std::cout << "\nFinal numbers of particles:\n           " <<  particle_counts[0]
            << " electrons\n           " << particle_counts[1]
            << " positive ions\n           " << particle_counts[2]
            << " negative ions\n";
  efile_metallic.close();
  efile_resistive.close();
  pfile.close();
  //app.Run(true);
}
