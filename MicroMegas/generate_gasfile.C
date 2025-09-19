#include <array>
#include <iostream>
#include <random>

#include "Garfield/Medium.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ComponentComsol.hh"

using namespace Garfield;


int main(int argc, char* argv[]) {

  // Check if the right number of arguments are provided
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <voltage> " << std::endl;
    return 1;
  }

  int i_dv = std::atoi(argv[1]);
  double dv = std::atof(argv[1]);
  std::string s_dv = std::to_string(i_dv);
  dv = -1*dv;
  s_dv = "_"+s_dv;
  
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
  
  // Create and setup the sensor
  Sensor sensor;
  sensor.AddComponent(&fm);

  // Add initial electrons

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(0.,0.3);

  double ex,ey,ez;
  double vx,vy,vz; // drift velocities
  double dl,dt; // diffusion coeffs
  Medium * m = nullptr;
  int stat; // status
  double interaction_distance,longitudinal_sd,transverse_sd,time_sd,total_distance,mean_time,drift_velocity;

  std::ofstream outfile;
  outfile.open("starting_positions.hh");

  sensor.ElectricField(0.,0.15 + posTopPlane,0.,ex,ey,ez,m,stat);
  double eMag = std::sqrt(ex*ex+ey*ey+ez*ez);
  gas.SetFieldGrid(eMag,eMag,1,false);  
  gas.GenerateGasTable(true);

  std::string gasfile_start = "gas";
  gas.WriteGasFile(gasfile_start + s_dv + path_end);
}