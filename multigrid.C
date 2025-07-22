#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <sstream>

#include <TCanvas.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <vector>
#include <iomanip>
#include <string>

#include <algorithm>
#include <cmath>

#include "Garfield/Plotting.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheGridSpaceCharge.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ComponentComsol.hh"
#include "Garfield/ComponentParallelPlate.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/Random.hh"


#define LOG(x) std::cout<<x<<std::endl

using namespace Garfield;

AvalancheMC drift;
AvalancheMicroscopic aval;
MediumMagboltz gas;
Medium * medi;

int N = 3;

AvalancheGridSpaceCharge avalsc1;
AvalancheGridSpaceCharge avalsc2;
AvalancheGridSpaceCharge avalsc3;
//AvalancheGridSpaceCharge avalsc4;
//std::vector<AvalancheGridSpaceCharge> scarray = {avalsc1,avalsc2,avalsc3,avalsc4};
std::vector<AvalancheGridSpaceCharge> scarray = {avalsc1,avalsc2,avalsc3};


// set time window parameters
double tmin = 0.;
double timestep = 0.05;

constexpr bool enable_penning = false;
double rPenning = 0.42; // from Penning transfer in argon-based gas mixtures, O Sahin et al (2010) 

std::vector<double> mean_pos = {0.,0.,0.};
int particle_count = 0;

void userHandleIonisation(double x, double y, double z, double t, int type, int level, Medium * m){
    drift.AddIon(x,y,z,tmin + timestep); // ion added after step.
}

void userHandleAttachment(double x, double y, double z, double t, int type, int level, Medium * m){
    drift.AddNegativeIon(x,y,z,tmin + timestep);
}

void add_electrons_to_grid(bool & electrons_remaining,const double & zmax){
    for (const auto& electron : aval.GetElectrons()) {
        // get electron positions
        double xf = electron.path.back().x;
        double yf = electron.path.back().y;
        double zf = electron.path.back().z;
        double tf = electron.path.back().t;
        int status = electron.status;
        // if still in a drift medium or at the end of timestep:
        if (status==0 || status == -17){
            electrons_remaining = true;
            if (yf < zmax){
                particle_count++;
                mean_pos[0] += xf;
                mean_pos[1] += yf;
                mean_pos[2] += zf;
                for (auto & a:scarray){
                    a.AddElectron(xf,yf,zf,tf,1);
                }
            }
        }
    }
}


void add_ions_to_grid(bool & ions_remaining,const double &zmax){
    for (const auto& ion : drift.GetIons()){
        double xf = ion.path.back().x;
        double yf = ion.path.back().y;
        double zf = ion.path.back().z;
        double tf = ion.path.back().t;
        double ti = ion.path.front().t;
        int status = ion.status;

        if (status==0 || status == -17){
            ions_remaining = true;
            if (yf < zmax){
                particle_count++;
                mean_pos[0] += xf;
                mean_pos[1] += yf;
                mean_pos[2] += zf;
                for (auto & a:scarray){
                    a.AddPositiveIon(xf,yf,zf,tf,1);
                }
            } 
        }
    }
}

void add_negative_ions_to_grid(bool & ions_remaining,const double &zmax){
    for (const auto& ion : drift.GetNegativeIons()){
        double xf = ion.path.back().x;
        double yf = ion.path.back().y;
        double zf = ion.path.back().z;
        double tf = ion.path.back().t;
        int status = ion.status;

        if (status==0 || status == -17){
            ions_remaining = true;
            if (yf < zmax) {
                for (auto & a:scarray){
                    a.AddNegativeIon(xf,yf,zf,tf,1);
                }
                particle_count++;
                mean_pos[0] += xf;
                mean_pos[1] += yf;
                mean_pos[2] += zf;
            }
        }
    }
}

int main(int argc, char * argv[]){

    std::ofstream outfile;
    outfile.open("outfile.txt");

    std::ofstream gridfile;
    gridfile.open("gridfile.txt");

    constexpr bool debug = false;
    constexpr bool enableSpaceCharge = true;
    
    // Defining the field map.
    // Load the field map.
    ComponentComsol fm;
    fm.Initialise("mesh.mphtxt","dielectrics.dat","Potential.txt", "mm");
    fm.EnablePeriodicityX();
    fm.EnablePeriodicityZ();
    fm.PrintRange();
   // fm.EnableConvergenceWarnings(false);

    
    // Setup the gas.
    
    gas.SetComposition("ar", 93, "co2", 7); // [%]
    gas.SetTemperature(293.15); // [K]
    gas.SetPressure(760); // [Torr]
    gas.Initialise(true);
    

    if (enable_penning) gas.EnablePenningTransfer();
    
    // Load the ion mobilities.
    const std::string path = std::getenv("GARFIELD_INSTALL");
    gas.LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");

    const unsigned int nMaterials = fm.GetNumberOfMaterials();
    for (unsigned int i = 0; i < nMaterials; ++i) {
        const double eps = fm.GetPermittivity(i);
        LOG("eps = "<< eps);
        if(eps==1) fm.SetMedium(i, &gas);
    }
    fm.PrintMaterials();

    // Create the sensor.
    Sensor s1;
    Sensor s2;
    Sensor s3;
    //Sensor s4;
    std::vector<Sensor *> sensarray = {&s1,&s2,&s3};//,&s4};
    for (auto &s:sensarray){
        s->AddComponent(&fm);
        if(debug) s->EnableDebugging();
    } 

    for (int i =0;i<N;i++){
        scarray[i].SetSensor(sensarray[i]);
        if (debug) scarray[i].EnableDebugging();
    }
    
    //AvalancheMicroscopic aval;
    aval.SetSensor(&s3);
    aval.EnableSignalCalculation();
    aval.SetUserHandleIonisation(userHandleIonisation);
    aval.SetUserHandleAttachment(userHandleAttachment);
    if(debug) aval.EnableDebugging();
    
    //AvalancheMC drift;
    drift.SetSensor(&s3);
    drift.SetTimeSteps(timestep);
    if(debug) drift.EnableDebugging();
    

    // We create a fake RPC which is just a gas gap with dV = 0
    // This fake RPC lives under the mesh of the micromegas.
    const double mesh_pos = 0.012;
    const double anode_pos = 0.;
    double d_gas = mesh_pos - anode_pos; // cm
    double d_top = 0.;
    double d_bottom = 0.;
    std::vector<double> layers = {d_gas,d_bottom};
    double y_mid = d_gas / 2;
    double e_gas = 1.;
    double e_top = 2.;
    double e_bottom = 2.;
    std::vector<double> eps = {e_gas,e_bottom};
    double voltage = 0.;

    ComponentParallelPlate cmp1;
    ComponentParallelPlate cmp2;
    ComponentParallelPlate cmp3;
    //ComponentParallelPlate cmp4;
    std::vector<ComponentParallelPlate> cmparray = {cmp1,cmp2,cmp3};//,cmp4};
  
  
    
    for (int i=0; i<N;i++){
        if (enableSpaceCharge) cmparray[i].EnableSpaceCharge();
        cmparray[i].SetMedium(&gas);
        cmparray[i].Setup(int(layers.size()), eps, layers, voltage, {});
        cmparray[i].SetAvalancheGridSpaceChargeObject(&scarray[i]);
        sensarray[i]->AddComponent(&cmparray[i]);
    }
    

  // Starting parameters for the electron
  const double x0 = 0.0;
  const double y0 = 0.013; // mesh is at y=0.012 more or less
  const double z0 = 0.;
  const double t0 = 0.;
  const double e0 = 0.1;

  // grid parameters
  const double zmax = d_gas;
  const double zmin = anode_pos; // < anode pos with a bufffer region
  const double rmax = std::max(std::sqrt(x0*x0 + z0*z0)*1.25,y0*1.25); // < with a buffer region


  // define grid spacings of simultaneous grids
  std::vector<int> steps_v = {50,75,100};


  std::vector<double> zstepsizes,rstepsizes;
  for (int i=0; i<N; i++){
    zstepsizes.push_back((zmax-zmin)/steps_v[i]);
    rstepsizes.push_back(rmax/steps_v[i]);
  }  

  gridfile << "zmin,zmax,zsteps,rmax,rsteps" << std::endl;
  aval.AddElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
  for (int i = 0; i < N; i++){
    scarray[i].Set2dGrid(zmin,zmax,steps_v[i],rmax,steps_v[i]);
    gridfile << zmin << "," << zmax << "," << steps_v[i] << "," <<  rmax << "," << steps_v[i] << std::endl;
    scarray[i].TransportOffGrid();
    scarray[i].EnableSpaceChargeEffect();
  }

  
  // set number of steps
  const int N_timesteps = 50;

  bool electrons_remaining;
  bool ions_remaining;

  // for each time step
  double x,y,z,ex,ey,ez;
  int status;
  for (int step = 0; step < N_timesteps; ++step){


    std::cout << "step " << step << "\n";
    electrons_remaining = false;
    ions_remaining = false;
    for (auto & a:scarray) a.ClearGrid();

    particle_count = 0;
    mean_pos = {0.,0.,0.};
    add_electrons_to_grid(electrons_remaining, zmax);
    add_ions_to_grid(ions_remaining,zmax);
    add_negative_ions_to_grid(ions_remaining,zmax);
    if(particle_count > 0){
        for (auto & a:scarray){
            a.UpdateMeanPosition(mean_pos[0]/particle_count,mean_pos[1]/particle_count,mean_pos[2]/particle_count);
        }
    }
    for (auto & a:scarray) a.UpdateFieldOnGrid();

    if (step % 5 == 0 && step>0){
        for (int i=0;i<N;i++){
            for (int j =-200; j <=200;j++){
                x = (mean_pos[0]/particle_count)+j/70000.;
                y = mean_pos[1]/particle_count;
                z = (mean_pos[2]/particle_count)+j/70000.;
                medi = cmparray[i].GetMedium(x,y,z);
                cmparray[i].ElectricField(x,y,z,ex,ey,ez,medi,status);
                outfile << step << "," << steps_v[i] << "," << x << "," << y << "," << z << "," << ex << "," << ey << "," << ez << std::endl;
            }
        }
    }

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
        if (ions_remaining){
            std::cout << "MC\n";
            drift.SetTimeWindow(tmin, tmin + timestep);
            drift.ResumeAvalanche();  

            tmin += timestep;
        }
        else{
            break;
        }
    }
  }
}
