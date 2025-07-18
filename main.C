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
AvalancheGridSpaceCharge avalsc;
MediumMagboltz gas;

// set time window parameters
double tmin = 0.;
double timestep = 0.05;

constexpr bool enable_penning = false;
double rPenning = 0.42; // from Penning transfer in argon-based gas mixtures, O Sahin et al (2010) 

std::vector<double> mean_pos = {0.,0.,0.};
int particle_count = 0;

void userHandleIonisation(double x, double y, double z, double t, int type, int level, Medium * m){
    if (type == 4 && enable_penning){
        drift.AddIon(x,y,z,tmin+timestep);
    }
    else{
        drift.AddIon(x,y,z,tmin + timestep); // ion added after step.
    }
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
                avalsc.AddElectron(xf,yf,zf,tf,1);
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
                avalsc.AddPositiveIon(xf,yf,zf,tf,1);
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
                avalsc.AddNegativeIon(xf,yf,zf,tmin,1);
                particle_count++;
                mean_pos[0] += xf;
                mean_pos[1] += yf;
                mean_pos[2] += zf;
            }
        }
    }
}

int main(int argc, char * argv[]){

    TApplication app("app", &argc, argv);

    std::ofstream ion_pos;
    ion_pos.open("ion_pos.txt");
    std::ofstream dfile;
    dfile.open("filenames_drift.txt");
    std::ofstream ffile;
    ffile.open("filenames_field.txt");
        
    plottingEngine.SetDefaultStyle();
    
    constexpr bool debug = false;
    constexpr bool plotDrift = true;
    constexpr bool plotSignal = false;
    constexpr bool plotField = true;
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
    LOG("GetNumberOfMaterials = "<< nMaterials);
    for (unsigned int i = 0; i < nMaterials; ++i) {
        const double eps = fm.GetPermittivity(i);
        LOG("eps = "<< eps);
        if(eps==1) fm.SetMedium(i, &gas);
    }
    fm.PrintMaterials();

    constexpr double pitch = 200.e-4; //cm
    constexpr double halfpitch = 0.5 * pitch;

    ComponentParallelPlate cmp;
    if (enableSpaceCharge) cmp.EnableSpaceCharge();
    cmp.SetMedium(&gas);

    // Create the sensor.
    Sensor sensor;
    sensor.AddComponent(&fm);
    // for(int i = 0;i<5;i++) sensor.AddElectrode(&wField, label[i]);


    TCanvas* cfield = new TCanvas("cfield", "", 600, 600);
    ViewField fieldView;
    if (plotField) {
        
        fieldView.SetComponent(&cmp);
        //fieldView.SetSensor(&sensor);
        // Set the normal vector of the viewing plane (xz plane).
        fieldView.SetPlane(0., 0., 1., 0., 0., 0);
        fieldView.SetArea(-pitch, 0, pitch, 0.02);        
        cfield->SetLeftMargin(0.16);
        fieldView.SetCanvas(cfield);
        fieldView.EnableAutoRange(false,false);
        fieldView.SetElectricFieldRange(1.e-3,12);
        //      fieldView.PlotContour();
    }
    
    if(debug) sensor.EnableDebugging();
    
    //AvalancheMicroscopic aval;
    aval.SetSensor(&sensor);
    aval.EnableSignalCalculation();
    aval.SetUserHandleIonisation(userHandleIonisation);
    aval.SetUserHandleAttachment(userHandleAttachment);
    if(debug) aval.EnableDebugging();
    
    //AvalancheMC drift;
    drift.SetSensor(&sensor);
    //drift.SetDistanceSteps(1e-4);
    drift.EnableSignalCalculation();
    drift.SetTimeSteps(timestep);
    if(debug) drift.EnableDebugging();

    //AvalancheGridSpaceCharge avalsc;
    avalsc.SetSensor(&sensor);
    if(debug) avalsc.EnableDebugging();
    
    if (enableSpaceCharge) avalsc.EnableSpaceChargeEffect();
    
    ViewDrift driftView;
    TCanvas* cd = new TCanvas();
    if (plotDrift) {
        aval.EnablePlotting(&driftView);
        drift.EnablePlotting(&driftView);
        driftView.SetCollisionMarkerSize(0.0001);
        driftView.SetColourIonisations(7);
        // track.EnablePlotting(&driftView);
        driftView.SetPlane(0, 0, 1, 0, 0, 0);
        driftView.SetArea(-pitch, 0, pitch, 0.02);
        driftView.SetCanvas(cd);
        constexpr bool twod = true;
    }


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
  cmp.Setup(int(layers.size()), eps, layers, voltage, {});
  
  cmp.SetAvalancheGridSpaceChargeObject(&avalsc);
  sensor.AddComponent(&cmp);
  
  

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
  int zsteps = 70;
  int rsteps = 70;

  aval.AddElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
  avalsc.Set2dGrid(zmin,zmax,zsteps,rmax,rsteps);
  avalsc.TransportOffGrid();

  
  // set number of steps
  const int N_timesteps = 50;

  bool electrons_remaining;
  bool ions_remaining;
  double xf,yf,zf,tf,xi,yi,zi,ti;
  int status;
  double tol = 1e-5;
  int ne,ni;
  int counter = 0;  

  /*
  Medium * medi = sensor.GetMedium(0.004,0.004,0.);
  int stat;
  double ex,ey,ez;*/
  
  

  // for each time step
  for (int step = 0; step < N_timesteps; ++step){
    if (step == N_timesteps - 1){
        for (const auto & ion:drift.GetIons()){
            ion_pos << ion.path.back().y << "\n";
        }
    }  

    std::cout << "step " << step << "\n";
    electrons_remaining = false;
    ions_remaining = false;
    avalsc.ClearGrid();

    particle_count = 0;
    mean_pos = {0.,0.,0.};
    add_electrons_to_grid(electrons_remaining, zmax);
    add_ions_to_grid(ions_remaining,zmax);
    add_negative_ions_to_grid(ions_remaining,zmax);
    if(particle_count > 0) avalsc.UpdateMeanPosition(mean_pos[0]/particle_count,mean_pos[1]/particle_count,mean_pos[2]/particle_count);

    avalsc.UpdateFieldOnGrid();

    if (electrons_remaining){
        std::cout << "Microscopic\n";
        aval.SetTimeWindow(tmin, tmin + timestep);
        aval.ResumeAvalanche();

        if (ions_remaining){
            std::cout << "MC\n";
            drift.SetTimeWindow(tmin, tmin + timestep);
            drift.ResumeAvalanche();
            
            
            if (plotDrift){
                driftView.Plot(true);
                cd->Update();
                char filename[50];
                snprintf(filename, 50, "frames/aval_%03d.png", counter);
                cd->SaveAs(filename);    
                dfile << filename << std::endl;
            }
            if (plotField){
                cfield->Clear();
                fieldView.Plot("E","zcol");
                char filename[50];
                snprintf(filename, 50, "frames/field_%03d.png", counter);
                cfield->SaveAs(filename);   
                ffile << filename << std::endl;
            }
            
            counter++; 
            tmin += timestep;
            continue;  
        }
        if (plotDrift){
            driftView.Plot(true);
            cd->Update();
            char filename[50];
            snprintf(filename, 50, "frames/aval_%03d.png", counter);
            cd->SaveAs(filename); 
            dfile << filename << std::endl;
        }
        if (plotField){
            cfield->Clear();
            fieldView.Plot("E","zcol");
            char filename[50];
            snprintf(filename, 50, "frames/field_%03d.png", counter);
            cfield->SaveAs(filename);  
            ffile << filename << std::endl;
        }
        tmin += timestep;
        counter++;
    }
    else{
        if (ions_remaining){
            std::cout << "MC\n";
            drift.SetTimeWindow(tmin, tmin + timestep);
            drift.ResumeAvalanche();  
            

            if (plotDrift){
                driftView.Plot(true);
                cd->Update();
                char filename[50];
                snprintf(filename, 50, "frames/aval_%03d.png", counter);
                cd->SaveAs(filename);  
                dfile << filename << std::endl; 
            }
            if (plotField){
                cfield->Clear();
                fieldView.Plot("E","zcol");
                char filename[50];
                snprintf(filename, 50, "frames/field_%03d.png", counter);
                cfield->SaveAs(filename); 
                ffile << filename << std::endl;  
            }
            tmin += timestep;
            counter++; 

        }
        else{
            break;
        }
    }
  }
}
