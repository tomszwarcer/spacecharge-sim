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
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ComponentComsol.hh"
#include "Garfield/ComponentParallelPlate.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/Random.hh"
#include "Garfield/ComponentAnalyticField.hh"


#define LOG(x) std::cout<<x<<std::endl

using namespace Garfield;

AvalancheMC drift;
AvalancheMicroscopic aval;
AvalancheGridSpaceCharge avalsc;
MediumMagboltz gas;

// set time window parameters
double tmin = 0.;
double timestep = 0.05;

std::vector<double> mean_pos = {0.,0.,0.};
int particle_count = 0;

void userHandleIonisation(double x, double y, double z, double t, int type, int level, Medium * m){
    drift.AddIon(x,y,z,tmin + timestep); // ion added after step.
}

void userHandleAttachment(double x, double y, double z, double t, int type, int level, Medium * m){
    drift.AddNegativeIon(x,y,z,tmin + timestep);
}

int add_electrons_to_grid(bool & electrons_remaining,const double & zmax){
    int ne = 0;
    for (const auto& electron : aval.GetElectrons()) {
        // get electron positions
        double xf = electron.path.back().x;
        double yf = electron.path.back().y;
        double zf = electron.path.back().z;
        double tf = electron.path.back().t;
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
                avalsc.AddElectron(xf,yf,zf,tf,1);
            }
        }
    }
    return ne;
}

int add_ions_to_grid(bool & ions_remaining,const double &zmax){
    int ni= 0;
    for (const auto& ion : drift.GetIons()){
        double xf = ion.path.back().x;
        double yf = ion.path.back().y;
        double zf = ion.path.back().z;
        double tf = ion.path.back().t;
        double ti = ion.path.front().t;
        int status = ion.status;

        if (status==0 || status == -17){
            ni++;
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
    return ni;
}

int add_negative_ions_to_grid(bool & ions_remaining,const double &zmax){
    int ni=0;
    for (const auto& ion : drift.GetNegativeIons()){
        double xf = ion.path.back().x;
        double yf = ion.path.back().y;
        double zf = ion.path.back().z;
        double tf = ion.path.back().t;
        int status = ion.status;

        if (status==0 || status == -17){
	    ni++;
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
    return ni;
}

int main(int argc, char * argv[]){

    constexpr bool plotDrift = false;
    constexpr bool plotField = false;

    TApplication app("app", &argc, argv);


    int run_number = std::atoi(argv[1]);
    std::string s_run_number = argv[1];
    int i_dv = std::atoi(argv[2]);
    double dv = std::atof(argv[2]);
    std::string s_dv = std::to_string(i_dv);
    dv = -1*dv;
    s_dv = "_"+s_dv;
    std::ofstream size;
    constexpr bool debug = false;
    constexpr bool plotSignal = false;
    constexpr bool enableSpaceCharge = false;
    constexpr bool enable_penning = true;

    
    // Setup the gas.
    
    gas.SetComposition("ar", 93, "co2", 7); // [%]
    gas.SetTemperature(293.15); // [K]
    gas.SetPressure(760); // [Torr]
    gas.Initialise(true);
    

    if (enable_penning) gas.EnablePenningTransfer();
    
    // Load the ion mobilities.
    gas.LoadIonMobility("/afs/cern.ch/work/t/tszwarce/garfieldpp/install/share/Garfield/Data/IonMobility_Ar+_Ar.txt");

    constexpr double pitch = 200.e-4; //cm
    constexpr double halfpitch = 0.5 * pitch;

    ComponentParallelPlate sc;
    std::string path_end = ".txt";
    if (enableSpaceCharge){
        sc.EnableSpaceCharge();
        std::string path_begin = "/afs/cern.ch/work/t/tszwarce/ppp/gain_output/on/size_sc_";
        
        size.open(path_begin+s_run_number+s_dv+path_end);
    }
    else{
        std::string path_begin = "/afs/cern.ch/work/t/tszwarce/ppp/gain_output/off/size_no_sc_";
        size.open(path_begin+s_run_number+s_dv+path_end);
    }   
    sc.SetMedium(&gas);

    

    // Create the sensor.
    Sensor sensor;


    if(debug) sensor.EnableDebugging();
    
    //AvalancheMicroscopic aval;
    aval.SetSensor(&sensor);
    aval.SetUserHandleIonisation(userHandleIonisation);
    aval.SetUserHandleAttachment(userHandleAttachment);
    if(debug) aval.EnableDebugging();
    
    //AvalancheMC drift;
    drift.SetSensor(&sensor);
    drift.EnableSignalCalculation();
    drift.SetTimeSteps(timestep);
    if(debug) drift.EnableDebugging();

    //AvalancheGridSpaceCharge avalsc;
    avalsc.SetSensor(&sensor);
    if(debug) avalsc.EnableDebugging();
    
    if (enableSpaceCharge) avalsc.EnableSpaceChargeEffect();
    
  // We create a fake RPC which is just a gas gap with dV = 0
  // This fake RPC lives under the mesh of the micromegas.
  const double mesh_pos = 0.0128;
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
  sc.Setup(int(layers.size()), eps, layers, voltage, {});
  
  sc.SetAvalancheGridSpaceChargeObject(&avalsc);
  sensor.AddComponent(&sc);

  ComponentAnalyticField pp;
  pp.SetMedium(&gas);
  pp.AddPlaneY(anode_pos, 0.);
  pp.AddPlaneY(mesh_pos, dv);
  sensor.AddComponent(&pp);


  // Starting parameters for the electron
  const double x0 = 0.0;
  const double y0 = 0.012; // mesh is at y=0.012 more or less
  const double z0 = 0.;
  const double t0 = 0.;
  const double e0 = 0.1;

  // grid parameters
  const double zmax = d_gas;
  const double zmin = anode_pos; // < anode pos
  const double rmax = std::max(std::sqrt(x0*x0 + z0*z0)*1.25,y0*1.25); // < with a buffer region
  int zsteps = 40;
  int rsteps = 40;

  aval.AddElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
  avalsc.Set2dGrid(zmin,zmax,zsteps,rmax,rsteps);
  avalsc.TransportOffGrid();

  
  // set number of steps
  const int N_timesteps = 40;

  bool electrons_remaining;
  bool ions_remaining;
  int ne,ni;

    /* TCanvas* cmsh = new TCanvas("cfield", "", 600, 600);
    ViewFEMesh meshView;
    meshView.SetPlane(0., 0., 1., 0., 0., 0);
    meshView.SetArea(-pitch, 0, pitch, 0.02);  
    meshView.SetComponent(&fm);
    meshView.SetCanvas(cmsh);
    meshView.SetFillMesh(true);
    meshView.EnableAxes();
    meshView.Plot();
    cmsh->SaveAs("test.png"); */

    TCanvas* cfield = new TCanvas("cfield", "", 600, 600);
    ViewField fieldView;
    if (plotField) {
        fieldView.SetComponent(&pp);
        //fieldView.SetSensor(&sensor);
        // Set the normal vector of the viewing plane (xz plane).
        fieldView.SetPlane(0., 0., 1., 0., 0., 0);
        fieldView.SetArea(-pitch, 0, pitch, 0.02);        
        cfield->SetLeftMargin(0.16);
        fieldView.SetCanvas(cfield);
        fieldView.EnableAutoRange(false,false);
        fieldView.SetElectricFieldRange(1.e-3,1.e2);
        //      fieldView.PlotContour();
    }

    ViewDrift driftView;
    TCanvas* cd = new TCanvas();
    if (plotDrift) {
        aval.EnablePlotting(&driftView);
        drift.EnablePlotting(&driftView);
        driftView.SetCollisionMarkerSize(0.0000001);
        driftView.SetColourIonisations(7);
        // track.EnablePlotting(&driftView);
        driftView.SetPlane(0, 0, 1, 0, 0, 0);
        driftView.SetArea(-pitch, 0, pitch, 0.02);
        driftView.SetCanvas(cd);
        constexpr bool twod = true;
    }

  // for each time step
  int counter = 0;
  for (int step = 0; step < N_timesteps; ++step){

    std::cout << "step " << step << "\n";
    electrons_remaining = false;
    ions_remaining = false;
    avalsc.ClearGrid();

    particle_count = 0;
    mean_pos = {0.,0.,0.};
    ne = add_electrons_to_grid(electrons_remaining, zmax);
    ni = add_ions_to_grid(ions_remaining,zmax);
    ni = ni + add_negative_ions_to_grid(ions_remaining,zmax);
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
                char filename[100];
                snprintf(filename, 100, "/afs/cern.ch/work/t/tszwarce/ppp/build/frames/aval_%03d.png", counter);
                cd->SaveAs(filename);    
            }
            if (plotField){
                cfield->Clear();
                fieldView.Plot("e","zcol");
                char filename[100];
                snprintf(filename, 100, "/afs/cern.ch/work/t/tszwarce/ppp/build/frames/field_%03d.png", counter);
                cfield->SaveAs(filename);   
            }
            
            counter++; 
            tmin += timestep;
            continue;  
        }
        if (plotDrift){
            driftView.Plot(true);
            cd->Update();
            char filename[100];
            snprintf(filename, 100, "/afs/cern.ch/work/t/tszwarce/ppp/build/frames/aval_%03d.png", counter);
            cd->SaveAs(filename); 
        }
        if (plotField){
            cfield->Clear();
            fieldView.Plot("e","zcol");
            char filename[100];
            snprintf(filename, 100, "/afs/cern.ch/work/t/tszwarce/ppp/build/frames/field_%03d.png", counter);
            cfield->SaveAs(filename);  
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
                char filename[100];
                snprintf(filename, 100, "/afs/cern.ch/work/t/tszwarce/ppp/build/frames/aval_%03d.png", counter);
                cd->SaveAs(filename);  
            }
            if (plotField){
                cfield->Clear();
                fieldView.Plot("e","zcol");
                char filename[100];
                snprintf(filename, 100, "/afs/cern.ch/work/t/tszwarce/ppp/build/frames/field_%03d.png", counter);
                cfield->SaveAs(filename); 
            }
            tmin += timestep;
            counter++; 

        }
        else{
            break;
        }
    }
  }
  size << i_dv << "," <<ni << "\n";
}
