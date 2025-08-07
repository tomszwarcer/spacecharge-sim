#include <TTree.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TFile.h>
#include <iostream>
#include <TGraph.h>
#include <TLatex.h>
#include <fstream>
#include <numeric>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentGrid.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Random.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ComponentComsol.hh"
#include "Garfield/ComponentUser.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewMedium.hh"
#include "Garfield/Random.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ComponentChargedRing.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/FundamentalConstants.hh"



#define LOG(x) std::cout<<x<<std::endl

using namespace Garfield;

// Create ions
AvalancheMC drift;

AvalancheMicroscopic aval;
// set time window parameters
double tmin = 0.;
double timestep = 0.025;

std::vector<double> mean_pos = {0.,0.,0.};
int particle_count = 0;



void userHandleIonisation(double x, double y, double z, double t, int type, int level, Medium * m){
    drift.AddIon(x,y,z,tmin + timestep); // ion added after step.
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
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <process_index>" << std::endl;
    return 1;
  }

  std::ofstream outfile("outfile.txt");
  
  bool plotting = false;
  // Convert the process index argument to an integer
  int processIndex = std::atoi(argv[1]);
  
  // Example logic: Print a message including the process index
  LOG("Script: Running job with process index: " << processIndex);
  
  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();
  
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
  const double voltage = -700; // V
  
  ComponentAnalyticField mm; // Micromegas field
  mm.SetMedium(&gas);
  mm.AddPlaneY(posBottomPlane, 0.);
  mm.AddPlaneY(posTopPlane, voltage);

  aval.SetUserHandleIonisation(userHandleIonisation);
  aval.SetUserHandleAttachment(userHandleAttachment);
  
  ComponentChargedRing rings;
  rings.SetArea(-0.02,posBottomPlane,-0.02,0.02,posTopPlane,0.02);
  rings.SetSpacingTolerance(0.00001);
  rings.SetMedium(&gas);
  
  // Create the sensor
  Sensor sensor;
  sensor.AddComponent(&mm);
  sensor.AddComponent(&rings);
  
  // Create electrons
  aval.SetSensor(&sensor);
  
  drift.SetSensor(&sensor);

  const bool enableDebug = false;
  const bool enableSpaceCharge = true;
  const bool plotField = (true && enableSpaceCharge);
  const bool plotDrift = true;

  if (enableDebug) { 
    rings.EnableDebugging();
    aval.EnableDebugging();
    drift.EnableDebugging();
  }
  
  // Plot electric field
  TCanvas* cfield = new TCanvas("cfield", "", 600, 600);  
  ViewField fieldView;
  if (plotField) {
    fieldView.SetComponent(&rings);
    //fieldView.SetSensor(&sensor);
    // Set the normal vector of the viewing plane (xz plane).
    fieldView.SetPlane(0., 0., 1., 0., 0., 0);
    //fieldView.SetArea(-0.001, 0.006, 0.0012, 0.0072);      
    fieldView.SetArea(-0.01, posBottomPlane, 0.01, posTopPlane);      
    cfield->SetLeftMargin(0.16);
    fieldView.SetCanvas(cfield);
    fieldView.EnableAutoRange(false,false);
    fieldView.SetElectricFieldRange(0.,75.);
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
      driftView.SetArea(-0.01, posBottomPlane, 0.01, posTopPlane);
      driftView.SetCanvas(cd);
      constexpr bool twod = true;
  }
    
  const double elCharge = 1.60217663e-19; // C
  
  double midpt = (posTopPlane - posBottomPlane)/2;
  const double x0 = 0.; // cm
  double y0 = midpt*1.2;
  const double z0 = 0.; // cm
  const double t0 = 0.; // ns
  const double e0 = 0.1; // eV

  const unsigned int nFrames = 40; 

  aval.AddElectron(x0,y0,z0,0.,0.); 

  std::array<int,3> particles;


  for (unsigned int iF = 0; iF < nFrames; ++iF) {

    std::cout << "Frame " << iF << "\n";
    
    rings.ClearActiveRings();

    double meanZ = 0;
    double meanX = 0;

    particle_count = 0;
    mean_pos = {0.,0.,0.};

    bool electrons_remaining = false;
    bool ions_remaining = false;

    particles = get_mean(electrons_remaining,ions_remaining,posTopPlane);
    std::cout << "e: " << particles[0] << ", +i: " << particles[1] << ", -i: " << particles[2] << "\n";

    meanX = mean_pos[0]/particle_count;
    meanZ = mean_pos[2]/particle_count;
    if(particle_count > 0) rings.UpdateCentre(meanX, meanZ);
    
    LOG("   meanZ = "<< meanZ <<" cm.");
    LOG("   meanX = "<< meanX <<" cm.");
        
    double x,y,z;

    if (enableSpaceCharge){

        for (auto & electron:aval.GetElectrons()){
        double xf = electron.path.back().x;
        double yf = electron.path.back().y;
        double zf = electron.path.back().z;
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
            
            
            if (plotDrift){
                driftView.Plot(true);
                cd->Update();
                char filename[100];
                snprintf(filename, 100, "/afs/cern.ch/work/t/tszwarce/nogrid/build/frames/aval_%03d.png", iF);
                cd->SaveAs(filename);    
            }
            if (plotField){
                cfield->Clear();
                fieldView.Plot("e","zcol");
                char filename[100];
                snprintf(filename, 100, "/afs/cern.ch/work/t/tszwarce/nogrid/build/frames/field_%03d.png", iF);
                cfield->SaveAs(filename);   
            }
            
            tmin += timestep;
            continue;  
        }
        if (plotDrift){
            driftView.Plot(true);
            cd->Update();
            char filename[100];
            snprintf(filename, 100, "/afs/cern.ch/work/t/tszwarce/nogrid/build/frames/aval_%03d.png", iF);
            cd->SaveAs(filename); 
        }
        if (plotField){
            cfield->Clear();
            fieldView.Plot("e","zcol");
            char filename[100];
            snprintf(filename, 100, "/afs/cern.ch/work/t/tszwarce/nogrid/build/frames/field_%03d.png", iF);
            cfield->SaveAs(filename);  
        }
        tmin += timestep;
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
                snprintf(filename, 100, "/afs/cern.ch/work/t/tszwarce/nogrid/build/frames/aval_%03d.png", iF);
                cd->SaveAs(filename);  
            }
            if (plotField){
                cfield->Clear();
                fieldView.Plot("e","zcol");
                char filename[100];
                snprintf(filename, 100, "/afs/cern.ch/work/t/tszwarce/nogrid/build/frames/field_%03d.png", iF);
                cfield->SaveAs(filename); 
            }
            tmin += timestep;

        }
        else{
            break;
        }
    }


  }
  outfile.close();
  
  LOG("Script: Done.\n");
}
