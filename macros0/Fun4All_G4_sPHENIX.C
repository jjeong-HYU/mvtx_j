#ifndef MACRO_FUN4ALLG4SPHENIX_C
#define MACRO_FUN4ALLG4SPHENIX_C

#include <GlobalVariables.C>

#include "./DisplayOn.C"
#include <G4Setup_sPHENIX.C>
#include <G4_Bbc.C>
#include <G4_CaloTrigger.C>
#include <G4_DSTReader.C>
#include <G4_Global.C>
#include <G4_HIJetReco.C>
#include <G4_Input.C>
#include <G4_Jets.C>
#include <G4_KFParticle.C>
#include <G4_ParticleFlow.C>
#include <G4_Production.C>
#include <G4_TopoClusterReco.C>
#include <G4_Tracking.C>
#include <G4_User.C>
#include <QA.C>

#include "FROG.h"

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <../writeHits/writeMVTXhits.cc>
//#include "/home/raise/Desktop/JungJae/Singularity/Workreal/MVTXalignment/install/physics/ReactionList.hh"

#include <g4centrality/PHG4CentralityReco.h>


#include <Geant4/G4MuonPlus.hh>
#include <Geant4/G4MuonMinus.hh>
#include <Geant4/G4RunManager.hh>
#include <Geant4/G4VUserPhysicsList.hh>
#include <Geant4/G4UImanager.hh>
R__ADD_LIBRARY_PATH("/cvmfs/sphenix.sdcc.bnl.gov/gcc-8.3/release/release_new/new.1/lib")
R__LOAD_LIBRARY(libfun4all.so)
R__ADD_LIBRARY_PATH("/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/release/release_new/new.8/lib")
R__LOAD_LIBRARY(libg4centrality.so.0.0.0)


#define PI 3.14159265

// For HepMC Hijing
// try inputFile = /sphenix/sim/sim01/sphnxpro/sHijing_HepMC/sHijing_0-12fm.dat

int Fun4All_G4_sPHENIX(
    const int nEvents = 1e0,
    const string &inputFile = "/sphenix/lustre01/sphnxpro/mdc2/shijing_hepmc/fm_0_20/g4hits/G4Hits_sHijing_0_20fm-0000000060-00000.root",
   // const string &inputFile = "/sphenix/lustre01/sphnxpro/mdc2/shijing_hepmc/fm_0_20/pileup/DST_BBC_G4HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000040-00000.root",
    const string &outputFile = "G4sPHENIX.root",
    const string &embed_input_file = "https://www.phenix.bnl.gov/WWW/publish/phnxbld/sPHENIX/files/sPHENIX_G4Hits_sHijing_9-11fm_00000_00010.root",
    const int skip = 0,
    const string &outdir = ".")
{
  bool runDisplay = false;

  FROG *fr = new FROG();

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  //Opt to print all random seed used for debugging reproducibility. Comment out to reduce stdout prints.
  PHRandomSeed::Verbosity(0);

  // just if we set some flags somewhere in this macro
  recoConsts *rc = recoConsts::instance();
  // By default every random number generator uses
  // PHRandomSeed() which reads /dev/urandom to get its seed
  // if the RANDOMSEED flag is set its value is taken as seed
  // You can either set this to a random value using PHRandomSeed()
  // which will make all seeds identical (not sure what the point of
  // this would be:
  //  rc->set_IntFlag("RANDOMSEED",PHRandomSeed());
  // or set it to a fixed value so you can debug your code
  //  rc->set_IntFlag("RANDOMSEED", 12345);
  
  //===============
  // Input options
  //===============
  // verbosity setting (applies to all input managers)
  Input::VERBOSITY = 0;
  // First enable the input generators
  // Either:
  // read previously generated g4-hits files, in this case it opens a DST and skips
  // the simulations step completely. The G4Setup macro is only loaded to get information
  // about the number of layers used for the cell reco code
  
  Input::READHITS = true;
  const string &inputFile0 = "/sphenix/lustre01/sphnxpro/mdc2/shijing_hepmc/fm_0_20/g4hits/G4Hits_sHijing_0_20fm-0000000060-00000.root";
  const string &inputFile1 = "/sphenix/lustre01/sphnxpro/mdc2/shijing_hepmc/fm_0_20/g4hits/G4Hits_sHijing_0_20fm-0000000060-00001.root";
  const string &inputFile2 = "/sphenix/lustre01/sphnxpro/mdc2/shijing_hepmc/fm_0_20/g4hits/G4Hits_sHijing_0_20fm-0000000060-00002.root";
  const string &inputFile3 = "/sphenix/lustre01/sphnxpro/mdc2/shijing_hepmc/fm_0_20/g4hits/G4Hits_sHijing_0_20fm-0000000060-00003.root";
  const string &inputFile4 = "/sphenix/lustre01/sphnxpro/mdc2/shijing_hepmc/fm_0_20/g4hits/G4Hits_sHijing_0_20fm-0000000060-00004.root";
  INPUTREADHITS::filename[0] = inputFile0;
  //INPUTREADHITS::filename[1] = inputFile1;
  //INPUTREADHITS::filename[2] = inputFile2;
  //INPUTREADHITS::filename[3] = inputFile3;
  //INPUTREADHITS::filename[4] = inputFile4;
  //INPUTREADHITS::filename[0] = inputFile2;
  // if you use a filelist
  // INPUTREADHITS::listfile[0] = inputFile;
  // Or:
  // Use particle generator
  // And
  // Further choose to embed newly simulated events to a previous simulation. Not compatible with `readhits = true`
  // In case embedding into a production output, please double check your G4Setup_sPHENIX.C and G4_*.C consistent with those in the production macro folder
  // E.g. /sphenix/sim//sim01/production/2016-07-21/single_particle/spacal2d/
  //  Input::EMBED = true;
  //INPUTEMBED::filename[0] = embed_input_file;
  // if you use a filelist
  //INPUTEMBED::listfile[0] = embed_input_file;
/*
   //Input::SIMPLE = true;
   Input::SIMPLE_NUMBER = 1; // if you need 2 of them
   Input::SIMPLE_VERBOSITY = 0;

   INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("mu-", 1);
   INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
                                                                             PHG4SimpleEventGenerator::Uniform,
                                                                             PHG4SimpleEventGenerator::Uniform);
   INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., 0.);
   INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0., 0., 0.);
   double eta_min=3.5*PI/180;
   double eta_max=3.5*PI/180;
   double phi_min=-90*PI/180;
   double phi_max=90*PI/180;
   double pt_min=20.0;
   double pt_max=20.0;
   INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(eta_min,eta_max);
   INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(phi_min,phi_max);
   INPUTGENERATOR::SimpleEventGenerator[0]->set_pt_range(pt_min, pt_max);
*/
  //  Input::PYTHIA6 = true;

  //Input::PYTHIA8 = true;
  //PYTHIA8::config_file = "/home/raise/Desktop/JungJae/Singularity_new/Singularity/cvmfs/sphenix.sdcc.bnl.gov/gcc-8.3/release/release_new/new.1/share/calibrations/Generators/HeavyFlavor_TG/phpythia8_d02kpi_MDC2.cfg";
  //  Input::GUN = true;
  //  Input::GUN_NUMBER = 3; // if you need 3 of them
  // Input::GUN_VERBOSITY = 1;

  //D0 generator
  //Input::DZERO = false;
  //Input::DZERO_VERBOSITY = 0;
  //Lambda_c generator //Not ready yet
  //Input::LAMBDAC = false;
  //Input::LAMBDAC_VERBOSITY = 0;
  // Upsilon generator
  //Input::UPSILON = true;
  //Input::UPSILON_NUMBER = 3; // if you need 3 of them
  //Input::UPSILON_VERBOSITY = 0;

//  Input::HEPMC = false;
//  INPUTHEPMC::filename = inputFile;

  // Event pile up simulation with collision rate in Hz MB collisions.
  //Input::PILEUPRATE = 100e3;
cout<< "hey hou"<< endl;
  //-----------------
  // Initialize the selected Input/Event generation
  //-----------------
  // This creates the input generator(s)
  InputInit();

  cout<< "say hey" << endl;

  //--------------
  // Set generator specific options
  //--------------
  // can only be set after InputInit() is called

  // Simple Input generator:
  // if you run more than one of these Input::SIMPLE_NUMBER > 1
  // add the settings for other with [1], next with [2]...

  ///////////////////////////////////////////////////////////////////////////////////////


/*/  //physics model fix_jjy
  ReactionList* fReactionList = ReactionList::Instance();
  fReactionList->SetReactionFlag("mu-","MultipleScattering",false);
  fReactionList->SetReactionFlag("mu-","CoulombScattering",false);
  fReactionList->SetReactionFlag("mu-","Bremsstrahlung",false);
  fReactionList->SetReactionFlag("mu-","PairProduction",true);
  fReactionList->SetReactionFlag("mu-","Ionisation",true);
  fReactionList->SetReactionFlag("All","eLossFluct",false);
  */
  //cout<<"Muon msc: "<<fReactionList->GetReactionFlag("mu-","MultipleScattering")<<endl;


  //--------------
  // Set Input Manager specific options
  //--------------
  // can only be set after InputInit() is called

  if (Input::HEPMC)
  {
    //! apply sPHENIX nominal beam parameter with 2mrad crossing as defined in sPH-TRG-2020-001
    Input::ApplysPHENIXBeamParameter(INPUTMANAGER::HepMCInputManager);

    // optional overriding beam parameters
    //INPUTMANAGER::HepMCInputManager->set_vertex_distribution_width(100e-4, 100e-4, 8, 0);  //optional collision smear in space, time
    //    INPUTMANAGER::HepMCInputManager->set_vertex_distribution_mean(0,0,0,0);//optional collision central position shift in space, time
    // //optional choice of vertex distribution function in space, time
    //INPUTMANAGER::HepMCInputManager->set_vertex_distribution_function(PHHepMCGenHelper::Gaus, PHHepMCGenHelper::Gaus, PHHepMCGenHelper::Gaus, PHHepMCGenHelper::Gaus);
    //! embedding ID for the event
    //! positive ID is the embedded event of interest, e.g. jetty event from pythia
    //! negative IDs are backgrounds, .e.g out of time pile up collisions
    //! Usually, ID = 0 means the primary Au+Au collision background
    //INPUTMANAGER::HepMCInputManager->set_embedding_id(Input::EmbedID);
    if (Input::PILEUPRATE > 0)
    {
      // Copy vertex settings from foreground hepmc input
      INPUTMANAGER::HepMCPileupInputManager->CopyHelperSettings(INPUTMANAGER::HepMCInputManager);
      // and then modify the ones you want to be different
      // INPUTMANAGER::HepMCPileupInputManager->set_vertex_distribution_width(100e-4,100e-4,8,0);
    }
  }
  if (Input::PILEUPRATE > 0)
  {
    //! apply sPHENIX nominal beam parameter with 2mrad crossing as defined in sPH-TRG-2020-001
    Input::ApplysPHENIXBeamParameter(INPUTMANAGER::HepMCPileupInputManager);
  }
  // register all input generators with Fun4All
  InputRegister();

  //======================
  // Write the DST
  //======================

  //Enable::DSTOUT = true;
  Enable::DSTOUT_COMPRESS = false;
  DstOut::OutputDir = outdir;
  DstOut::OutputFile = outputFile;

  // turn the display on (default off)
  Enable::DISPLAY = runDisplay ? true : false;

  //======================
  // What to run
  //======================

  // Global options (enabled for all enables subsystems - if implemented)
  //  Enable::ABSORBER = true;
  //  Enable::OVERLAPCHECK = true;
  //  Enable::VERBOSITY = 1;

  Enable::BBC = true;
  //Enable::BBCFAKE = true;  // Smeared vtx and t0, use if you don't want real BBC in simulation

  Enable::PIPE = runDisplay ? false : true;
  Enable::PIPE_ABSORBER = true;

  // central tracking
  Enable::MVTX = true;
  Enable::MVTX_CELL = Enable::MVTX && true;
  Enable::MVTX_CLUSTER = Enable::MVTX_CELL && true;
  Enable::MVTX_SERVICE = Enable::MVTX && runDisplay ? false : true;

  Enable::INTT = runDisplay ? false : true;
  Enable::INTT_CELL = Enable::INTT && true;
  Enable::INTT_CLUSTER = Enable::INTT_CELL && true;

  Enable::TPC = runDisplay ? false : true;
  Enable::TPC_ABSORBER = false;
  Enable::TPC_CELL = Enable::TPC && true;
  Enable::TPC_CLUSTER = Enable::TPC_CELL && true;

  Enable::MICROMEGAS=true;

  Enable::TRACKING_TRACK = runDisplay ? false : true;

  Enable::CEMC = false;
  Enable::HCALIN = false;
  Enable::MAGNET = false;
  Enable::MAGNET_ABSORBER = true;
  Enable::HCALOUT = false;
  Enable::EPD = false;
  Enable::PLUGDOOR_ABSORBER = false;

  Enable::GLOBAL_RECO = false;
  //Enable::GLOBAL_FASTSIM = true;
  Enable::JETS = false;

  // new settings using Enable namespace in GlobalVariables.C
  Enable::BLACKHOLE = true;

  //Enable::USER = true;

  //---------------
  // Magnet Settings
  //---------------

   // const string magfield = "1.5"; // alternatively to specify a constant magnetic field, give a float number, which will be translated to solenoidal field in T, if string use as fieldmap name (including path)
   // G4MAGNET::magfield = string(getenv("CALIBRATIONROOT")) + string("/Field/Map/sPHENIX.2d.root");  // default map from the calibration database
   // G4MAGNET::magfield_rescale = 1.5;//-1.4 / 1.5;  // make consistent with expected Babar field strength of 1.4T

  // Initialize the selected subsystems
  G4Init();

  //---------------------
  // GEANT4 Detector description
  //---------------------
  if (!Input::READHITS)
  {
    G4Setup();
  }
/*/
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  //UImanager->ApplyCommand("/run/setCut 1 km");
  UImanager->ApplyCommand("/process/eLoss/fluct false");
  UImanager->ApplyCommand("/process/inactivate msc");
//  UImanager->ApplyCommand("/process/inactivate hadElastic");
//  //UImanager->ApplyCommand("/process/msc/LateralDisplacement false");
//  //UImanager->ApplyCommand("/process/msc/MuHadLateralDisplacement false");
//  //UImanager->ApplyCommand("/process/msc/LateralDisplacement false");
//  //UImanager->ApplyCommand("/process/msc/Pixe false");
  UImanager->ApplyCommand("/process/inactivate CoulombScat");
  UImanager->ApplyCommand("/process/inactivate muBrems");
  UImanager->ApplyCommand("/process/inactivate eBrems");
  UImanager->ApplyCommand("/process/inactivate muPairProd");
//  UImanager->ApplyCommand("/process/verbose 5");
//  UImanager->ApplyCommand("/process/eLoss/verbose 5");

*/
  //jjy
  PHG4CentralityReco *cent = new PHG4CentralityReco();
  cent->Verbosity(0);
  cent->GetCalibrationParameters().ReadFromFile("centrality", "xml", 0, 0, string(getenv("CALIBRATIONROOT")) + string("/Centrality/"));
  se->registerSubsystem( cent );

  //------------------
  // Detector Division
  //------------------

  if (Enable::BBC || Enable::BBCFAKE) Bbc_Reco();

  if (Enable::MVTX_CELL) Mvtx_Cells();
  if (Enable::INTT_CELL) Intt_Cells();
  if (Enable::TPC_CELL) TPC_Cells();
  if (Enable::MICROMEGAS_CELL) Micromegas_Cells();

  //--------------
  // SVTX tracking
  //--------------
  if(Enable::TRACKING_TRACK)
  {
    TrackingInit();
  }
  if (Enable::MVTX_CLUSTER) Mvtx_Clustering();
  if (Enable::INTT_CLUSTER) Intt_Clustering();
  if (Enable::TPC_CLUSTER) TPC_Clustering();
  if (Enable::MICROMEGAS_CLUSTER) Micromegas_Clustering();

  if (Enable::TRACKING_TRACK)
  {
    Tracking_Reco();
  }
  //-----------------
  // Global Vertexing
  //-----------------

  if (Enable::GLOBAL_RECO && Enable::GLOBAL_FASTSIM)
  {
    cout << "You can only enable Enable::GLOBAL_RECO or Enable::GLOBAL_FASTSIM, not both" << endl;
    gSystem->Exit(1);
  }
  if (Enable::GLOBAL_RECO)
  {
    Global_Reco();
  }
  else if (Enable::GLOBAL_FASTSIM)
  {
    Global_FastSim();
  }

   writeMVTXhits *myMVTXhits = new writeMVTXhits("myMVTXhits");
   myMVTXhits->Verbosity(0);
   se->registerSubsystem(myMVTXhits);

  //--------------
  // Set up Input Managers
  //--------------

  InputManagers();

  if (Enable::PRODUCTION)
  {
    Production_CreateOutputDir();
  }

  if (Enable::DSTOUT)
  {
    string FullOutFile = DstOut::OutputDir + "/" + DstOut::OutputFile;
    Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", FullOutFile);
    if (Enable::DSTOUT_COMPRESS)
    {
      ShowerCompress();
      DstCompress(out);
    }
    se->registerOutputManager(out);
  }
  //-----------------
  // Event processing
  //-----------------
  if (Enable::DISPLAY)
  {
    DisplayOn();

    gROOT->ProcessLine("Fun4AllServer *se = Fun4AllServer::instance();");
    gROOT->ProcessLine("PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");");

    cout << "-------------------------------------------------" << endl;
    cout << "You are in event display mode. Run one event with" << endl;
    cout << "se->run(1)" << endl;
    cout << "Run Geant4 command with following examples" << endl;
    gROOT->ProcessLine("displaycmd()");

    return 0;
  }

  // if we use a negative number of events we go back to the command line here
  if (nEvents < 0)
  {
    return 0;
  }
  // if we run the particle generator and use 0 it'll run forever
  if (nEvents == 0 && !Input::HEPMC && !Input::READHITS)
  {
    cout << "using 0 for number of events is a bad idea when using particle generators" << endl;
    cout << "it will run forever, so I just return without running anything" << endl;
    return 0;
  }

  se->skip(skip);
  se->run(nEvents);


/*
  for(Double_t angl=-178;angl<=180;angl+=2)
  {
    for(double h=-0.9;h<=1.0;h+=0.1)
    {
      if (Input::SIMPLE)
      {
        INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("mu-", 1);
          INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
                                                                                    PHG4SimpleEventGenerator::Uniform,
                                                                                    PHG4SimpleEventGenerator::Uniform);
          INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., h);
          INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0., 0., 0.);
        double eta_min=-0.002;
        double eta_max=-0.002;
        double phi_min=angl*PI/180;
        double phi_max=angl*PI/180;
        double pt_min=2.0;
        double pt_max=2.0;
        INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(eta_min,eta_max);
        INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(phi_min,phi_max);
        INPUTGENERATOR::SimpleEventGenerator[0]->set_pt_range(pt_min, pt_max);
      }
      se->run(nEvents);
    }
  }

for(Double_t angl=-178;angl<=180;angl+=2)
  {
    for(double h=-0.9;h<=1.0;h+=0.1)
    {
      if (Input::SIMPLE)
      {
        INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("mu-", 1);
          INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
                                                                                    PHG4SimpleEventGenerator::Uniform,
                                                                                    PHG4SimpleEventGenerator::Uniform);
          INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., h);
          INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0., 0., 0.);
        double eta_min=-0.001;
        double eta_max=-0.001;
        double phi_min=angl*PI/180;
        double phi_max=angl*PI/180;
        double pt_min=2.0;
        double pt_max=2.0;
        INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(eta_min,eta_max);
        INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(phi_min,phi_max);
        INPUTGENERATOR::SimpleEventGenerator[0]->set_pt_range(pt_min, pt_max);
      }
      se->run(nEvents);
    }
  }

    for(Double_t angl=-178;angl<=180;angl+=2)
  {
    for(double h=-0.9;h<=1.0;h+=0.1)
    {
      if (Input::SIMPLE)
      {
        INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("mu-", 1);
          INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
                                                                                    PHG4SimpleEventGenerator::Uniform,
                                                                                    PHG4SimpleEventGenerator::Uniform);
          INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., h);
          INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0., 0., 0.);
        double eta_min=0.001;
        double eta_max=0.001;
        double phi_min=angl*PI/180;
        double phi_max=angl*PI/180;
        double pt_min=2.0;
        double pt_max=2.0;
        INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(eta_min,eta_max);
        INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(phi_min,phi_max);
        INPUTGENERATOR::SimpleEventGenerator[0]->set_pt_range(pt_min, pt_max);
      }
      se->run(nEvents);
    }
  }

  for(Double_t angl=-178;angl<=180;angl+=2)
  {
    for(double h=-0.9;h<=1.0;h+=0.1)
    {
      if (Input::SIMPLE)
      {
        INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("mu-", 1);
          INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
                                                                                    PHG4SimpleEventGenerator::Uniform,
                                                                                    PHG4SimpleEventGenerator::Uniform);
          INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., h);
          INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0., 0., 0.);
        double eta_min=0.002;
        double eta_max=0.002;
        double phi_min=angl*PI/180;
        double phi_max=angl*PI/180;
        double pt_min=2.0;
        double pt_max=2.0;
        INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(eta_min,eta_max);
        INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(phi_min,phi_max);
        INPUTGENERATOR::SimpleEventGenerator[0]->set_pt_range(pt_min, pt_max);
      }
      se->run(nEvents);
    }
  }*/

  //-----
  // Exit
  //-----

  se->End();
  std::cout << "All done" << std::endl;
  delete se;

  gSystem->Exit(0);
  return 0;
}
#endif
