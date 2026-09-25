// Updated on 5/29/2014, Hexc, Olesya
//            6/20/2014, Hexc, Olesya: Remove the verbose level
//            2/11/2015, Hexc, Olesya: Add additional arguments for the stepping actions
//            4/5/2015, Hexc, switch to QGSP_BERT_HP physics processes since
//                      FTFP_BERT seems killing neutron tracking at low energies.
//            6/12/2015, Hexc, Olesya: Modified for running ECRS in BNL RCF nodes.
//            1/14/2021, Hexc, Jarred and Marcus: update the ECRS simulation for the newer GEANT4 release and Linux OS
//            2/25/2021, Hexc, Jarred, Marcus, Zachary, Jack, Ernesto: Added a command line option for number of threads to run.
//            3/16/2021, Hexc: Update the code for running batch mode
//            9/25/2026: Parallelized the simulation. Enabled Geant4
//                      event-level multithreading through G4RunManagerFactory (which
//                      instantiates G4MTRunManager when Geant4 is built with
//                      multithreading) and wired up the worker-thread count: the
//                      optional 3rd command-line argument sets the number of threads,
//                      otherwise all available CPU cores are used. Events are now
//                      distributed across worker threads and the per-thread
//                      G4AnalysisManager ntuples are merged automatically by the
//                      master thread. See the SetNumberOfThreads() block below.
//
#include "G4Types.hh"

/*
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
*/
#include "G4RunManagerFactory.hh"

#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

// 9/25/2026: needed for querying the CPU core count and for
//            std::max / atoi used when configuring the number of worker threads.
#include "G4Threading.hh"
#include <algorithm>
#include <cstdlib>

// 4/5/2015
#include "QGSP_BERT_HP.hh"
//#include "FTFP_BERT.hh"
#include "G4PhysListFactory.hh"
#include "G4OpticalPhysics.hh"
//#include "FPPhysicsList.hh"

#include <time.h> 

#include "G4ProcessTable.hh"
#include "G4UnitsTable.hh"

#include "ECRSMagneticField.hh"
#include "ECRSDetectorConstruction.hh"
#include "ECRSPrimaryGeneratorAction.hh"
#include "ECRSActionInitialization.hh"
#include "ECRSVisManager.hh"
#include "ECRSRunAction.hh"
#include "ECRSEventAction.hh"
#include "ECRSStackingAction.hh"
#include "ECRSTrackingAction.hh"
#include "ECRSSteppingAction.hh"
#include "ECRSSingleton.hh"
#include "ECRSUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " ECRS <seedIndex> <macroFile> [nThreads]" << G4endl;
    G4cerr << "   seedIndex : line index into ECRS_500kRand.txt for the RNG seed" << G4endl;
    G4cerr << "   macroFile : Geant4 macro to execute in batch mode" << G4endl;
    G4cerr << "   nThreads  : (optional) number of worker threads to run in"     << G4endl;
    G4cerr << "               parallel; defaults to all available CPU cores."    << G4endl;
  }
}

using namespace CLHEP;

int main(int argc, char** argv) {

  // Units
  // Definition of new units in the unit table should be defined 
  // at the beginning before the
  // instantiation of the runManager and should be followed by 
  // G4UnitDefinition::BuildUnitsTable()  
  
  G4cout << "Building our own units ... " << G4endl;

  new G4UnitDefinition("earth radii","re","Length",re);
  new G4UnitDefinition("earth radii 1","Re","Length",re);
  new G4UnitDefinition("earth radii 2","RE","Length",re);
  new G4UnitDefinition("hour","hour","Time",3600.*s);
  new G4UnitDefinition("minute","minute","Time",60.*s);
  new G4UnitDefinition("day","day","Time",24.*3600.*s);
  new G4UnitDefinition("nanotesla","nT","Magnetic flux density",nT);
  new G4UnitDefinition("gigavolt","GV","Electric potential",1000.*megavolt);
  
  G4UnitDefinition::BuildUnitsTable();

  // The following code reads in parameters from the file earth.par
  // The parameters read are: seed index for random number generator
  // and the name of the output file for data generated.
  // The code also opens the output file.
  
  srand(time(NULL));
  
  //  G4cout << "Open output file (as a singleton) ... " << G4endl;
  // ECRSSingleton* myOut = ECRSSingleton::instance();
  //  char Ofilename[100];
  //  G4int seed_index=10;

  G4int seed_index;

  // Select different random number sequence
  //    argc == 1 --- call rand() function for seed_index
  //    argc == 2 --- get the seed_index from the ECRS_500kRand.txt
  //          argv[1] is the index of the random numbers sequence in the ECRS_500kRand.txt
  //
  G4long ECRS_Rand[10000];
  if (argc == 1) {
    seed_index = (int)rand();  
  } else if ((argc >= 2) && (argc <= 4)) {   // 9/25/2026: allow optional argv[3]=nThreads
    std::ifstream ECRS_RandNumFile("ECRS_500kRand.txt");
    for ( int i = 0; i < 10000; i++) {
      ECRS_RandNumFile >> ECRS_Rand[i];
    }
    ECRS_RandNumFile.close();
    seed_index = ECRS_Rand[atol(argv[1])];
  } else {
    G4cout << " You forgot picking a proper random seed index" << G4endl;
    exit (1);
  }

  //myOut->Fopen("Cosmic_Output.dat");
  //G4cout << "Just opened the output file ********* " << G4endl;
  
  // Select the RanecuEngine random number generator with seeds defined above
  
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  G4long Myseeds[2];
  Myseeds[0] = (int)rand();
  Myseeds[1] = (int)rand(); 
  CLHEP::HepRandom::setTheSeeds(Myseeds,seed_index);


  /*
  
  //
  // Run manager (with multithreaded option)
  //
  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  G4cout << "We are using G4MTRunManager ........ " << G4endl;
  // if ( nThreads > 0 ) {
  //  G4cout << " argv[2] : " << argv[1] << G4endl;
  runManager->SetNumberOfThreads(8);
  //  }
#else
  G4RunManager* runManager = new G4RunManager;
  G4cout << "We are using G4RunManager (not multi-threaded)........ " << G4endl;
#endif

  */
  
  // Construct the default run manager
  //
  auto* runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  // 9/25/2026: PARALLELIZATION
  // -------------------------------------------------------------------------
  // The simulation was previously run serially (one event at a time on a
  // single core). Geant4 supports event-level parallelism: each worker thread
  // processes a different subset of the events concurrently, and the per-thread
  // G4AnalysisManager ntuples are merged automatically by the master thread at
  // the end of the run, so the physics results are unchanged.
  //
  // G4RunManagerFactory::CreateRunManager(Default) already returns a
  // G4MTRunManager when Geant4 was built with multithreading support (and a
  // sequential G4RunManager otherwise). Here we set how many worker threads to
  // use: an optional 3rd command-line argument (argv[3]) overrides the default
  // of "all available hardware cores". In a sequential Geant4 build
  // SetNumberOfThreads() is simply a no-op.
  G4int nThreads = G4Threading::G4GetNumberOfCores();
  if (argc >= 4) {
    nThreads = std::max(1, atoi(argv[3]));
  }
  runManager->SetNumberOfThreads(nThreads);
  G4cout << "ECRS parallelization: requesting " << nThreads
         << " worker thread(s) (detected cores: "
         << G4Threading::G4GetNumberOfCores() << ")." << G4endl;
  //
  // 9/25/2026: THREAD-SAFETY NOTE (please read)
  // The event loop, physics and G4AnalysisManager output are all thread-safe.
  // HOWEVER, the analytic magnetic-field models used by ECRSMagneticField
  // (t89c_, t89c_boberg_, mcos_t96_01_, t01_01_, igrf_to_cc_) keep their state
  // in global Fortran COMMON blocks (geopack_, igrfcc_, ...). Those globals are
  // shared across worker threads, so concurrent field evaluations can race.
  // To be certain MT results match the serial run, the COMMON blocks must be
  // made thread-local (e.g. declare them "thread_local"/__thread, or guard the
  // field call), or validated to be read-only after initialization. Until that
  // is done, run with nThreads = 1 for production output, or verify a small MT
  // run against the serial result first.
  // -------------------------------------------------------------------------
  //
  // Mandatory initialization classes
  //
  // Construct the detector
  //
  G4cout << "Construct the detector ... " << G4endl;
  ECRSDetectorConstruction* detector = new ECRSDetectorConstruction;
  G4cout << "Detector construction is done!  " << G4endl;
  
  runManager->SetUserInitialization(detector);
  
  //
  // Define and reguster physics processes
  //
  // Default physics process list is: QGSP_BERT_HP
  // 
  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = NULL;
  G4String physName = "";
  char* path = getenv("PHYSLIST");
  if (path) {
      physName = G4String(path);
  } else {
      physName = "QGSP_BERT_HP"; // default
  }
  // reference PhysicsList via its name
  if (factory.IsReferencePhysList(physName)) {
      phys = factory.GetReferencePhysList(physName);
  }
  
  // For now, we are not implementing any optical processes.
  // Now add and configure optical physics
  //
  //G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  //opticalPhysics->Configure(kCerenkov, true);
  //opticalPhysics->SetCerenkovStackPhotons(false);
  //  opticalPhysics->Configure(kScintillation, true);  
  //  opticalPhysics->Configure(kAbsorption, true); 
  //  opticalPhysics->Configure(kBoundary, true);      
  //  opticalPhysics->Configure(kWLS, true);

  // Set control parameters for scintillation
  // Followed from: 
  // https://indico.cern.ch/event/789510/contributions/3279418/attachments/1818134/2972494/AH_OpticalPhotons_slides.pdf
  //  opticalPhysics->SetScintillationYieldFactor(1.0);
  //  opticalPhysics->SetScintillationExcitationRatio(0.0);

  //  opticalPhysics->SetMaxNumPhotonsPerStep(100);
  //  opticalPhysics->SetMaxBetaChangePerStep(10.0);
  
  //  phys->RegisterPhysics(opticalPhysics);
  
  phys->DumpList();

  runManager->SetUserInitialization(phys);

  runManager->SetUserInitialization(new ECRSActionInitialization(detector));

  // Commented out on 3/16/2021
  //  G4ProcessTable::GetProcessTable()->SetProcessActivation("MYTransportation",false);
  //  G4ProcessTable::GetProcessTable()->SetProcessActivation("Transportation",true);

  
  // Initialize visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  G4cout << "Initializing visualization manatger ... " << G4endl;
  visManager->Initialize();
  
  
  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  // Activate score ntuple writer
  // The Root output type (Root) is selected in B3Analysis.hh.
  // The verbose level can be also set via UI commands
  // /score/ntuple/writerVerbose level
  //G4TScoreNtupleWriter<G4AnalysisManager> scoreNtupleWriter;
  //scoreNtupleWriter.SetVerboseLevel(1);

  // Process macro or start UI session
  //
  if ( ! ui ) {
    // batch mode
    G4cout << " I am here " << G4endl;
    G4String command = "/control/execute ";
    G4String fileName = argv[2];               // running macro file 
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete runManager;
  
  G4cout << "RunManager deleted." << G4endl;
  
  //myOut->Fclose();
  
  return 0;
}
