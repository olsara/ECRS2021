// 6/2/2014, Modified by hexc, Olesya
// 7/21/2014, Hexc, Olesya
// Added code to write out timing information
// 10/1/2014: Hexc, Olesya - Adding code to write root output for data
//    analysis.
// 12/1/2014: Hexc, Olesya - Add "(aTrack->GetParentID() != 0)" condition for not writing out
//    the track info of the incident (i.e. primary) particle
// 2/12/2015: Hexc, Olesya - Add event ID into the output ntuple. Removed altitude calculation.
// 9/20/2023: Hexc - replaced g4root.hh with G4AnalysisManager.hh

#include "ECRSTrackingAction.hh"
#include "ECRSTrackingMessenger.hh"
#include "ECRSParticleFlag.hh"   // 9/28/2026: shared particle-name -> flag lookup

#include "G4UserTrackingAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4UnitsTable.hh"
#include "G4AnalysisManager.hh"

//#include "g4root.hh"

using namespace CLHEP;

ECRSTrackingAction::ECRSTrackingAction(ECRSRunAction* run_action_in, ECRSEventAction* evt_action_in)
{
  run_action = run_action_in;
  evt_action = evt_action_in;
  fileOut = ECRSSingleton::instance();
  trackMessenger = new ECRSTrackingMessenger(this);
}

ECRSTrackingAction::~ECRSTrackingAction()
{
  delete trackMessenger;
}

void ECRSTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // Here we print the Particle Data Group particle id number as well as the
  // energy when a new particle is created for tracking
  if (aTrack->GetParentID() != 0) {
    G4ParticleDefinition* particle = aTrack->GetDefinition();
    G4int PID = particle->GetPDGEncoding();
    G4double KE = aTrack->GetKineticEnergy();
    G4double Global_time = aTrack->GetGlobalTime();
    G4double Local_time = aTrack->GetLocalTime();
    G4double Proper_time = aTrack->GetProperTime();
    G4ThreeVector P = aTrack->GetMomentum();
    
    G4String particleName = particle->GetParticleName();
    G4ThreeVector position = aTrack->GetPosition();
    
    // find and set event ID
    G4int evtID = evt_action->getCurrentEventID();
    
    G4int flagParticle = ECRSParticleFlag(particleName);   // 9/28/2026: shared lookup (was a 17-way if/else chain)

    /* We comment this text output since the data is stored in the tracking ntuple
    fileOut->fout << 4  << "  "  << particle->GetParticleName() <<"  "<<PID<<"  " << "  " <<  KE/MeV  << "  "
		  << position.getX()/m  << "  " << position.getY()/m  << "  " <<  position.getZ()/m  << "   " 
		  << Global_time/ms << "   "  << Local_time/ms << "   "  << Proper_time/ms << G4endl;
    */    

    // get analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    /*
    G4double altitude = sqrt( position.getX()/m * position.getX()/m + 
			      position.getY()/m * position.getY()/m + 
			      position.getZ()/m * position.getZ()/m) - 6371200;
    */

    // fill ntuple
    analysisManager->FillNtupleDColumn(1, run_action->GetNtColID(0), flagParticle);
    analysisManager->FillNtupleDColumn(1, run_action->GetNtColID(1), KE/MeV);
    analysisManager->FillNtupleDColumn(1, run_action->GetNtColID(2), position.getX()/m);
    analysisManager->FillNtupleDColumn(1, run_action->GetNtColID(3), position.getY()/m);
    analysisManager->FillNtupleDColumn(1, run_action->GetNtColID(4), position.getZ()/m);
    //    analysisManager->FillNtupleDColumn(1, run_action->GetNtColID(5), altitude);
    analysisManager->FillNtupleDColumn(1, run_action->GetNtColID(5), P.getX()/MeV);
    analysisManager->FillNtupleDColumn(1, run_action->GetNtColID(6), P.getY()/MeV);
    analysisManager->FillNtupleDColumn(1, run_action->GetNtColID(7), P.getZ()/MeV);
    analysisManager->FillNtupleDColumn(1, run_action->GetNtColID(8), Global_time/ms);
    analysisManager->FillNtupleDColumn(1, run_action->GetNtColID(9), Local_time/ms);
    analysisManager->FillNtupleDColumn(1, run_action->GetNtColID(10), Proper_time/ms);
    analysisManager->FillNtupleDColumn(1, run_action->GetNtColID(11), evtID);
    analysisManager->AddNtupleRow(1);
    
  }
  
}

void ECRSTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
}
