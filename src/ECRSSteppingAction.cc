// 7/21/2014: Hexc, Olesya:
// Modified the output format to write out hits info when the hits are close to 
// the surface of the Earth.
// 2/11/2015: Hexc, Olesya - Add eventAction and the primaryGeneratorAction in order to 
//              output the event number and the primary particle info.
// 4/28/2015: Hexc, Olesya - Add aTrack->SetTrackStatus(fStopAndKill); when the hit is below the surface.
// 10/31/2017: Hexc - Stop tracking a particle after 10 seconds.
// 9/20/2023: Hexc - replaced g4root.hh with G4AnalysisManager.hh
//
#include "ECRSSteppingAction.hh"
#include "ECRSSteppingMessenger.hh"
#include "ECRSDetectorConstruction.hh"
#include "ECRSSingleton.hh"
#include "ECRSParticleFlag.hh"   // 9/28/2026: shared particle-name -> flag lookup

#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4AnalysisManager.hh"

//#include "g4root.hh"

using namespace CLHEP;

ECRSSteppingAction::ECRSSteppingAction(ECRSDetectorConstruction* det, 
				       ECRSRunAction* run_action_in, 
				       ECRSEventAction* evt_action_in)
  :detector(det),run_action(run_action_in),evt_action(evt_action_in),
   aTrack(NULL),particle(NULL),
   currentVolume(NULL),nextVolume(NULL),point(NULL)
{
  fileOut = ECRSSingleton::instance();
  KE = 0;
  particleName = "";
  theMessenger = new ECRSSteppingMessenger(this);
}

ECRSSteppingAction::~ECRSSteppingAction()
{
  delete theMessenger;
}

void ECRSSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  aTrack = aStep->GetTrack();
  currentVolume = aTrack->GetVolume();
  nextVolume = aTrack->GetNextVolume();
  G4double Global_time = aTrack->GetGlobalTime();

  // Stop tracking when the global time exceeds 10 seconds;
  if (Global_time > 10.0*s){
    aTrack->SetTrackStatus(fStopAndKill);
    return;
  }
  
  if ((currentVolume != detector->GetECRS()) &&
      (nextVolume == detector->GetECRS()))         // output the hits info while crossing the earth surface boundary
    {
      //  G4int PID = particle->GetPDGEncoding();
 
      G4double Local_time = aTrack->GetLocalTime();
      G4double Proper_time = aTrack->GetProperTime();
      KE = aTrack->GetKineticEnergy();
      particle = aTrack->GetDefinition();
      particleName = particle->GetParticleName();   // 9/28/2026: removed now-unused PID (only used by the deleted fout write)

      // find and set event ID
      G4int evtID = evt_action->getCurrentEventID();

      G4ThreeVector position = aTrack->GetPosition();
      G4ThreeVector mom = aTrack->GetMomentum();
      
      // 9/28/2026: Removed a per-hit write to the shared ECRSSingleton::fout
      // stream. That stream is never opened, so the write produced no output,
      // and being a shared global it was a data race under multithreading.
      // The same information is already stored in the hits ntuple below.
      // (Was: fileOut->fout << 0 << "  " << particle->GetParticleName() << ... )


      G4int flagParticle = ECRSParticleFlag(particleName);   // 9/28/2026: shared lookup (was a 17-way if/else chain)

      // get analysis manager
      G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
      
      // fill ntuple (ntup[le ID = 2)
      analysisManager->FillNtupleDColumn(2, run_action->GetNtColID(12), flagParticle);
      analysisManager->FillNtupleDColumn(2, run_action->GetNtColID(13), KE/MeV);
      analysisManager->FillNtupleDColumn(2, run_action->GetNtColID(14), position.getX()/m);
      analysisManager->FillNtupleDColumn(2, run_action->GetNtColID(15), position.getY()/m);
      analysisManager->FillNtupleDColumn(2, run_action->GetNtColID(16), position.getZ()/m);
      analysisManager->FillNtupleDColumn(2, run_action->GetNtColID(17), mom.getX()/MeV);
      analysisManager->FillNtupleDColumn(2, run_action->GetNtColID(18), mom.getY()/MeV);
      analysisManager->FillNtupleDColumn(2, run_action->GetNtColID(19), mom.getZ()/MeV);
      analysisManager->FillNtupleDColumn(2, run_action->GetNtColID(20), Global_time/ms);
      analysisManager->FillNtupleDColumn(2, run_action->GetNtColID(21), Local_time/ms);  
      analysisManager->FillNtupleDColumn(2, run_action->GetNtColID(22), evtID);  // added on 2/10/2015

      analysisManager->AddNtupleRow(2);

      // A quick solution for stopping tracks below the surface
      // 
      aTrack->SetTrackStatus(fStopAndKill);  
      
    } 
}
