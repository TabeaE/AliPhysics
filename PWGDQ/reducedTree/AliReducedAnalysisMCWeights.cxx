// Analysis task for calculating PCC weights from reduced trees
// Creation date: 2025/09/14
// Authors: Tabea Eder (tabea.maria.eder@cern.ch)

#include "AliReducedAnalysisMCWeights.h"

#include <iostream>
using std::cout;
using std::endl;

#include <TClonesArray.h>
#include <TRandom3.h>
#include <TCanvas.h>

#include "AliHistogramManager.h"
#include "AliReducedVarManager.h"
#include "AliReducedEventInfo.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "AliMCSpectraWeights.h"

ClassImp(AliReducedAnalysisMCWeights);


//_____________________________________________________________________________
AliReducedAnalysisMCWeights::AliReducedAnalysisMCWeights() :
AliReducedAnalysisTaskSE(),
fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
fEventCuts(),
fTrackCuts(),
fComputeMult(kTRUE),
fMeasMultTrackCuts(),
fWeightsTrackCuts(),
fTrueMultTrackCuts(),
fOptionRunOverMC(kFALSE),
fReweightPC(kFALSE),
fMCPCWeights(0x0),
fMCSpectraWeights(0x0)
{
  //
  // default constructor
  //
}


//_____________________________________________________________________________
AliReducedAnalysisMCWeights::AliReducedAnalysisMCWeights(const Char_t* name, const Char_t* title) :
AliReducedAnalysisTaskSE(name, title),
fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
fEventCuts(),
fTrackCuts(),
fComputeMult(kTRUE),
fMeasMultTrackCuts(),
fWeightsTrackCuts(),
fTrueMultTrackCuts(),
fOptionRunOverMC(kFALSE),
fReweightPC(kFALSE),
fMCPCWeights(0x0),
fMCSpectraWeights(0x0)
{
  //
  // named constructor
  //
  fEventCuts.SetOwner(kTRUE);
  fTrackCuts.SetOwner(kTRUE);
  fMeasMultTrackCuts.SetOwner(kTRUE);
  fWeightsTrackCuts.SetOwner(kTRUE);
  fTrueMultTrackCuts.SetOwner(kTRUE);
  fRand.reset(new TRandom3());
}


//_____________________________________________________________________________
// Destructor
AliReducedAnalysisMCWeights::~AliReducedAnalysisMCWeights()
{
  fEventCuts.Clear("C");
  fTrackCuts.Clear("C");
  fMeasMultTrackCuts.Clear("C");
  fWeightsTrackCuts.Clear("C");
  fTrueMultTrackCuts.Clear("C");
  if(fHistosManager) delete fHistosManager;
}


//_____________________________________________________________________________
// Initialize stuff
void AliReducedAnalysisMCWeights::Init()
{
  AliReducedVarManager::SetDefaultVarNames();

  AliReducedVarManager::SetUseVariable(AliReducedVarManager::kDeltaVtxZMC);

  fHistosManager->SetUseDefaultVariableNames(kTRUE);
  fHistosManager->SetDefaultVarNames(AliReducedVarManager::fgVariableNames,
                                     AliReducedVarManager::fgVariableUnits);

}


//_____________________________________________________________________________
// Process the current event
void AliReducedAnalysisMCWeights::Process()
{
  if(!fEvent) return;
  if(!(fEvent->IsA() == AliReducedEventInfo::Class())) {
    cout << "ERROR: AliReducedAnalysisMCWeights::Process() needs AliReducedEventInfo events" << endl;
    return;
  }

  if(fEventCounter%10000 == 0) cout << "Event no. " << fEventCounter << endl;
  fEventCounter++;

  AliReducedVarManager::SetEvent(fEvent);
  // Reset the values array, keep only the run wise data (LHC and ALICE GRP information)
  // NOTE: The run wise data will be updated automatically in the VarManager in case a run number change
  //       is detected.
  for(Int_t i=AliReducedVarManager::kNRunWiseVariables; i<AliReducedVarManager::kNVars; ++i)
    fValues[i] = -9999.;

  // Fill event information before applying event cuts
  AliReducedVarManager::FillEventInfo(fEvent, fValues);
  if(fComputeMult) {
    CountNch05();
    FillMultiplicity();
  }
  Int_t nGlobalEstimators = (fComputeMult ? GetNMeasMultCuts() : 0);
  // Assuming that we have selected 2% unbiased events only for data and not for MC
  Bool_t isEventUnbiased = fEvent->TestEventTag(14) || GetRunOverMC();
  if(isEventUnbiased) {
    fHistosManager->FillHistClass("Event_BeforeCuts", fValues);
    for(UShort_t ibit=0; ibit<64; ++ibit) {
      AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
      fHistosManager->FillHistClass("EventTag_BeforeCuts", fValues);
      AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
      fHistosManager->FillHistClass("EventTriggers_BeforeCuts", fValues);
    }
  } else {
    fHistosManager->FillHistClass("Event_noTag14_BeforeCuts", fValues);
    for(UShort_t ibit=0; ibit<64; ++ibit) {
      AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
      fHistosManager->FillHistClass("EventTag_noTag14_BeforeCuts", fValues);
      AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
      fHistosManager->FillHistClass("EventTriggers_noTag14_BeforeCuts", fValues);
    }
  }

  // Histograms for event selection efficiencies
  if(isEventUnbiased) {
    for(Int_t iev=0; iev<fEventCuts.GetEntries(); ++iev) {
      AliReducedInfoCut* evCut = (AliReducedInfoCut*)fEventCuts.At(iev);
      if(evCut->IsSelected(fEvent, fValues)) {
        fHistosManager->FillHistClass(Form("Event_%s",GetEventCutName(iev)), fValues);
        for(Int_t iest=0; iest<nGlobalEstimators; iest++)
          fHistosManager->FillHistClass(Form("EventMult_%s_%s",GetMeasMultcutName(iest),
                                             GetEventCutName(iev)), fValues);
      }
    }
  }

  // Apply event selection
  if(!IsEventSelected(fEvent, fValues)) return;

  // Fill event info histograms after cuts
  if(isEventUnbiased) {
    fHistosManager->FillHistClass("Event_AfterCuts", fValues);
    for(Int_t icut=0; icut<nGlobalEstimators; icut++)
      fHistosManager->FillHistClass(Form("EventMult_%s",GetMeasMultcutName(icut)), fValues);
    // Correlations between different multiplicity estimators
    fHistosManager->FillHistClass("CorrelMult", fValues);
    for(UShort_t ibit=0; ibit<64; ++ibit) {
      AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
      fHistosManager->FillHistClass("EventTag_AfterCuts", fValues);
      AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
      fHistosManager->FillHistClass("EventTriggers_AfterCuts", fValues);
    }
  } else {
    fHistosManager->FillHistClass("Event_noTag14_AfterCuts", fValues);
    for(UShort_t ibit=0; ibit<64; ++ibit) {
      AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
      fHistosManager->FillHistClass("EventTag_noTag14_AfterCuts", fValues);
      AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
      fHistosManager->FillHistClass("EventTriggers_noTag14_AfterCuts", fValues);
    }
  }

  if(fMCSpectraWeights->GetTaskStatus() < AliMCSpectraWeights::TaskState::kMCSpectraObtained) {
    FillMCSpectra();
  } /*else {
    fMCSpectraWeights->SetCurrentEvent(fMCEvent);
    fMCSpectraWeights->StartNewEvent();
    TString fStoredObjectName = "fMCSpectraWeights";
    // Add to AliVEvent
    auto tmpObject = static_cast<AliMCSpectraWeightsHandler*>(fEvent->FindListObject(fStoredObjectName.Data()));
    if(!tmpObject) {
      AliMCSpectraWeightsHandler* handler = new AliMCSpectraWeightsHandler(fMCSpectraWeights,
                                                                           fStoredObjectName.Data());
      fEvent->AddObject(handler);
    }
  }*/

}

//_____________________________________________________________________________
// Fill MC spectra for first train run
void AliReducedAnalysisMCWeights::FillMCSpectra() {
  if(fMCSpectraWeights->GetTaskStatus() >= AliMCSpectraWeights::TaskState::kMCSpectraObtained)
    return;

  AliReducedTrackInfo* particle = nullptr;
  for(Int_t iArray=1; iArray<=2; iArray++) {
    TClonesArray* tracklist = (iArray==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
    TIter nextTrack(tracklist);
    for(Int_t ipart=0; ipart<tracklist->GetEntries(); ++ipart) {
      particle = (AliReducedTrackInfo*)nextTrack();
      for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i)
        fValues[i] = -9999.;
      AliReducedVarManager::FillTrackInfo(particle, fValues);
      if(!particle->IsMCKineParticle()) continue;

      // Get rapidity and apply cut
      Float_t maxRap = 0.5;  // hard coded max eta; in all papers 0.5
      fValues[AliReducedVarManager::kPCCPartRap] = particle->Rapidity(GetParticleMass(particle));
      if(TMath::Abs(fValues[AliReducedVarManager::kPCCPartRap]) > maxRap)
        continue;

      // Get particle type
      fValues[AliReducedVarManager::kPCCPartType] = GetParticleType(particle);
      if(fValues[AliReducedVarManager::kPCCPartType] < 0)
        continue;

      fHistosManager->FillHistClass("fHistMCGenPrimTrackParticle", fValues);
    }
  }
}


//_____________________________________________________________________________
// Fill multiplicity values
void AliReducedAnalysisMCWeights::FillMultiplicity()
{

  // Fill global tracks (both signal and MC and MC truth number of Jpsi)
  for(Int_t icut=0; icut<GetNMeasMultCuts(); icut++) fValues[AliReducedVarManager::kNGlobalTracks+icut] = 0.;
  fValues[AliReducedVarManager::kMCNchWoPileup] = 0.;
  fValues[AliReducedVarManager::kMCNch09]       = 0.;
  fValues[AliReducedVarManager::kMCNJpsi]       = 0.;

  // Loop over both arrays
  Int_t trackID = 0;  // use unique track ID as seed for GetNRepetitions
  for(Int_t iArray=1; iArray<=2; iArray++) {
    AliReducedTrackInfo* track;
    TClonesArray* tracklist = (iArray==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
    TIter nextTrack(tracklist);
    // Loop over tracks
    for(Int_t it=0; it<tracklist->GetEntries(); ++it) {
      track = (AliReducedTrackInfo*)nextTrack();
      for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i)
        fValues[i] = -9999.;
      AliReducedVarManager::FillTrackInfo(track, fValues);

      // TODO kMCNchWoPileup is filled here, because kMCNch is defined such that it includes pileup tracks if
      //      one has MC with pileup. The same thing could also be achieved by introducing a TrueMultTrackCut
      //      for |eta|<1, but atm this code is not safe for multiple TrueMultTrackCuts.
      // TODO kMCNch excludes tracks from jpsi daughters per default. Is this done properly here? Since I test
      //      this on MC w/o pileup, I would think kMCNch should not include pileup and be "correct".
      if(track->IsMCTruth() && track->Charge() && abs(track->Eta()) < 1.)
        fValues[AliReducedVarManager::kMCNchWoPileup] ++;

      // Get measured multiplicity (track has to be reconstructed & selected by at least one meas mult cut)
      if(!track->IsMCTruth() && IsTrackMeasMultSelected(track,fValues)) {
        Float_t weightPCC = 1.;
        Int_t   nRepsPCC  = 1;
        if(fOptionRunOverMC && fReweightPC) {
          weightPCC = GetParticleWeight(track);
          nRepsPCC  = GetNRepetitions(weightPCC, trackID);
          fValues[AliReducedVarManager::kPCCWeight] = weightPCC;
        }

        for(Int_t icut=0; icut<GetNMeasMultCuts(); icut++) {
          if(track->TestMultFlag(icut)) {
            TH2F*   hWeightsTrack = (TH2F*) fWeightsTrackCuts.At(icut);
            Float_t weightTotal   = weightPCC * hWeightsTrack->GetBinContent(hWeightsTrack->FindBin(
                                                  fValues[AliReducedVarManager::kRunNo],track->Pt()));
            Int_t   nRepsTotal    = GetNRepetitions(weightTotal, trackID);

            fValues[AliReducedVarManager::kNGlobalTracks+icut] += nRepsTotal;
            fValues[AliReducedVarManager::kPCCnRepetitions]     = nRepsPCC;

            const char* cutName = GetMeasMultcutName(icut);
            fHistosManager->FillHistClass(Form("TrackMult_%s_MeasMult",cutName), fValues);
          }
        }
      }

      // Get true multiplicity (track has to be true MC (not rec) & selected by at least one true mult cut)
      if(track->IsMCTruth() && IsTrackTrueMultSelected(track, fValues)) {
        Int_t nReps = 1;
        if(fOptionRunOverMC && fReweightPC) {
          Float_t weightPCC = GetParticleWeight(track);
          nReps = GetNRepetitions(weightPCC, trackID);
          fValues[AliReducedVarManager::kPCCWeight] = weightPCC;
        }

        fValues[AliReducedVarManager::kMCNch09]         += nReps;
        fValues[AliReducedVarManager::kPCCnRepetitions]  = nReps;
        for(UShort_t flag=0; flag<(8*sizeof(UInt_t)); ++flag) {
          AliReducedVarManager::FillTrackMCFlag(track, flag, fValues);
          fHistosManager->FillHistClass("TrackMultMCFlag", fValues);
        }
        fHistosManager->FillHistClass(Form("TrackMult_TrueMult"), fValues);
      }

      trackID++;
    }  // end loop over tracks
  }  // end loop over arrays

  // Fill kMCNch09+1 for MC truth accepted events only
  // (necessary to compute contamination of events with |vtx_z|>10cm and Nch=0)
  if(fOptionRunOverMC) {
    fValues[AliReducedVarManager::kMCNch09+1] = fValues[AliReducedVarManager::kMCNch09];
    if(fValues[AliReducedVarManager::kMCNch09] < 1 || abs(fValues[AliReducedVarManager::kVtxZMC]) > 10.) {
      fValues[AliReducedVarManager::kMCNch09+1] = -9999.;
    }
  }
}


//_____________________________________________________________________________
// Apply event cuts
Bool_t AliReducedAnalysisMCWeights::IsEventSelected(AliReducedBaseEvent* event, Float_t* values/*=0x0*/)
{
  if(fEventCuts.GetEntries() == 0) return kTRUE;
  // Loop over all the cuts and make a logical and between all cuts in the list
  for(Int_t i=0; i<fEventCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fEventCuts.At(i);
    if(values) { if(!cut->IsSelected(event,values)) return kFALSE; }
    else       { if(!cut->IsSelected(event))        return kFALSE; }
  }
  return kTRUE;
}


//_____________________________________________________________________________
// Apply cuts for determining measured multiplicity
Bool_t AliReducedAnalysisMCWeights::IsTrackMeasMultSelected(AliReducedBaseTrack* track,
                                                            Float_t* values/*=0x0*/)
{
  if(fMeasMultTrackCuts.GetEntries() == 0) return kTRUE;
  if(track->IsMCTruth())                   return kFALSE;
  track->ResetMultFlags();

  for(Int_t i=0; i<fMeasMultTrackCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*) fMeasMultTrackCuts.At(i);
    if(values) { if(cut->IsSelected(track, values)) track->SetMultFlag(i); }
    else       { if(cut->IsSelected(track))         track->SetMultFlag(i); }
  }
  return (track->GetMultFlags()>0 ? kTRUE : kFALSE);
}


//_____________________________________________________________________________
// Apply cuts for determining MC true multiplicity
Bool_t AliReducedAnalysisMCWeights::IsTrackTrueMultSelected(AliReducedBaseTrack* track,
                                                            Float_t* values/*=0x0*/)
{
  if(fTrueMultTrackCuts.GetEntries() == 0) return kTRUE;
  if(!track->IsMCTruth())                  return kFALSE;
  track->ResetMultFlags();

  for(Int_t i=0; i<fTrueMultTrackCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*) fTrueMultTrackCuts.At(i);
    if(values) { if(cut->IsSelected(track, values)) track->SetMultFlag(i); }
    else       { if(cut->IsSelected(track))         track->SetMultFlag(i); }
  }
  return (track->GetMultFlags()>0 ? kTRUE : kFALSE);
}


//_____________________________________________________________________________
// Count Nch in |eta|<0.5 with cuts to match expert input for PCC.
void AliReducedAnalysisMCWeights::CountNch05()
{
  const Double_t lowPt = 0.05;
  const Double_t eta   = 0.5;

  fValues[AliReducedVarManager::kMCNch05] = 0.;

  // Loop over both arrays
  for(Int_t iArray=1; iArray<=2; iArray++) {
    AliReducedTrackInfo* track;
    TClonesArray* tracklist = (iArray==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
    TIter nextTrack(tracklist);
    // Loop over tracks
    for(Int_t it=0; it<tracklist->GetEntries(); ++it) {
      track = (AliReducedTrackInfo*)nextTrack();
      if(!track->IsMCTruth())                 continue;
      if(!track->TestMCFlag(0))               continue;  // secondary rejection  // TODO hard-coded
      if(TMath::Abs(track->Charge()) < 0.01)  continue;  // neutral rejection
      if(TMath::Abs(track->Eta())    > eta)   continue;  // acceptance cut
      if(track->Pt()                 < lowPt) continue;  // low pT cut
      fValues[AliReducedVarManager::kMCNch05]++;
    }  // end track loop
  }  // end array loop
}


//_____________________________________________________________________________
// Get PCC weight
Float_t AliReducedAnalysisMCWeights::GetParticleWeight(AliReducedTrackInfo* track)
{
  // The weight depends on pt, species and multiplicity class
  // It is used while reweighting particle composition in MC in order to match it to data
  // Particle type: pion = 0, proton = 1, kaon = 2, sigma- = 3, sigma+ = 4, rest = 5, lambda = 6
  //                weight_lambda = weight_xi = weight_sigma+
  // See $ALICE_PHYSICS/PWG/Tools/AliMCSpectraWeights.h and https://alice-notes.web.cern.ch/node/1311

  Float_t ptMC = track->Pt();
  Int_t   type = 5;
  // TODO For now flags are hard-coded
  if      (track->TestMCFlag(1)) type = 0;  // pion   primary
  else if (track->TestMCFlag(2)) type = 1;  // proton primary
  else if (track->TestMCFlag(3)) type = 2;  // kaon   primary
  else if (track->TestMCFlag(4)) type = 3;  // sigma- primary
  else if (track->TestMCFlag(5)) type = 4;  // sigma+ primary
  else if (track->TestMCFlag(6)) type = 4;  // lambda primary

  // Check if particle is a detected secondary from weak decay
  // (as opposed to secondaries from detector material interaction).
  // In that case the weight of the mother has to be applied for all daughters.
  // In that case the mother particle was put in array 1.
  if(type == 5 && !track->IsMCTruth() && !track->TestMCFlag(0)) {
    Int_t pdgmother = abs(track->MCPdg(1));
    if      (pdgmother == 3222) type = 4;  // mother is sigma+  -> use sigma+ weight
    else if (pdgmother == 3112) type = 3;  // mother is sigma-  -> use sigma- weight
    else if (pdgmother == 310)  type = 2;  // mother is K0short -> use kaon   weight
    else if (pdgmother == 130)  type = 2;  // mother is K0long  -> use kaon   weight
    else if (pdgmother == 311)  type = 2;  // mother is K0short -> use kaon   weight
    else if (pdgmother == 321)  type = 2;  // mother is kaon    -> use kaon   weight
    else if (pdgmother == 211)  type = 0;  // mother is pion    -> use pion   weight
    else if (pdgmother == 3122) type = 6;  // mother is lambda  -> use sigma+ weight

    if(type != 5) {
      // The mother is strange primary - we look for its pt in array 1
      Int_t labelmother = track->MCLabel(1);
      // loop over the first track array
      TClonesArray* trackList = fEvent->GetTracks();
      TIter nextTrack(trackList);
      AliReducedTrackInfo* mother;
      for(Int_t it=0; it<trackList->GetEntries(); ++it) {
        mother = (AliReducedTrackInfo*)nextTrack();
        if(mother->IsMCTruth() && (mother->MCLabel(0)==labelmother)) {
          // If mother is a primary, get the pt of the mother to determine weight factor
          if(mother->TestMCFlag(0)) ptMC = mother->PtMC();
          else                      type = 5;
          break;
        }
        // We did not find the mother: do nothing
        if(it == trackList->GetEntries()-1) type = 5;
      }
    }
  }
  if(type == 5) return 1;

  // Get Nch (|eta|<0.5) for multiplicity class
  Float_t mult = fValues[AliReducedVarManager::kMCNch05];

  Int_t iptMC = fMCPCWeights->GetXaxis()->FindBin(ptMC);
  Int_t imult = fMCPCWeights->GetYaxis()->FindBin(mult);
  Int_t itype = fMCPCWeights->GetZaxis()->FindBin(type);

  // In case out of the range
  if(ptMC < fMCPCWeights->GetXaxis()->GetBinLowEdge(1) ||
     ptMC > fMCPCWeights->GetXaxis()->GetBinUpEdge(fMCPCWeights->GetNbinsX()))
  { std::cout << "PCC: pT out of range!!!!!!!!" << std::endl; }

  if(mult < fMCPCWeights->GetYaxis()->GetBinLowEdge(1) ||
     mult > fMCPCWeights->GetYaxis()->GetBinUpEdge(fMCPCWeights->GetNbinsY()))
  { std::cout << "PCC: mult out of range!!!!!!!!" << std::endl; }

  Float_t weight = fMCPCWeights->GetBinContent(iptMC, imult, itype);

  return weight;
}


//_____________________________________________________________________________
// Get particle type
Int_t AliReducedAnalysisMCWeights::GetParticleType(AliReducedTrackInfo* track)
{
  // Particle type: pion = 0, proton = 1, kaon = 2, sigma- = 3, sigma+ = 4, rest = 5, lambda = 6

  Int_t   type = 5;
  // TODO For now MC flags are hard-coded
  if      (track->TestMCFlag(1)) type = AliMCSpectraWeights::ParticleType::kPion;        // pion   primary
  else if (track->TestMCFlag(2)) type = AliMCSpectraWeights::ParticleType::kProtons;     // proton primary
  else if (track->TestMCFlag(3)) type = AliMCSpectraWeights::ParticleType::kKaon;        // kaon   primary
  else if (track->TestMCFlag(4)) type = AliMCSpectraWeights::ParticleType::kSigmaMinus;  // sigma- primary
  else if (track->TestMCFlag(5)) type = AliMCSpectraWeights::ParticleType::kSigmaPlus;   // sigma+ primary
  else if (track->TestMCFlag(6)) type = AliMCSpectraWeights::ParticleType::kLambda;      // lambda primary
  else                           type = AliMCSpectraWeights::ParticleType::kRest;        // other (rest)

  return type;
}


//_____________________________________________________________________________
// Get particle true mass
Float_t AliReducedAnalysisMCWeights::GetParticleMass(AliReducedTrackInfo* track)
{
  Float_t mass = -999.;
  // TODO For now flags and masses are hard-coded
  if      (track->TestMCFlag(1))  mass = 0.13957039;     // pion
  else if (track->TestMCFlag(2))  mass = 0.93827208816;  // proton
  else if (track->TestMCFlag(3))  mass = 0.493677;       // kaon
  else if (track->TestMCFlag(4))  mass = 1.197449;       // sigma-
  else if (track->TestMCFlag(5))  mass = 1.18937;        // sigma+
  else if (track->TestMCFlag(6))  mass = 1.115683;       // lambda
  else if (track->TestMCFlag(7))  mass = 0.497611;       // K0S  TODO AliReducedVarManager::fgkParticleMass: 614
  else if (track->TestMCFlag(8))  mass = 0.497611;       // K0L
  else if (track->TestMCFlag(9))  mass = 0.493677;       // kaon
  else if (track->TestMCFlag(10)) mass = 1.115683;       // lambda

  return mass;
}


//_____________________________________________________________________________
// Decide how often to repeat particle in MC to match data.
Int_t AliReducedAnalysisMCWeights::GetNRepetitions(Float_t scalingFactor, Int_t part)
{
  Int_t   nRepetitions = (Int_t)scalingFactor;
  Float_t rest         = scalingFactor - nRepetitions;

  fRand->SetSeed(GetSeed(part));
  nRepetitions += (fRand->Rndm()<=rest) ? 1 : 0;
  return nRepetitions;
}


//_____________________________________________________________________________
// Define random (but reproducable) seed
unsigned long AliReducedAnalysisMCWeights::GetSeed(Int_t part)
{
  unsigned long seed = (unsigned long)fValues[AliReducedVarManager::kEventNumberInFile];
  seed <<= 5;
  seed += (int)fValues[AliReducedVarManager::kRunNo];
  seed <<= 2;
  seed += part;
  seed <<= 3;
  seed += (unsigned int)fValues[AliReducedVarManager::kTimeStamp];
  return seed;
}


//_____________________________________________________________________________
// Run stuff after the event loop
void AliReducedAnalysisMCWeights::Finish() {}


//_____________________________________________________________________________
// Search the track list for pure MC track with label and return the track pointer
AliReducedTrackInfo* AliReducedAnalysisMCWeights::FindTrackByLabel(Int_t label, Bool_t isTruth)
{
  AliReducedTrackInfo* track     = nullptr;
  TClonesArray*        trackList = fEvent->GetTracks();
  TIter nextTrack(trackList);
  for(Int_t i=0; i<trackList->GetEntries(); ++i) {
    track = (AliReducedTrackInfo*)nextTrack();
    if(isTruth  && (!track->IsMCKineParticle())) continue;
    if(!isTruth &&   track->IsMCKineParticle())  continue;
    if(track->MCLabel(0) == label) return track;
  }
  return nullptr;
}

