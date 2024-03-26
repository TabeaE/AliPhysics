//
// Creation date: 2017/08/09
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#include "AliReducedAnalysisFilterTrees.h"

#include <iostream>
using std::cout;
using std::endl;

#include <TClonesArray.h>
#include <TIterator.h>
#include <TRandom.h>
#include <TF1.h>

#include "AliReducedVarManager.h"
#include "AliReducedEventInfo.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedPairInfo.h"
#include "AliHistogramManager.h"

ClassImp(AliReducedAnalysisFilterTrees);


//___________________________________________________________________________
AliReducedAnalysisFilterTrees::AliReducedAnalysisFilterTrees() :
AliReducedAnalysisTaskSE(),
fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
fMixingHandler(new AliMixingHandler("J/psi signal extraction","", AliMixingHandler::kMixResonanceLegs)),
fMixingHandlerMult(),
fMultBinsMixing(),
fNMultBinsMixing(0),
fEventCuts(),
fTrackCuts(),
fWriteFilteredTracks(kTRUE),
fWriteFilteredTracksCandidatesOnly(kFALSE),
fFillTrackV0Histograms(kFALSE),
fOptionRunMixing(kTRUE),
fOptionRunMixingMult(kFALSE),
fComputeMult(kTRUE),
fPairCuts(),
fWriteFilteredPairs(kTRUE),
fRejectEmptyEvents(kFALSE),
fBuildCandidatePairs(kFALSE),
fBuildCandidateLikePairs(kFALSE),
fMCTruthJpsi2eeOnly(kFALSE),
fRegionsToMCTruth(kFALSE),
fDefaultRandomPhi(kTRUE),
fMinPtLeading(0.0),
fReweightCut(-1),
fCandidateType(AliReducedPairInfo::kJpsiToEE),
fLeg1Cuts(),
fLeg2Cuts(),
fCandidatePairCuts(),
fRunCandidatePrefilter(kFALSE),
fRunCandidatePrefilterOnSameCharge(kFALSE),
fJpsiMassDist(),
fMeasuredMultTrackCuts(),
fTrueMultTrackCuts(),
fLeg1PrefilterCuts(),
fLeg2PrefilterCuts(),
fLeg1PairPrefilterCuts(),
fLeg2PairPrefilterCuts(),
fLeg1Tracks(),
fLeg2Tracks(),
fLeg1PrefilteredTracks(),
fLeg2PrefilteredTracks(),
fOptionRunOverMC(kFALSE),
fLegCandidatesMCcuts(),
fLegCandidatesMCcuts_RequestSameMother(),
fJpsiMotherMCcuts(),
fMCJpsiPtWeights(0x0),
fSkipMCEvent(kFALSE),
fJpsiElectronMCcuts()
{
  //
  // default constructor
  //
}


//___________________________________________________________________________
AliReducedAnalysisFilterTrees::AliReducedAnalysisFilterTrees(const Char_t* name, const Char_t* title) :
AliReducedAnalysisTaskSE(name,title),
fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
fMixingHandler(new AliMixingHandler("J/psi signal extraction","", AliMixingHandler::kMixResonanceLegs)),
fMixingHandlerMult(),  
fMultBinsMixing(),
fNMultBinsMixing(0),  
fEventCuts(),
fTrackCuts(),
fWriteFilteredTracks(kTRUE),
fWriteFilteredTracksCandidatesOnly(kFALSE),
fFillTrackV0Histograms(kFALSE),
fRejectEmptyEvents(kFALSE),
fOptionRunMixing(kTRUE),
fOptionRunMixingMult(kFALSE),
fComputeMult(kTRUE),
fPairCuts(),
fWriteFilteredPairs(kTRUE),
fBuildCandidatePairs(kFALSE),
fBuildCandidateLikePairs(kFALSE),
fMCTruthJpsi2eeOnly(kFALSE),
fRegionsToMCTruth(kFALSE),
fDefaultRandomPhi(kTRUE),
fMinPtLeading(0.0),
fReweightCut(-1),
fCandidateType(AliReducedPairInfo::kJpsiToEE),
fLeg1Cuts(),
fLeg2Cuts(),
fCandidatePairCuts(),
fRunCandidatePrefilter(kFALSE),
fRunCandidatePrefilterOnSameCharge(kFALSE),
fJpsiMassDist(),
fMeasuredMultTrackCuts(),
fTrueMultTrackCuts(),
fLeg1PrefilterCuts(),
fLeg2PrefilterCuts(),
fLeg1PairPrefilterCuts(),
fLeg2PairPrefilterCuts(),
fLeg1Tracks(),
fLeg2Tracks(),
fLeg1PrefilteredTracks(),
fLeg2PrefilteredTracks(),
fOptionRunOverMC(kFALSE),
fLegCandidatesMCcuts(),
fLegCandidatesMCcuts_RequestSameMother(),
fJpsiMotherMCcuts(),
fMCJpsiPtWeights(0x0),
fSkipMCEvent(kFALSE),
fJpsiElectronMCcuts()  
{
  //
  // named constructor
  //

  fMixingHandlerMult.SetOwner(kTRUE);
  fEventCuts.SetOwner(kTRUE);
  fTrackCuts.SetOwner(kTRUE);
  fPairCuts.SetOwner(kTRUE);
  fLeg1Cuts.SetOwner(kTRUE);
  fLeg2Cuts.SetOwner(kTRUE);
  fCandidatePairCuts.SetOwner(kTRUE);
  fMeasuredMultTrackCuts.SetOwner(kTRUE);
  fTrueMultTrackCuts.SetOwner(kTRUE);
  fLeg1PrefilterCuts.SetOwner(kTRUE);
  fLeg2PrefilterCuts.SetOwner(kTRUE);
  fLeg1PairPrefilterCuts.SetOwner(kTRUE);
  fLeg2PairPrefilterCuts.SetOwner(kTRUE);
  fLeg1Tracks.SetOwner(kFALSE);
  fLeg2Tracks.SetOwner(kFALSE);
  fLeg1PrefilteredTracks.SetOwner(kFALSE);
  fLeg2PrefilteredTracks.SetOwner(kFALSE);
  fLegCandidatesMCcuts.SetOwner(kTRUE);
  fJpsiMotherMCcuts.SetOwner(kTRUE);
  fJpsiElectronMCcuts.SetOwner(kTRUE);
  for(Int_t i=0; i<32; ++i) fLegCandidatesMCcuts_RequestSameMother[i] = kTRUE;
}


//___________________________________________________________________________
AliReducedAnalysisFilterTrees::~AliReducedAnalysisFilterTrees() 
{
  //
  // destructor
  //

  fEventCuts.Clear("C");             fTrackCuts.Clear("C");
  fPairCuts.Clear("C");              fLeg1Cuts.Clear("C");
  fLeg2Cuts.Clear("C");              fCandidatePairCuts.Clear("C");
  fMeasuredMultTrackCuts.Clear("C"); fTrueMultTrackCuts.Clear("C");
  fLeg1PrefilterCuts.Clear("C");     fLeg1PairPrefilterCuts.Clear("C");
  fLeg2PrefilterCuts.Clear("C");     fLeg2PairPrefilterCuts.Clear("C");
  fLeg1Tracks.Clear("C");            fLeg1PrefilteredTracks.Clear("C");
  fLeg2Tracks.Clear("C");            fLeg2PrefilteredTracks.Clear("C");
  if(fHistosManager) delete fHistosManager;
  if(fMixingHandler) delete fMixingHandler;
  if(fJpsiMassDist)  delete fJpsiMassDist;
  fMixingHandlerMult.Clear("C");
}


//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::Init() {
  //
  // initialize stuff
  //

  AliReducedVarManager::SetDefaultVarNames();

  AliReducedVarManager::SetUseVariable(AliReducedVarManager::kPseudoProperDecayTime);
  AliReducedVarManager::SetUseVariable(AliReducedVarManager::kPseudoProperDecayTimeMC);
  AliReducedVarManager::SetUseVariable(AliReducedVarManager::kPseudoProperTimeError);
  AliReducedVarManager::SetUseVariable(AliReducedVarManager::kTriggerPseudoProperDecayTime);
  AliReducedVarManager::SetUseVariable(AliReducedVarManager::kTripletPseudoProperDecayTime);

  fHistosManager->SetUseDefaultVariableNames(kTRUE);
  fHistosManager->SetDefaultVarNames(AliReducedVarManager::fgVariableNames,
                                     AliReducedVarManager::fgVariableUnits);
  fMixingHandler->SetHistogramManager(fHistosManager);
}


//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::Process() {
  //
  // process the current event
  //  
  
  if(!fEvent) return;
  if(!(fEvent->IsA() == AliReducedEventInfo::Class())) {
    cout << "ERROR: AliReducedAnalysisFilterTrees::Process() needs AliReducedEventInfo events" << endl;
    return;
  }

  if(fEventCounter%10000 == 0) cout << "Event no. " << fEventCounter << endl;
  fEventCounter++;
  
  AliReducedVarManager::SetEvent(fEvent);
  // reset the values array, keep only the run wise data (LHC and ALICE GRP information)
  // NOTE: The run wise data will be updated automatically in the VarManager in case a run number change
  //       is detected.
  for(Int_t i=AliReducedVarManager::kNRunWiseVariables; i<AliReducedVarManager::kNVars; ++i)
    fValues[i] = -9999.;
  
  // fill event information before event cuts
  AliReducedVarManager::FillEventInfo(fEvent, fValues);
  
  if(fComputeMult) FillMultiplicity(kFALSE);
  Int_t nGlobalEstimators = (fComputeMult ? GetNMeasMultCuts() : 0);
  
  // Assuming that we have selected 2% unbiased events only for data and not for MC
  Bool_t isEventUnbiased = fEvent->TestEventTag(14) || GetRunOverMC();
  
  if(isEventUnbiased) {
    fHistosManager->FillHistClass("Event_MB_BeforeCuts", fValues);
    for(UShort_t ibit=0; ibit<64; ++ibit) {
      AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
      fHistosManager->FillHistClass("EventTag_MB_BeforeCuts", fValues);
      AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
      fHistosManager->FillHistClass("EventTriggers_MB_BeforeCuts", fValues);
    }
  }
  else {
    fHistosManager->FillHistClass("Event_noTag14_BeforeCuts", fValues);
    for(UShort_t ibit=0; ibit<64; ++ibit) {
      AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
      fHistosManager->FillHistClass("EventTag_noTag14_BeforeCuts", fValues);
      AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
      fHistosManager->FillHistClass("EventTriggers_noTag14_BeforeCuts", fValues);
    }
  }
  
  
  Bool_t isEventSelected = IsEventSelected(fEvent);
  
  // Trigger efficiency
  if(isEventUnbiased && isEventSelected && GetRunOverMC() && (fValues[AliReducedVarManager::kMCNch]>0)) {
    fHistosManager->FillHistClass("Event_INELGT0", fValues);
    for(Int_t icut=0; icut<nGlobalEstimators; icut++)
      fHistosManager->FillHistClass(Form("EventMult_%s_INELGT0",GetMeasMultcutName(icut)), fValues);
  }
  if(GetRunOverMC() && !(fValues[AliReducedVarManager::kINT7Triggered])) isEventSelected = kFALSE;
  
  //Vtx reconstruction
  if(isEventUnbiased && isEventSelected) {
    fHistosManager->FillHistClass("Event_NoVtxRec", fValues);
    for(Int_t icut=0; icut<nGlobalEstimators; icut++)
      fHistosManager->FillHistClass(Form("EventMult_%s_NoVtxRec",GetMeasMultcutName(icut)), fValues);
  }
  if(fValues[AliReducedVarManager::kNVtxContributors] < 0.1) isEventSelected = kFALSE;
  
  //Vtx z cut
  if(isEventUnbiased && isEventSelected) {
    fHistosManager->FillHistClass("Event_NoVtxzCut", fValues);
    for(Int_t icut=0; icut<nGlobalEstimators; icut++)
      fHistosManager->FillHistClass(Form("EventMult_%s_NoVtxzCut",GetMeasMultcutName(icut)), fValues);
  }
  if(abs(fValues[AliReducedVarManager::kVtxZ]) > 10) isEventSelected = kFALSE;
  
  // Reset the multiplicity estimators if event was not selected
  if(!isEventSelected) {
    for(Int_t iest=AliReducedVarManager::kMultiplicity;
        iest<AliReducedVarManager::kMultiplicity+AliReducedVarManager::kNMultiplicityEstimators; iest++)
      fValues[iest] = -9999.;
    fValues[AliReducedVarManager::kMCNch09+2] = -9999.;
  }
  

  if(isEventUnbiased) {
    // For multiplicity unfolding
    for(Int_t cutMode=0; cutMode<2*nGlobalEstimators; cutMode=cutMode+2) {
      // For MC, only the smearing matrix is important (supposed independent of the trigger)
      fHistosManager->FillHistClass(Form("pPb_5TeV_Data_cutMode_%d", cutMode+100), fValues);
      fHistosManager->FillHistClass(Form("pPb_5TeV_MC_cutMode_%d",   cutMode+100), fValues);
      
      // Fill for each jpsi counts
      for(Int_t j=0; j<fValues[AliReducedVarManager::kMCNJpsi]; j++)
        fHistosManager->FillHistClass(Form("pPb_5TeV_MC_cutMode_%d", cutMode+101), fValues);
    }

    // TODO
    /*// Fill histograms for tracks used for measured multiplicity estimation
    for(Int_t iArray=1; iArray<=2; iArray++) {
      AliReducedTrackInfo* track;
      TClonesArray* tracklist = (iArray==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
      TIter nextTrack(tracklist);
      // Loop over tracks
      for(Int_t it=0; it<tracklist->GetEntries(); ++it) {
        track = (AliReducedTrackInfo*)nextTrack();
        if(!(track->IsMCTruth()) && IsTrackMeasuredMultSelected(track,fValues)) {
          for(Int_t icut=0; icut<GetNMeasMultCuts(); icut++) {
            if(track->TestFlag(icut)) {
              fHistosManager->FillHistClass("TrackMult_AfterCuts", fValues);
            }
          }
        }
      }
    }*/
  }
  
  // apply event selection
  if(!isEventSelected) return;
  
  if(fOptionRunOverMC) {
    fSkipMCEvent = kFALSE;
    FillMCTruthHistograms();
    if(fSkipMCEvent) return;
  }
  
  // fill event info histograms after cuts
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
  }
  else {
    fHistosManager->FillHistClass("Event_noTag14_AfterCuts", fValues);
    for(Int_t icut=0; icut<nGlobalEstimators; icut++)
      fHistosManager->FillHistClass(Form("EventMult_noTag14_%s",GetMeasMultcutName(icut)), fValues);

    // Correlations between different multiplicity estimators
    fHistosManager->FillHistClass("CorrelMult_noTag_14", fValues);

    for(UShort_t ibit=0; ibit<64; ++ibit) {
      AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
      fHistosManager->FillHistClass("EventTag_noTag14_AfterCuts", fValues);
      AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
      fHistosManager->FillHistClass("EventTriggers_noTag14_AfterCuts", fValues);
    }
  }
  
  
  CreateFilteredEvent();
  if(fRejectEmptyEvents && (fFilteredEvent->NPairs()+fFilteredEvent->NTracks1()+fFilteredEvent->NTracks2())==0)
    return;
  fFilteredTree->Fill();
  
}


//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::CreateFilteredEvent() {
  //
  // create the filtered event
  //
  // NOTE: The following information is filtered from the input event
  //     1) Event header (either base event or full event header)
  //     2) selected V0 candidates -> see WriteFilteredPairs()
  //     3) selected tracks and legs of selected V0 candidates -> see WriteFilteredTracks()
  //     4) create candidates of the types specified in AliReducedPairInfo -> see BuildCandidatePairs()

  if(fFilteredTreeWritingOption==kFullEventsWithBaseTracks ||
     fFilteredTreeWritingOption==kFullEventsWithFullTracks)
    ((AliReducedEventInfo*)fFilteredEvent)->CopyEventHeader((AliReducedEventInfo*)fEvent);
  else
    fFilteredEvent->CopyEventHeader(fEvent);
  
  if(fWriteFilteredPairs)  WriteFilteredPairs();
  if(fBuildCandidatePairs) BuildCandidatePairs();
  if(fWriteFilteredTracks) {
    WriteFilteredTracks();
    if(!fWriteFilteredTracksCandidatesOnly) WriteFilteredTracks(2);
  }
  
  if(fComputeMult) {
    FillMultiplicity(kTRUE);
    for(Int_t icut=0; icut<GetNMeasMultCuts(); icut++) {
      
      // If event is unbiased
      if(fEvent->TestEventTag(14) || GetRunOverMC()) {
        // Regions relative to jpsi
        fHistosManager->FillHistClass(Form("EventMultRegions_%s", GetMeasMultcutName(icut)), fValues);
        // Regions relative to leading particle
        if(fValues[AliReducedVarManager::kPtLeading+icut] > fMinPtLeading)
          fHistosManager->FillHistClass(Form("EventMultRegions2Leading_%s", GetMeasMultcutName(icut)), fValues);
      }
      
      fFilteredEvent->SetNGlobalTracks(fValues[AliReducedVarManager::kNGlobalTracks+icut], icut);
      // Fill mult regions relative to jpsi
      fFilteredEvent->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksToward+icut],     0, kTRUE,
                                        icut);
      fFilteredEvent->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksTransverse+icut], 1, kTRUE,
                                        icut);
      fFilteredEvent->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksAway+icut],       2, kTRUE,
                                        icut);
      // Fill mult regions relative to leading pt
      fFilteredEvent->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksToward +
                                        AliReducedVarManager::kNMaxCutsGlobalTracks+icut], 0, kFALSE, icut);
      fFilteredEvent->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksTransverse +
                                        AliReducedVarManager::kNMaxCutsGlobalTracks+icut], 1, kFALSE, icut);
      fFilteredEvent->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksAway +
                                        AliReducedVarManager::kNMaxCutsGlobalTracks+icut], 2, kFALSE, icut);
      
      fFilteredEvent->SetLeadingParticle(fValues[AliReducedVarManager::kPtLeading +icut],
                                         fValues[AliReducedVarManager::kPhiLeading+icut],
                                         fValues[AliReducedVarManager::kEtaLeading+icut], icut);
    }  // end loop over multiplicity cuts
  }  // end if fComputeMult
}


//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::WriteFilteredPairs() {
  //
  // select and add filtered pair candidates to the filtered event
  //

  // loop over the pair list in the unfiltered event and evaluate all the pair cuts
  AliReducedPairInfo* pair = 0x0;
  TClonesArray* pairList = fEvent->GetPairs();
  if(!pairList) return;
  TIter nextPair(pairList);
  for(Int_t ip=0; ip<fEvent->NPairs(); ++ip) {
    pair = (AliReducedPairInfo*)nextPair();
    AliReducedVarManager::FillPairInfo(pair, fValues);
    fHistosManager->FillHistClass("Pair_BeforeCuts", fValues);
    for(UShort_t iflag=0; iflag<32; ++iflag) {
      AliReducedVarManager::FillPairQualityFlag(pair, iflag, fValues);
      fHistosManager->FillHistClass("PairQualityFlags_BeforeCuts", fValues);
    }

    TString pairTypeStr = "";
    if(pair->CandidateId() >= 0) {
      fCandidateType = pair->CandidateId();
      if      (pair->PairType()==0) pairTypeStr = "Offline";
      else if (pair->PairType()==1) pairTypeStr = "OnTheFly";
    }
    
    if(IsPairSelected(pair,fValues)) {
      for(Int_t icut=0; icut<fPairCuts.GetEntries(); ++icut) {
        if(pair->TestFlag(icut)) {
//           fHistosManager->FillHistClass(Form("Pair_%s", fPairCuts.At(icut)->GetName()), fValues);
          for(UShort_t iflag=0; iflag<32; ++iflag) {
            AliReducedVarManager::FillPairQualityFlag(pair, iflag, fValues);
            fHistosManager->FillHistClass(Form("PairQualityFlags_%s", fPairCuts.At(icut)->GetName()), fValues);
          }
          switch(fCandidateType) {
            case AliReducedPairInfo::kJpsiToEE :
              fHistosManager->FillHistClass(  Form("Pair_%s", fPairCuts.At(icut)->GetName()), fValues);
              break;
            case AliReducedPairInfo::kGammaConv :
              fHistosManager->FillHistClass(  Form("Pair_%s_%sGamma",
                                                   fPairCuts.At(icut)->GetName(),pairTypeStr.Data()), fValues);
              if(pair->IsPureV0Gamma())
                fHistosManager->FillHistClass(Form("Pair_%s_%sPureGamma",
                                                   fPairCuts.At(icut)->GetName(),pairTypeStr.Data()), fValues);
              break;
            case AliReducedPairInfo::kK0sToPiPi :
              fHistosManager->FillHistClass(  Form("Pair_%s_%sK0s",
                                                   fPairCuts.At(icut)->GetName(),pairTypeStr.Data()), fValues);
              if(pair->IsPureV0K0s())
                fHistosManager->FillHistClass(Form("Pair_%s_%sPureK0s",
                                                   fPairCuts.At(icut)->GetName(),pairTypeStr.Data()), fValues);
              break;
            case AliReducedPairInfo::kLambda0ToPPi :
              fHistosManager->FillHistClass(  Form("Pair_%s_%sLambda",
                                                   fPairCuts.At(icut)->GetName(),pairTypeStr.Data()), fValues);
              if(pair->IsPureV0Lambda())
                fHistosManager->FillHistClass(Form("Pair_%s_%sPureLambda",
                                                   fPairCuts.At(icut)->GetName(),pairTypeStr.Data()), fValues);
              break;
            case AliReducedPairInfo::kALambda0ToPPi :
              fHistosManager->FillHistClass(  Form("Pair_%s_%sALambda",
                                                   fPairCuts.At(icut)->GetName(),pairTypeStr.Data()), fValues);
              if(pair->IsPureV0ALambda())
                fHistosManager->FillHistClass(Form("Pair_%s_%sPureALambda",
                                                   fPairCuts.At(icut)->GetName(),pairTypeStr.Data()), fValues);
              break;
          };
        }
      }
      TClonesArray& pairs = *(fFilteredEvent->fCandidates);
      AliReducedPairInfo* filteredPair = NULL;
      filteredPair = new(pairs[fFilteredEvent->fNV0candidates[1]]) AliReducedPairInfo(*pair);
      fFilteredEvent->fNV0candidates[1] += 1;
    }
  }  // end loop over pairs
}


//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::WriteFilteredTracks(Int_t array /*=1*/) {
  //
  // select and add filtered tracks to the filtered event
  //
  
  // loop over the track list and evaluate all the track cuts
  AliReducedBaseTrack* track = 0x0;
  TClonesArray* trackList = (array==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
  if(!trackList) return;
  // Comment following line for comparison with AliReducedAnalysisJpsi2ee Track_BeforeCuts histogram
  if(fWriteFilteredTracksCandidatesOnly && (fFilteredEvent->NV0Candidates() == 0)) return;
  TIter nextTrack(trackList);
  for(Int_t it=0; it<trackList->GetEntries(); ++it) {
    track = (AliReducedBaseTrack*) nextTrack();
    AliReducedVarManager::FillTrackInfo(track, fValues);
    AliReducedVarManager::FillClusterMatchedTrackInfo(track, fValues);
    fHistosManager->FillHistClass("Track_WriteFilteredTracks_BeforeCuts", fValues);
    
    Bool_t writeTrack = IsTrackSelected(track, fValues);
    writeTrack |= (fWriteFilteredPairs && TrackIsCandidateLeg(track));
    
    if(fWriteFilteredTracksCandidatesOnly)
      // Comment following line for comparison with AliReducedAnalysisJpsi2ee Track_BeforePrefilter_standard
      //   histogram.
      writeTrack = TrackIsCandidateLeg(track) && !(track->IsMCTruth());
    
    if(writeTrack) {
      for(Int_t icut=0; icut<fTrackCuts.GetEntries(); ++icut) {
        if(track->TestFlag(icut)) {
          for(UInt_t iflag=0; iflag<64; ++iflag) {
            AliReducedVarManager::FillTrackQualityFlag(track, iflag, fValues);
            fHistosManager->FillHistClass(Form("TrackQualityFlags_WriteFilteredTracks_%s",
                                               fTrackCuts.At(icut)->GetName()), fValues);
          }
          if(!fFillTrackV0Histograms)
            fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%s",
                                               fTrackCuts.At(icut)->GetName()), fValues);
          else {
            if(track->IsGammaLeg())
              fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sGammaLeg",
                                                 fTrackCuts.At(icut)->GetName()), fValues);
            if(track->IsPureGammaLeg())
              fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sPureGammaLeg",
                                                 fTrackCuts.At(icut)->GetName()), fValues);
            if(track->IsK0sLeg())
              fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sK0sLeg",
                                                 fTrackCuts.At(icut)->GetName()), fValues);
            if(track->IsPureK0sLeg())
              fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sPureK0sLeg",
                                                 fTrackCuts.At(icut)->GetName()), fValues);
            if(track->IsLambdaLeg()) {
              if(track->Charge()>0)
                fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sLambdaPosLeg",
                                                   fTrackCuts.At(icut)->GetName()), fValues);
              else
                fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sLambdaNegLeg",
                                                   fTrackCuts.At(icut)->GetName()), fValues);
            }
            if(track->IsPureLambdaLeg()) {
              if(track->Charge()>0)
                fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sPureLambdaPosLeg",
                                                   fTrackCuts.At(icut)->GetName()), fValues);
              else
                fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sPureLambdaNegLeg",
                                                   fTrackCuts.At(icut)->GetName()), fValues);
            }
            if(track->IsALambdaLeg()) {
              if(track->Charge()>0)
                fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sALambdaPosLeg",
                                                   fTrackCuts.At(icut)->GetName()), fValues);
              else
                fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sALambdaNegLeg",
                                                   fTrackCuts.At(icut)->GetName()), fValues);
            }
            if(track->IsPureALambdaLeg()) {
              if(track->Charge()>0)
                fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sPureALambdaPosLeg",
                                                   fTrackCuts.At(icut)->GetName()), fValues);
              else
                fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sPureALambdaNegLeg",
                                                   fTrackCuts.At(icut)->GetName()), fValues);
            }
          }
        }
      }
      TClonesArray& tracks = (array==1 ? *(fFilteredEvent->fTracks) : *(fFilteredEvent->fTracks2));
      
      
      if(track->IsA() == AliReducedBaseTrack::Class()) {
        AliReducedBaseTrack* filteredParticle = NULL;
        filteredParticle = new(tracks[tracks.GetEntries()]) AliReducedBaseTrack(*track);
      }
      if(track->IsA() == AliReducedTrackInfo::Class()) {
        AliReducedTrackInfo* filteredParticle = NULL;
        AliReducedTrackInfo* tempTrack = dynamic_cast<AliReducedTrackInfo*>(track);
        //filteredParticle = (AliReducedTrackInfo*) tracks.ConstructedAt(tracks.GetEntries());
        new(tracks[tracks.GetEntries()]) AliReducedTrackInfo(*tempTrack);
        //filteredParticle = new(tracks[tracks.GetEntries()]) AliReducedTrackInfo(*tempTrack);
      }
      fFilteredEvent->fNtracks[1] += 1;
    }
  }  // end loop over tracks
}


//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::BuildCandidatePairs() {
  //
  // Build candidate pairs and add them to the filtered event
  //

  // clear the track arrays
  fLeg1Tracks.Clear("C"); fLeg1PrefilteredTracks.Clear("C");
  fLeg2Tracks.Clear("C"); fLeg2PrefilteredTracks.Clear("C");
  RunCandidateLegsSelection(1);
  //RunCandidateLegsSelection(2);
  
  if(fRunCandidatePrefilter) {
    RunCandidateLegsPrefilter(1);
    RunCandidateLegsPrefilter(2);
  }
  
  // Feed the selected tracks to the event mixing handler 
  if(fOptionRunMixing) {
    fMixingHandler->FillEvent(&fLeg1Tracks, &fLeg2Tracks, fValues, AliReducedPairInfo::kJpsiToEE);
    // Event mixing in multiplicity bins
    if(fOptionRunMixingMult) {
      AliMixingHandler* handler;
      TIter nextHandler(&fMixingHandlerMult); 
      for(Int_t i=0; i<fMixingHandlerMult.GetEntries(); i++) {
        handler = (AliMixingHandler*) nextHandler();
        if(fValues[AliReducedVarManager::kNGlobalTracks] >= fMultBinsMixing[i]  &&
           fValues[AliReducedVarManager::kNGlobalTracks] <  fMultBinsMixing[i+1]) {
          handler->FillEvent(&fLeg1Tracks, &fLeg2Tracks, fValues, AliReducedPairInfo::kJpsiToEE);
        }
      }
    }
  }
  
  if(fLeg1Tracks.GetEntries()+fLeg2Tracks.GetEntries() > 1)
    RunSameEventPairing();
}


//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::RunCandidateLegsSelection(Int_t arrayOption /*=1*/) {
  //
  // select leg candidates and prefilter tracks
  // NOTE: In the case of symmetric decay channels, the candidates are separated by charge in fLeg1Tracks
  //       and fLeg2Tracks. For asymmetric decays, the candidates for each of the LEG1 and LEG2 are not a
  //       priori separated by charge. It is the responsability of the analyzer to setup the track cuts such
  //       that these are separated.
  //

  Bool_t isAsymmetricDecayChannel = IsAsymmetricDecayChannel();
//   Bool_t mcDecision = kTRUE;
//   UInt_t mcDecisionMap = 1;
  
  // loop over the track list and evaluate all the track cuts
  AliReducedBaseTrack* track = 0x0;
  TClonesArray* trackList = (arrayOption==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
  if(!trackList) return;
  TIter nextTrack(trackList);
  for(Int_t it=0; it<trackList->GetEntries(); ++it) {
    track = (AliReducedTrackInfo*)nextTrack();
    if(fOptionRunOverMC && track->IsMCKineParticle()) continue;
    for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i)
      fValues[i] = -9999.;
    AliReducedVarManager::FillTrackInfo(track, fValues);
    AliReducedVarManager::FillClusterMatchedTrackInfo(track, fValues);
    fHistosManager->FillHistClass("Track_BeforeCuts", fValues);
    // NOTE: mcDecisionMap is not implemented properly for asymmetric decay channels.
//     if(fOptionRunOverMC && (fLegCandidatesMCcuts.GetEntries()>0))
//       mcDecisionMap = CheckReconstructedLegMCTruth(track);
//     if(fOptionRunOverMC && fLegCandidatesMCcuts) mcDecision = CheckReconstructedLegMCTruth(track);
    if(isAsymmetricDecayChannel) {
      if(IsCandidateLegSelected(track,fValues,1)) {
        fLeg1Tracks.Add(track);
        FillCandidateLegHistograms("Track_LEG1_BeforePrefilter", track, 1, isAsymmetricDecayChannel);
      }
      if(IsCandidateLegSelected(track,fValues,2)) {
        fLeg2Tracks.Add(track);
        FillCandidateLegHistograms("Track_LEG2_BeforePrefilter", track, 2, isAsymmetricDecayChannel);
      }
    }
    else {
      // mcDecisionMap is by default true, it can be false only if running on MC and the track fails the test
      if(IsCandidateLegSelected(track,fValues,1)) {// && mcDecisionMap) {
        if(track->Charge() > 0) {
          fLeg1Tracks.Add(track);
          FillCandidateLegHistograms("Track_LEG1_BeforePrefilter", track, 1, isAsymmetricDecayChannel);
        }
        else if(track->Charge() < 0) {
          fLeg2Tracks.Add(track);
          FillCandidateLegHistograms("Track_LEG2_BeforePrefilter", track, 2, isAsymmetricDecayChannel);
        }
      }
    }
    
    if(!fRunCandidatePrefilter) continue;
    
    if(isAsymmetricDecayChannel) {
      if(IsCandidateLegPrefilterSelected(track, fValues, 1)) {
        fLeg1PrefilteredTracks.Add(track);
        fHistosManager->FillHistClass("Track_LEG1_PrefilterTrack", fValues);
      }
      if(IsCandidateLegPrefilterSelected(track, fValues, 2)) {
        fLeg2PrefilteredTracks.Add(track);
        fHistosManager->FillHistClass("Track_LEG2_PrefilterTrack", fValues);
      }
    }
    else {
      if(IsCandidateLegPrefilterSelected(track, fValues)) {
        if(track->Charge() > 0) {
          fLeg1PrefilteredTracks.Add(track);
          fHistosManager->FillHistClass("Track_LEG1_PrefilterTrack", fValues);
        }
        else if(track->Charge() < 0) {
          fLeg2PrefilteredTracks.Add(track);
          fHistosManager->FillHistClass("Track_LEG2_PrefilterTrack", fValues);
        }
      }
    }
  }  // end loop over tracks
  
}


//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::RunCandidateLegsPrefilter(Int_t leg) {
  //
  // Run the prefilter selection
  // At this point it is assumed that the track lists are filled
  //

  Bool_t isAsymmetricDecayChannel = IsAsymmetricDecayChannel();
  
  // initialize iterators
  TIter iterLeg((leg==1 ? &fLeg1Tracks : &fLeg2Tracks));
  TIter iterPrefLeg1(&fLeg1PrefilteredTracks);
  TIter iterPrefLeg2(&fLeg2PrefilteredTracks);
  
  // Pair the LEG candidates with the prefilter selected tracks
  AliReducedBaseTrack* track     = 0;
  AliReducedBaseTrack* prefTrack = 0;
  for(Int_t it=0; it<(leg==1?fLeg1Tracks.GetEntries():fLeg2Tracks.GetEntries()); ++it) {
    track = (AliReducedBaseTrack*)iterLeg();
    
    if((leg==2 && !isAsymmetricDecayChannel) ||
       (leg==1 && (isAsymmetricDecayChannel ||
                  (!isAsymmetricDecayChannel&&fRunCandidatePrefilterOnSameCharge))))
    {
      iterPrefLeg1.Reset();
      for(Int_t it2=0; it2<fLeg1PrefilteredTracks.GetEntries(); ++it2) {
        prefTrack = (AliReducedBaseTrack*)iterPrefLeg1();
        
        if(track->TrackId()==prefTrack->TrackId()) continue;  // avoid self-pairing
        AliReducedVarManager::FillPairInfo(track, prefTrack, fCandidateType, fValues);
        if(!IsCandidateLegPairPrefilterSelected(fValues,1)) {
          track->ResetFlags();
          break;
        }
      }  // end loop over prefiltered selected leg1 candidates tracks
    }  // end if (if symmetric decay channel: selected leg2 candidates)
    else if((leg==1 && !isAsymmetricDecayChannel) ||
            (leg==2 && (isAsymmetricDecayChannel ||
                       (!isAsymmetricDecayChannel&&fRunCandidatePrefilterOnSameCharge))))
    {
      iterPrefLeg2.Reset();
      for(Int_t it2 = 0; it2<fLeg2PrefilteredTracks.GetEntries(); ++it2) {
        prefTrack = (AliReducedBaseTrack*)iterPrefLeg2();

        if(track->TrackId()==prefTrack->TrackId()) continue;  // avoid self-pairing
        AliReducedVarManager::FillPairInfo(track, prefTrack, fCandidateType, fValues);
        if(!IsCandidateLegPairPrefilterSelected(fValues,2)) {
          track->ResetFlags();
          break;
        }
      }  // end loop over prefiltered selected leg2 candidates tracks
    }  // end if (if symmetric decay channel: selected leg1 candidates)
  }  // end loop over selected leg1/2 candidates tracks

  // remove tracks
  iterLeg.Reset();
  for(Int_t it=(leg==1?fLeg1Tracks.GetEntries():fLeg2Tracks.GetEntries())-1; it>=0; --it) {
    track = (AliReducedBaseTrack*)iterLeg();
    if(!track->GetFlags()) {
      if(leg==1) fLeg1Tracks.Remove(track);
      if(leg==2) fLeg2Tracks.Remove(track);
    }
  }
  
  // fill histograms after the prefilter
  iterLeg.Reset();
  for(Int_t it=0; it<(leg==1?fLeg1Tracks.GetEntries():fLeg2Tracks.GetEntries()); ++it) {
    track = (AliReducedBaseTrack*)iterLeg();
    AliReducedVarManager::FillTrackInfo(track, fValues);
    AliReducedVarManager::FillClusterMatchedTrackInfo(track, fValues);
    FillCandidateLegHistograms(Form("Track_LEG%d_AfterPrefilter",leg), track,
                               (leg==2 && isAsymmetricDecayChannel ? 2 : 1), isAsymmetricDecayChannel);

    if(track->IsA() != AliReducedTrackInfo::Class()) continue;
    AliReducedTrackInfo* trackInfo = dynamic_cast<AliReducedTrackInfo*>(track);
    if(!trackInfo) continue;
    for(Int_t icut=0; icut<fLeg1Cuts.GetEntries(); ++icut) {
      if(track->TestFlag(leg==2 && isAsymmetricDecayChannel ? icut+32 : icut)) {
        TString legCutName = (leg==2 && isAsymmetricDecayChannel? fLeg2Cuts.At(icut)->GetName() :
                                                                  fLeg1Cuts.At(icut)->GetName());
        for(UInt_t iflag=0; iflag<AliReducedVarManager::kNTrackingFlags; ++iflag) {
          AliReducedVarManager::FillTrackingFlag(trackInfo, iflag, fValues);
          fHistosManager->FillHistClass(Form("TrackStatusFlags_%s",legCutName.Data()), fValues);
        }
        for(UInt_t iflag=0; iflag<64; ++iflag) {
          AliReducedVarManager::FillTrackQualityFlag(trackInfo, iflag, fValues);
          fHistosManager->FillHistClass(Form("TrackQualityFlags_%s",legCutName.Data()), fValues);
        }
        for(Int_t iLayer=0; iLayer<6; ++iLayer) {
          AliReducedVarManager::FillITSlayerFlag(trackInfo, iLayer, fValues);
          fHistosManager->FillHistClass(Form("TrackITSclusterMap_%s",legCutName.Data()), fValues);
          AliReducedVarManager::FillITSsharedLayerFlag(trackInfo, iLayer, fValues);
          fHistosManager->FillHistClass(Form("TrackITSsharedClusterMap_%s",legCutName.Data()), fValues);
        }
        for(Int_t iLayer=0; iLayer<8; ++iLayer) {
          AliReducedVarManager::FillTPCclusterBitFlag(trackInfo, iLayer, fValues);
          fHistosManager->FillHistClass(Form("TrackTPCclusterMap_%s",legCutName.Data()), fValues);
        }
      }
    }  // end loop over leg cuts
  }  // end loop over legs
}


//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::RunSameEventPairing() {
  //
  // Run the same event pairing
  //   
  
  Bool_t isAsymmetricDecayChannel = IsAsymmetricDecayChannel();
  TIter iterLeg1(&fLeg1Tracks);
  TIter iterLeg2(&fLeg2Tracks);
  
  AliReducedTrackInfo* leg1Track   = 0;
  AliReducedTrackInfo* leg2Track   = 0;
  AliReducedTrackInfo* leg1Track_2 = 0;
  
  for(Int_t it1=0; it1<fLeg1Tracks.GetEntries(); ++it1) {
    leg1Track = (AliReducedTrackInfo*)iterLeg1();

    iterLeg2.Reset();
    for(Int_t it2=0; it2<fLeg2Tracks.GetEntries(); ++it2) {
      leg2Track = (AliReducedTrackInfo*)iterLeg2();

      // verify that the two current tracks have at least 1 common bit
      ULong_t compatibilityMask = CheckTrackCompatibility(leg1Track, leg2Track, isAsymmetricDecayChannel);
      if(!compatibilityMask) continue;
      AliReducedVarManager::FillPairInfo(leg1Track, leg2Track, fCandidateType, fValues);
      
      if(fOptionRunOverMC) {
        Bool_t isJpsi = abs(leg1Track->MCPdg(0))==11 && abs(leg2Track->MCPdg(0))==11 &&
                        leg1Track->MCLabel(1)==leg2Track->MCLabel(1) && leg1Track->MCPdg(1)==443;
        Bool_t isJpsiFromB = isJpsi && ((abs(leg1Track->MCPdg(2))>500  && abs(leg1Track->MCPdg(2))<599) ||
                                        (abs(leg1Track->MCPdg(2))>5000 && abs(leg1Track->MCPdg(2))<5999));
        fValues[AliReducedVarManager::kPairMCMap] = isJpsi+2*isJpsiFromB;  // TODO implement with fMCMap|=(UShort_t(1)<<i
      }

      ULong_t pairCutMask = IsCandidatePairSelected(fValues);
      if(!pairCutMask) continue;
      FillCandidatePairHistograms(compatibilityMask, pairCutMask, 1, "Pair_Candidate", isAsymmetricDecayChannel,
                                  (fOptionRunOverMC?CheckReconstructedLegMCTruth(leg1Track,leg2Track):0));

      if(fOptionRunOverMC && (fLegCandidatesMCcuts.GetEntries()>0) &&
         (!CheckReconstructedLegMCTruth(leg1Track,leg2Track)))
        continue;

      TClonesArray& pairs = *(fFilteredEvent->fCandidates);
      AliReducedPairInfo* candidatePair = new(pairs[fFilteredEvent->fNV0candidates[1]]) AliReducedPairInfo();
      candidatePair->SetLegIds(leg1Track->TrackId(), leg2Track->TrackId());
      candidatePair->SetFlags(compatibilityMask);
      SetupPair(candidatePair, fValues);
      fFilteredEvent->fNV0candidates[1] += 1;
//       FillCandidatePairHistograms("Pair_Candidate12", candidatePair, fValues, isAsymmetricDecayChannel);
    }  // end loop over leg2 tracks

    if(fBuildCandidateLikePairs) {
      for(Int_t it1_2=it1+1; it1_2<fLeg1Tracks.GetEntries(); ++it1_2) {
        leg1Track_2 = (AliReducedTrackInfo*)fLeg1Tracks.At(it1_2);
        
        // verify that the two current tracks have at least 1 common bit
        ULong_t compatibilityMask = CheckTrackCompatibility(leg1Track, leg1Track_2, isAsymmetricDecayChannel);
        if(!compatibilityMask) continue;
        AliReducedVarManager::FillPairInfo(leg1Track, leg1Track_2, fCandidateType, fValues);
        fValues[AliReducedVarManager::kPairMCMap] = 0;
        
        ULong_t pairCutMask = IsCandidatePairSelected(fValues);
        if(!pairCutMask) continue;
        TClonesArray& pairs = *(fFilteredEvent->fCandidates);
        AliReducedPairInfo* candidatePair = new(pairs[fFilteredEvent->fNV0candidates[1]]) AliReducedPairInfo();
        candidatePair->SetLegIds(leg1Track->TrackId(), leg1Track_2->TrackId());
        candidatePair->SetFlags(compatibilityMask);
        SetupPair(candidatePair, fValues);
        fFilteredEvent->fNV0candidates[1] += 1;            
//         FillCandidatePairHistograms("Pair_Candidate11", candidatePair, fValues, isAsymmetricDecayChannel);
        FillCandidatePairHistograms(compatibilityMask, pairCutMask, 0, "Pair_Candidate",
                                    isAsymmetricDecayChannel,
                                    (fOptionRunOverMC?CheckReconstructedLegMCTruth(leg1Track,leg1Track_2):0));
      }  // end loop over leg1 tracks
    }  // end if fBuildCandidateLikePairs
  }  // end loop over leg1 tracks

  if(fBuildCandidateLikePairs) {
    AliReducedBaseTrack* leg2Track_2 = 0;
    iterLeg2.Reset();
    for(Int_t it2=0; it2<fLeg2Tracks.GetEntries(); ++it2) {
      leg2Track = (AliReducedTrackInfo*)iterLeg2();
      
      for(Int_t it2_2=it2+1; it2_2<fLeg2Tracks.GetEntries(); ++it2_2) {
        leg2Track_2 = (AliReducedTrackInfo*)fLeg2Tracks.At(it2_2);
        
        // verify that the two current tracks have at least 1 common bit
        ULong_t compatibilityMask = CheckTrackCompatibility(leg2Track, leg2Track_2, isAsymmetricDecayChannel);
        if(!compatibilityMask) continue;
        AliReducedVarManager::FillPairInfo(leg2Track, leg2Track_2, fCandidateType, fValues);
        fValues[AliReducedVarManager::kPairMCMap] = 0;

        ULong_t pairCutMask = IsCandidatePairSelected(fValues);
        if(!pairCutMask) continue;
        TClonesArray& pairs = *(fFilteredEvent->fCandidates);
        AliReducedPairInfo* candidatePair = new(pairs[fFilteredEvent->fNV0candidates[1]]) AliReducedPairInfo();
        candidatePair->SetLegIds(leg2Track->TrackId(), leg2Track_2->TrackId());
        candidatePair->SetFlags(compatibilityMask);
        SetupPair(candidatePair, fValues);
        fFilteredEvent->fNV0candidates[1] += 1;           
//         FillCandidatePairHistograms("Pair_Candidate22", candidatePair, fValues, isAsymmetricDecayChannel);
        FillCandidatePairHistograms(compatibilityMask, pairCutMask, 2, "Pair_Candidate",
                                    isAsymmetricDecayChannel,
                                    (fOptionRunOverMC?CheckReconstructedLegMCTruth(leg2Track,leg2Track_2):0));
      }  // end loop over leg2 tracks
    }  // end loop over leg2 tracks
  }  // end if fBuildCandidateLikePairs
}


//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::FillMultiplicity(Bool_t regions /*= kFALSE*/) {
  //
  // Fill multiplicity values (global if regions = false; in the regions if regions = true).
  //
  //
  
  // Fill global tracks (both signal and MC and MC truth number of Jpsi)
  if(!regions) {
    for(Int_t icut=0; icut<GetNMeasMultCuts(); icut++)
      fValues[AliReducedVarManager::kNGlobalTracks+icut] = 0.;
    fValues[AliReducedVarManager::kMCNch09] = 0.;
    fValues[AliReducedVarManager::kMCNJpsi] = 0.;
  }
  Float_t phi = 0;
  if(regions) {
    // Define the phi reference
    if(fDefaultRandomPhi) phi = TMath::TwoPi() * gRandom->Rndm();
    else                  phi = fValues[AliReducedVarManager::kPhiLeading];
    
    if(!fRegionsToMCTruth) {
      TClonesArray* pairs = fFilteredEvent->fCandidates;
      TIter nextPair(pairs);
      AliReducedPairInfo* jpsiPair;
      Float_t maxProbJpsi = 0.;
      
      // Loop over jpsi pair candidates
      for(Int_t i=0; i<pairs->GetEntries(); i++) {
        jpsiPair = (AliReducedPairInfo*) nextPair();
        if(jpsiPair->PairType() != 1) continue;  // reject like-sign pairs
        if(fJpsiMassDist) {
          // Take the jpsi with the maximal probability to be a real Jpsi (given its mass, if possible)
          Float_t probJpsi = fJpsiMassDist->Eval(jpsiPair->Mass(0));
          if(probJpsi >= maxProbJpsi) {
            maxProbJpsi = probJpsi;
            phi = jpsiPair->Phi();
          }
        }
        else phi = jpsiPair->Phi();
      }
      if((fValues[AliReducedVarManager::kMCNJpsi]==1) && (pairs->GetEntries()>0.)) {
        fValues[AliReducedVarManager::kPhiJpsiMCTruth] = abs(fValues[AliReducedVarManager::kPhiJpsiMCTruth]
                                                             - phi);
        fHistosManager->FillHistClass("DeltaPhi_JpsiTruth_JpsiCandidate", fValues);
      }
    }
    else if(fRegionsToMCTruth && (fValues[AliReducedVarManager::kMCNJpsi]==1))
      phi = fValues[AliReducedVarManager::kPhiJpsiMCTruth];
    
    for(Int_t icut=0; icut<GetNMeasMultCuts(); icut++) {
      fValues[AliReducedVarManager::kNGlobalTracksToward+icut]     = 0.;
      fValues[AliReducedVarManager::kNGlobalTracksTransverse+icut] = 0.;
      fValues[AliReducedVarManager::kNGlobalTracksAway+icut]       = 0.;
      fValues[AliReducedVarManager::kNGlobalTracksToward+
              AliReducedVarManager::kNMaxCutsGlobalTracks+icut] = 0.;
      fValues[AliReducedVarManager::kNGlobalTracksTransverse+
              AliReducedVarManager::kNMaxCutsGlobalTracks+icut] = 0.;
      fValues[AliReducedVarManager::kNGlobalTracksAway+
              AliReducedVarManager::kNMaxCutsGlobalTracks+icut] = 0.;
    }
    fValues[AliReducedVarManager::kMCNch09Toward]       = 0.;
    fValues[AliReducedVarManager::kMCNch09Transverse]   = 0.;
    fValues[AliReducedVarManager::kMCNch09Away]         = 0.;
    fValues[AliReducedVarManager::kMCNch09Toward+1]     = 0.;
    fValues[AliReducedVarManager::kMCNch09Transverse+1] = 0.;
    fValues[AliReducedVarManager::kMCNch09Away+1]       = 0.;
  }
  

  // Loop over both arrays
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

      // Get measured multiplicity (track has to be reconstructed & selected by at least one meas mult cut)
      if(!(track->IsMCTruth()) && IsTrackMeasuredMultSelected(track,fValues)) {
        for(Int_t icut=0; icut<GetNMeasMultCuts(); icut++) {
          if(track->TestFlag(icut)) {
            if(!regions) {
              fValues[AliReducedVarManager::kNGlobalTracks+icut] += 1;
              // Look for leading particle
              if((track->Pt()>fValues[AliReducedVarManager::kPtLeading+icut]) && !fRegionsToMCTruth) {
                fValues[AliReducedVarManager::kPtLeading +icut] = track->Pt();
                fValues[AliReducedVarManager::kPhiLeading+icut] = track->Phi();
                fValues[AliReducedVarManager::kEtaLeading+icut] = track->Eta();
              }
            }
            else {
              // Regions relative to jpsi
              float phiIcut   = fValues[AliReducedVarManager::kPhiLeading+icut];
              float delta_phi = abs(track->Phi()-phi);
              if(phi == fValues[AliReducedVarManager::kPhiLeading])
                delta_phi = abs(track->Phi()-phiIcut);
              if(delta_phi<M_PI/3. || delta_phi>5*M_PI/3.)
                fValues[AliReducedVarManager::kNGlobalTracksToward+icut] += 1.;
              else if((delta_phi>  M_PI/3. && delta_phi<2*M_PI/3.) ||
                      (delta_phi>4*M_PI/3. && delta_phi<5*M_PI/3.))
                fValues[AliReducedVarManager::kNGlobalTracksTransverse+icut] += 1.;
              else if(delta_phi>2*M_PI/3. && delta_phi<4*M_PI/3.)
                fValues[AliReducedVarManager::kNGlobalTracksAway+icut] += 1.;

              // Regions relative to leading particle
              float delta_phi_leading = abs(track->Phi()-phiIcut);
              if(delta_phi_leading<M_PI/3. || delta_phi_leading>5*M_PI/3.) 
                fValues[AliReducedVarManager::kNGlobalTracksToward+
                        AliReducedVarManager::kNMaxCutsGlobalTracks+icut] += 1.;
              else if((delta_phi_leading>  M_PI/3. && delta_phi_leading<2*M_PI/3.) ||
                      (delta_phi_leading>4*M_PI/3. && delta_phi_leading<5*M_PI/3.))
                fValues[AliReducedVarManager::kNGlobalTracksTransverse+
                        AliReducedVarManager::kNMaxCutsGlobalTracks+icut] += 1.;
              else if(delta_phi_leading>2*M_PI/3. && delta_phi_leading<4*M_PI/3.)
                fValues[AliReducedVarManager::kNGlobalTracksAway+
                        AliReducedVarManager::kNMaxCutsGlobalTracks+icut] += 1.;
            }
          }
        }
      }

      // Get true multiplicity (track has to be true MC (not rec) & selected by at least one true mult cut)
      if(track->IsMCTruth() && IsTrackTrueMultSelected(track, fValues)) {
        if(!regions) {
          fValues[AliReducedVarManager::kMCNch09] += 1;
          // Look for leading particle
          if((track->Pt()>fValues[AliReducedVarManager::kPtLeading]) && fRegionsToMCTruth) {
            for(Int_t icut=0; icut<GetNMeasMultCuts(); icut++) {
              fValues[AliReducedVarManager::kPtLeading +icut] = track->Pt();
              fValues[AliReducedVarManager::kPhiLeading+icut] = track->Phi();
              fValues[AliReducedVarManager::kEtaLeading+icut] = track->Eta();
            }
          }
        }
        else {
          // Regions relative to jpsi
          float delta_phi = abs(track->Phi()-phi);
          if(delta_phi<M_PI/3. || delta_phi>5*M_PI/3.) 
            fValues[AliReducedVarManager::kMCNch09Toward] += 1.;
          else if((delta_phi>  M_PI/3. && delta_phi<2*M_PI/3.) ||
                  (delta_phi>4*M_PI/3. && delta_phi<5*M_PI/3.))
            fValues[AliReducedVarManager::kMCNch09Transverse] += 1.;
          else if(delta_phi>2*M_PI/3. && delta_phi<4*M_PI/3.)
            fValues[AliReducedVarManager::kMCNch09Away] += 1.;

          // Regions relative to leading particle
          float delta_phi_leading = abs(track->Phi()-fValues[AliReducedVarManager::kPhiLeading]);
          if(delta_phi_leading<M_PI/3. || delta_phi_leading>5*M_PI/3.) 
            fValues[AliReducedVarManager::kMCNch09Toward+1] += 1.;
          else if((delta_phi_leading>  M_PI/3. && delta_phi_leading<2*M_PI/3.) ||
                  (delta_phi_leading>4*M_PI/3. && delta_phi_leading<5*M_PI/3.))
            fValues[AliReducedVarManager::kMCNch09Transverse+1] += 1.;
          else if(delta_phi_leading>2*M_PI/3. && delta_phi_leading<4*M_PI/3.)
            fValues[AliReducedVarManager::kMCNch09Away+1] += 1.;
        }
      }

      // Get MC truth number of jpsi and their azimuthal angle to construct the regions
      if(!regions && (iArray==1) && track->IsMCTruth()) {
        if((track->MCPdg(0)==443) && CheckMotherMCTruth(track)) {
          // find the jpsi daughters (to check if dielectron)
          Int_t daughter1Label = 0;
          Int_t daughter2Label = 0;
          FindJpsiTruthLegs(track, daughter1Label, daughter2Label);
          // Exclude if jpsi does not decay into dielectrons
          if(fMCTruthJpsi2eeOnly && ((daughter1Label==0) || (daughter2Label==0))) continue;
          fValues[AliReducedVarManager::kMCNJpsi] += 1;
          fValues[AliReducedVarManager::kPhiJpsiMCTruth] = track->Phi();
        }
      }
    }  // end loop over tracks
  }  // end loop over arrays

  // some quantities necessary to compute efficiency/contamination TODO is this implemented properly?
  if(!regions) {
    fValues[AliReducedVarManager::kMCNch09+1] = fValues[AliReducedVarManager::kMCNch09];
    fValues[AliReducedVarManager::kMCNch09+2] = fValues[AliReducedVarManager::kMCNch09];
    
    Bool_t isTriggered = kTRUE;
    if(!isTriggered) {
      fValues[AliReducedVarManager::kMCNch09+1] = -9999.;
      fValues[AliReducedVarManager::kMCNch09+2] = -9999.;
    }

    Bool_t isMCAccepted = (fValues[AliReducedVarManager::kMCNch]>0) &&
                          (abs(fValues[AliReducedVarManager::kVtxZMC])<=10.);
    if(!isMCAccepted) {
      fValues[AliReducedVarManager::kMCNch09]   = -9999.;
      fValues[AliReducedVarManager::kMCNch09+1] = -9999.;
      fValues[AliReducedVarManager::kMCNch09+2] = -9999.;
    }
  }
}


//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsEventSelected(AliReducedBaseEvent* event, Float_t* values/*=0x0*/) {
  //
  // apply event cuts
  //

  if(fEventCuts.GetEntries() == 0) return kTRUE;
  // loop over all the cuts and make a logical and between all cuts in the list
  for(Int_t i=0; i<fEventCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fEventCuts.At(i);
    if(values) { if(!cut->IsSelected(event,values)) return kFALSE; }
    else       { if(!cut->IsSelected(event))        return kFALSE; }
  }
  return kTRUE;
}


//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsTrackSelected(AliReducedBaseTrack* track, Float_t* values/*=0x0*/) {
  //
  // apply track cuts
  //

  if(fTrackCuts.GetEntries() == 0) return kTRUE;
  track->ResetFlags();
  
  // loop over all the cuts and toggle a filter bit if the track passes the corresponding cut
  for(Int_t i=0; i<fTrackCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fTrackCuts.At(i);
    if(values) { if(cut->IsSelected(track,values)) track->SetFlag(i); }
    else       { if(cut->IsSelected(track))        track->SetFlag(i); }
  }
  return (track->GetFlags()>0 ? kTRUE : kFALSE);
}


//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::TrackIsCandidateLeg(AliReducedBaseTrack* track) {
  //
  // if the track is a leg belonging to at least one of the candidate pairs written, then return kTRUE
  //

  AliReducedPairInfo* pair = 0x0;
  TClonesArray* pairList = fFilteredEvent->GetPairs();
  TIter nextPair(pairList);
  for(Int_t ip=0; ip<fFilteredEvent->NPairs(); ++ip) {
    pair = (AliReducedPairInfo*)nextPair();
    if(track->TrackId() == pair->LegId(0)) return kTRUE;
    if(track->TrackId() == pair->LegId(1)) return kTRUE;
  }
  return kFALSE;
}


//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsPairSelected(AliReducedPairInfo* pair, Float_t* values /*=0x0*/) {
  //
  // apply pair cuts
  // NOTE: Multiple cut sets are supported. The decisions are encoded in the fFlags inherited from
  //       AliReducedBaseTrack.

  if(fPairCuts.GetEntries()==0) return kTRUE;
  pair->ResetFlags();
  
  // Loop over all the cuts and toggle a filter bit if the pair passes the corresponding cut
  for(Int_t i=0; i<fPairCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fPairCuts.At(i);
    if(values) { if(cut->IsSelected(pair, values)) pair->SetFlag(i); }
    else       { if(cut->IsSelected(pair))         pair->SetFlag(i); }
  }
  return (pair->GetFlags()>0 ? kTRUE : kFALSE);
}


//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsTrackMeasuredMultSelected(AliReducedBaseTrack* track,
                                                                  Float_t* values /*=0x0*/) {
  //
  // apply cuts for determining measured multiplicity
  //
  
  if(fMeasuredMultTrackCuts.GetEntries() == 0) return kTRUE;
  track->ResetFlags();  
  
  for(Int_t i=0; i<fMeasuredMultTrackCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*) fMeasuredMultTrackCuts.At(i);
    if(values) { if(cut->IsSelected(track, values)) track->SetFlag(i); }
    else       { if(cut->IsSelected(track))         track->SetFlag(i); }
  }
  return (track->GetFlags()>0 ? kTRUE : kFALSE);
}


//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsTrackTrueMultSelected(AliReducedBaseTrack* track,
                                                              Float_t* values /*=0x0*/) {
  //
  // apply cuts for determining MC true multiplicity
  //

  if(fTrueMultTrackCuts.GetEntries()==0) return kTRUE;
  track->ResetFlags();  
  
  for(Int_t i=0; i<fTrueMultTrackCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*) fTrueMultTrackCuts.At(i) ;
    if(values) { if(cut->IsSelected(track, values)) track->SetFlag(i); }
    else       { if(cut->IsSelected(track))         track->SetFlag(i); }
  }
  return (track->GetFlags()>0 ? kTRUE : kFALSE);
}


//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsCandidateLegSelected(AliReducedBaseTrack* track,
                                                             Float_t* values /*=0x0*/, Int_t whichLeg /*=1*/) {
  //
  // apply cuts to the candidate leg
  //

  Bool_t isAsymmetricDecayChannel = IsAsymmetricDecayChannel();
  
  if(fLeg1Cuts.GetEntries() == 0) return kTRUE;
  
  // reset the flags for the track only if these are the cuts on LEG1
  // For LEG2, we use the same bit map as for LEG1, which it is assumed was already evaluated so the
  // ResetFlags() should not be called.
  // IMPORTANT: In the case of asymmetric decay channels, the LEG1 cuts have to be always evaluated before
  //            LEG2 cuts.
  if(whichLeg == 1) track->ResetFlags();
  
  for(Int_t i=0; i<fLeg1Cuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)(whichLeg==2 && isAsymmetricDecayChannel ? fLeg2Cuts.At(i)
                                                                                          : fLeg1Cuts.At(i));
    if(values) {
      if(cut->IsSelected(track,values)) track->SetFlag(whichLeg==2 && isAsymmetricDecayChannel ? 32+i : i);
    } else {
      if(cut->IsSelected(track))        track->SetFlag(whichLeg==2 && isAsymmetricDecayChannel ? 32+i : i);
    }
  }
  return (track->GetFlags()>0 ? kTRUE : kFALSE);
}


//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsCandidateLegPrefilterSelected(AliReducedBaseTrack* track,
    Float_t* values /*=0x0*/, Int_t whichLeg /*=1*/) {
  //
  // apply prefilter cuts on the candidate legs
  //

  if(fLeg1PrefilterCuts.GetEntries() == 0) return kTRUE;
  Bool_t isAsymmetricDecayChannel = IsAsymmetricDecayChannel();
  
  for(Int_t i=0; i<fLeg1PrefilterCuts.GetEntries(); ++i) {
    // if there are more cuts specified, we apply an AND on all of them
    // NOTE: The analysis task could also be configured such that there is a prefilter track cut
    //       corresponding to each track cut in which case one needs to make sure the number of prefilter
    //       cuts is the same as the number of track cuts.
    AliReducedInfoCut* cut = (AliReducedInfoCut*)((whichLeg==2 && isAsymmetricDecayChannel) ?
                                                  fLeg2PrefilterCuts.At(i) : fLeg1PrefilterCuts.At(i));
    if(values) { if(!cut->IsSelected(track, values)) return kFALSE; }
    else       { if(!cut->IsSelected(track))         return kFALSE; }
  }
  return kTRUE;
}


//___________________________________________________________________________
ULong_t AliReducedAnalysisFilterTrees::IsCandidatePairSelected(Float_t* values) {
  //
  // apply pair cuts
  //

  if(fCandidatePairCuts.GetEntries() == 0) return 1;

  ULong_t mask = 0;
  for(Int_t i=0; i<fCandidatePairCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fCandidatePairCuts.At(i);
    if(cut->IsSelected(values)) mask|=(ULong_t(1)<<i);
  }
  return mask;
}


// //___________________________________________________________________________
// Bool_t AliReducedAnalysisFilterTrees::IsCandidatePairSelected(Float_t* values) {
//   //
//   // apply pair cuts to the built candidates
//   //
//
//   if(fCandidatePairCuts.GetEntries()==0) return kTRUE;
//   // loop over all the cuts and make a logical and between all cuts in the list
//   for(Int_t i=0; i<fCandidatePairCuts.GetEntries(); ++i) {
//     AliReducedInfoCut* cut = (AliReducedInfoCut*)fCandidatePairCuts.At(i);
//     if(!cut->IsSelected(values)) return kFALSE;
//   }
//   return kTRUE;
// }


//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsCandidateLegPairPrefilterSelected(Float_t* values,
                                                                          Int_t whichLeg /*=1*/) {
  //
  // apply the prefilter pair cuts
  //

  if(fLeg1PairPrefilterCuts.GetEntries()==0) return kTRUE;
  Bool_t isAsymmetricDecayChannel = IsAsymmetricDecayChannel();
  
  // loop over all the cuts and make a logical OR between all cuts in the list
  for(Int_t i=0; i<fLeg1PairPrefilterCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)(whichLeg==2 && isAsymmetricDecayChannel?
                                                  fLeg2PairPrefilterCuts.At(i):fLeg1PairPrefilterCuts.At(i));
    if(cut->IsSelected(values)) return kTRUE;
  }
  return kFALSE;
}


//___________________________________________________________________________
ULong_t AliReducedAnalysisFilterTrees::CheckTrackCompatibility(AliReducedBaseTrack* leg1,
    AliReducedBaseTrack* leg2, Bool_t isAsymmetricDecayChannel) {
  //
  // check whether the 2 tracks fulfill at least one set of common cuts
  // NOTE: In the case of asymmetric decay channels, a track can fulfill simultaneously both the cuts of
  //       LEG1 and LEG2 and since we work with just one instance of the track object in memory, the cut
  //       flags map is shared between the two LEG cut sets: [0-31) cuts for LEG1, [32-63) cuts for LEG2.
  //       We should pair just tracks fulfilling cuts from the same "doublet": e.g. 0 - 32, 1-33, 2-34, ...
  //       In the case of symmetric decays, there is only one list of cuts for which the tracks are
  //       evaluated so here there is no problem. The flags of the tracks forming a pair are evaluated by a
  //       simple bitwise AND.
  if(!isAsymmetricDecayChannel) return (leg1->GetFlags() & leg2->GetFlags());
  
  ULong_t mask = 0;
  for(Int_t i=0; i<32; ++i) {
    if(leg1->TestFlag(i) && leg2->TestFlag(32+i))
      mask |= (ULong_t(1)<<i);
  }
  return mask;
}


//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::SetupPair(AliReducedPairInfo* pair, Float_t* values) {
  //
  // setup pair information
  //

  pair->Pt               (values[AliReducedVarManager::kPt]);
  pair->Phi              (values[AliReducedVarManager::kPhi]);
  pair->Eta              (values[AliReducedVarManager::kEta]);
  pair->Charge           (0);
  pair->CandidateId      (fCandidateType);
  pair->PairType         (values[AliReducedVarManager::kPairType]);
  pair->PairTypeSPD      (values[AliReducedVarManager::kPairTypeSPD]);
  pair->SetMass          (values[AliReducedVarManager::kMass]);
  pair->SetLxy           (values[AliReducedVarManager::kPairLxy]);
  pair->SetPseudoProper  (values[AliReducedVarManager::kPseudoProperDecayTime]);
  pair->SetPointingAngle (values[AliReducedVarManager::kPairPointingAngle]);
  pair->SetChisquare     (values[AliReducedVarManager::kPairChisquare]);
  pair->SetMCMap         (values[AliReducedVarManager::kPairMCMap]);
}


//___________________________________________________________________________
Bool_t AliReducedAnalysisFilterTrees::IsAsymmetricDecayChannel() {
  //
  // Check if the studied decay channel is assymetric
  // NOTE: If the decay channel studied is asymmetric, it is assumed that the cuts on the two legs are
  //       different and that the 2 lists of leg cuts hold the same number of cuts.
  if(fCandidateType==AliReducedPairInfo::kLambda0ToPPi)         return kTRUE;
  if(fCandidateType==AliReducedPairInfo::kALambda0ToPPi)        return kTRUE;
  if(fCandidateType==AliReducedPairInfo::kDplusToK0sPiplus)     return kTRUE;
  if(fCandidateType==AliReducedPairInfo::kDplusToK0sKplus)      return kTRUE;
  if(fCandidateType==AliReducedPairInfo::kDplusToPhiPiplus)     return kTRUE;
  if(fCandidateType==AliReducedPairInfo::kDminusToK0sPiminus)   return kTRUE;
  if(fCandidateType==AliReducedPairInfo::kDminusToK0sKminus)    return kTRUE;
  if(fCandidateType==AliReducedPairInfo::kDminusToPhiPiminus)   return kTRUE;
  if(fCandidateType==AliReducedPairInfo::kDzeroToKminusPiplus)  return kTRUE;
  if(fCandidateType==AliReducedPairInfo::kADzeroToKplusPiminus) return kTRUE;
  if(fCandidateType==AliReducedPairInfo::kDsplusToK0sKplus)     return kTRUE;
  if(fCandidateType==AliReducedPairInfo::kDsminusToK0sKminus)   return kTRUE;
  return kFALSE;
}


//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::FillCandidateLegHistograms(TString histClass, AliReducedBaseTrack* track,
                                                               Int_t leg, Bool_t isAsymmetricDecayChannel) {
  //
  // fill track histogram lists according to the track flags 
  //

  for(Int_t icut=0; icut<fLeg1Cuts.GetEntries(); ++icut) {
    if(track->TestFlag(leg==2 && isAsymmetricDecayChannel ? icut+32 : icut)) {
      fHistosManager->FillHistClass(Form("%s_%s",histClass.Data(),
          (leg==2 && isAsymmetricDecayChannel?fLeg2Cuts.At(icut)->GetName():fLeg1Cuts.At(icut)->GetName())),
          fValues);
    }
  }
}


/*//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::FillCandidatePairHistograms(TString histClass, AliReducedPairInfo* pair, Bool_t isAsymmetricDecayChannel) {
  //
  // fill track histogram lists according to the track flags 
  //

  for(Int_t icut=0; icut<fLeg1Cuts.GetEntries(); ++icut) {
    if(pair->TestFlag(icut)) {
      fHistosManager->FillHistClass(Form("%s_%s%s", histClass.Data(), fLeg1Cuts.At(icut)->GetName(),
          (isAsymmetricDecayChannel?Form("_%s",fLeg2Cuts.At(icut)->GetName()):"")), fValues);
    }
  }
}*/

//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::FillCandidatePairHistograms(ULong_t trackMask, ULong_t pairMask,
    Int_t pairType, TString pairClass /*="Pair_Candidate"*/, Bool_t isAsymmetricDecayChannel,
    UInt_t mcDecisions /* = 0*/) {

  //
  // Fill track histogram lists according to the track flags
  // PairType can be 0,1 or 2 corresponding to leg1&leg1, leg1&leg2 or leg2&leg2 pairs
  //

  TString typeStr[3] = {"11", "12", "22"};
  if(fPairCuts.GetEntries() > 1) {
    // loop over leg cuts (for asymmetric decay channels there should always be a pair of cuts for leg1 and
    // leg2)
    for(Int_t iTrackCut=0; iTrackCut<fLeg1Cuts.GetEntries(); ++iTrackCut) {
      // check if the legs of the pair fulfill the respective track cuts
      if(trackMask & (ULong_t(1)<<iTrackCut)) {
        // loop over pair cuts
        for(Int_t iPairCut=0; iPairCut<fPairCuts.GetEntries(); ++iPairCut) {
          // check if the pair fulfills the pair cut
          if(pairMask & (ULong_t(1)<<iPairCut)) {
            fHistosManager->FillHistClass(Form("%s%s_%s%s_%s", pairClass.Data(), typeStr[pairType].Data(),
                fLeg1Cuts.At(iTrackCut)->GetName(),
                (isAsymmetricDecayChannel?Form("_%s",fLeg2Cuts.At(iTrackCut)->GetName()):""),
                fPairCuts.At(iPairCut)->GetName()), fValues);
            if(mcDecisions && pairType==1) {
              for(Int_t iMC=0; iMC<=fLegCandidatesMCcuts.GetEntries(); ++iMC) {
                if(mcDecisions & (UInt_t(1)<<iMC))
                  fHistosManager->FillHistClass(Form("%s%s_%s%s_%s_%s", pairClass.Data(),
                      typeStr[pairType].Data(), fLeg1Cuts.At(iTrackCut)->GetName(),
                      (isAsymmetricDecayChannel?Form("_%s",fLeg2Cuts.At(iTrackCut)->GetName()):""),
                      fPairCuts.At(iPairCut)->GetName(), fLegCandidatesMCcuts.At(iMC)->GetName()), fValues);
              }
            }
          }
        }  // end loop over pair cuts
      }
    }  // end loop over leg cuts
  }
  else {
    for(Int_t iTrackCut=0; iTrackCut<fLeg1Cuts.GetEntries(); ++iTrackCut) {
      if(trackMask & (ULong_t(1)<<iTrackCut)) {
        fHistosManager->FillHistClass(Form("%s%s_%s%s", pairClass.Data(), typeStr[pairType].Data(),
            fLeg1Cuts.At(iTrackCut)->GetName(),
            (isAsymmetricDecayChannel?Form("_%s",fLeg2Cuts.At(iTrackCut)->GetName()):"")), fValues);
        if(mcDecisions && pairType==1) {
          for(Int_t iMC=0; iMC<=fLegCandidatesMCcuts.GetEntries(); ++iMC) {
            if(mcDecisions & (UInt_t(1)<<iMC))
              fHistosManager->FillHistClass(Form("%s%s_%s%s_%s", pairClass.Data(), typeStr[pairType].Data(),
                  fLeg1Cuts.At(iTrackCut)->GetName(),
                  (isAsymmetricDecayChannel?Form("_%s",fLeg2Cuts.At(iTrackCut)->GetName()):""),
                  fLegCandidatesMCcuts.At(iMC)->GetName()), fValues);
          }
        }
      }
    }  // end loop over leg cuts
  }
}


//___________________________________________________________________________
const Char_t* AliReducedAnalysisFilterTrees::GetCandidateLegCutName(Int_t i, Int_t leg) {
  //
  // get candidate leg cut name
  //
  if(leg==2 && IsAsymmetricDecayChannel())
    return (i<fLeg2Cuts.GetEntries() ? fLeg2Cuts.At(i)->GetName() : "");
  return (i<fLeg1Cuts.GetEntries() ? fLeg1Cuts.At(i)->GetName() : "");
}


//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::Finish() {
  //
  // run stuff after the event loop
  //
  
  if(fOptionRunMixing && !fOptionRunOverMC) {
    fMixingHandler->RunLeftoverMixing(AliReducedPairInfo::kJpsiToEE);
    if(fOptionRunMixingMult) {  // Mixed event in multiplicity bins
      AliMixingHandler* handler;
      TIter nextHandler(&fMixingHandlerMult); 
      for(Int_t i=0; i<fMixingHandlerMult.GetEntries(); i++) {
        handler = (AliMixingHandler*) nextHandler();
        handler->RunLeftoverMixing(AliReducedPairInfo::kJpsiToEE);
      }
    }
  }
  
}


//___________________________________________________________________________
UInt_t AliReducedAnalysisFilterTrees::CheckReconstructedLegMCTruth(AliReducedBaseTrack* track) {
// Bool_t AliReducedAnalysisFilterTrees::CheckReconstructedLegMCTruth(AliReducedBaseTrack* track) {
  //
  // Check a reconstructed track against all the specified MC truth cuts
  //
  // TODO: In the fLegCandidatesMCcuts one can also add AliSignalMC objects which can then be tested
  //       using the AliReducedTrackInfo::fMCPdg[]
  //

  if(fLegCandidatesMCcuts.GetEntries()==0) return 0;

  UInt_t decisionMap = 0;
  for(Int_t i=0; i<fLegCandidatesMCcuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fLegCandidatesMCcuts.At(i);
    if(cut->IsSelected(track)) decisionMap |= (UInt_t(1)<<i);
  }

  return decisionMap;

//   if(!fLegCandidatesMCcuts) return 0;
//
//   Bool_t decision = fLegCandidatesMCcuts->IsSelected(track);
//   //for(Int_t i=0; i<fLegCandidatesMCcuts.GetEntries(); ++i) {
//   //AliReducedInfoCut* cut = (AliReducedInfoCut*)fLegCandidatesMCcuts.At(i);
//   // if(->IsSelected(track))
//   //   decisionMap |= (UInt_t(1)<<i);
//   //}
//
//   return decision;
}


//___________________________________________________________________________
UInt_t AliReducedAnalysisFilterTrees::CheckReconstructedLegMCTruth(AliReducedBaseTrack* ptrack,
                                                                   AliReducedBaseTrack* ntrack) {
// Bool_t AliReducedAnalysisFilterTrees::CheckReconstructedLegMCTruth(AliReducedBaseTrack* ptrack, AliReducedBaseTrack* ntrack) {
  //
  // check the pair of tracks to see if they match the defined MC cuts and in addition
  // that they have the same mother
  // NOTE: The condition for the 2 tracks to have the same mother requires information on the MC label,
  //       which is available just in the full track information (AliReducedTrackInfo::fMCLabels[]).
  //       The consequence is that for the jpsi2ee analysis, the reconstructed tracks need to be always
  //       written as full tracks.
  //

  // check that both tracks are full tracks
  if(ptrack->IsA() != AliReducedTrackInfo::Class()) return 0;
  if(ntrack->IsA() != AliReducedTrackInfo::Class()) return 0;

  // check the MC requirements on each of the leg and their logical intersection
  if(fLegCandidatesMCcuts.GetEntries()==0) return 0;
  UInt_t pTrackDecisions = CheckReconstructedLegMCTruth(ptrack);
  if(!pTrackDecisions) return 0;
  UInt_t nTrackDecisions = CheckReconstructedLegMCTruth(ntrack);

  UInt_t decisions = 0;
  for(Int_t i=0; i<fLegCandidatesMCcuts.GetEntries(); ++i) {
    Bool_t pDecision = (pTrackDecisions & (UInt_t(1)<<i));
    Bool_t nDecision = (nTrackDecisions & (UInt_t(1)<<i));

    // if needed, check that the tracks have the same mother
    Bool_t sameMotherDecision = kTRUE;
    if(fLegCandidatesMCcuts_RequestSameMother[i])
      sameMotherDecision = (TMath::Abs(((AliReducedTrackInfo*)ptrack)->MCLabel(1)) ==
                            TMath::Abs(((AliReducedTrackInfo*)ntrack)->MCLabel(1)));
    if(sameMotherDecision && pDecision && nDecision) decisions |= (UInt_t(1)<<i);
  }

  return decisions;

  /*// check that both tracks are full tracks
  if(ptrack->IsA() != AliReducedTrackInfo::Class()) return 0;
  if(ntrack->IsA() != AliReducedTrackInfo::Class()) return 0;
  
  // check that the tracks have the same mother
  if(TMath::Abs(((AliReducedTrackInfo*)ptrack)->MCLabel(1)) !=
     TMath::Abs(((AliReducedTrackInfo*)ntrack)->MCLabel(1))) return 0;
  
  // check the MC requirements on each of the leg and their logical intersection
  if(!fLegCandidatesMCcuts) return 0;
  UInt_t pTrackDecisions = CheckReconstructedLegMCTruth(ptrack);
  if(!pTrackDecisions) return 0;
  UInt_t nTrackDecisions = CheckReconstructedLegMCTruth(ntrack);
  return (pTrackDecisions & nTrackDecisions);*/
}


//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::FillMCTruthHistograms() {
  //
  // fill histograms with pure signal
  //   
  // loop over the first track array
  LoopOverMCTracks(1);
  // and over the second
  // NOTE: In the current model, handling the MC truth info requires the labels, which are properties of the
  //       full track, so there is no point in looping over the second track array which, if it exists,
  //       contains just base tracks.
  //LoopOverMCTracks(2);
}


//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::LoopOverMCTracks(Int_t trackArray /*=1*/) {
  //
  // loop over the track array and check the pure MC tracks against the defined MC selections
  //   
  AliReducedTrackInfo* mother    = 0x0;
  AliReducedTrackInfo* daughter1 = 0x0;
  AliReducedTrackInfo* daughter2 = 0x0;
  
  TClonesArray* trackList = (trackArray==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
  if(!trackList) return;
  TIter nextTrack(trackList);

  // if the pt dependent weights were set, check the weight and reject randomly the event
  if(fMCJpsiPtWeights) {
    for(Int_t it=0; it<trackList->GetEntries(); ++it) {
      mother = (AliReducedTrackInfo*)nextTrack();
      if(!mother->IsMCKineParticle()) continue;
      if(!(mother->MCPdg(0)==443) || mother->MCPdg(1)==443) continue;
      
      // apply selections on the jpsi mother
      UInt_t motherDecisions = CheckMotherMCTruth(mother);
      if(!motherDecisions) continue;
      // Apply only selections to one of the MCJpsi cuts. TODO Why needed?
      if((fReweightCut>=0) && (fReweightCut<GetNJpsiMotherMCCuts()) &&
         !(motherDecisions&UInt_t(1)<<fReweightCut)) continue;
      
      Double_t pt = mother->Pt();
      if(pt>fMCJpsiPtWeights->GetXaxis()->GetXmax()) pt = fMCJpsiPtWeights->GetXaxis()->GetXmax();
      Double_t weight = fMCJpsiPtWeights->GetBinContent(fMCJpsiPtWeights->FindBin(pt));
      if(weight>1.0) weight = 1.0;
      Double_t rnd = gRandom->Rndm();
      if(weight<rnd) {
        fSkipMCEvent = kTRUE; 
        return;
      }
    }
  }
  
  nextTrack.Reset();
  for(Int_t it=0; it<trackList->GetEntries(); ++it) {
    mother = (AliReducedTrackInfo*)nextTrack();
    if(!mother->IsMCKineParticle()) continue;
    // Next line reduces a lot running time when all primaries are stored
    if(!(mother->MCPdg(0)==443) || mother->MCPdg(1)==443) continue;
    
    // apply selections on the jpsi mother
    UInt_t motherDecisions = CheckMotherMCTruth(mother);
    if(!motherDecisions) continue;
    
    // find the jpsi daughters (needed to compute 2-track properties like the polarization, etc.)
    Int_t daughter1Label = 0;
    Int_t daughter2Label = 0;
    FindJpsiTruthLegs(mother, daughter1Label, daughter2Label);
    daughter1 = FindTrackByLabel(daughter1Label, kTRUE);
    daughter2 = FindTrackByLabel(daughter2Label, kTRUE);
    
    
    // reset track variables and fill info
    for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i)
      fValues[i] = -9999.;
    AliReducedVarManager::FillMCTruthInfo(mother, fValues, daughter1, daughter2);

    // loop over jpsi mother selections and fill histograms before the kine cuts on electrons
    for(Int_t iCut=0; iCut<fJpsiMotherMCcuts.GetEntries(); ++iCut) {
      if(!(motherDecisions & (UInt_t(1)<<iCut))) continue;
      fHistosManager->FillHistClass(Form("PureMCTRUTH_BeforeSelection_%s",
                                         fJpsiMotherMCcuts.At(iCut)->GetName()), fValues);
    }

    if(!daughter1) continue;
    if(!daughter2) continue;
    
    // apply selections on pure MC daughter electrons (kine cuts)
    UInt_t daughter1Decisions = CheckDaughterMCTruth(daughter1);
    if(!daughter1Decisions) continue;
    UInt_t daughtersDecisions = daughter1Decisions & CheckDaughterMCTruth(daughter2);
    if(!daughtersDecisions) continue;
    for(Int_t iCut=0; iCut<fJpsiMotherMCcuts.GetEntries(); ++iCut) {
      // Fill histogram for MCTruth jpsi after selection on daughters
      // reset track variables and fill info. Necessary because of FillPairInfo of detected daughters in loop
      for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i)
        fValues[i] = -9999.;
      AliReducedVarManager::FillMCTruthInfo(mother, fValues, daughter1, daughter2);
      if(!(motherDecisions    & (UInt_t(1)<<iCut))) continue;
      if(!(daughtersDecisions & (UInt_t(1)<<iCut))) continue;
      fHistosManager->FillHistClass(Form("PureMCTRUTH_AfterSelection_%s",
                                         fJpsiMotherMCcuts.At(iCut)->GetName()), fValues);
      for(Int_t jcut=0; jcut<GetNMeasMultCuts(); jcut++)
        fHistosManager->FillHistClass(Form("JpsiPtMultCorrel_%s_%s", GetMeasMultcutName(jcut),
                                           fJpsiMotherMCcuts.At(iCut)->GetName()), fValues);

      // Fill histograms for reconstructed MCTruth jpsi daughters
      AliReducedTrackInfo* daughter1Det = FindTrackByLabel(daughter1Label, kFALSE);
      if(!daughter1Det) continue;
      if(!IsCandidateLegSelected(daughter1Det)) continue;        
      AliReducedTrackInfo* daughter2Det = FindTrackByLabel(daughter2Label, kFALSE);
      if(!daughter2Det) continue;
      if(!IsCandidateLegSelected(daughter2Det)) continue;

      AliReducedVarManager::FillPairInfo(daughter1Det, daughter2Det, fCandidateType, fValues);
      fHistosManager->FillHistClass(Form("PureMCTRUTH_DetectedDaughters_%s",
                                         fJpsiMotherMCcuts.At(iCut)->GetName()), fValues);
      for(Int_t jcut=0; jcut<GetNMeasMultCuts(); jcut++)
        fHistosManager->FillHistClass(Form("JpsiPtMultCorrelMeas_%s_%s", GetMeasMultcutName(jcut),
                                           fJpsiMotherMCcuts.At(iCut)->GetName()), fValues);
    }  // end loop over fJpsiMotherMCcuts
  }  // end loop over tracks

  return;
}


//___________________________________________________________________________
UInt_t AliReducedAnalysisFilterTrees::CheckMotherMCTruth(AliReducedTrackInfo* mother) {
  //
  // Check the mother pure MC truth against all defined selections and return a bit map with all decisions
  //
  if(fJpsiMotherMCcuts.GetEntries()==0) return 0;
  
  UInt_t decisionMap = 0;
  for(Int_t i=0; i<fJpsiMotherMCcuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fJpsiMotherMCcuts.At(i);
    if(cut->IsSelected(mother)) decisionMap |= (UInt_t(1)<<i);
  }
  
  return decisionMap;
}


//___________________________________________________________________________
UInt_t AliReducedAnalysisFilterTrees::CheckDaughterMCTruth(AliReducedTrackInfo* daughter) {
  //
  // Check the daughter pure MC truth against all defined selections and return a bit map with all decisions
  //
  if(fJpsiElectronMCcuts.GetEntries()==0) return 0;
  
  UInt_t decisionMap = 0;
  for(Int_t i=0; i<fJpsiElectronMCcuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fJpsiElectronMCcuts.At(i);
    if(cut->IsSelected(daughter)) decisionMap |= (UInt_t(1)<<i);
  }
  
  return decisionMap;
}


//___________________________________________________________________________
AliReducedTrackInfo* AliReducedAnalysisFilterTrees::FindTrackByLabel(Int_t label, Bool_t isTruth) {
  //
  // search the track list for pure MC track with label and return the track pointer
  //
  AliReducedTrackInfo* track = 0x0;
  TClonesArray* trackList = fEvent->GetTracks();
  TIter nextTrack(trackList);
  for(Int_t i=0; i<trackList->GetEntries(); ++i) {
    track = (AliReducedTrackInfo*)nextTrack();
    if(isTruth  && (!track->IsMCKineParticle())) continue;
    if(!isTruth &&   track->IsMCKineParticle())    continue;
    if(track->MCLabel(0)==label) return track;
  }
  return 0x0;
}


//___________________________________________________________________________
void AliReducedAnalysisFilterTrees::FindJpsiTruthLegs(AliReducedTrackInfo* mother, Int_t& leg1Label,
                                                      Int_t& leg2Label) {
  //
  // find the jpsi legs in the list of pure MC truth particles
  //
  Int_t mLabel = mother->MCLabel(0);
  AliReducedTrackInfo* track = 0x0;
  Int_t legsFound = 0;
  
  // loop over the first track array
  TClonesArray* trackList = fEvent->GetTracks();
  TIter nextTrack(trackList);
  for(Int_t i=0; i<trackList->GetEntries(); ++i) {
    if(legsFound==2) return;
    track = (AliReducedTrackInfo*)nextTrack();
    if(!track->IsMCKineParticle()) continue;
    if(track->MCLabel(1)==mLabel && TMath::Abs(track->MCPdg(0))==11) {
      legsFound += 1;
      if(legsFound==1) leg1Label = track->MCLabel(0);
      if(legsFound==2) leg2Label = track->MCLabel(0);
    }
  }
  return;
}

