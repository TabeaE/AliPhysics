//
// Creation date: 2017/08/09
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#include "AliReducedAnalysisFilterTrees.h"

#include <iostream>
using std::cout;
using std::endl;

#include <TClonesArray.h>
#include <TIterator.h>
#include <TRandom3.h>
#include <TF1.h>

#include "AliHistogramManager.h"
#include "AliReducedVarManager.h"
#include "AliReducedEventInfo.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedPairInfo.h"

ClassImp(AliReducedAnalysisFilterTrees);


//_____________________________________________________________________________
AliReducedAnalysisFilterTrees::AliReducedAnalysisFilterTrees() :
AliReducedAnalysisTaskSE(),
fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
fMixingHandler(new AliMixingHandler("J/psi signal extraction", "", 
                   AliMixingHandler::kMixResonanceLegs)),
fMixingHandlerMult                  (),
fMultBinsMixing                     (),
fNMultBinsMixing                    (0),
fEventCuts                          (),
fTrackCuts                          (),
fWriteFilteredTracks                (kTRUE),
fWriteFilteredTracksCandidatesOnly  (kFALSE),
fFillTrackV0Histograms              (kFALSE),
fOptionRunMixing                    (kTRUE),
fOptionRunMixingMult                (kFALSE),
fComputeMult                        (kTRUE),
fPairCuts                           (),
fWriteFilteredPairs                 (kTRUE),
fRejectEmptyEvents                  (kFALSE),
fBuildCandidatePairs                (kFALSE),
fBuildCandidateLikePairs            (kFALSE),
fMCTruthJpsi2eeOnly                 (kFALSE),
fRegionsToMCTruth                   (kFALSE),
fDefaultRandomPhi                   (kTRUE),
fMinPtLeading                       (0.0),
fReweightCut                        (-1),
fCandidateType                      (AliReducedPairInfo::kJpsiToEE),
fLeg1Cuts                           (),
fLeg2Cuts                           (),
fCandidatePairCuts                  (),
fRunCandidatePrefilter              (kFALSE),
fRunCandidatePrefilterOnSameCharge  (kFALSE),
fJpsiMassDist                       (),
fMeasMultTrackCuts                  (),
fWeightsTrackCuts                   (),
fTrueMultTrackCuts                  (),
fLeg1PrefilterCuts                  (),
fLeg2PrefilterCuts                  (),
fLeg1PairPrefilterCuts              (),
fLeg2PairPrefilterCuts              (),
fLeg1Tracks                         (),
fLeg2Tracks                         (),
fLeg1PrefilteredTracks              (),
fLeg2PrefilteredTracks              (),
fOptionRunOverMC                    (kFALSE),
fLegCandidatesMCcuts                (),
fLegCandidatesMCcuts_RequestSameMother (),
fJpsiMotherMCcuts                   (),
fMCJpsiPtWeights                    (0x0),
fSkipMCEvent                        (kFALSE),
fJpsiElectronMCcuts                 (),
fReweightPC                         (kFALSE),
fMCPCWeights                        (0x0)
{
  //
  // default constructor
  //
}


//_____________________________________________________________________________
AliReducedAnalysisFilterTrees::AliReducedAnalysisFilterTrees(const Char_t* name, const Char_t* title) :
AliReducedAnalysisTaskSE(name, title),
fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
fMixingHandler(new AliMixingHandler("J/psi signal extraction", "", 
                   AliMixingHandler::kMixResonanceLegs)),
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
fMeasMultTrackCuts(),
fWeightsTrackCuts(),
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
fJpsiElectronMCcuts(),
fReweightPC(kFALSE),
fMCPCWeights(0x0)
{
  //
  // named constructor
  //
  fMixingHandlerMult     .SetOwner(kTRUE);
  fEventCuts             .SetOwner(kTRUE);
  fTrackCuts             .SetOwner(kTRUE);
  fPairCuts              .SetOwner(kTRUE);
  fLeg1Cuts              .SetOwner(kTRUE);
  fLeg2Cuts              .SetOwner(kTRUE);
  fCandidatePairCuts     .SetOwner(kTRUE);
  fMeasMultTrackCuts     .SetOwner(kTRUE);
  fWeightsTrackCuts      .SetOwner(kTRUE);
  fTrueMultTrackCuts     .SetOwner(kTRUE);
  fLeg1PrefilterCuts     .SetOwner(kTRUE);
  fLeg2PrefilterCuts     .SetOwner(kTRUE);
  fLeg1PairPrefilterCuts .SetOwner(kTRUE);
  fLeg2PairPrefilterCuts .SetOwner(kTRUE);
  fLeg1Tracks            .SetOwner(kFALSE);
  fLeg2Tracks            .SetOwner(kFALSE);
  fLeg1PrefilteredTracks .SetOwner(kFALSE);
  fLeg2PrefilteredTracks .SetOwner(kFALSE);
  fLegCandidatesMCcuts   .SetOwner(kTRUE);
  fJpsiMotherMCcuts      .SetOwner(kTRUE);
  fJpsiElectronMCcuts    .SetOwner(kTRUE);
  for(Int_t i=0; i<32; ++i) fLegCandidatesMCcuts_RequestSameMother[i] = kTRUE;
  fRand.reset(new TRandom3());
}


//_____________________________________________________________________________
// Destructor
AliReducedAnalysisFilterTrees::~AliReducedAnalysisFilterTrees() 
{
  fEventCuts         .Clear("C"); 
  fTrackCuts         .Clear("C");  fPairCuts              .Clear("C");
  fLeg1Cuts          .Clear("C");  fLeg2Cuts              .Clear("C");
  fCandidatePairCuts .Clear("C");  fMeasMultTrackCuts     .Clear("C");
  fWeightsTrackCuts  .Clear("C");  fTrueMultTrackCuts     .Clear("C");
  fLeg1PrefilterCuts .Clear("C");  fLeg1PairPrefilterCuts .Clear("C");
  fLeg2PrefilterCuts .Clear("C");  fLeg2PairPrefilterCuts .Clear("C");
  fLeg1Tracks        .Clear("C");  fLeg1PrefilteredTracks .Clear("C");
  fLeg2Tracks        .Clear("C");  fLeg2PrefilteredTracks .Clear("C");
  if(fHistosManager) delete fHistosManager;
  if(fMixingHandler) delete fMixingHandler;
  if(fJpsiMassDist)  delete fJpsiMassDist;
  fMixingHandlerMult.Clear("C");
}


//_____________________________________________________________________________
// Initialize stuff
void AliReducedAnalysisFilterTrees::Init()
{
  AliReducedVarManager::SetDefaultVarNames();

  AliReducedVarManager::SetUseVariable(AliReducedVarManager::kDeltaVtxZMC);
  AliReducedVarManager::SetUseVariable(AliReducedVarManager::kPseudoProperDecayTime);
  AliReducedVarManager::SetUseVariable(AliReducedVarManager::kPseudoProperDecayTimeMC);
  AliReducedVarManager::SetUseVariable(AliReducedVarManager::kPseudoProperDecayTimeError);
  AliReducedVarManager::SetUseVariable(AliReducedVarManager::kTriggerPseudoProperDecayTime);
  AliReducedVarManager::SetUseVariable(AliReducedVarManager::kTripletPseudoProperDecayTime);

  fHistosManager->SetUseDefaultVariableNames(kTRUE);
  fHistosManager->SetDefaultVarNames(AliReducedVarManager::fgVariableNames,
                                     AliReducedVarManager::fgVariableUnits);
  fMixingHandler->SetHistogramManager(fHistosManager);

}


//_____________________________________________________________________________
// Process the current event
void AliReducedAnalysisFilterTrees::Process()
{

  if(!fEvent) return;
  if(!(fEvent->IsA() == AliReducedEventInfo::Class())) {
    cout << "ERROR: AliReducedAnalysisFilterTrees::Process() needs AliReducedEventInfo events" << endl;
    return;
  }

  if(fEventCounter%10000 == 0) cout << "Event no. " << fEventCounter << endl;
  fEventCounter++;
  
  AliReducedVarManager::SetEvent(fEvent);
  // Reset the values array, keep only the run wise data (LHC and ALICE GRP information).
  // NOTE: The run wise data will be updated automatically in the VarManager in case a run number change
  //       is detected.
  for(Int_t i=AliReducedVarManager::kNRunWiseVariables; i<AliReducedVarManager::kNVars; ++i)
    fValues[i] = -9999.;
  
  // Fill event information before applying event cuts
  AliReducedVarManager::FillEventInfo(fEvent, fValues);
  
  if(fComputeMult) {
    CountNch05();
    FillMultiplicity(kFALSE);
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
  
  // Fill histograms for multiplicity unfolding
  if(isEventUnbiased) {
    // For MC, only the smearing matrix is important (supposed to be independent of the trigger)
    for(Int_t cutMode=0; cutMode<4*nGlobalEstimators; cutMode=cutMode+4) {
      // MB triggered
      if(fValues[AliReducedVarManager::kINT7Triggered]) {
        fHistosManager->FillHistClass(Form("pPb_5TeV_Data_cutMode_%d",cutMode+100), fValues);
        fHistosManager->FillHistClass(Form("pPb_5TeV_MC_cutMode_%d",  cutMode+100), fValues);
      }
      // Fill for each jpsi counts
      for(Int_t j=0; j<fValues[AliReducedVarManager::kMCNJpsi]; j++) {
        fHistosManager->FillHistClass(Form("pPb_5TeV_MC_cutMode_%d",cutMode+103), fValues);
      }
    }
  }
  
  if(fOptionRunOverMC) {
    fSkipMCEvent = kFALSE;
    FillMCTruthHistograms();
    if(fSkipMCEvent) return;
  }
  
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
  
  CreateFilteredEvent();

  if(fRejectEmptyEvents && (fFilteredEvent->NPairs()+fFilteredEvent->NTracks1()+fFilteredEvent->NTracks2())==0)
    return;
  fFilteredTree->Fill();
  
}


//_____________________________________________________________________________
// Create the filtered event
void AliReducedAnalysisFilterTrees::CreateFilteredEvent()
{
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
        fHistosManager->FillHistClass(Form("EventMultRegions_%s",GetMeasMultcutName(icut)), fValues);
        // Regions relative to leading particle
        if(fValues[AliReducedVarManager::kPtLeading+icut] > fMinPtLeading) {
          fHistosManager->FillHistClass(Form("EventMultRegions2Leading_%s",GetMeasMultcutName(icut)),
                                        fValues);
        }
      }
      
      fFilteredEvent->SetNGlobalTracks(fValues[AliReducedVarManager::kNGlobalTracks+icut], icut);
      // Fill mult regions relative to jpsi
      fFilteredEvent->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksToward+icut],
                                        0, kTRUE, icut);
      fFilteredEvent->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksTransverse+icut],
                                        1, kTRUE, icut);
      fFilteredEvent->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksAway+icut],
                                        2, kTRUE, icut);
      // Fill mult regions relative to leading pt
      fFilteredEvent->SetNTracksRegions(fValues[(int)AliReducedVarManager::kNGlobalTracksToward +
                                        (int)AliReducedVarManager::kNMaxCutsGlobalTracks+icut],
                                        0, kFALSE, icut);
      fFilteredEvent->SetNTracksRegions(fValues[(int)AliReducedVarManager::kNGlobalTracksTransverse +
                                        (int)AliReducedVarManager::kNMaxCutsGlobalTracks+icut],
                                        1, kFALSE, icut);
      fFilteredEvent->SetNTracksRegions(fValues[(int)AliReducedVarManager::kNGlobalTracksAway +
                                        (int)AliReducedVarManager::kNMaxCutsGlobalTracks+icut],
                                        2, kFALSE, icut);
      
      fFilteredEvent->SetLeadingParticle(fValues[AliReducedVarManager::kPtLeading +icut],
                                         fValues[AliReducedVarManager::kPhiLeading+icut],
                                         fValues[AliReducedVarManager::kEtaLeading+icut], icut);
    }  // end loop over multiplicity cuts
  }  // end if fComputeMult
}


//_____________________________________________________________________________
// Select and add filtered pair candidates to the filtered event
void AliReducedAnalysisFilterTrees::WriteFilteredPairs()
{

  // Loop over the pair list in the unfiltered event and evaluate all the pair cuts
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
      if      (pair->PairType() == 0) pairTypeStr = "Offline";
      else if (pair->PairType() == 1) pairTypeStr = "OnTheFly";
    }
    
    if(IsPairSelected(pair,fValues)) {
      for(Int_t icut=0; icut<fPairCuts.GetEntries(); ++icut) {
        if(pair->TestFlag(icut)) {
//           fHistosManager->FillHistClass(Form("Pair_%s", fPairCuts.At(icut)->GetName()), fValues);
          for(UShort_t iflag=0; iflag<32; ++iflag) {
            AliReducedVarManager::FillPairQualityFlag(pair, iflag, fValues);
            fHistosManager->FillHistClass(Form("PairQualityFlags_%s",fPairCuts.At(icut)->GetName()), fValues);
          }
          switch(fCandidateType) {
            case AliReducedPairInfo::kJpsiToEE :
              fHistosManager->FillHistClass(  Form("Pair_%s",fPairCuts.At(icut)->GetName()), fValues);
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
        }  // end if pair cut selected
      }  // end loop over pair cuts
      TClonesArray& pairs = *(fFilteredEvent->fCandidates);
      AliReducedPairInfo* filteredPair = NULL;
      filteredPair = new(pairs[fFilteredEvent->fNV0candidates[1]]) AliReducedPairInfo(*pair);
      fFilteredEvent->fNV0candidates[1] += 1;
    }  // end if IsPairSelected
  }  // end loop over pairs
}


//_____________________________________________________________________________
// Select and add filtered tracks to the filtered event
void AliReducedAnalysisFilterTrees::WriteFilteredTracks(Int_t array/*=1*/)
{
  // Loop over the track list and evaluate all the track cuts
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
          if(!fFillTrackV0Histograms) {
            fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%s",
                                               fTrackCuts.At(icut)->GetName()), fValues);
          } else {
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
              if(track->Charge() > 0)
                fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sLambdaPosLeg",
                                                   fTrackCuts.At(icut)->GetName()), fValues);
              else
                fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sLambdaNegLeg",
                                                   fTrackCuts.At(icut)->GetName()), fValues);
            }
            if(track->IsPureLambdaLeg()) {
              if(track->Charge() > 0)
                fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sPureLambdaPosLeg",
                                                   fTrackCuts.At(icut)->GetName()), fValues);
              else
                fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sPureLambdaNegLeg",
                                                   fTrackCuts.At(icut)->GetName()), fValues);
            }
            if(track->IsALambdaLeg()) {
              if(track->Charge() > 0)
                fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sALambdaPosLeg",
                                                   fTrackCuts.At(icut)->GetName()), fValues);
              else
                fHistosManager->FillHistClass(Form("Track_WriteFilteredTracks_%sALambdaNegLeg",
                                                   fTrackCuts.At(icut)->GetName()), fValues);
            }
            if(track->IsPureALambdaLeg()) {
              if(track->Charge() > 0)
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


//_____________________________________________________________________________
// Build candidate pairs and add them to the filtered event
void AliReducedAnalysisFilterTrees::BuildCandidatePairs()
{
  // Clear the track arrays
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


//_____________________________________________________________________________
// Select leg candidates and prefilter tracks
void AliReducedAnalysisFilterTrees::RunCandidateLegsSelection(Int_t arrayOption/*=1*/)
{
  // NOTE: In the case of symmetric decay channels, the candidates are separated by charge in 
  //       fLeg1Tracks and fLeg2Tracks. For asymmetric decays, the candidates for each of the LEG1 and 
  //       LEG2 are not a priori separated by charge. It is the responsability of the analyzer to setup 
  //       the track cuts such that these are separated.

  Bool_t isAsymmetricDecayChannel = IsAsymmetricDecayChannel();
  UInt_t mcDecisionMap = 1;
  
  // Loop over the track list and evaluate all the track cuts
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
    if(fOptionRunOverMC && (fLegCandidatesMCcuts.GetEntries()>0))
      mcDecisionMap = CheckReconstructedLegMCTruth(track);
    if(isAsymmetricDecayChannel) {
      if(IsCandidateLegSelected(track,fValues,1)) {
        fLeg1Tracks.Add(track);
        FillCandidateLegHistograms("Track_LEG1_BeforePrefilter", track, 1, isAsymmetricDecayChannel);
      }
      if(IsCandidateLegSelected(track,fValues,2)) {
        fLeg2Tracks.Add(track);
        FillCandidateLegHistograms("Track_LEG2_BeforePrefilter", track, 2, isAsymmetricDecayChannel);
      }
    } else {
      // mcDecisionMap is by default true, it can be false only if running on MC and the track fails 
      // the test
      if(IsCandidateLegSelected(track,fValues,1) && mcDecisionMap) {
        if(track->Charge() > 0) {
          fLeg1Tracks.Add(track);
          FillCandidateLegHistograms("Track_LEG1_BeforePrefilter", track, 1, isAsymmetricDecayChannel);
        } else if(track->Charge() < 0) {
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
    } else {
      if(IsCandidateLegPrefilterSelected(track, fValues)) {
        if(track->Charge() > 0) {
          fLeg1PrefilteredTracks.Add(track);
          fHistosManager->FillHistClass("Track_LEG1_PrefilterTrack", fValues);
        } else if(track->Charge() < 0) {
          fLeg2PrefilteredTracks.Add(track);
          fHistosManager->FillHistClass("Track_LEG2_PrefilterTrack", fValues);
        }
      }
    }
  }  // end loop over tracks

  // TEST
  TIter iterLeg1(&fLeg1Tracks);
  TIter iterLeg2(&fLeg2Tracks);

  AliReducedTrackInfo* leg1Track = 0;
  AliReducedTrackInfo* leg2Track = 0;

  for(Int_t it1=0; it1<fLeg1Tracks.GetEntries(); ++it1) {
    leg1Track = (AliReducedTrackInfo*)iterLeg1();

    iterLeg2.Reset();
    for(Int_t it2=0; it2<fLeg2Tracks.GetEntries(); ++it2) {
      leg2Track = (AliReducedTrackInfo*)iterLeg2();

      // Verify that the two current tracks have at least 1 common bit
      ULong_t compatibilityMask = CheckTrackCompatibility(leg1Track, leg2Track, isAsymmetricDecayChannel);
      if(!compatibilityMask) continue;
      AliReducedVarManager::FillPairInfo(leg1Track, leg2Track, fCandidateType, fValues);

      Bool_t isJpsiFromB = kFALSE;
      if(fOptionRunOverMC) {
        Bool_t isJpsi = abs(leg1Track->MCPdg(0))==11 && abs(leg2Track->MCPdg(0))==11 &&
                        leg1Track->MCLabel(1)==leg2Track->MCLabel(1) && leg1Track->MCPdg(1)==443;
        isJpsiFromB = isJpsi && ((abs(leg1Track->MCPdg(2))>500  && abs(leg1Track->MCPdg(2))<599) ||
                                 (abs(leg1Track->MCPdg(2))>5000 && abs(leg1Track->MCPdg(2))<5999));
        // TODO implement with fMCMap|=(UShort_t(1)<<i
        fValues[AliReducedVarManager::kPairMCMap] = isJpsi+2*isJpsiFromB;
      }

      FillCandidatePairHistograms(compatibilityMask, 0, 1, "Pair_Candidate_AfterTrackCut",
                                  isAsymmetricDecayChannel,
                                  (fOptionRunOverMC?CheckReconstructedLegMCTruth(leg1Track,leg2Track):0));
    }
  }
}


//_____________________________________________________________________________
// Run the prefilter selection
void AliReducedAnalysisFilterTrees::RunCandidateLegsPrefilter(Int_t leg)
{
  // At this point it is assumed that the track lists are filled

  Bool_t isAsymmetricDecayChannel = IsAsymmetricDecayChannel();
  
  // Initialize iterators
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
        
        if(track->TrackId() == prefTrack->TrackId()) continue;  // avoid self-pairing
        AliReducedVarManager::FillPairInfo(track, prefTrack, fCandidateType, fValues);
        if(!IsCandidateLegPairPrefilterSelected(fValues,1)) {
          track->ResetFlags();
          break;
        }
      }  // end loop over prefiltered selected leg1 candidates tracks
    }  // end if (if symmetric decay channel: selected leg2 candidates)
    else if((leg==1 && !isAsymmetricDecayChannel) ||
            (leg==2 && (isAsymmetricDecayChannel  ||
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

  // Remove tracks
  iterLeg.Reset();
  for(Int_t it=(leg==1?fLeg1Tracks.GetEntries():fLeg2Tracks.GetEntries())-1; it>=0; --it) {
    track = (AliReducedBaseTrack*)iterLeg();
    if(!track->GetFlags()) {
      if(leg == 1) fLeg1Tracks.Remove(track);
      if(leg == 2) fLeg2Tracks.Remove(track);
    }
  }
  
  // Fill histograms after the prefilter
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


//_____________________________________________________________________________
// Run the same event pairing
void AliReducedAnalysisFilterTrees::RunSameEventPairing()
{
  
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

      // Verify that the two current tracks have at least 1 common bit
      ULong_t compatibilityMask = CheckTrackCompatibility(leg1Track, leg2Track, isAsymmetricDecayChannel);
      if(!compatibilityMask) continue;
      AliReducedVarManager::FillPairInfo(leg1Track, leg2Track, fCandidateType, fValues);
      
      Bool_t isJpsiFromB = kFALSE;
      if(fOptionRunOverMC) {
        Bool_t isJpsi = abs(leg1Track->MCPdg(0))==11 && abs(leg2Track->MCPdg(0))==11 &&
                        leg1Track->MCLabel(1)==leg2Track->MCLabel(1) && leg1Track->MCPdg(1)==443;
        isJpsiFromB = isJpsi && ((abs(leg1Track->MCPdg(2))>500  && abs(leg1Track->MCPdg(2))<599) ||
                                 (abs(leg1Track->MCPdg(2))>5000 && abs(leg1Track->MCPdg(2))<5999));
        // TODO implement with fMCMap|=(UShort_t(1)<<i
        fValues[AliReducedVarManager::kPairMCMap] = isJpsi+2*isJpsiFromB;
      }

      UInt_t recLegMCTruthMask = fOptionRunOverMC ? CheckReconstructedLegMCTruth(leg1Track,leg2Track) : 0;
      FillCandidatePairHistograms(compatibilityMask, 0, 1, "Pair_Candidate_AfterPrefilter",
                                  isAsymmetricDecayChannel, recLegMCTruthMask);
      ULong_t pairCutMask = IsCandidatePairSelected(fValues);
      if(!pairCutMask) continue;
      FillCandidatePairHistograms(compatibilityMask, pairCutMask, 1, "Pair_Candidate",
                                  isAsymmetricDecayChannel, recLegMCTruthMask);

      if(fOptionRunOverMC && fLegCandidatesMCcuts.GetEntries()>0 && !recLegMCTruthMask) continue;

      TClonesArray& pairs = *(fFilteredEvent->fCandidates);
      AliReducedPairInfo* candidatePair = new(pairs[fFilteredEvent->fNV0candidates[1]]) AliReducedPairInfo();

      if(fOptionRunOverMC && isJpsiFromB) {
        AliReducedBaseTrack* jpsimother = FindTrackByLabel(leg1Track->MCLabel(2), kTRUE);
        if(jpsimother) {
          candidatePair->PtMother (jpsimother->Pt());
          candidatePair->PhiMother(jpsimother->Phi());
          candidatePair->EtaMother(jpsimother->Eta());
        }
      } else {
        candidatePair->PtMother (0.);
        candidatePair->PhiMother(0.);
        candidatePair->EtaMother(0.);
      }

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
        FillCandidatePairHistograms(compatibilityMask, pairCutMask, 0, "Pair_Candidate",
                                    isAsymmetricDecayChannel,
                                    (fOptionRunOverMC?CheckReconstructedLegMCTruth(leg1Track,leg1Track_2):0));
        SetupPair(candidatePair, fValues);
        fFilteredEvent->fNV0candidates[1] += 1;
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
        
        // Verify that the two current tracks have at least 1 common bit
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
        FillCandidatePairHistograms(compatibilityMask, pairCutMask, 2, "Pair_Candidate",
                                    isAsymmetricDecayChannel,
                                    (fOptionRunOverMC?CheckReconstructedLegMCTruth(leg2Track,leg2Track_2):0));
        SetupPair(candidatePair, fValues);
        fFilteredEvent->fNV0candidates[1] += 1;
      }  // end loop over leg2 tracks
    }  // end loop over leg2 tracks
  }  // end if fBuildCandidateLikePairs
}


//_____________________________________________________________________________
// Fill multiplicity values (regions=false: global; regions=true: in regions).
void AliReducedAnalysisFilterTrees::FillMultiplicity(Bool_t regions/*=kFALSE*/)
{
  // Fill global tracks (both signal and MC and MC truth number of Jpsi)
  Float_t phi = 0;
  if(!regions) {
    for(Int_t icut=0; icut<GetNMeasMultCuts(); icut++) fValues[AliReducedVarManager::kNGlobalTracks+icut] = 0.;
    fValues[AliReducedVarManager::kMCNchWoPileup] = 0.;
    fValues[AliReducedVarManager::kMCNch09]       = 0.;
    fValues[AliReducedVarManager::kMCNJpsi]       = 0.;
  } else {
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
            phi         = jpsiPair->Phi();
          }
        } else {
          phi = jpsiPair->Phi();
        }
      }
      if((fValues[AliReducedVarManager::kMCNJpsi]==1) && (pairs->GetEntries()>0.)) {
        fValues[AliReducedVarManager::kPhiJpsiMCTruth] = abs(fValues[AliReducedVarManager::kPhiJpsiMCTruth] -
                                                             phi);
        fHistosManager->FillHistClass("DeltaPhi_JpsiTruth_JpsiCandidate", fValues);
      }
    } else if(fRegionsToMCTruth && (fValues[AliReducedVarManager::kMCNJpsi]==1)) {
      phi = fValues[AliReducedVarManager::kPhiJpsiMCTruth];
    }
    
    for(Int_t icut=0; icut<GetNMeasMultCuts(); icut++) {
      fValues[AliReducedVarManager::kNGlobalTracksToward     +icut] = 0.;
      fValues[AliReducedVarManager::kNGlobalTracksTransverse +icut] = 0.;
      fValues[AliReducedVarManager::kNGlobalTracksAway       +icut] = 0.;
      fValues[(int)AliReducedVarManager::kNGlobalTracksToward+
              (int)AliReducedVarManager::kNMaxCutsGlobalTracks+icut] = 0.;
      fValues[(int)AliReducedVarManager::kNGlobalTracksTransverse+
              (int)AliReducedVarManager::kNMaxCutsGlobalTracks+icut] = 0.;
      fValues[(int)AliReducedVarManager::kNGlobalTracksAway+
              (int)AliReducedVarManager::kNMaxCutsGlobalTracks+icut] = 0.;
    }
    fValues[AliReducedVarManager::kMCNch09Toward]       = 0.;
    fValues[AliReducedVarManager::kMCNch09Toward+1]     = 0.;
    fValues[AliReducedVarManager::kMCNch09Transverse]   = 0.;
    fValues[AliReducedVarManager::kMCNch09Transverse+1] = 0.;
    fValues[AliReducedVarManager::kMCNch09Away]         = 0.;
    fValues[AliReducedVarManager::kMCNch09Away+1]       = 0.;
  }
  
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

      // TODO: kMCNchWoPileup is filled here, because kMCNch is defined such that it includes pileup 
      //       tracks if one has MC with pileup. The same thing could also be achieved by introducing a 
      //       TrueMultTrackCut for |eta|<1, but atm this code is not safe for multiple 
      //       TrueMultTrackCuts.
      // TODO: kMCNch excludes tracks from jpsi daughters per default. Is this done properly here? 
      //       Since I test this on MC w/o pileup, I would think kMCNch should not include pileup and 
      //       be "correct".
      if(!regions && track->IsMCTruth() && track->Charge() && abs(track->Eta())<1.)
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

            if(!regions) {
              fValues[AliReducedVarManager::kNGlobalTracks+icut] += nRepsTotal;
              fValues[AliReducedVarManager::kPCCnRepetitions]     = nRepsPCC;

              const char* cutName = GetMeasMultcutName(icut);
              fHistosManager->FillHistClass(Form("TrackMult_%s_MeasMult",cutName), fValues);

              // Look for leading particle
              if((track->Pt()>fValues[AliReducedVarManager::kPtLeading+icut]) && !fRegionsToMCTruth) {
                fValues[AliReducedVarManager::kPtLeading +icut] = track->Pt();
                fValues[AliReducedVarManager::kPhiLeading+icut] = track->Phi();
                fValues[AliReducedVarManager::kEtaLeading+icut] = track->Eta();
              }
            } else {
              // Regions relative to jpsi
              Float_t phiIcut   = fValues[AliReducedVarManager::kPhiLeading+icut];
              Float_t delta_phi = abs(track->Phi()-phi);
              if(phi == fValues[AliReducedVarManager::kPhiLeading]) delta_phi = abs(track->Phi()-phiIcut);
              if(delta_phi<M_PI/3. || delta_phi>5*M_PI/3.)
              { fValues[AliReducedVarManager::kNGlobalTracksToward+icut]     += nRepsTotal; }
              else if((delta_phi>  M_PI/3. && delta_phi<2*M_PI/3.) ||
                      (delta_phi>4*M_PI/3. && delta_phi<5*M_PI/3.))
              { fValues[AliReducedVarManager::kNGlobalTracksTransverse+icut] += nRepsTotal; }
              else if(delta_phi>2*M_PI/3. && delta_phi<4*M_PI/3.)
              { fValues[AliReducedVarManager::kNGlobalTracksAway+icut]       += nRepsTotal; }

              // Regions relative to leading particle
              Float_t delta_phi_leading = abs(track->Phi()-phiIcut);
              if(delta_phi_leading<M_PI/3. || delta_phi_leading>5*M_PI/3.)
              { fValues[(int)AliReducedVarManager::kNGlobalTracksToward+
                        (int)AliReducedVarManager::kNMaxCutsGlobalTracks+icut] += nRepsTotal;
              } else if((delta_phi_leading>  M_PI/3. && delta_phi_leading<2*M_PI/3.) ||
                      (delta_phi_leading>4*M_PI/3. && delta_phi_leading<5*M_PI/3.))
              { fValues[(int)AliReducedVarManager::kNGlobalTracksTransverse+
                        (int)AliReducedVarManager::kNMaxCutsGlobalTracks+icut] += nRepsTotal;
              } else if(delta_phi_leading>2*M_PI/3. && delta_phi_leading<4*M_PI/3.)
              { fValues[(int)AliReducedVarManager::kNGlobalTracksAway+
                        (int)AliReducedVarManager::kNMaxCutsGlobalTracks+icut] += nRepsTotal;
              }
            }
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

        if(!regions) {
          fValues[AliReducedVarManager::kMCNch09]         += nReps;
          fValues[AliReducedVarManager::kPCCnRepetitions]  = nReps;
          for(UShort_t flag=0; flag<(8*sizeof(UInt_t)); ++flag) {
            AliReducedVarManager::FillTrackMCFlag(track, flag, fValues);
            fHistosManager->FillHistClass("TrackMultMCFlag", fValues);
          }
          fHistosManager->FillHistClass(Form("TrackMult_TrueMult"), fValues);

          // Look for leading particle
          if((track->Pt()>fValues[AliReducedVarManager::kPtLeading]) && fRegionsToMCTruth) {
            for(Int_t icut=0; icut<GetNMeasMultCuts(); icut++) {
              fValues[AliReducedVarManager::kPtLeading +icut] = track->Pt();
              fValues[AliReducedVarManager::kPhiLeading+icut] = track->Phi();
              fValues[AliReducedVarManager::kEtaLeading+icut] = track->Eta();
            }
          }
        } else {
          // Regions relative to jpsi
          Float_t delta_phi = abs(track->Phi()-phi);
          if(delta_phi<M_PI/3. || delta_phi>5*M_PI/3.)
          { fValues[AliReducedVarManager::kMCNch09Toward]     += nReps; }
          else if((delta_phi>  M_PI/3. && delta_phi<2*M_PI/3.) || (delta_phi>4*M_PI/3. && delta_phi<5*M_PI/3.))
          { fValues[AliReducedVarManager::kMCNch09Transverse] += nReps; }
          else if(delta_phi>2*M_PI/3. && delta_phi<4*M_PI/3.)
          { fValues[AliReducedVarManager::kMCNch09Away]       += nReps; }

          // Regions relative to leading particle
          Float_t delta_phi_leading = abs(track->Phi()-fValues[AliReducedVarManager::kPhiLeading]);
          if(delta_phi_leading<M_PI/3. || delta_phi_leading>5*M_PI/3.)
          { fValues[AliReducedVarManager::kMCNch09Toward+1]     += nReps; }
          else if((delta_phi_leading>  M_PI/3. && delta_phi_leading<2*M_PI/3.) ||
                  (delta_phi_leading>4*M_PI/3. && delta_phi_leading<5*M_PI/3.))
          { fValues[AliReducedVarManager::kMCNch09Transverse+1] += nReps; }
          else if(delta_phi_leading>2*M_PI/3. && delta_phi_leading<4*M_PI/3.)
          { fValues[AliReducedVarManager::kMCNch09Away+1]       += nReps; }
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
    trackID++;
    }  // end loop over tracks
  }  // end loop over arrays

  // Fill kMCNch09+1 for MC truth accepted events only
  // (necessary to compute contamination of events with |vtx_z|>10cm and Nch=0).
  if(fOptionRunOverMC && !regions) {
    fValues[AliReducedVarManager::kMCNch09+1] = fValues[AliReducedVarManager::kMCNch09];
    if(fValues[AliReducedVarManager::kMCNch09] < 1 || abs(fValues[AliReducedVarManager::kVtxZMC]) > 10.) {
      fValues[AliReducedVarManager::kMCNch09+1] = -9999.;
    }
  }
  
}



//_____________________________________________________________________________
// Apply event cuts
Bool_t AliReducedAnalysisFilterTrees::IsEventSelected(AliReducedBaseEvent* event, Float_t* values/*=0x0*/)
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
// Apply track cuts
Bool_t AliReducedAnalysisFilterTrees::IsTrackSelected(AliReducedBaseTrack* track, Float_t* values/*=0x0*/)
{
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


//_____________________________________________________________________________
// Check if track is a leg of at least one candidate pair written to filtered tree
Bool_t AliReducedAnalysisFilterTrees::TrackIsCandidateLeg(AliReducedBaseTrack* track)
{
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


//_____________________________________________________________________________
// Apply pair cuts
Bool_t AliReducedAnalysisFilterTrees::IsPairSelected(AliReducedPairInfo* pair, Float_t* values/*=0x0*/)
{
  // NOTE: Multiple cut sets are supported. The decisions are encoded in the fFlags inherited from
  //       AliReducedBaseTrack.

  if(fPairCuts.GetEntries() == 0) return kTRUE;
  pair->ResetFlags();
  
  // Loop over all the cuts and toggle a filter bit if the pair passes the corresponding cut
  for(Int_t i=0; i<fPairCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fPairCuts.At(i);
    if(values) { if(cut->IsSelected(pair, values)) pair->SetFlag(i); }
    else       { if(cut->IsSelected(pair))         pair->SetFlag(i); }
  }
  return (pair->GetFlags()>0 ? kTRUE : kFALSE);
}


//_____________________________________________________________________________
// Apply cuts for determining measured multiplicity
Bool_t AliReducedAnalysisFilterTrees::IsTrackMeasMultSelected(AliReducedBaseTrack* track,
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
Bool_t AliReducedAnalysisFilterTrees::IsTrackTrueMultSelected(AliReducedBaseTrack* track,
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
// Apply cuts to the candidate leg
Bool_t AliReducedAnalysisFilterTrees::IsCandidateLegSelected(AliReducedBaseTrack* track,
                                                             Float_t* values/*=0x0*/, Int_t whichLeg/*=1*/)
{
  Bool_t isAsymmetricDecayChannel = IsAsymmetricDecayChannel();
  
  if(fLeg1Cuts.GetEntries() == 0) return kTRUE;
  
  // Reset the flags for the track only if these are the cuts on LEG1.
  // For LEG2, we use the same bit map as for LEG1, which it is assumed was already evaluated so the
  // ResetFlags() should not be called.
  // IMPORTANT: In the case of asymmetric decay channels, the LEG1 cuts have to be always evaluated 
  //            before LEG2 cuts.
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


//_____________________________________________________________________________
// Apply prefilter cuts on the candidate legs
Bool_t AliReducedAnalysisFilterTrees::IsCandidateLegPrefilterSelected(AliReducedBaseTrack* track,
    Float_t* values/*=0x0*/, Int_t whichLeg/*=1*/)
{
  if(fLeg1PrefilterCuts.GetEntries() == 0) return kTRUE;
  Bool_t isAsymmetricDecayChannel = IsAsymmetricDecayChannel();
  
  for(Int_t i=0; i<fLeg1PrefilterCuts.GetEntries(); ++i) {
    // If there are more cuts specified, we apply an AND on all of them.
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


//_____________________________________________________________________________
// Apply pair cuts
ULong_t AliReducedAnalysisFilterTrees::IsCandidatePairSelected(Float_t* values)
{
  if(fCandidatePairCuts.GetEntries() == 0) return 1;

  ULong_t mask = 0;
  for(Int_t i=0; i<fCandidatePairCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fCandidatePairCuts.At(i);
    if(cut->IsSelected(values)) mask|=(ULong_t(1)<<i);
  }
  return mask;
}


//_____________________________________________________________________________
// Apply the prefilter pair cuts
Bool_t AliReducedAnalysisFilterTrees::IsCandidateLegPairPrefilterSelected(Float_t* values, Int_t whichLeg/*=1*/)
{
  if(fLeg1PairPrefilterCuts.GetEntries() == 0) return kTRUE;
  Bool_t isAsymmetricDecayChannel = IsAsymmetricDecayChannel();
  
  // loop over all the cuts and make a logical OR between all cuts in the list
  for(Int_t i=0; i<fLeg1PairPrefilterCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)(whichLeg==2 && isAsymmetricDecayChannel ?
                                                  fLeg2PairPrefilterCuts.At(i) : fLeg1PairPrefilterCuts.At(i));
    if(cut->IsSelected(values)) return kTRUE;
  }
  return kFALSE;
}


//_____________________________________________________________________________
// Check whether the 2 tracks fulfill at least one set of common cuts
ULong_t AliReducedAnalysisFilterTrees::CheckTrackCompatibility(AliReducedBaseTrack* leg1,
  AliReducedBaseTrack* leg2, Bool_t isAsymmetricDecayChannel)
{
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
    if(leg1->TestFlag(i) && leg2->TestFlag(32+i)) {
      mask |= (ULong_t(1)<<i);
    }
  }
  return mask;
}


//_____________________________________________________________________________
// Setup pair information
void AliReducedAnalysisFilterTrees::SetupPair(AliReducedPairInfo* pair, Float_t* values)
{
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
  pair->SetChisquare     (values[AliReducedVarManager::kPairChisquare]);  // TODO: Gauthier: kPairChi2prNDOF
  pair->SetMCMap         (values[AliReducedVarManager::kPairMCMap]);

  pair->SetPairTopology(values[AliReducedVarManager::kPairLxy],  0);
  pair->SetPairTopology(values[AliReducedVarManager::kPairLxyz], 1);
  //pair->SetPairTopology(values[AliReducedVarManager::kPairLxy], 2);
  //pair->SetPairTopology(values[AliReducedVarManager::kPairLxy], 3);
  pair->SetPairTopology(values[AliReducedVarManager::kPseudoProperDecayTime],    4);
  pair->SetPairTopology(values[AliReducedVarManager::kPseudoProperDecayTimeXYZ], 5);
  if(values[AliReducedVarManager::kPseudoProperDecayTimeError] > 0.) {
    pair->SetPairTopology(values[AliReducedVarManager::kPseudoProperDecayTime] /
                          values[AliReducedVarManager::kPseudoProperDecayTimeError], 6);
  }
  if(values[AliReducedVarManager::kPseudoProperDecayTimeXYZError] > 0.) {
    pair->SetPairTopology(values[AliReducedVarManager::kPseudoProperDecayTimeXYZ] /
                          values[AliReducedVarManager::kPseudoProperDecayTimeXYZError], 7);
  }
  pair->SetPairTopology(values[AliReducedVarManager::kPairCosPointingAngle],   8);
  pair->SetPairTopology(values[AliReducedVarManager::kPairCosPointingAngleXY], 9);
  pair->SetPairTopology(values[AliReducedVarManager::kPairDCAXY], 10);
  pair->SetPairTopology(values[AliReducedVarManager::kPairDCAZ],  11);
/* TODO: uncommenting this leads to different nof entries in Pair_Candidate and Track_WriteFilteredTracks
 * histograms. I guess it has to do with flags being overriden in FillMultiplicity. Possible solution would be
 * to either
 * - fix the using of the same flags for mult and electron cuts or
 * - to use different flags.*/
  if(fComputeMult) {
    FillMultiplicity(kTRUE);
//     // Fill mult regions relative to jpsi
//     for(Int_t icut=0; icut<GetNMeasMultCuts(); icut++) {
//       pair->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksToward+icut],     0, icut);
//       pair->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksTransverse+icut], 1, icut);
//       pair->SetNTracksRegions(fValues[AliReducedVarManager::kNGlobalTracksAway+icut],       2, icut);
//     }
  }
}


//_____________________________________________________________________________
// Check if the studied decay channel is assymetric
Bool_t AliReducedAnalysisFilterTrees::IsAsymmetricDecayChannel()
{
  // NOTE: If the decay channel studied is asymmetric, it is assumed that the cuts on the two legs are
  //       different and that the two lists of leg cuts hold the same number of cuts.
  
  if(fCandidateType == AliReducedPairInfo::kLambda0ToPPi         ) return kTRUE;
  if(fCandidateType == AliReducedPairInfo::kALambda0ToPPi        ) return kTRUE;
  if(fCandidateType == AliReducedPairInfo::kDplusToK0sPiplus     ) return kTRUE;
  if(fCandidateType == AliReducedPairInfo::kDplusToK0sKplus      ) return kTRUE;
  if(fCandidateType == AliReducedPairInfo::kDplusToPhiPiplus     ) return kTRUE;
  if(fCandidateType == AliReducedPairInfo::kDminusToK0sPiminus   ) return kTRUE;
  if(fCandidateType == AliReducedPairInfo::kDminusToK0sKminus    ) return kTRUE;
  if(fCandidateType == AliReducedPairInfo::kDminusToPhiPiminus   ) return kTRUE;
  if(fCandidateType == AliReducedPairInfo::kDzeroToKminusPiplus  ) return kTRUE;
  if(fCandidateType == AliReducedPairInfo::kADzeroToKplusPiminus ) return kTRUE;
  if(fCandidateType == AliReducedPairInfo::kDsplusToK0sKplus     ) return kTRUE;
  if(fCandidateType == AliReducedPairInfo::kDsminusToK0sKminus   ) return kTRUE;
  return kFALSE;
}


//_____________________________________________________________________________
// Fill track histogram lists according to the track flags
void AliReducedAnalysisFilterTrees::FillCandidateLegHistograms(TString histClass,
  AliReducedBaseTrack* track, Int_t leg, Bool_t isAsymmetricDecayChannel)
{
  for(Int_t icut=0; icut<fLeg1Cuts.GetEntries(); ++icut) {
    if(track->TestFlag(leg==2 && isAsymmetricDecayChannel ? icut+32 : icut)) {
      fHistosManager->FillHistClass(Form("%s_%s",histClass.Data(),
        (leg==2&&isAsymmetricDecayChannel ? fLeg2Cuts.At(icut)->GetName() : fLeg1Cuts.At(icut)->GetName())),
        fValues);
    }
  }
}


//_____________________________________________________________________________
// Fill track histogram lists according to the track flags
void AliReducedAnalysisFilterTrees::FillCandidatePairHistograms(ULong_t trackMask, ULong_t pairMask,
  Int_t pairType, TString pairClass/*="Pair_Candidate"*/, Bool_t isAsymmetricDecayChannel,
  UInt_t mcDecisions/*=0*/)
{
  // PairType can be 0, 1 or 2 corresponding to leg1&leg1, leg1&leg2 or leg2&leg2 pairs
  
  TString typeStr[3] = {"11", "12", "22"};
  if(fPairCuts.GetEntries() > 1) {
    // Loop over leg cuts
    // (for asymmetric decay channels there should always be a pair of cuts for leg1 and leg2).
    for(Int_t iTrackCut=0; iTrackCut<fLeg1Cuts.GetEntries(); ++iTrackCut) {
      // Check if the legs of the pair fulfill the respective track cuts
      if(trackMask & (ULong_t(1)<<iTrackCut)) {
        // Loop over pair cuts
        for(Int_t iPairCut=0; iPairCut<fPairCuts.GetEntries(); ++iPairCut) {
          // Check if the pair fulfills the pair cut
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
  } else {
    for(Int_t iTrackCut=0; iTrackCut<fLeg1Cuts.GetEntries(); ++iTrackCut) {
      if(trackMask & (ULong_t(1)<<iTrackCut)) {
        fHistosManager->FillHistClass(Form("%s%s_%s%s", pairClass.Data(), typeStr[pairType].Data(),
          fLeg1Cuts.At(iTrackCut)->GetName(),
          (isAsymmetricDecayChannel?Form("_%s",fLeg2Cuts.At(iTrackCut)->GetName()):"")), fValues);
        if(mcDecisions && pairType==1) {
          for(Int_t iMC=0; iMC<=fLegCandidatesMCcuts.GetEntries(); ++iMC) {
            if(mcDecisions & (UInt_t(1)<<iMC)) {
              fHistosManager->FillHistClass(Form("%s%s_%s%s_%s", pairClass.Data(), typeStr[pairType].Data(),
                fLeg1Cuts.At(iTrackCut)->GetName(),
                (isAsymmetricDecayChannel?Form("_%s",fLeg2Cuts.At(iTrackCut)->GetName()):""),
                fLegCandidatesMCcuts.At(iMC)->GetName()), fValues);
            }
          }  // end loop over MC leg cuts
        }  // end if MC truth accepted
      }  // end if track cut accepted
    }  // end loop over leg cuts
  }  // end if symmetric decay channel
}


//_____________________________________________________________________________
// Get candidate leg cut name
const Char_t* AliReducedAnalysisFilterTrees::GetCandidateLegCutName(Int_t i, Int_t leg)
{
  if(leg==2 && IsAsymmetricDecayChannel())
    return (i<fLeg2Cuts.GetEntries() ? fLeg2Cuts.At(i)->GetName() : "");
  return (i<fLeg1Cuts.GetEntries() ? fLeg1Cuts.At(i)->GetName() : "");
}


//_____________________________________________________________________________
// Count Nch in |eta|<0.5 with cuts to match expert input for PCC.
void AliReducedAnalysisFilterTrees::CountNch05()
{
  const Double_t lowPt = 0.05;
  const Double_t eta   = 0.5;

  fValues[AliReducedVarManager::kMCNch05] = 0.;
  // TODO: need to FillTrackInfo()?

  // Loop over both arrays
  for(Int_t iArray=1; iArray<=2; iArray++) {
    AliReducedTrackInfo* track;
    TClonesArray* tracklist = (iArray==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
    TIter nextTrack(tracklist);
    // Loop over tracks
    for(Int_t it=0; it<tracklist->GetEntries(); ++it) {
      track = (AliReducedTrackInfo*)nextTrack();
      if(!track->IsMCTruth())                 continue;
      // TODO: Hard-coded flag
      if(!track->TestMCFlag(6))               continue;  // secondary rejection
      if(TMath::Abs(track->Charge()) < 0.01)  continue;  // neutral rejection
      if(TMath::Abs(track->Eta())    > eta)   continue;  // acceptance cut
      if(track->Pt()                 < lowPt) continue;  // low pT cut
      fValues[AliReducedVarManager::kMCNch05]++;
    }  // end track loop
  }  // end array loop
}


//_____________________________________________________________________________
// Get PCC weight
Float_t AliReducedAnalysisFilterTrees::GetParticleWeight(AliReducedTrackInfo* track)
{
  // The weight depends on pt, species and multiplicity class.
  // It is used while reweighting particle composition in MC in order to match it to data
  // Particle type: pion = 0, proton = 1, kaon = 2, sigma- = 3, sigma+ = 4, rest = 5,
  //                weight_lambda = weight_xi = weight_sigma+
  // See $ALICE_PHYSICS/PWG/Tools/AliMCSpectraWeights.h and https://alice-notes.web.cern.ch/node/1311

  Float_t ptMC = track->Pt();
  Int_t   type = 5;
  // TODO: Hard-coded flags
  if      (track->TestMCFlag(7))  type = 0;  // pion   primary
  else if (track->TestMCFlag(8))  type = 2;  // kaon   primary
  else if (track->TestMCFlag(9))  type = 1;  // proton primary
  else if (track->TestMCFlag(10)) type = 4;  // sigma+ primary
  else if (track->TestMCFlag(11)) type = 3;  // sigma- primary

  // Check if particle is a detected secondary from weak decay
  // (as opposed to secondaries from detector material interaction).
  // In that case the weight of the mother has to be applied for all daughters.
  // In that case the mother particle was put in array 1.
  if(type == 5 && !track->IsMCTruth() && !track->TestMCFlag(6)) {
    Int_t pdgmother = abs(track->MCPdg(1));
    if      (pdgmother == 3222) type = 4;  // mother is sigma+  -> use sigma+ weight
    else if (pdgmother == 3112) type = 3;  // mother is sigma-  -> use sigma- weight
    else if (pdgmother == 310)  type = 2;  // mother is K0short -> use kaon   weight
    else if (pdgmother == 211)  type = 0;  // mother is pion    -> use pion   weight
    else if (pdgmother == 3122) type = 4;  // mother is lambda  -> use sigma+ weight

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
          if(mother->TestMCFlag(6)) ptMC = mother->PtMC();
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
// Decide how often to repeat particle in MC to match data.
Int_t AliReducedAnalysisFilterTrees::GetNRepetitions(Float_t scalingFactor, Int_t part)
{
  Int_t   nRepetitions = (Int_t)scalingFactor;
  Float_t rest         = scalingFactor - nRepetitions;

  fRand->SetSeed(GetSeed(part));
  nRepetitions += (fRand->Rndm()<=rest) ? 1 : 0;
  return nRepetitions;
}


//_____________________________________________________________________________
// Define random (but reproducable) seed
unsigned long AliReducedAnalysisFilterTrees::GetSeed(Int_t part)
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
void AliReducedAnalysisFilterTrees::Finish()
{
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


//_____________________________________________________________________________
// Check a reconstructed track against all the specified MC truth cuts
UInt_t AliReducedAnalysisFilterTrees::CheckReconstructedLegMCTruth(AliReducedBaseTrack* track)
{
  // TODO: In the fLegCandidatesMCcuts one can also add AliSignalMC objects which can then be tested
  //       using the AliReducedTrackInfo::fMCPdg[].

  if(fLegCandidatesMCcuts.GetEntries() == 0) return 0;

  UInt_t decisionMap = 0;
  for(Int_t i=0; i<fLegCandidatesMCcuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fLegCandidatesMCcuts.At(i);
    if(cut->IsSelected(track)) decisionMap |= (UInt_t(1)<<i);
  }

  return decisionMap;
}


//_____________________________________________________________________________
// Check the pair of tracks to see if they match the defined MC cuts and in addition
// that they have the same mother, if requested.
UInt_t AliReducedAnalysisFilterTrees::CheckReconstructedLegMCTruth(AliReducedBaseTrack* ptrack,
                                                                   AliReducedBaseTrack* ntrack)
{
  // NOTE: The condition for the two tracks to have the same mother requires information on the MC 
  //       label, which is available just in the full track information
  //       (AliReducedTrackInfo::fMCLabels[]). The consequence is that for the jpsi2ee analysis, the 
  //       reconstructed tracks need to be always written as full tracks.

  // Check that both tracks are full tracks
  if(ptrack->IsA() != AliReducedTrackInfo::Class()) return 0;
  if(ntrack->IsA() != AliReducedTrackInfo::Class()) return 0;

  // Check the MC requirements on each of the leg and their logical intersection
  if(fLegCandidatesMCcuts.GetEntries() == 0) return 0;
  UInt_t pTrackDecisions = CheckReconstructedLegMCTruth(ptrack);
  if(!pTrackDecisions) return 0;
  UInt_t nTrackDecisions = CheckReconstructedLegMCTruth(ntrack);

  // If the tracks have the same mother, check that the mother fullfills any MC truth requirement
  if(TMath::Abs(((AliReducedTrackInfo*)ptrack)->MCLabel(1)) ==
     TMath::Abs(((AliReducedTrackInfo*)ntrack)->MCLabel(1)))
  {
    AliReducedTrackInfo* mother = FindTrackByLabel(((AliReducedTrackInfo*)ptrack)->MCLabel(1), kTRUE);
    UInt_t motherDecisions = CheckMotherMCTruth(mother);
    if(!motherDecisions) return 0;
  }

  // Check the tracks agains the leg MCtruth cuts
  UInt_t decisions = 0;
  for(Int_t i=0; i<fLegCandidatesMCcuts.GetEntries(); ++i) {
    Bool_t pDecision = (pTrackDecisions & (UInt_t(1)<<i));
    Bool_t nDecision = (nTrackDecisions & (UInt_t(1)<<i));
    // If requested, check that the tracks have the same mother
    Bool_t sameMotherDecision = kTRUE;
    if(fLegCandidatesMCcuts_RequestSameMother[i])
      sameMotherDecision = (TMath::Abs(((AliReducedTrackInfo*)ptrack)->MCLabel(1)) ==
                            TMath::Abs(((AliReducedTrackInfo*)ntrack)->MCLabel(1)));
    if(sameMotherDecision && pDecision && nDecision) decisions |= (UInt_t(1)<<i);
  }

  return decisions;
}


//_____________________________________________________________________________
// Fill histograms with pure MC signal
void AliReducedAnalysisFilterTrees::FillMCTruthHistograms()
{
  // Loop over the first track array
  LoopOverMCTracks(1);
  // and over the second
  // NOTE: In the current model, handling the MC truth info requires the labels, which are properties of the
  //       full track, so there is no point in looping over the second track array which, if it exists,
  //       contains just base tracks.
  //LoopOverMCTracks(2);
}


//_____________________________________________________________________________
// Loop over the track array and check the pure MC tracks against the defined MC selections
void AliReducedAnalysisFilterTrees::LoopOverMCTracks(Int_t trackArray/*=1*/)
{
  AliReducedTrackInfo* mother    = 0x0;
  AliReducedTrackInfo* daughter1 = 0x0;
  AliReducedTrackInfo* daughter2 = 0x0;
  
  TClonesArray* trackList = (trackArray==1 ? fEvent->GetTracks() : fEvent->GetTracks2());
  if(!trackList) return;
  TIter nextTrack(trackList);

  // If the pt dependent weights were set, check the weight and reject randomly the event
  if(fMCJpsiPtWeights) {
    for(Int_t it=0; it<trackList->GetEntries(); ++it) {
      mother = (AliReducedTrackInfo*)nextTrack();
      if(!mother->IsMCKineParticle()) continue;
      if(!(mother->MCPdg(0)==443) || mother->MCPdg(1)==443) continue;

      // Apply selections on the Jpsi mother
      UInt_t motherDecisions = CheckMotherMCTruth(mother);
      if(!motherDecisions) continue;
      // Apply only selections to one of the MCJpsi cuts. TODO Why needed?
      if((fReweightCut>=0) && (fReweightCut<GetNJpsiMotherMCCuts()) &&
         !(motherDecisions&UInt_t(1)<<fReweightCut)) continue;

      Double_t pt = mother->Pt();
      if(pt > fMCJpsiPtWeights->GetXaxis()->GetXmax()) pt = fMCJpsiPtWeights->GetXaxis()->GetXmax();
      Double_t weight = fMCJpsiPtWeights->GetBinContent(fMCJpsiPtWeights->FindBin(pt));
      if(weight > 1.0) weight = 1.0;
      Double_t rnd = gRandom->Rndm();
      if(weight < rnd) {
        fSkipMCEvent = kTRUE;
        return;
      }
    }
  }
  
  nextTrack.Reset();
  for(Int_t it=0; it<trackList->GetEntries(); ++it) {
    mother = (AliReducedTrackInfo*)nextTrack();
    if(!mother->IsMCKineParticle()) continue;
    // The next line reduces the run time when all primaries are stored
    if(!(mother->MCPdg(0)==443) || mother->MCPdg(1)==443) continue;

    // Apply selections on the jpsi mother
    UInt_t motherDecisions = CheckMotherMCTruth(mother);
    if(!motherDecisions) continue;
    
    // Find Jpsi daughters (needed to compute 2-track properties like the polarization, etc.)
    Int_t daughter1Label = 0;
    Int_t daughter2Label = 0;
    FindJpsiTruthLegs(mother, daughter1Label, daughter2Label);
    daughter1 = FindTrackByLabel(daughter1Label, kTRUE);
    daughter2 = FindTrackByLabel(daughter2Label, kTRUE);

    // Reset track variables and fill info
    for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i)
      fValues[i] = -9999.;
    AliReducedVarManager::FillMCTruthInfo(mother, fValues, daughter1, daughter2);

    // Loop over jpsi mother selections and fill histograms before the kine cuts on electrons
    for(Int_t iCut=0; iCut<fJpsiMotherMCcuts.GetEntries(); ++iCut) {
      if(!(motherDecisions & (UInt_t(1)<<iCut))) continue;
      fHistosManager->FillHistClass(Form("PureMCTRUTH_BeforeSelection_%s",
                                         fJpsiMotherMCcuts.At(iCut)->GetName()), fValues);
    }

    if(!daughter1) continue;
    if(!daughter2) continue;

    // Apply selections on pure MC daughter electrons (kine cuts)
    UInt_t daughter1Decisions = CheckDaughterMCTruth(daughter1);
    if(!daughter1Decisions) continue;
    UInt_t daughtersDecisions = daughter1Decisions & CheckDaughterMCTruth(daughter2);
    if(!daughtersDecisions) continue;
    for(Int_t iCut=0; iCut<fJpsiMotherMCcuts.GetEntries(); ++iCut) {
      // Fill histogram for MCTruth Jpsi after selection on daughters.
      // Reset track variables and fill info.
      // Necessary because of FillPairInfo of detected daughters in loop.
      for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kNTrackVars; ++i)
        fValues[i] = -9999.;
      AliReducedVarManager::FillMCTruthInfo(mother, fValues, daughter1, daughter2);
      if(!(motherDecisions    & (UInt_t(1)<<iCut))) continue;
      if(!(daughtersDecisions & (UInt_t(1)<<iCut))) continue;
      fHistosManager->FillHistClass(Form("PureMCTRUTH_AfterSelection_%s",
                                         fJpsiMotherMCcuts.At(iCut)->GetName()), fValues);
      for(Int_t jcut=0; jcut<GetNMeasMultCuts(); jcut++) {
        fHistosManager->FillHistClass(Form("JpsiPtMultCorrel_%s_%s", GetMeasMultcutName(jcut),
                                           fJpsiMotherMCcuts.At(iCut)->GetName()), fValues);
      }
    }  // end loop over fJpsiMotherMCcuts
  }  // end loop over tracks

  return;
}


//_____________________________________________________________________________
// Check the mother pure MC truth against all defined selections and return a bit map with all decisions
UInt_t AliReducedAnalysisFilterTrees::CheckMotherMCTruth(AliReducedTrackInfo* mother)
{
  if(fJpsiMotherMCcuts.GetEntries() == 0) return 0;
  
  UInt_t decisionMap = 0;
  for(Int_t i=0; i<fJpsiMotherMCcuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fJpsiMotherMCcuts.At(i);
    if(cut->IsSelected(mother)) decisionMap |= (UInt_t(1)<<i);
  }
  
  return decisionMap;
}


//_____________________________________________________________________________
// Check the daughter pure MC truth against all defined selections and return a bit map with all 
// decisions
UInt_t AliReducedAnalysisFilterTrees::CheckDaughterMCTruth(AliReducedTrackInfo* daughter)
{
  if(fJpsiElectronMCcuts.GetEntries() == 0) return 0;
  
  UInt_t decisionMap = 0;
  for(Int_t i=0; i<fJpsiElectronMCcuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fJpsiElectronMCcuts.At(i);
    if(cut->IsSelected(daughter)) decisionMap |= (UInt_t(1)<<i);
  }
  
  return decisionMap;
}


//_____________________________________________________________________________
// Search the track list for pure MC track with label and return the track pointer
AliReducedTrackInfo* AliReducedAnalysisFilterTrees::FindTrackByLabel(Int_t label, Bool_t isTruth)
{
  AliReducedTrackInfo* track     = 0x0;
  TClonesArray*        trackList = fEvent->GetTracks();
  TIter nextTrack(trackList);
  for(Int_t i=0; i<trackList->GetEntries(); ++i) {
    track = (AliReducedTrackInfo*)nextTrack();
    if(isTruth  && (!track->IsMCKineParticle())) continue;
    if(!isTruth &&   track->IsMCKineParticle())  continue;
    if(TMath::Abs(track->MCLabel(0)) == label) return track;
  }
  return 0x0;
}


//_____________________________________________________________________________
// Find the Jpsi legs in the list of pure MC truth particles
void AliReducedAnalysisFilterTrees::FindJpsiTruthLegs(AliReducedTrackInfo* mother, Int_t& leg1Label,
                                                      Int_t& leg2Label)
{
  Int_t mLabel    = mother->MCLabel(0);
  Int_t legsFound = 0;
  AliReducedTrackInfo* track = 0x0;
  
  // loop over the first track array
  TClonesArray* trackList = fEvent->GetTracks();
  TIter nextTrack(trackList);
  for(Int_t i=0; i<trackList->GetEntries(); ++i) {
    if(legsFound == 2) return;
    track = (AliReducedTrackInfo*)nextTrack();
    if(!track->IsMCKineParticle()) continue;
    if(track->MCLabel(1)==mLabel && TMath::Abs(track->MCPdg(0))==11) {
         legsFound += 1;
      if(legsFound == 1) leg1Label = track->MCLabel(0);
      if(legsFound == 2) leg2Label = track->MCLabel(0);
    }
  }
  return;
}
