// Analysis task for calculating PCC weights from reduced trees
// Creation date: 2025/09/14
// Authors: Tabea Eder (tabea.maria.eder@cern.ch)

#ifndef ALIREDUCEDANALYSISMCWEIGHTS_H
#define ALIREDUCEDANALYSISMCWEIGHTS_H

#include <TList.h>
#include <TH3F.h>

#include "AliHistogramManager.h"
#include "AliReducedAnalysisTaskSE.h"
#include "AliReducedInfoCut.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedBaseTrack.h"
#include "AliReducedTrackInfo.h"
#include "/data/t_eder04/local/alice/AliPhysics/PWG/Tools/AliMCSpectraWeights.h"

class TRandom3;

//________________________________________________________________
class AliReducedAnalysisMCWeights : public AliReducedAnalysisTaskSE {

public:
  AliReducedAnalysisMCWeights();
  AliReducedAnalysisMCWeights(const Char_t* name, const Char_t* title);
  virtual ~AliReducedAnalysisMCWeights();

  virtual void Init();
  virtual void Process();
  virtual void Finish();


  // -------------------- Setters ------------------------------------------------------------ //

  void AddTrackCut         (AliReducedInfoCut* cut) {fTrackCuts.Add(cut);}
  void AddEventCut         (AliReducedInfoCut* cut) {fEventCuts.Add(cut);}
  void AddTrueMultTrackCut (AliReducedInfoCut* cut) {fTrueMultTrackCuts.Add(cut);}
  void AddMeasMultTrackCut (AliReducedInfoCut* cut, TH2F* hWeights=nullptr) {
    fMeasMultTrackCuts.Add(cut);
    if(!hWeights) {
      hWeights = new TH2F(Form("hWeightsTrackCuts_%s",cut->GetName()), "weights", 1, 0, 3e5, 1, 0, 1e3);
      hWeights->SetBinContent(1, 1, 1.);
    }
    fWeightsTrackCuts.Add(hWeights);
  }

  void SetComputeMult                 (Bool_t option)                {fComputeMult      = option;}
  void SetRunOverMC                   (Bool_t option)                {fOptionRunOverMC  = option;}
  void SetReweightPC                  (Bool_t option)                {fReweightPC       = option;}
  void SetMCPCWeights                 (TH3F* weights)                {fMCPCWeights      = weights;}
  void SetMCSpectraWeightObject       (AliMCSpectraWeights* weights) {fMCSpectraWeights = weights;}


  // -------------------- Getters ------------------------------------------------------------ //

  virtual AliHistogramManager* GetHistogramManager () const {return fHistosManager;}
  Int_t         GetNEventCuts              ()         const {return fEventCuts.GetEntries();}
  const Char_t* GetEventCutName            (Int_t i)  const {return (i<fEventCuts.GetEntries() ?
                                                             fEventCuts.At(i)->GetName() : "");}
  Int_t         GetNTrackCuts              ()         const {return fTrackCuts.GetEntries();}
  const Char_t* GetTrackCutName            (Int_t i)  const {return (i<fTrackCuts.GetEntries() ?
                                                             fTrackCuts.At(i)->GetName() : "");}
  Bool_t        GetComputeMult             ()         const {return fComputeMult;};
  Bool_t        GetRunOverMC               ()         const {return fOptionRunOverMC;};
  Int_t         GetNMeasMultCuts           ()         const {return fMeasMultTrackCuts.GetEntries();}
  const Char_t* GetMeasMultcutName         (Int_t i)  const {return (i<fMeasMultTrackCuts.GetEntries() ?
                                                             fMeasMultTrackCuts.At(i)->GetName() : "");}
  Bool_t        GetReweightPC              ()         const {return fReweightPC;}

  void FillMCSpectra();

protected:
  Bool_t fReweightPC;   // Whether to re-weight the particle composition
  TH3F*  fMCPCWeights;  // The weights used to correct the particle composition
  std::unique_ptr<TRandom3>fRand{};  // random generator

  AliHistogramManager* fHistosManager;  // Histogram manager

  TList  fEventCuts;            // array of event cuts used for filtering
  TList  fTrackCuts;            // array of track cuts used for filtering

  Bool_t fComputeMult;          // if true, count the tracks to compute true and measured multiplicity
  TList  fMeasMultTrackCuts;    // cuts for tracks for determination of measured multiplicity
  TList  fWeightsTrackCuts;     // histograms giving weights applied to tracks, to smear the efficiency,
                                // vs pt vs run number
  TList  fTrueMultTrackCuts;    // cuts for MC tracks for determination of true multiplicity

  Bool_t fOptionRunOverMC;      // true: trees contain MC info -> fill histos to compute efficiencies,
                                // false: run normally as on data

  Bool_t IsEventSelected(AliReducedBaseEvent* event, Float_t* values=0x0);

  AliReducedTrackInfo* FindTrackByLabel (Int_t label, Bool_t isTruth=true);

  void   FillMultiplicity               ();
  Bool_t IsTrackMeasMultSelected        (AliReducedBaseTrack* track, Float_t* values=0x0);
  Bool_t IsTrackTrueMultSelected        (AliReducedBaseTrack* track, Float_t* values=0x0);

  void          CountNch05            ();
  Float_t       GetParticleWeight     (AliReducedTrackInfo* track);
  Int_t         GetParticleType       (AliReducedTrackInfo* track);
  Float_t       GetParticleMass       (AliReducedTrackInfo* track);
  Int_t         GetNRepetitions       (Float_t scalingFactor, Int_t part);
  unsigned long GetSeed               (Int_t part);

  AliMCSpectraWeights* fMCSpectraWeights;  //-> object to determine efficiency scaling

  ClassDef(AliReducedAnalysisMCWeights, 1);
};

#endif
