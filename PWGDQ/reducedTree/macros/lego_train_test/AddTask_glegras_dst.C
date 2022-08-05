//===============================================================================================================
// addtask to create trees for J/psi - hadron correlation analysis in pp 13TeV (last updated: 26/04/2019)
//===============================================================================================================

/*R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGGA/GammaConv/macros/AddTask_V0Reader.C>*/
/*
#include "$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskImproveITS.C"
#include "$ALICE_PHYSICS/PWGHF/vertexingHF/AliAnalysisTaskSEImproveITS.h"
#include "$ALICE_PHYSICS/PWGHF/vertexingHF/AliAnalysisTaskSEImproveITS.cxx"*/

AliAnalysisCuts* EventFilter(Bool_t isAOD);
AliAnalysisCuts* CaloClusterFilter(Bool_t isAOD);
AliAnalysisCuts* JPsiElectronFilterLowMom(Bool_t isAOD);
AliAnalysisCuts* JPsiElectronFilterHighMom(Bool_t isAOD);
AliAnalysisCuts* JPsiElectronPreFilter(Bool_t isAOD);
AliAnalysisCuts* AssociatedHadronFilter(Bool_t isAOD, Int_t spec);
AliAnalysisCuts* AssociatedHadronFilterNoDCA(Bool_t isAOD, Int_t spec);
AliAnalysisCuts* AssociatedPionFilter(Bool_t isAOD, Int_t spec);
AliAnalysisCuts* TightDCACut(Bool_t isAOD);
AliAnalysisCuts* TPCChi2Cut(Bool_t isAOD, Int_t spec);
AliAnalysisCuts* QualityCut(Bool_t isAOD, Int_t spec);

void AddMCSignals(AliAnalysisTaskReducedTreeMaker* task);
void SetInactiveBranches(AliAnalysisTaskReducedTreeMaker *task);


//_______________________________________________________________________________________________________________
AliAnalysisTask *AddTask_glegras_dst(Int_t reducedEventType=-1, Bool_t writeTree=kTRUE, TString  prod="") {

    //get current analysis manager
    //-----------------------------------------------------------------------------------------------------------
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) { Error("AddTask_glegras_dst", "No analysis manager found."); return 0; }

    // query MC handler and AOD
    //-----------------------------------------------------------------------------------------------------------
    Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
    Bool_t isAOD = mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

    //create task
    //-----------------------------------------------------------------------------------------------------------
    AliAnalysisTaskReducedTreeMaker* task = new AliAnalysisTaskReducedTreeMaker("DSTTreeMaker", kTRUE);
    /*
    if(hasMC){ 
      ///add improver task before the tree maker 
      TGrid::Connect("alien://"); // if not connected input files are not loaded by the improver task
      gSystem->Load("libPWGHFvertexingHF.so"); // load the needed library
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskImproveITS.C");
      AliAnalysisTaskSEImproveITS *itsImpr = AddTaskImproveITS(kFALSE, prod,"central");
      itsImpr->SetMimicData(kTRUE);
      itsImpr->SetAOD(isAOD);
    }
    else {
      printf("WARNING: AddTaskImproveITS() should be added on MC only ! \n");
    }*/
  
    // select trigger (according to production)
    //-----------------------------------------------------------------------------------------------------------
    Int_t triggerChoice = 4;
    
    cout << "reducedEventType: " << reducedEventType << endl;
    
    
    printf("AddTask_glegras__dst() trigger choice set to %d (%s)\n", triggerChoice, prod.Data());
    if (triggerChoice==0)       // MB (INT7)
      task->SetTriggerMask(AliVEvent::kINT7);
    else if (triggerChoice==1)  // high mult. (SPD, V0)
      task->SetTriggerMask(AliVEvent::kHighMultSPD+AliVEvent::kHighMultV0);
    else if (triggerChoice==2)  // EMCal (L0, L1)
      task->SetTriggerMask(AliVEvent::kEMC7+AliVEvent::kEMCEGA);
    else if (triggerChoice==3)  // TRD
      task->SetTriggerMask(AliVEvent::kTRD);
    else if (triggerChoice==4)  // MB + high mult.
      task->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kHighMultSPD+AliVEvent::kHighMultV0);
    else if (triggerChoice==5)  // MB + EMCal
      task->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kEMC7+AliVEvent::kEMCEGA);
    else if (triggerChoice==6)  // MB + TRD
      task->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kTRD);
    else if (triggerChoice==9)  // MB + high mult. + EMCal + TRD
      task->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kHighMultSPD+AliVEvent::kHighMultV0+AliVEvent::kEMC7+AliVEvent::kEMCEGA+AliVEvent::kTRD);
    else if (triggerChoice==10)  // MB + high mult. + EMCal
      task->SetTriggerMask(AliVEvent::kINT7+AliVEvent::kHighMultSPD+AliVEvent::kHighMultV0+AliVEvent::kEMCEGA);
    else {
      printf("WARNING: In AddTask_glegras__dst(), no trigger specified, or not supported! Using kINT7.\n");
      task->SetTriggerMask(AliVEvent::kINT7);
    }

    // pile-up, physics selection and analysis utils
    //-----------------------------------------------------------------------------------------------------------
    task->SetRejectPileup(kTRUE);    // done at analysis level   
    task->UsePhysicsSelection(kTRUE);
    task->SetUseAnalysisUtils(kTRUE);

    // toggle filling of branches of tree
    //-----------------------------------------------------------------------------------------------------------
    task->SetFillTrackInfo(kTRUE);
    task->SetFillCaloClusterInfo(kTRUE);
    task->SetFillV0Info(kFALSE);
    task->SetFillGammaConversions(kFALSE);
    task->SetFillK0s(kFALSE);
    task->SetFillLambda(kFALSE);
    task->SetFillALambda(kFALSE);
    task->SetFillFMDInfo(kFALSE);
    task->SetFillEventPlaneInfo(kFALSE);
    task->SetFillTRDMatchedTracks(kFALSE);
    if (hasMC) {
      task->SetFillMCInfo(kTRUE);
      AddMCSignals(task);
    }

    // event selection
    //-----------------------------------------------------------------------------------------------------------
    task->SetEventFilter(EventFilter(isAOD));

    // cluster selection
    //-----------------------------------------------------------------------------------------------------------
    task->AddCaloClusterFilter(CaloClusterFilter(isAOD));

    // Only keep 2% of unbiased events - accessible with event filter bit 14
    task->SetWriteUnbiasedEvents(0.02);

  
    // track selection
    //-----------------------------------------------------------------------------------------------------------
    // electron tracks
    task->AddTrackFilter(JPsiElectronFilterLowMom(isAOD),   kFALSE, 2);  // full track info, fQualityFlag bit 32 - all momentum
    // Written biased events are events with at least 2 tracks passing this filter
    //task->AddTrackFilter(JPsiElectronFilterHighMom(isAOD),  kFALSE);  // full track info, fQualityFlag bit 33 - high momentum
    task->AddTrackFilter(JPsiElectronPreFilter(isAOD),      kFALSE);  // full track info, fQualityFlag bit 33 - prefilter
    // associated unidentified hadrons
    task->AddTrackFilter(AssociatedHadronFilter(isAOD, 0),  kTRUE);   // base track info, fQualityFlag bit 34 - unid. hadron, TPC refit
    task->AddTrackFilter(AssociatedHadronFilter(isAOD, 1),  kTRUE);   // base track info, fQualityFlag bit 35 - unid. hadron, ITS + TPC refit
    task->AddTrackFilter(AssociatedHadronFilter(isAOD, 2),  kTRUE);   // base track info, fQualityFlag bit 36 - unid. hadron, ITS + TPC refit, SPD any
    task->AddTrackFilter(AssociatedHadronFilter(isAOD, 3),  kTRUE);   // base track info, fQualityFlag bit 37 - unid. hadron, TPC refit, no kinks
    task->AddTrackFilter(AssociatedHadronFilter(isAOD, 4),  kTRUE);   // base track info, fQualityFlag bit 38 - unid. hadron, ITS + TPC refit, no kinks
    task->AddTrackFilter(AssociatedHadronFilter(isAOD, 5),  kTRUE);   // base track info, fQualityFlag bit 39 - unid. hadron, ITS + TPC refit, SPD any, no kinks
    // associated pions
    task->AddTrackFilter(AssociatedPionFilter(isAOD, 0),    kTRUE);   // base track info, fQualityFlag bit 40 - pions, TPC refit
    task->AddTrackFilter(AssociatedPionFilter(isAOD, 1),    kTRUE);   // base track info, fQualityFlag bit 41 - pions, ITS + TPC refit
    task->AddTrackFilter(AssociatedPionFilter(isAOD, 2),    kTRUE);   // base track info, fQualityFlag bit 42 - pions, ITS + TPC refit, SPD any
    task->AddTrackFilter(AssociatedPionFilter(isAOD, 3),    kTRUE);   // base track info, fQualityFlag bit 43 - pions, TPC refit, no kinks
    task->AddTrackFilter(AssociatedPionFilter(isAOD, 4),    kTRUE);   // base track info, fQualityFlag bit 44 - pions, ITS + TPC refit, no kinks
    task->AddTrackFilter(AssociatedPionFilter(isAOD, 5),    kTRUE);   // base track info, fQualityFlag bit 45 - pions, ITS + TPC refit, SPD any, no kinks
    // additional cuts
    task->AddTrackFilter(TightDCACut(isAOD),                kTRUE);   // base track info, fQualityFlag bit 46 - unid. hadrons, TPC refit + tight DCA cut
    task->AddTrackFilter(TPCChi2Cut(isAOD, 0),                 kTRUE);   // base track info, fQualityFlag bit 47 - unid. hadrons, TPC refit + TPC chi2 cut
    task->AddTrackFilter(TPCChi2Cut(isAOD, 1),                 kTRUE);   // base track info, fQualityFlag bit 48 - unid. hadrons, TPC refit + TPC Ncluster cut

    task->AddTrackFilter(QualityCut(isAOD, 0),                  kTRUE);   // base track info, fQualityFlag bit 49 - unid. hadrons, cuts on TPC comparable to PWGLF cut
    task->AddTrackFilter(QualityCut(isAOD, 1),                  kTRUE);   // base track info, fQualityFlag bit 50 - unid. hadrons, cuts on TPC comparable to PWGLF cut
    task->AddTrackFilter(QualityCut(isAOD, 2),                  kTRUE);   // base track info, fQualityFlag bit 51 - unid. hadrons, cuts on TPC comparable to PWGLF cut
    task->AddTrackFilter(QualityCut(isAOD, 3),                  kTRUE);   // base track info, fQualityFlag bit 52 - unid. hadrons, cuts on TPC comparable to PWGLF cut
    task->AddTrackFilter(QualityCut(isAOD, 4),                  kTRUE);   // base track info, fQualityFlag bit 53 - unid. hadrons, cuts on TPC comparable to PWGLF cut
    task->AddTrackFilter(QualityCut(isAOD, 5),                  kTRUE);   // base track info, fQualityFlag bit 54 - unid. hadrons, cuts on TPC comparable to PWGLF cut
  

    // associated unidentified hadrons without DCA (because DCA does not make any sense when we do not have vertex)
    //task->AddTrackFilter(AssociatedHadronFilterNoDCA(isAOD, 0),  kTRUE);   // base track info, fQualityFlag bit 53 - unid. hadron, TPC refit
    task->AddTrackFilter(AssociatedHadronFilterNoDCA(isAOD, 1),  kTRUE);   // base track info, fQualityFlag bit 55 - unid. hadron, ITS + TPC refit
    task->AddTrackFilter(AssociatedHadronFilterNoDCA(isAOD, 2),  kTRUE);   // base track info, fQualityFlag bit 56 - unid. hadron, ITS + TPC refit, SPD any
    //task->AddTrackFilter(AssociatedHadronFilterNoDCA(isAOD, 3),  kTRUE);   // base track info, fQualityFlag bit 56 - unid. hadron, TPC refit, no kinks
    task->AddTrackFilter(AssociatedHadronFilterNoDCA(isAOD, 4),  kTRUE);   // base track info, fQualityFlag bit 57 - unid. hadron, ITS + TPC refit, no kinks
    task->AddTrackFilter(AssociatedHadronFilterNoDCA(isAOD, 5),  kTRUE);   // base track info, fQualityFlag bit 58 - unid. hadron, ITS + TPC refit, SPD any, no kinks



    // active/inactive branches of tree
    //-----------------------------------------------------------------------------------------------------------
    SetInactiveBranches(task);

    // set writing options
    //  AliAnalysisTaskReducedTreeMaker::kBaseEventsWithBaseTracks = 0
    //  AliAnalysisTaskReducedTreeMaker::kBaseEventsWithFullTracks = 1
    //  AliAnalysisTaskReducedTreeMaker::kFullEventsWithBaseTracks = 2
    //  AliAnalysisTaskReducedTreeMaker::kFullEventsWithFullTracks = 3
    //-----------------------------------------------------------------------------------------------------------
    if(reducedEventType!=-1) task->SetTreeWritingOption(reducedEventType);
    task->SetWriteTree(writeTree);
    //task->SetEventWritingRequirement(2,1,0.02);

    // add task to manager
    //-----------------------------------------------------------------------------------------------------------
    mgr->AddTask(task);

    // connect output containers
    //-----------------------------------------------------------------------------------------------------------
    AliAnalysisDataContainer* cReducedEvent   = mgr->CreateContainer("ReducedEventDQ", AliReducedBaseEvent::Class(), AliAnalysisManager::kExchangeContainer, "reducedEvent");
    AliAnalysisDataContainer* cReducedTree    = 0x0;
    AliAnalysisDataContainer* cEventStatsInfo = 0x0;
    if(task->WriteTree()) {
      cReducedTree    = mgr->CreateContainer("dstTree",     TTree::Class(), AliAnalysisManager::kOutputContainer, "dstTree.root");
      cEventStatsInfo = mgr->CreateContainer("EventStats",  TList::Class(), AliAnalysisManager::kOutputContainer, "dstTree.root");
    }
    mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, cReducedEvent);
    if(task->WriteTree()) {
      mgr->ConnectOutput(task, 2, cReducedTree);
      mgr->ConnectOutput(task, 3, cEventStatsInfo);
    }

    // done
    //-----------------------------------------------------------------------------------------------------------
    return task;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* EventFilter(Bool_t isAOD) {
  //
  // event cuts
  //
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  //if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  //eventCuts->SetRequireVertex();
  //eventCuts->SetMinVtxContributors(1);
  //eventCuts->SetVertexZ(-10.,10.);
  return eventCuts;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* CaloClusterFilter(Bool_t isAOD) {
  //
  // cluster cuts
  //
  AliDielectronClusterCuts* clusterCuts = new AliDielectronClusterCuts("clusters", "clusters");
  clusterCuts->SetCaloType(AliDielectronClusterCuts::kEMCal);
  clusterCuts->SetNCellsMinCut(2);
  clusterCuts->SetRejectExotics();
  clusterCuts->SetM02Cut(0.1,   0.7);
  clusterCuts->SetM20Cut(0.01,  0.7);
  clusterCuts->SetRequireTrackMatch(kTRUE);
  return clusterCuts;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* JPsiElectronFilterLowMom(Bool_t isAOD) {
  //
  // electron track cuts for p < 5 GeV
  //
  AliDielectronCutGroup*  jpsiElectrons = new AliDielectronCutGroup("jpsiElectrons", "J/psi candidate electrons (low mom. PID)", AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts*   trackCuts     = new AliDielectronVarCuts("trackCuts", "track cuts");
  trackCuts->AddCut(AliDielectronVarManager::kImpactParXY,  -1.0, 1.0);
  trackCuts->AddCut(AliDielectronVarManager::kImpactParZ,   -3.0, 3.0);
  trackCuts->AddCut(AliDielectronVarManager::kEta,          -0.9, 0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,           1.0, 1.0e+30);
  trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,      70.0, 161.0);
  //trackCuts->AddCut(AliDielectronVarManager::kPIn,          0.0, 5.0);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaEle, -4.0, 4.0);
  /*trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPro, -2000.0, 2.0, kTRUE);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPio, -2000.0, 2.0, kTRUE);*/
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPio, -2000.0, -2.0, kTRUE);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPro, -2000.0, -2.0, kTRUE);

  jpsiElectrons->AddCut(trackCuts);
  AliDielectronTrackCuts* trackCuts1 = new AliDielectronTrackCuts("trackCuts1", "track cuts 1");
  trackCuts1->SetRequireITSRefit(kTRUE);
  trackCuts1->SetRequireTPCRefit(kTRUE);
  if (isAOD) trackCuts1->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  jpsiElectrons->AddCut(trackCuts1);
  return jpsiElectrons;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* JPsiElectronFilterHighMom(Bool_t isAOD) {
  //
  // electron track cuts for p > 5 GeV
  //
  AliDielectronCutGroup*  jpsiElectrons = new AliDielectronCutGroup("jpsiElectronsHighMom", "J/psi candidate electrons (high mom. PID)", AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts*   trackCuts     = new AliDielectronVarCuts("trackCuts", "track cuts");
  trackCuts->AddCut(AliDielectronVarManager::kImpactParXY,  -1.0, 1.0);
  trackCuts->AddCut(AliDielectronVarManager::kImpactParZ,   -3.0, 3.0);
  trackCuts->AddCut(AliDielectronVarManager::kEta,          -0.9, 0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,           1.0, 1.0e+30);
  trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,      70.0, 161.0);
  trackCuts->AddCut(AliDielectronVarManager::kPIn,          5.0, 1.0e+30);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaEle, -4.0, 4.0);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPio, -2000.0, -2.0, kTRUE);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPro, -2000.0, -2.0, kTRUE);
  jpsiElectrons->AddCut(trackCuts);
  AliDielectronTrackCuts* trackCuts1 = new AliDielectronTrackCuts("trackCuts1", "track cuts 1");
  trackCuts1->SetRequireITSRefit(kTRUE);
  trackCuts1->SetRequireTPCRefit(kTRUE);
  if (isAOD) trackCuts1->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  jpsiElectrons->AddCut(trackCuts1);
  return jpsiElectrons;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* JPsiElectronPreFilter(Bool_t isAOD) {
  //
  // electron prefilter track cuts (full information stored in first array)
  //
  AliDielectronCutGroup*  jpsiElectrons = new AliDielectronCutGroup("jpsiElectronsPrefilterCut", "J/psi candidate electrons(prefilter cut)", AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts*   trackCuts     = new AliDielectronVarCuts("trackCuts", "track cuts");
  trackCuts->AddCut(AliDielectronVarManager::kImpactParXY,  -3.0, 3.0);
  trackCuts->AddCut(AliDielectronVarManager::kImpactParZ,   -10.0, 10.0);
  trackCuts->AddCut(AliDielectronVarManager::kEta,          -0.9, 0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,           0.2, 1.0e+30);
  trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,      50.0, 161.0);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaEle, -4.0, 4.0);
  //trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPro, -2000.0, 2.0, kTRUE);
  //trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPio, -2000.0, 2.0, kTRUE);
  jpsiElectrons->AddCut(trackCuts);
  AliDielectronTrackCuts* trackCuts1 = new AliDielectronTrackCuts("trackCuts1","track cuts 1");
  trackCuts1->SetRequireTPCRefit(kTRUE);
  if (isAOD) trackCuts1->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  jpsiElectrons->AddCut(trackCuts1);
  return jpsiElectrons;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* AssociatedHadronFilter(Bool_t isAOD, Int_t spec) {
  //
  // associated hadron track cuts
  // 0 - TPC refit
  // 1 - ITS + TPC refit
  // 2 - ITS + TPC refit, SPD any
  // 3 - TPC refit, no kinks
  // 4 - ITS + TPC refit, no kinks
  // 5 - ITS + TPC refit, SPD any, no kinks
  //

  AliDielectronCutGroup*  assocHadr = new AliDielectronCutGroup(Form("assocHadron%d", spec), Form("Associated hadrons %d", spec), AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts*   trackCuts = new AliDielectronVarCuts("trackCuts", "track cuts");
  trackCuts->AddCut(AliDielectronVarManager::kImpactParXY,              -1.0, 1.0);
  trackCuts->AddCut(AliDielectronVarManager::kImpactParZ,               -3.0, 3.0);
  trackCuts->AddCut(AliDielectronVarManager::kEta,                      -0.9, 0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,                        0.1, 1.0e+30);
  //trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,                  70.0, 161.0); //Moved to TPC cuts
  if (spec>=3) trackCuts->AddCut(AliDielectronVarManager::kKinkIndex0,  -0.5, 0.5);
  assocHadr->AddCut(trackCuts);
  AliDielectronTrackCuts* trackCuts1 = new AliDielectronTrackCuts("trackCuts1", "track cuts 1");
  if (spec==0 || spec==3) {
    trackCuts1->SetRequireTPCRefit(kTRUE);
  } else if (spec==1 || spec==4) {
    trackCuts1->SetRequireITSRefit(kTRUE);
    trackCuts1->SetRequireTPCRefit(kTRUE);
  } else if (spec==2 || spec==5) {
    trackCuts1->SetRequireITSRefit(kTRUE);
    trackCuts1->SetRequireTPCRefit(kTRUE);
    trackCuts1->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kAny);
  } else {
    printf("AddTask_ailecVertexcorr_dst(): Associated hadron cut spec %d not recognized! Will default to ITS+TPC refit! \n", spec);
    trackCuts1->SetRequireITSRefit(kTRUE);
    trackCuts1->SetRequireTPCRefit(kTRUE);
  }
  if (isAOD) trackCuts1->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  assocHadr->AddCut(trackCuts1);
  return assocHadr;
}

//_______________________________________________________________________________________________________________
AliAnalysisCuts* AssociatedHadronFilterNoDCA(Bool_t isAOD, Int_t spec) {
  //
  // associated hadron track cuts
  // 0 - TPC refit
  // 1 - ITS + TPC refit
  // 2 - ITS + TPC refit, SPD any
  // 3 - TPC refit, no kinks
  // 4 - ITS + TPC refit, no kinks
  // 5 - ITS + TPC refit, SPD any, no kinks
  //

  AliDielectronCutGroup*  assocHadr = new AliDielectronCutGroup(Form("assocHadronNoDCA%d", spec), Form("Associated hadrons no DCA Cut %d", spec), AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts*   trackCuts = new AliDielectronVarCuts("trackCuts", "track cuts");
  trackCuts->AddCut(AliDielectronVarManager::kEta,                      -0.9, 0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,                        0.1, 1.0e+30);
  //trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,                  70.0, 161.0);
  if (spec>=3) trackCuts->AddCut(AliDielectronVarManager::kKinkIndex0,  -0.5, 0.5);
  assocHadr->AddCut(trackCuts);
  AliDielectronTrackCuts* trackCuts1 = new AliDielectronTrackCuts("trackCuts1", "track cuts 1");
  if (spec==0 || spec==3) {
    trackCuts1->SetRequireTPCRefit(kTRUE);
  } else if (spec==1 || spec==4) {
    trackCuts1->SetRequireITSRefit(kTRUE);
    trackCuts1->SetRequireTPCRefit(kTRUE);
  } else if (spec==2 || spec==5) {
    trackCuts1->SetRequireITSRefit(kTRUE);
    trackCuts1->SetRequireTPCRefit(kTRUE);
    trackCuts1->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kAny);
  } else {
    printf("AddTask_ailecVertexcorr_dst(): Associated hadron cut spec %d not recognized! Will default to ITS+TPC refit! \n", spec);
    trackCuts1->SetRequireITSRefit(kTRUE);
    trackCuts1->SetRequireTPCRefit(kTRUE);
  }
  if (isAOD) trackCuts1->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  assocHadr->AddCut(trackCuts1);
  return assocHadr;
}

//_____________________________________________________________________________________________
AliAnalysisCuts* TightDCACut(Bool_t isAOD) {
  //
  //Tighter DCA cut
  //

  AliDielectronCutGroup*  tightDCACut = new AliDielectronCutGroup("tightDCACut", "tightDCACut", AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts*   trackCuts = new AliDielectronVarCuts("trackCuts", "track cuts");
  /*TF1* cutIPXY = new TF1("cutIPXY","0.0182+0.0350/pow(x,1.01)", 0.15, 1e3);
  TF1* cutIPXYNeg = new TF1("cutIPXYNeg","-0.0182-0.0350/pow(x,1.01)", 0.15, 1e3);
  trackCuts->AddCut(AliDielectronVarManager::kImpactParXY, cutIPXYNeg, cutIPXY, kFALSE, AliDielectronVarManager::kPt);*/
  //trackCuts->AddCut(AliDielectronVarManager::kImpactParXYsigma,         -7.0, 7.0);
  trackCuts->AddCut(AliDielectronVarManager::kImpactParZ,               -2.0, 2.0);
  trackCuts->AddCut(AliDielectronVarManager::kEta,                      -0.9, 0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,                        0.1, 1.0e+30);
  //trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,                  70.0, 161.0);
  tightDCACut->AddCut(trackCuts);

  AliDielectronTrackCuts* trackCuts1 = new AliDielectronTrackCuts("trackCuts1", "track cuts 1");
  trackCuts1->SetRequireTPCRefit(kTRUE);
  if (isAOD) trackCuts1->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  tightDCACut->AddCut(trackCuts1);

  AliESDtrackCuts* fTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
  fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  if (!isAOD) {tightDCACut->AddCut(fTrackCuts);}
  return tightDCACut;
}

//_________________________________________________________________
AliAnalysisCuts* TPCChi2Cut(Bool_t isAOD, int spec) {
  //
  // TPC chi2 cut and/or TPC Ncls
  //

  AliDielectronCutGroup*  tpcChi2Cut = new AliDielectronCutGroup(Form("TPCChi2Cut%d", spec), Form("TPCChi2Cut%d", spec), AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts*   trackCuts = new AliDielectronVarCuts("trackCuts", "track cuts");
  trackCuts->AddCut(AliDielectronVarManager::kImpactParXY,              -1.0, 1.0);
  trackCuts->AddCut(AliDielectronVarManager::kImpactParZ,               -3.0, 3.0);
  trackCuts->AddCut(AliDielectronVarManager::kEta,                      -0.9, 0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,                        0.1, 1.0e+30);
  if(spec == 0) trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,     70.0, 161.0);
  if(spec == 1) trackCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,   0., 4.0);
  tpcChi2Cut->AddCut(trackCuts);
  AliDielectronTrackCuts* trackCuts1 = new AliDielectronTrackCuts("trackCuts1", "track cuts 1");
  trackCuts1->SetRequireTPCRefit(kTRUE);
  if (isAOD) trackCuts1->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  tpcChi2Cut->AddCut(trackCuts1);
  return tpcChi2Cut;
}


//___________________________________________________________________________
AliAnalysisCuts* QualityCut(Bool_t isAOD, Int_t spec) {
  //
  // cut for charged particle multiplicity as in 
  // spec = 0 or 1: alice-notes.web.cern.ch/node/1193
  // spec = 2 or 3 or 4 or 5: PWGLF/SPECTRA/ChargedHadrons/MultDepSpec/AliMultDepSpecAnalysisTask.cxx

  AliDielectronCutGroup*  qualityCut = new AliDielectronCutGroup(Form("qualityCut%d",spec), Form("qualityCut%d",spec), AliDielectronCutGroup::kCompAND);

  AliESDtrackCuts* fTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");

  if (spec == 0 || spec == 1) {
    // same as alice-notes.web.cern.ch/node/1193 

    fTrackCuts->SetMinNCrossedRowsTPC(70);
    fTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    fTrackCuts->SetMaxChi2PerClusterTPC(4);
    fTrackCuts->SetAcceptKinkDaughters(kFALSE);
    fTrackCuts->SetRequireTPCRefit(kTRUE);
    fTrackCuts->SetRequireITSRefit(kTRUE);
    if(spec == 0) fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny); 
    if(spec == 1) fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
    fTrackCuts->SetMaxDCAToVertexZ(2);
    fTrackCuts->SetDCAToVertex2D(kFALSE);
    fTrackCuts->SetRequireSigmaToVertex(kFALSE);
    fTrackCuts->SetMaxChi2PerClusterITS(36);

  }
  if (spec == 2 || spec == 3 || spec == 4 || spec == 5) {
    // these are the default track selections (cutMode == 100)
    fTrackCuts->SetRequireTPCRefit(true);
    fTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    if (spec == 2 || spec == 4) {
      // all reconstructions since mid 2020 require a tighter cut
      fTrackCuts->SetMaxChi2PerClusterTPC(2.5);
    } 
    if (spec == 3 || spec == 5) {
      fTrackCuts->SetMaxChi2PerClusterTPC(4.0);
    }
    fTrackCuts->SetRequireITSRefit(true);
    fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    fTrackCuts->SetMaxChi2PerClusterITS(36.);
    fTrackCuts->SetDCAToVertex2D(false);
    fTrackCuts->SetRequireSigmaToVertex(false);
    fTrackCuts->SetMaxDCAToVertexZ(2.0);
    fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01"); // 7 sigma cut
    fTrackCuts->SetAcceptKinkDaughters(false);
    fTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
    if (spec == 2 || spec == 3) fTrackCuts->SetMaxChi2TPCConstrainedGlobal(36.);     // golden chi2 cut - might pose problems ConstrainTPCInner failed
    fTrackCuts->SetCutGeoNcrNcl(3, 130, 1.5, 0.85, 0.7); // geometrical length cut - this is a strong cut
  }

  AliDielectronVarCuts*   trackCuts = new AliDielectronVarCuts("trackCuts", "track cuts");
  trackCuts->AddCut(AliDielectronVarManager::kEta,                      -0.9, 0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,                        0.1, 1.0e+30);

  qualityCut->AddCut(fTrackCuts);
  qualityCut->AddCut(trackCuts);
  return qualityCut;
}


//_______________________________________________________________________________________________________________
AliAnalysisCuts* AssociatedPionFilter(Bool_t isAOD, Int_t spec) {
  //
  // associated pion track cuts
  // 0 - TPC refit
  // 1 - ITS + TPC refit
  // 2 - ITS + TPC refit, SPD any
  // 3 - TPC refit, no kinks
  // 4 - ITS + TPC refit, no kinks
  // 5 - ITS + TPC refit, SPD any, no kinks
  //
  AliDielectronCutGroup*  assocPi   = new AliDielectronCutGroup(Form("assocPion%d", spec), Form("Associated pions %d", spec), AliDielectronCutGroup::kCompAND);
  AliDielectronVarCuts*   trackCuts = new AliDielectronVarCuts("trackCuts", "track cuts");
  trackCuts->AddCut(AliDielectronVarManager::kImpactParXY,              -1.0, 1.0);
  trackCuts->AddCut(AliDielectronVarManager::kImpactParZ,               -3.0, 3.0);
  trackCuts->AddCut(AliDielectronVarManager::kEta,                      -0.9, 0.9);
  trackCuts->AddCut(AliDielectronVarManager::kPt,                        0.1, 1.0e+30);
  trackCuts->AddCut(AliDielectronVarManager::kNclsTPC,                  70.0, 161.0);
  if (spec>=3) trackCuts->AddCut(AliDielectronVarManager::kKinkIndex0,  -0.5, 0.5);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaPio,             -3.0, 3.0);
  trackCuts->AddCut(AliDielectronVarManager::kTPCnSigmaEle,             -2.0, 2000., kTRUE);
  assocPi->AddCut(trackCuts);
  AliDielectronTrackCuts* trackCuts1 = new AliDielectronTrackCuts("trackCuts1", "track cuts 1");
  if (spec==0 || spec==3) {
    trackCuts1->SetRequireTPCRefit(kTRUE);
  } else if (spec==1 || spec==4) {
    trackCuts1->SetRequireITSRefit(kTRUE);
    trackCuts1->SetRequireTPCRefit(kTRUE);
  } else if (spec==2 || spec==5) {
    trackCuts1->SetRequireITSRefit(kTRUE);
    trackCuts1->SetRequireTPCRefit(kTRUE);
    trackCuts1->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kAny);
  } else {
    printf("AddTask_ailecVertexcorr_dst(): Associated pion cut spec %d not recognized! Will default to ITS+TPC refit! \n", spec);
    trackCuts1->SetRequireITSRefit(kTRUE);
    trackCuts1->SetRequireTPCRefit(kTRUE);
  }
  if (isAOD) trackCuts1->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  assocPi->AddCut(trackCuts1);
  return assocPi;
}

//_______________________________________________________________________________________________________________
void AddMCSignals(AliAnalysisTaskReducedTreeMaker* task) {
  //
  // MC signals
  //

  // 0 = kJpsiInclusive
  AliSignalMC* jpsiInclusive = new AliSignalMC("JpsiInclusive", "",1,1);
  jpsiInclusive->SetPDGcode(0, 0, 443, kFALSE);
  task->AddMCsignal(jpsiInclusive, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 1 = kJpsiNonPrompt
  AliSignalMC* jpsiFromB = new AliSignalMC("JpsiNonPrompt","",1,2);
  jpsiFromB->SetPDGcode(0, 0, 443, kFALSE);
  jpsiFromB->SetPDGcode(0, 1, 503, kTRUE);
  task->AddMCsignal(jpsiFromB, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 2 = kJpsiPrompt
  AliSignalMC* jpsiPrompt = new AliSignalMC("JpsiPrompt","",1,2);
  jpsiPrompt->SetPDGcode(0, 0, 443, kFALSE);
  jpsiPrompt->SetPDGcode(0, 1, 503, kTRUE, kTRUE);
  task->AddMCsignal(jpsiPrompt, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 3 = kJpsiRadiative
  AliSignalMC* jpsiInclusiveRadiative = new AliSignalMC("JpsiInclusiveRadiative","",1,1);
  jpsiInclusiveRadiative->SetPDGcode(0, 0, 443, kFALSE);
  jpsiInclusiveRadiative->SetSourceBit(0, 0, AliSignalMC::kRadiativeDecay);
  task->AddMCsignal(jpsiInclusiveRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 4 = kJpsiNonRadiative
  AliSignalMC* jpsiInclusiveNonRadiative = new AliSignalMC("JpsiInclusiveNonRadiative","",1,1);
  jpsiInclusiveNonRadiative->SetPDGcode(0, 0, 443, kFALSE);
  jpsiInclusiveNonRadiative->SetSourceBit(0, 0, AliSignalMC::kRadiativeDecay, kTRUE);
  task->AddMCsignal(jpsiInclusiveNonRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 5 = kJpsiNonPromptRadiative
  AliSignalMC* jpsiFromBRadiative = new AliSignalMC("JpsiFromBRadiative","",1,2);
  jpsiFromBRadiative->SetPDGcode(0, 0, 443, kFALSE);
  jpsiFromBRadiative->SetPDGcode(0, 1, 503, kTRUE);
  jpsiFromBRadiative->SetSourceBit(0, 0, AliSignalMC::kRadiativeDecay, kFALSE);
  task->AddMCsignal(jpsiFromBRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 6 = kJpsiNonPromptNonRadiative
  AliSignalMC* jpsiFromBNonRadiative = new AliSignalMC("JpsiFromBNonRadiative","",1,2);
  jpsiFromBNonRadiative->SetPDGcode(0, 0, 443, kFALSE);
  jpsiFromBNonRadiative->SetPDGcode(0, 1, 503, kTRUE);
  jpsiFromBNonRadiative->SetSourceBit(0, 0, AliSignalMC::kRadiativeDecay, kTRUE);
  task->AddMCsignal(jpsiFromBNonRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 7 = kJpsiDecayElectron
  AliSignalMC* electronFromJpsi = new AliSignalMC("electronFromJpsiInclusive","",1,2);
  electronFromJpsi->SetPDGcode(0, 0, 11, kTRUE);
  electronFromJpsi->SetPDGcode(0, 1, 443);
  task->AddMCsignal(electronFromJpsi, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 8 = kJpsiNonPromptDecayElectron
  AliSignalMC* electronFromJpsiNonPrompt = new AliSignalMC("electronFromJpsiNonPrompt","",1,3);
  electronFromJpsiNonPrompt->SetPDGcode(0, 0, 11, kTRUE);
  electronFromJpsiNonPrompt->SetPDGcode(0, 1, 443);
  electronFromJpsiNonPrompt->SetPDGcode(0, 2, 503, kTRUE);
  task->AddMCsignal(electronFromJpsiNonPrompt, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 9 = kJpsiPromptDecayElectron
  AliSignalMC* electronFromJpsiPrompt = new AliSignalMC("electronFromJpsiPrompt","", 1,3);
  electronFromJpsiPrompt->SetPDGcode(0, 0, 11, kTRUE);
  electronFromJpsiPrompt->SetPDGcode(0, 1, 443);
  electronFromJpsiPrompt->SetPDGcode(0, 2, 503, kTRUE, kTRUE);
  task->AddMCsignal(electronFromJpsiPrompt, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 10 = kJpsiRadiativeDecayElectron
  AliSignalMC* electronFromJpsiRadiative = new AliSignalMC("electronFromJpsiRadiative","",1,2);
  electronFromJpsiRadiative->SetPDGcode(0, 0, 11, kTRUE);
  electronFromJpsiRadiative->SetPDGcode(0, 1, 443);
  electronFromJpsiRadiative->SetSourceBit(0, 1, AliSignalMC::kRadiativeDecay);
  task->AddMCsignal(electronFromJpsiRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);
  
  // 11 = kJpsiNonRadiativeDecayElectron
  AliSignalMC* electronFromJpsiNonRadiative = new AliSignalMC("electronFromJpsiNonRadiative","",1,2);
  electronFromJpsiNonRadiative->SetPDGcode(0, 0, 11, kTRUE);
  electronFromJpsiNonRadiative->SetPDGcode(0, 1, 443);
  electronFromJpsiNonRadiative->SetSourceBit(0, 1, AliSignalMC::kRadiativeDecay, kTRUE);
  task->AddMCsignal(electronFromJpsiNonRadiative, AliAnalysisTaskReducedTreeMaker::kFullTrack);

  // 12 = kJpsiDecayPhoton
  AliSignalMC* photonFromJpsiDecay = new AliSignalMC("photonFromJpsiDecay","",1,2);
  photonFromJpsiDecay->SetPDGcode(0, 0, 22);
  photonFromJpsiDecay->SetPDGcode(0, 1, 443);
  task->AddMCsignal(photonFromJpsiDecay, AliAnalysisTaskReducedTreeMaker::kFullTrack);
  
  // 13 = physical primaries
  AliSignalMC* physicalPrimary=new AliSignalMC("physicalPrimary","",1,1);
  physicalPrimary->SetSourceBit(0, 0, AliSignalMC::kPhysicalPrimary);
  task->AddMCsignal(physicalPrimary, AliAnalysisTaskReducedTreeMaker::kFullTrack);
}

//_______________________________________________________________________________________________________________
void SetInactiveBranches(AliAnalysisTaskReducedTreeMaker *task) {
  //
  // set inactive branches for treee
  //

  // event
  task->SetTreeInactiveBranch("fCentrality*");
  task->SetTreeInactiveBranch("fCentQuality");
  task->SetTreeInactiveBranch("fNV0candidates*");
  //task->SetTreeInactiveBranch("fCandidates.*");
  task->SetTreeInactiveBranch("fEventNumberInFile");
  task->SetTreeInactiveBranch("fL0TriggerInputs");
  task->SetTreeInactiveBranch("fL1TriggerInputs");
  task->SetTreeInactiveBranch("fL2TriggerInputs");
  task->SetTreeInactiveBranch("fTimeStamp");
  task->SetTreeInactiveBranch("fEventType");
 // task->SetTreeInactiveBranch("fMultiplicityEstimators*");
  task->SetTreeInactiveBranch("fMultiplicityEstimatorPercentiles*");
  task->SetTreeInactiveBranch("fIRIntClosestIntMap*");
 // task->SetTreeInactiveBranch("fVtxTPC*");
  task->SetTreeInactiveBranch("fNVtxTPCContributors");
  //task->SetTreeInactiveBranch("fVtxSPD*");
  task->SetTreeInactiveBranch("fNVtxSPDContributors");
  task->SetTreeInactiveBranch("fNPMDtracks");
  task->SetTreeInactiveBranch("fNTRDtracks");
  task->SetTreeInactiveBranch("fNTRDtracklets");
  //task->SetTreeInactiveBranch("fSPDntrackletsEta*");
//  task->SetTreeInactiveBranch("fVZEROMult*");
  task->SetTreeInactiveBranch("fZDCnEnergy*");
  task->SetTreeInactiveBranch("fZDCpEnergy*");
  task->SetTreeInactiveBranch("fZDCnTotalEnergy*");
  task->SetTreeInactiveBranch("fZDCpTotalEnergy*");
  task->SetTreeInactiveBranch("fT0amplitude*");
  task->SetTreeInactiveBranch("fT0zVertex");
  task->SetTreeInactiveBranch("fT0sattelite");
  task->SetTreeInactiveBranch("fFMD.*");
  task->SetTreeInactiveBranch("fEventPlane.*");
  task->SetTreeInactiveBranch("fTRDfired");
  task->SetTreeInactiveBranch("fBC");
  task->SetTreeInactiveBranch("fNtracksPerTrackingFlag*");
  task->SetTreeInactiveBranch("fT0TOF*");
  task->SetTreeInactiveBranch("fT0start");
  task->SetTreeInactiveBranch("fT0pileup");
  //task->SetTreeInactiveBranch("fITSClusters*");
  //task->SetTreeInactiveBranch("fSPDnSingle");

  // tracks
  //task->SetTreeInactiveBranch("fTracks.fTPCPhi");
  //task->SetTreeInactiveBranch("fTracks.fTPCPt");
 // task->SetTreeInactiveBranch("fTracks.fTPCEta");
 // task->SetTreeInactiveBranch("fTracks.fTPCDCA*");
  task->SetTreeInactiveBranch("fTracks.fTrackLength");
  task->SetTreeInactiveBranch("fTracks.fMassForTracking");
  task->SetTreeInactiveBranch("fTracks.fHelixCenter*");
  task->SetTreeInactiveBranch("fTracks.fHelixRadius");
  //task->SetTreeInactiveBranch("fTracks.fITSsignal");
  //task->SetTreeInactiveBranch("fTracks.fITSnSig*");
  task->SetTreeInactiveBranch("fTracks.fTPCActiveLength");
  task->SetTreeInactiveBranch("fTracks.fTPCGeomLength");
  task->SetTreeInactiveBranch("fTracks.fTPCdEdxInfoQmax*");
  task->SetTreeInactiveBranch("fTracks.fTPCdEdxInfoQtot*");
  task->SetTreeInactiveBranch("fTracks.fTRDntracklets*");
  //task->SetTreeInactiveBranch("fTracks.fTRDpid*");
  //task->SetTreeInactiveBranch("fTracks.fTRDpidLQ2D*");
  task->SetTreeInactiveBranch("fTracks.fTRDGTUtracklets");
  task->SetTreeInactiveBranch("fTracks.fTRDGTUlayermask");
  task->SetTreeInactiveBranch("fTracks.fTRDGTUpt");
  task->SetTreeInactiveBranch("fTracks.fTRDGTUsagitta");
  task->SetTreeInactiveBranch("fTracks.fTRDGTUPID");
  task->SetTreeInactiveBranch("fTracks.fTOFbeta");
  task->SetTreeInactiveBranch("fTracks.fTOFtime");
  task->SetTreeInactiveBranch("fTracks.fTOFdx");
  task->SetTreeInactiveBranch("fTracks.fTOFdz");
  task->SetTreeInactiveBranch("fTracks.fTOFmismatchProbab");
  task->SetTreeInactiveBranch("fTracks.fTOFchi2");
  //task->SetTreeInactiveBranch("fTracks.fTOFnSig*");
  task->SetTreeInactiveBranch("fTracks.fTOFdeltaBC");
}
