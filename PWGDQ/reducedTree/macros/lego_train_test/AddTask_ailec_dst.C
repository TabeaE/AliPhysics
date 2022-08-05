//===============================================================================================================
// addtask to create trees for ChiC detection//
//===============================================================================================================


AliAnalysisCuts* EventFilter(Bool_t isAOD);
AliAnalysisCuts* CaloClusterFilter(Bool_t isAOD);
AliAnalysisCuts* GammaConvElectronCuts(Bool_t isAOD);
AliESDv0KineCuts* V0StrongCuts(Int_t mode, Int_t type);

void AddMCSignals(AliAnalysisTaskReducedTreeMaker* task);
void SetInactiveBranches(AliAnalysisTaskReducedTreeMaker *task);
void SetActiveBranches(AliAnalysisTaskReducedTreeMaker *task);


//_______________________________________________________________________________________________________________
AliAnalysisTask *AddTask_ailec_dst(Int_t reducedEventType=-1, Bool_t writeTree=kTRUE, TString  prod="") {

    //get current analysis manager
    //-----------------------------------------------------------------------------------------------------------
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) { Error("AddTask_ailec_dst", "No analysis manager found."); return 0; }

    // query MC handler and AOD
    //-----------------------------------------------------------------------------------------------------------
    Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
    Bool_t isAOD = mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

    //create task
    //-----------------------------------------------------------------------------------------------------------
    AliAnalysisTaskReducedTreeMaker* task = new AliAnalysisTaskReducedTreeMaker("DSTTreeMaker", kTRUE);
    
    
  
    // select trigger (according to production)
    //-----------------------------------------------------------------------------------------------------------
    Int_t triggerChoice = 4;
    
    cout << "reducedEventType: " << reducedEventType << endl;
    
    
    printf("AddTask_ailec_dst() trigger choice set to %d (%s)\n", triggerChoice, prod.Data());
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
      printf("WARNING: In AddTask_ailec_dst(), no trigger specified, or not supported! Using kINT7.\n");
      task->SetTriggerMask(AliVEvent::kINT7);
    }

    // pile-up, physics selection and analysis utils
    //-----------------------------------------------------------------------------------------------------------
    task->SetRejectPileup(kFALSE);    
    task->UsePhysicsSelection(kTRUE); 
    task->SetUseAnalysisUtils(kTRUE);

    // toggle filling of branches of tree
    //-----------------------------------------------------------------------------------------------------------
    
    task->SetFillV0Info(kTRUE);
    task->SetFillGammaConversions(kTRUE);
    task->SetFillTrackInfo(kTRUE);
    task->SetFillCaloClusterInfo(kTRUE);
    task->SetFillK0s(kTRUE);
    task->SetFillLambda(kTRUE);
    task->SetFillALambda(kTRUE);
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
  
    // track selection
    //-----------------------------------------------------------------------------------------------------------
    // electron tracks
    task->SetK0sMassRange(0.44,0.55);
    task->SetLambdaMassRange(1.090,1.14);
    task->SetGammaConvMassRange(0.0,0.1);
    task->SetV0StrongCuts(V0StrongCuts(AliESDv0KineCuts::kPurity, AliESDv0KineCuts::kPP));
    
    
    // active/inactive branches of tree
    //-----------------------------------------------------------------------------------------------------------
    SetInactiveBranches(task);
    //SetActiveBranches(task);

    
    task->SetTreeWritingOption(AliAnalysisTaskReducedTreeMaker::kFullEventsWithBaseTracks);
  task->SetWriteTree(writeTree);

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
  if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
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

//______________________________________________________________________________________
AliAnalysisCuts* GammaConvElectronCuts(Bool_t isAOD) {
  //
  // Cuts for the selection of electrons from gamma conversions
  //
  if(!isAOD) {
    AliESDtrackCuts* electronCuts = new AliESDtrackCuts;
    electronCuts->SetPtRange(0.05,1.0e+30);
    electronCuts->SetEtaRange(-0.9,0.9);
    electronCuts->SetRequireTPCRefit(kTRUE);
    electronCuts->SetMinNClustersTPC(60);
    return electronCuts;
  }
}

//______________________________________________________________________________________
AliESDv0KineCuts* V0StrongCuts(Int_t mode, Int_t type) {
  //
  // cuts for the selection of gamma conversions
  //
  AliESDv0KineCuts* cuts = new AliESDv0KineCuts();
  cuts->SetMode(mode, type);
  
  // leg cuts
  cuts->SetNTPCclusters(50);
  cuts->SetTPCrefit(kTRUE);
  cuts->SetTPCchi2perCls(4.0);
  cuts->SetTPCclusterratio(0.6);
  cuts->SetNoKinks(kTRUE);
  // gamma cuts                                                                                                                      
  cuts->SetGammaCutChi2NDF(10.0);
  Float_t cosPoint[2] = {0.0, 0.02};
  cuts->SetGammaCutCosPoint(cosPoint);
  Float_t cutDCA[2] = {0.0, 0.25};
  cuts->SetGammaCutDCA(cutDCA);
  Float_t vtxR[2] = {3.0, 90.0};
  cuts->SetGammaCutVertexR(vtxR);
  Float_t psiPairCut[2]={0.0,0.05};
  cuts->SetGammaCutPsiPair(psiPairCut);
  cuts->SetGammaCutInvMass(0.05);

  
  return cuts;
}

//_______________________________________________________________________________________________________________
void AddMCSignals(AliAnalysisTaskReducedTreeMaker* task) {
  //
  // MC signals
  //

  // 0 = kJpsiInclusive
  AliSignalMC* jpsiInclusive = new AliSignalMC("JpsiInclusive", "",1,1);
  jpsiInclusive->SetPDGcode(0, 0, 443, kFALSE);
  task->AddMCsignal(jpsiInclusive, AliAnalysisTaskReducedTreeMaker::kBaseTrack);

  // 1 = kJpsiNonPrompt
  AliSignalMC* jpsiFromB = new AliSignalMC("JpsiNonPrompt","",1,2);
  jpsiFromB->SetPDGcode(0, 0, 443, kFALSE);
  jpsiFromB->SetPDGcode(0, 1, 503, kTRUE);
  task->AddMCsignal(jpsiFromB, AliAnalysisTaskReducedTreeMaker::kBaseTrack);

  // 2 = kJpsiPrompt
  AliSignalMC* jpsiPrompt = new AliSignalMC("JpsiPrompt","",1,2);
  jpsiPrompt->SetPDGcode(0, 0, 443, kFALSE);
  jpsiPrompt->SetPDGcode(0, 1, 503, kTRUE, kTRUE);
  task->AddMCsignal(jpsiPrompt, AliAnalysisTaskReducedTreeMaker::kBaseTrack);
  
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
  AliSignalMC* primaryParticles=new AliSignalMC("primary","",1,1);
  primaryParticles->SetSourceBit(0, 0, AliSignalMC::kPhysicalPrimary);
  task->AddMCsignal(primaryParticles, AliAnalysisTaskReducedTreeMaker::kBaseTrack);
}

//_______________________________________________________________________________________________________________
void SetInactiveBranches(AliAnalysisTaskReducedTreeMaker *task) {
  //
  // set inactive branches for tree
  //

  // event
  task->SetTreeInactiveBranch("fCentrality*");
  task->SetTreeInactiveBranch("fCentQuality");
  task->SetTreeInactiveBranch("fEventNumberInFile");
  task->SetTreeInactiveBranch("fL0TriggerInputs");
  task->SetTreeInactiveBranch("fL1TriggerInputs");
  task->SetTreeInactiveBranch("fL2TriggerInputs");
  task->SetTreeInactiveBranch("fTimeStamp");
  task->SetTreeInactiveBranch("fEventType");
  task->SetTreeInactiveBranch("fMultiplicityEstimators*");
  task->SetTreeInactiveBranch("fMultiplicityEstimatorPercentiles*");
  task->SetTreeInactiveBranch("fIRIntClosestIntMap*");
  task->SetTreeInactiveBranch("fVtxSPD*");
  task->SetTreeInactiveBranch("fNVtxSPDContributors");
  task->SetTreeInactiveBranch("fNPMDtracks");
  task->SetTreeInactiveBranch("fNTRDtracks");
  task->SetTreeInactiveBranch("fNTRDtracklets");
  task->SetTreeInactiveBranch("fSPDntrackletsEta*");
  task->SetTreeInactiveBranch("fVZEROMult*");
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
  task->SetTreeInactiveBranch("fITSClusters*");
  task->SetTreeInactiveBranch("fSPDnSingle");

  // tracks
  task->SetTreeInactiveBranch("fTracks.fTPCPhi");
  task->SetTreeInactiveBranch("fTracks.fTPCPt");
  task->SetTreeInactiveBranch("fTracks.fTPCEta");
  task->SetTreeInactiveBranch("fTracks.fTPCDCA*");
  task->SetTreeInactiveBranch("fTracks.fTrackLength");
  task->SetTreeInactiveBranch("fTracks.fMassForTracking");
  task->SetTreeInactiveBranch("fTracks.fHelixCenter*");
  task->SetTreeInactiveBranch("fTracks.fHelixRadius");
  task->SetTreeInactiveBranch("fTracks.fITSsignal");
  task->SetTreeInactiveBranch("fTracks.fITSnSig*");
  task->SetTreeInactiveBranch("fTracks.fTPCActiveLength");
  task->SetTreeInactiveBranch("fTracks.fTPCGeomLength");
  task->SetTreeInactiveBranch("fTracks.fTPCdEdxInfoQmax*");
  task->SetTreeInactiveBranch("fTracks.fTPCdEdxInfoQtot*");
  task->SetTreeInactiveBranch("fTracks.fTRDntracklets*");
  task->SetTreeInactiveBranch("fTracks.fTRDpid*");
  task->SetTreeInactiveBranch("fTracks.fTRDpidLQ2D*");
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
  task->SetTreeInactiveBranch("fTracks.fTOFnSig*");
  task->SetTreeInactiveBranch("fTracks.fTOFdeltaBC");

  // V0s
  /*task->SetTreeInactiveBranch("fCandidates.fP*");
  task->SetTreeInactiveBranch("fCandidates.fIsCartesian");
  task->SetTreeInactiveBranch("fCandidates.fCharge");
  task->SetTreeInactiveBranch("fCandidates.fFlags");
  task->SetTreeInactiveBranch("fCandidates.fQualityFlags");
  task->SetTreeInactiveBranch("fCandidates.fCandidateId");
  task->SetTreeInactiveBranch("fCandidates.fPairType");
  task->SetTreeInactiveBranch("fCandidates.fPairTypeSPD");
  task->SetTreeInactiveBranch("fCandidates.fLegIds*");
  task->SetTreeInactiveBranch("fCandidates.fMass*");
  task->SetTreeInactiveBranch("fCandidates.fLxy");
  task->SetTreeInactiveBranch("fCandidates.fPsProper");
  task->SetTreeInactiveBranch("fCandidates.fPointingAngle");
  task->SetTreeInactiveBranch("fCandidates.fChisquare");*/
}

//_______________________________________________________________________________________________________________
void SetActiveBranches(AliAnalysisTaskReducedTreeMaker *task) {
  //
  // set active branches for tree
  //

  // V0s
  task->SetTreeActiveBranch("fCandidates.fP*");
  task->SetTreeActiveBranch("fCandidates.fIsCartesian");
  task->SetTreeActiveBranch("fCandidates.fCharge");
  task->SetTreeActiveBranch("fCandidates.fFlags");
  task->SetTreeActiveBranch("fCandidates.fQualityFlags");
  task->SetTreeActiveBranch("fCandidates.fCandidateId");
  task->SetTreeActiveBranch("fCandidates.fPairType");
  task->SetTreeActiveBranch("fCandidates.fPairTypeSPD");
  task->SetTreeActiveBranch("fCandidates.fLegIds*");
  task->SetTreeActiveBranch("fCandidates.fMass*");
  task->SetTreeActiveBranch("fCandidates.fLxy");
  task->SetTreeActiveBranch("fCandidates.fPsProper");
  task->SetTreeActiveBranch("fCandidates.fPointingAngle");
  task->SetTreeActiveBranch("fCandidates.fChisquare");
  task->SetTreeActiveBranch("fCandidates.fPolarAngleTheta");
}
