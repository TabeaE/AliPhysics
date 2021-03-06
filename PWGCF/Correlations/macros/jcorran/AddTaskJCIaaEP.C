//_____________________________________________________________________
AliAnalysisTask *AddTaskJCIaaEP(TString taskName, Bool_t enableEP, int EPdetID, TString cardName, TString jtrigg, TString jassoc, TString cardSetting, TString inclusFileName=""){
	// Load Custom Configuration and parameters
	// override values with parameters

	cout<<"### DEBUG Input is "<< cardName <<"\t"<<jtrigg<<"\t"<<jassoc<<"\t"<<inclusFileName<<"\t"<<"#########"<<endl;
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	//==== Set up di-hadron correlation jT task ====
	AliJCIaaEPTask *myTask = new AliJCIaaEPTask(taskName.Data(),"JOD");
	myTask->SetDebugLevel(5);
  	myTask->SetJFlowBaseTaskName("JFlowBaseTask");  // AliJCatalystTask has this name hard coded
	myTask->SetEPDector( EPdetID );
	cout << myTask->GetName() << endl;


	// === Set up JCard ====
	AliJCard *card = new AliJCard(cardName.Data());
	card->PrintOut();
	card->ReadLine( cardSetting.Data() );
	card->ReCompile();
	card->PrintOut();

  // === Create analysis object ===
  
	AliJIaaAna *fAna;
	fAna = new AliJIaaAna( kFALSE );

	fAna->SetCard( card );
	fAna->SetTrigger( jtrigg.Data() );
	fAna->SetAssoc( jassoc.Data() );
	fAna->SetEnableEP( enableEP );
	if( inclusFileName ) fAna->SetInclusiveFile(inclusFileName.Data());

	myTask->SetAnalysis( fAna );

	mgr->AddTask((AliAnalysisTask*) myTask);

	// Create containers for input/output
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();


	// Connect input/output
	mgr->ConnectInput(myTask, 0, cinput);
	AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",myTask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), myTask->GetName()));
	mgr->ConnectOutput(myTask, 1, jHist );

	return myTask;
}

