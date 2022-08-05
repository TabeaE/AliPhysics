TChain* CreateChain(const char *xmlfile, const char *type="ESD");

const char *anatype = "ESD";

void lego_train()
{
// Analysis using ESD data
// Automatically generated analysis steering macro executed in grid subjobs

   TStopwatch timer;
   timer.Start();

// Set temporary merging directory to current one
   gSystem->Setenv("TMPDIR", gSystem->pwd());

// Set temporary compilation directory to current one
   gSystem->SetBuildDir(gSystem->pwd(), kTRUE);

// Reset existing include path and add current directory first in the search
   gSystem->SetIncludePath("-I.");
// load base root libraries
   gSystem->Load("libTree");
   gSystem->Load("libGeom");
   gSystem->Load("libVMC");
   gSystem->Load("libPhysics");

   gSystem->Load("libMinuit");

// Add aditional libraries
   gSystem->Load("libVMC");
   gSystem->Load("libPhysics");
   gSystem->Load("libTree");
   gSystem->Load("libMinuit");
   gSystem->Load("libProof");
   gSystem->Load("libSTEERBase");
   gSystem->Load("libESD");
   gSystem->Load("libAOD");
// Load analysis framework libraries
   gSystem->Load("libSTEERBase");
   gSystem->Load("libESD");
   gSystem->Load("libAOD");
   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");
   gSystem->Load("libOADB");
   gSystem->Load("libCORRFW");

// include path
   TString intPath = gInterpreter->GetIncludePath();
   TObjArray *listpaths = intPath.Tokenize(" ");
   TIter nextpath(listpaths);
   TObjString *pname;
   while ((pname=(TObjString*)nextpath())) {
      TString current = pname->GetName();
      if (current.Contains("AliRoot") || current.Contains("ALICE_ROOT")) continue;
      gSystem->AddIncludePath(current);
   }
   if (listpaths) delete listpaths;
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   printf("Include path: %s\n", gSystem->GetIncludePath());


// analysis source to be compiled at runtime (if any)

// read the analysis manager from file
   AliAnalysisManager *mgr = AliAnalysisAlien::LoadAnalysisManager("lego_train.root");
   if (!mgr) return;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   plugin->SetRunMode("test");
   plugin->SetFileForTestMode("./__alice__data__2018__LHC18c__000285957__pass2_CENT_root_archive_AliESDs_3.txt");
   plugin->SetNtestFiles(3);
   mgr->SetGridHandler(plugin);
   mgr->SetDebugLevel(10);
   mgr->SetNSysInfo(100);
   mgr->PrintStatus();
   gEnv->SetValue("XNet.Debug", "1");
   mgr->StartAnalysis("localfile", 123456789, 0);
   timer.Stop();
   timer.Print();
}

