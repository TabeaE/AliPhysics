#ifdef __CLING__
// Tell  ROOT where to find AliRoot headers
R__ADD_INCLUDE_PATH($ALICE_ROOT)
#include <ANALYSIS/macros/train/AddESDHandler.C>

R__ADD_INCLUDE_PATH($ALICE_ROOT)
#include <ANALYSIS/macros/train/AddMCHandler.C>

R__ADD_INCLUDE_PATH($ALICE_ROOT)
#include <ANALYSIS/macros/train/AddAODHandler.C>

#endif

void handlers()
{
  {
   // AddAODHandler();
    AddESDHandler();
     // AddMCHandler();
      
   /* gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C"));*/
    AliVEventHandler* handler = AddESDHandler();
   //   AliVEventHandler* handler = AddAODHandler();
   //gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/train/AddMCHandler.C"));
   //AliMCEventHandler* handlerMC = AddMCHandler(kFALSE);
		

//  handlerMC->SetReadTR(1);
//  handlerMC->SetPreReadMode(AliMCEventHandler::kLmPreRead);*/
  }
}
