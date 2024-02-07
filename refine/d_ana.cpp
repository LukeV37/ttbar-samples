#include "d.h"

#include <string>
using std::string;

int thfff = -1;

int main(int argc, char* argv[])
{
  if (argc<2) return 1;

  const char* chtree = "bTag_AntiKt4EMTopoJets";
  TChain* tt = new TChain(chtree);
  for (int i1 = 1; i1<argc; ++i1) {
    if (strlen(argv[i1])>2 && argv[i1][0]=='-' && argv[i1][1]=='f') { thfff = argv[i1][2]-'0'; continue; }
    tt->Add(argv[i1]);
  }

  TTree* mytree = (TTree*)gROOT->FindObject(chtree);
  Long64_t nentries = mytree->GetEntries();
  d* t = new d(mytree); // all branches are set

  TFile* fout = new TFile("a1.root","RECREATE");
  t->Loop();
  fout->Write();
  delete fout;
}
