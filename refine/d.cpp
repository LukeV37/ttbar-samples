#define d_cxx
#include "d.h"

#include <TH2.h>

#include <iostream>
#include <cmath>
#include <map>
#include <set>
#include <algorithm>

using namespace std;

void d::Loop()
{
  if (fChain == 0) return;

  // book histograms
  TH1F* h_dz = new TH1F("dz","",50,-0.5,0.5);
  TH1F* h_z_all = new TH1F("z_all","",100,-500.,500.);
  TH1F* h_z_good = (TH1F*)h_z_all->Clone("z_good");

  TH1* h_njet = new TH1F("njet","",100,0.,100.);
  //TH1* h_njpv = new TH1F("num_jpv","",100,0.,100.);  // Lacking jet labels
  TH1* h_ntpv = new TH1F("num_tpv","",100,0.,100.);

  int count[4] = {0};

  // Initialize output tree
  TTree *t1 = new TTree("Data", "Preprocessed Features and Labels");

  // Initialize feature vars, labels, and other info
  int num_features = 19;
  vector<float> features(num_features); 
  int jet_label;
  int trk_label;
  int jet_ID=0;
  int Event_ID=0;
  int trk_ID=0;

  // Initialize tree branches
  t1->Branch("features", &features);
  t1->Branch("jet_label", &jet_label);
  t1->Branch("trk_label", &trk_label);
  t1->Branch("trk_ID", &trk_ID);
  t1->Branch("jet_ID", &jet_ID);
  t1->Branch("Event_ID", &Event_ID);

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries = 10000;

  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (jentry && jentry%100==0) cout << jentry << " / " << nentries << '\r'; cout.flush();
    fChain->GetEntry(jentry);

    // skip events with badly reconstructed PV
    double dz = PVz-truth_PVz;
    h_dz->Fill(dz);
    h_z_all->Fill(truth_PVz);
    if (fabs(dz)>0.1) continue;
    h_z_good->Fill(truth_PVz);

    int njet = 0; if (jet_pt) njet = jet_pt->size();
    h_njet->Fill(njet);

    if (njet==0) continue;

    vector<int> vertex_per_event;
    // Loop through jets
    for (int ijet = 0; ijet<njet; ++ijet) {
        // Store jet info
        features[0] = (*jet_pt_orig)[ijet];
        features[1] = (*jet_eta)[ijet];
        features[2] = (*jet_phi)[ijet];

        // Loop through tracks
        int ntr = (*jet_trackAssoc_index)[ijet].size();         // Number of tracks associated with current jet. Note some jets do not have associated tracks.
        for (int itr = 0; itr<ntr; ++itr) {
            int jtr = (*jet_trackAssoc_index)[ijet][itr];       // jtr is the index value that maps jet_var[jet_idx][trk_idx] -> trk_var[jtr]

            // Store track info
            features[3] = (*trk_pt)[jtr] / (*trk_charge)[jtr];
            features[4] = (*trk_eta)[jtr];
            features[5] = (*trk_phi)[jtr];
            features[6] = (*trk_d0)[jtr];
            features[7] = (*trk_d0_err)[jtr];
            features[8] = (*trk_z0SinTheta)[jtr];
            features[9] = (*trk_z0SinTheta_err)[jtr];
            features[10] = (*trk_nPixHits)[jtr];
            features[11] = (*trk_nsharedPixHits)[jtr];
            features[12] = (*trk_nsplitPixHits)[jtr];
            features[13] = (*trk_nSCTHits)[jtr];
            features[14] = (*trk_nSCTSharedHits)[jtr];
            features[15] = (*trk_nInnermostPixelLayerHits)[jtr];
            features[16] = (*trk_nInnermostPixelLayerSharedHits)[jtr];
            features[17] = (*trk_nInnermostPixelLayerSplitHits)[jtr];
            features[18] = (*trk_nNextToInnermostPixelLayerHits)[jtr];

            // Store label info
            int imc = (*trk_mc_index)[jtr];  // I have no idea what imc is, but I trust sasha
            if(imc<0){                       // Treat "fake tracks" as independent class per sashas suggestion
                jet_label=-999;
                trk_label=-999;
                vertex_per_event.push_back(trk_label);
            }
            else{
                jet_label = (*mc_puevent)[imc];  // jet label that defines jet vertex. -1=HS ; else=PU
                trk_label = (*mc_putype)[imc];   // trk label that defines PU. 1=PU ; 0 = Not PU
                vertex_per_event.push_back(trk_label);
            }

            trk_ID++;       // Increment trk_ID after each jet

            t1->Fill();
        }
        jet_ID++;        // Increment jet_ID after each jet
    }


    int Print = 0;
    if(Print){
    cout << "* Event " << jentry << endl;
    for (int ijet = 0; ijet<njet; ++ijet) {
        cout << "jet " << ijet << ", pt_orig=" << (*jet_pt_orig)[ijet] << ", pt=" << (*jet_pt)[ijet] << ", eta=" << (*jet_eta)[ijet] << ", phi=" << (*jet_phi)[ijet] << endl;
        int ntr = (*jet_trackAssoc_index)[ijet].size();
        cout << "* jet_trackAssoc_index size" << ntr << endl;
        for (int itr = 0; itr<ntr; ++itr) {
            int jtr = (*jet_trackAssoc_index)[ijet][itr];
            cout << " track " << jtr << ", pt=" << (*trk_pt)[jtr] << ", eta=" << (*trk_eta)[jtr] << ", phi=" << (*trk_phi)[jtr];
            int imc = (*trk_mc_index)[jtr];
            cout << " -> mc " << imc << ", pt=" << (*mc_pt)[imc] << ", eta=" << (*mc_eta)[imc] << ", phi=" << (*mc_phi)[imc];
            cout << ", puevent=" << (*mc_puevent)[imc] << ", putype=" << (*mc_putype)[imc];
            cout << endl;
        }
      }   
    }
  set<int> vertices(vertex_per_event.begin(),vertex_per_event.end());
  vector<int> unique_vertices;
  unique_vertices.assign(vertices.begin(),vertices.end());
  for (int i=0;i<unique_vertices.size();i++){
    h_ntpv->Fill(std::count(vertex_per_event.begin(),vertex_per_event.end(),unique_vertices[i]));
  }
  vertex_per_event.clear();

  Event_ID++;
  }
  h_ntpv->Write();
}
