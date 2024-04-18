#include <vector>
#include <iostream>
#include <numeric>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <set>
#include <fstream>

using namespace std;

void edges(int type){
    // Declare pointer to the root file
    TFile *file = new TFile("output.root");
    cout << "Type: " << type << endl;
    // Declare pointer to the TTree
    TString run_type;
    if (type==0) run_type = "train";
    if (type==1) run_type = "test";
    TTree *tree = (TTree*)file->Get(run_type);

    // Declare local vars and have branches point to these local vars
    int num_feat = 19;
    vector<float> *features = 0;              // Used to store training variables
    int j_label;                            // Pileup label for jet
    int t_label;                            // Pileup label for trk
    int t_ID;                               // Unique trk ID used for edges
    int j_ID;                               // Unique trk ID used for balancing
    int E_ID;                             // Event ID used for splitting test train
    tree->SetBranchAddress("features",&features);
    tree->SetBranchAddress("jet_label",&j_label);
    tree->SetBranchAddress("trk_label",&t_label);
    tree->SetBranchAddress("trk_ID",&t_ID);
    tree->SetBranchAddress("jet_ID",&j_ID);
    tree->SetBranchAddress("Event_ID",&E_ID);

    // Declare c++ vectors
    vector<vector<float> > feats;
    vector<float> trk_feat;
    vector<int> jet_label;
    vector<int> trk_label;
    vector<int> trk_ID;
    vector<int> jet_ID;
    vector<int> Event_ID;
    vector<int> isKept;
    vector<int> isTrain;
    vector<int> isTest;

    // Read entries and store in c++ vectors for easy indexing
    Long64_t nentries = tree->GetEntries();
    for (int i=0;i<nentries;i++){
        tree->GetEntry(i);
        jet_label.push_back(j_label);
        trk_label.push_back(t_label);
        trk_ID.push_back(t_ID);
        jet_ID.push_back(j_ID);
        Event_ID.push_back(E_ID);
        isKept.push_back(1);
        isTrain.push_back(0);
        isTest.push_back(0);
        //cout << "size" << features->size() << endl;
        float *dat = features->data(); // MUST be done for every tree entry
        trk_feat.clear();              // Make sure to clean up after each track
        for(int j=0;j<num_feat;j++){
            trk_feat.push_back(dat[j]);
            //cout << dat[j] << endl;
        }
        feats.push_back(trk_feat);
    }
    // Now we are done reading data from root file and now we can use c++ vectors! :)
    
  //create edges and labels vectors
  struct Edge {
    int source;
	int neighbor;
  };

  struct Edge edge;
  vector<Edge> edges;
  vector<int> labels;

  set<int> J;
  vector<int> Jets;

	// Loop through Events

    // Before looping through Dataset, first determine the total number of Events
    set<int> E;
    map<int, float> SoftKiller;
    for (size_t i=0;i<trk_label.size();i++){
        E.insert(Event_ID[i]);
    }
    int num_Events = E.size();

    int E_start_idx=0;
    int E_end_idx=0;
	int E_current_ID=0;
    int middle_idx;

    //num_Events=2;
    for(int i=0;i<num_Events;i++){
        E_start_idx = E_end_idx;
        E_current_ID = Event_ID[E_start_idx];

        J.clear();
        Jets.clear();

        // Figure out the Event_ending_idx using a while loop
        while(E_current_ID == Event_ID[E_end_idx]){
            E_end_idx++;
        }

        // Figure out how many jets per event
		for(int j=E_start_idx;j<E_end_idx;j++){
            J.insert(jet_ID[j]);
        }
        Jets.assign(J.begin(),J.end());           
        //cout << "Num of Jets in Event: " << Jets.size() << endl;
        for(int k=0;k<int(Jets.size());k++){       
            vector<float> pT_vector;
            for(int j=E_start_idx;j<E_end_idx;j++){
                if(jet_ID[j]==Jets[k]){            
                    pT_vector.push_back(abs(feats[j][3]));
                }
            }
            sort(pT_vector.begin(), pT_vector.end());
            middle_idx = int(pT_vector.size() / 2.0);
            float median_pT = pT_vector[middle_idx];
            SoftKiller[Jets[k]] = median_pT;
            pT_vector.clear();
        }

        // Figure out edges
        for(int j=E_start_idx;j<E_end_idx;j++){
            for(int k=j+1;k<E_end_idx;k++){
                //if((jet_ID[j]==jet_ID[k])&&(trk_label[j]==trk_label[k])){
                if(jet_ID[j]==jet_ID[k]){
                    int jet_ID_m = jet_ID[j];
                    float trk_pt_1 = abs(feats[j][3]);
                    float trk_pt_2 = abs(feats[k][3]);
                    //cout << jet_ID_m << " " <<SoftKiller[jet_ID_m] << endl;
                    if(trk_pt_1>SoftKiller[jet_ID_m] && trk_pt_2>SoftKiller[jet_ID_m]){
                        edge.source = trk_ID[j];
                        edge.neighbor = trk_ID[k];
                        edges.push_back(edge);
                        edge.source = trk_ID[k];
                        edge.neighbor = trk_ID[j];
                        edges.push_back(edge);
                    }
                    if(trk_pt_1<SoftKiller[jet_ID_m] && trk_pt_2<SoftKiller[jet_ID_m]){
                        edge.source = trk_ID[j];
                        edge.neighbor = trk_ID[k];
                        edges.push_back(edge);
                        edge.source = trk_ID[k];
                        edge.neighbor = trk_ID[j];
                        edges.push_back(edge);
                    }
                }
            }
        }
    }

    SoftKiller.clear();

  ofstream data_file(run_type+"_data.txt");
  ofstream edges_file(run_type+"_edges.txt");

  for (size_t i=0;i<edges.size();i++){
    edges_file << edges[i].source <<" "<< edges[i].neighbor << endl;
  }

  
  for (size_t i=0;i<trk_label.size();i++){
    for (int j=0;j<num_feat;j++){
      data_file << feats[i][j] << " ";
    }
    data_file << trk_label[i] << " " << trk_ID[i] << endl;
  }

  cout << "Success: "+run_type+"_data.txt and "+run_type+"_edges.txt have been written to disk." << endl;
  cout << "Mean of Labels: " << accumulate(trk_label.begin(), trk_label.end(), 0.0) / trk_label.size() << endl;
}
