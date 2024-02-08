#include <vector>
#include <iostream>
#include <set>
#include <time.h>
#include <cstdlib>

void preprocess(){
    // Declare pointer to the root file
    TFile *file = new TFile("../refine/refined_ttbar.root");

    // Declare pointer to the TTree
    TTree *tree = (TTree*)file->Get("Refined_TTree");

    // Declare local vars and have branches point to these local vars
    int num_feat = 19;
    vector<float> *features = 0;              // Used to store training variables
    int j_label;                            // Pileup label for jet
    int t_label;                            // Pileup label for trk
    int t_ID;                               // Unique trk ID used for edges
    int E_ID;                             // Event ID used for balancing
    tree->SetBranchAddress("features",&features);
    tree->SetBranchAddress("jet_label",&j_label);
    tree->SetBranchAddress("trk_label",&t_label);
    tree->SetBranchAddress("trk_ID",&t_ID);
    tree->SetBranchAddress("Event_ID",&E_ID);

    // Declare c++ vectors
    vector<vector<float> > feats;
    vector<float> trk_feat;
    vector<int> jet_label;
    vector<int> trk_label;
    vector<int> trk_ID;
    vector<int> Event_ID;

    // Read entries and store in c++ vectors for easy indexing
    for (int i=0;i<10;i++){
        tree->GetEntry(i);
        jet_label.push_back(j_label);
        trk_label.push_back(t_label);
        trk_ID.push_back(t_ID);
        Event_ID.push_back(E_ID);
        //cout << "size" << features->size() << endl;
        float *dat = features->data(); // MUST be done for every tree entry
        trk_feat.clear();              // Make sure to clean up after each track
        for(int j=0;j<num_feat;j++){
            trk_feat.push_back(*dat);
        }
        feats.push_back(trk_feat);
    }

    // Balance Events by cutting subgraphs
    // loop through each event and store vertex labels of jets in that event
    // generate random number to see if background vertex is cut
    srand( time( 0 ) );         // Set the seed of srand() with the current time
    float threshold = 0.7;      // Threshold which defines rejection of background
    float rng = 0;              // A (psuedo)random number
    set<int> s;                 // The set of all vertex labels in each event
    vector<int> labels;         // An iterable of all vertex labels in each event
    int unique_trk_ID = 1;      // Will be incremented to give tracks unique IDs

    int E_current;  // Current Event Unique ID
    int E_next;     // Next Event Unique ID
    
    E_current=Event_ID[0]; // Initilize current event ID

    for (int i=0;i<trk_label.size();i++){
        
    }

}

// TO-DO:
// -Then balance with vectors
// -Then split
// -Then convert to zero based indices
// -Then create edges
// -Then write features and edges to root file
// -Then cout info into a text file
