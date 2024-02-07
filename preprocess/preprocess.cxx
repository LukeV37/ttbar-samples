#include <vector>

void preprocess(){
    // Declare pointer to the root file
    TFile *file = new TFile("../refine/refined_ttbar.root");

    // Declare pointer to the TTree
    TTree *tree = (TTree*)file->Get("Refined_TTree");

    // Declare local variables and have branches point to these local vars
    int num_feat = 19;
    vector<float> *features = 0;     // Used to store training variables
    int jet_label;                            // Pileup label for jet
    int trk_label;                            // Pileup label for trk
    int trk_ID;                               // Unique trk ID used for edges
    int Event_ID;                             // Event ID used for balancing
    tree->SetBranchAddress("features",&features);
    tree->SetBranchAddress("jet_label",&jet_label);
    tree->SetBranchAddress("trk_label",&trk_label);
    tree->SetBranchAddress("trk_ID",&trk_ID);
    tree->SetBranchAddress("Event_ID",&Event_ID);

    // Read entries
    for (int i=0;i<10;i++){
        tree->GetEntry(i);
        cout << "jet label" << jet_label << endl;
        
        //cout << "size" << features->size() << endl;
        float *feat = features->data(); // MUST be done for every tree entry
        for(int j=0;j<num_feat;j++){
            cout << feat[j] << endl;
        }
    }
}

// TO-DO:
// -Write function to read each element and return vector
// -For example, read jet_label reads only jet label and returns 1D vector
// -Read features reads only features and returns 2D
// -Then balance with vectors
// -Then split
// -Then convert to zero based indices
// -Then create edges
// -Then write features and edges to root file
// -Then cout info into a text file
