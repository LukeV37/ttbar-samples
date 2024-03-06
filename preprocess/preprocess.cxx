#include <vector>
#include <iostream>
#include <set>
#include <time.h>
#include <cstdlib>
#include <algorithm>
#include <cmath>

void preprocess(){
    // Declare pointer to the root file
    TFile *file = new TFile("../refine/refined_ttbar.root");
    TH1* Jet_Score_h = new TH1D("Jet_Score","",20,0.,1.);
    TH1* Jet_Score_AFTER_h = new TH1D("Jet_Score_AFTER","",20,0.,1.);
    TH1* Balance_Score_h = new TH1D("Balance_Score","",20,0.,1.);

    // Declare pointer to the TTree
    TTree *tree = (TTree*)file->Get("Data");

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
        }
        feats.push_back(trk_feat);
    }
    // Now we are done reading data from root file and now we can use c++ vectors! :)

    // Balance Events by cutting subgraphs
    // loop through each event and store vertex labels of jets in that event
    // generate random number to see if background vertex is cut
    srand( time( 0 ) );         // Set the seed of srand() with the current time
    float balance_threshold = 0.91;      // Threshold which defines rejection of background
    float HS_threshold = 0.5;
    float rng = 0;              // A (psuedo)random number generated by rand() function
    set<int> s;                 // The set of all vertex labels in each event
    set<int> E;                 // The set of all event labels
    vector<int> labels;         // An iterable of all vertex labels in each event
    vector<int> cut_labels;     // A list of which verticies to cut
    set<int> J;                 // The set of all jet IDs in each event
    vector<int> Jets;           // An iterable of all jet IDs in each event
    vector<float> Jet_Score;      // HS Jet will have higher score PU will have lower score
    int Score;                    // HS trks add to score PU trks do not add to score


    // Before looping through Dataset, first determine the total number of Events
    for (int i=0;i<trk_label.size();i++){
        E.insert(Event_ID[i]);
    }
    int num_Events = E.size();
    
    int E_current_ID;           // Set the starting seed for event ID
    int E_start_idx=0;          // Used to save idx of when event starts
    int E_end_idx=0;            // Used to save idx of when event ends
    int num_trk;                // Used to save the number of tracks per jet


    for(int i=0;i<num_Events;i++){
        E_start_idx = E_end_idx;                // At the beginning of the event pick up where the last event left off.
        E_current_ID = Event_ID[E_start_idx];   // Reset the current event ID after while loop breaks
        cut_labels.clear();                     // Reset the cut labels at the beginning of each event
        J.clear();      // Clear set of jets from event
        Jets.clear();   // Clear vecotr of jets from event
        Jet_Score.clear(); // Clear vector of jet score from each event

        // Figure out the Event_ending_idx using a while loop
        while(E_current_ID == Event_ID[E_end_idx]){
            E_end_idx++;
        }

        // Figure out how mant jets per event
        for(int j=E_start_idx;j<E_end_idx;j++){
            J.insert(jet_ID[j]);
        }
        Jets.assign(J.begin(),J.end());                 // Convert set to vector
        for(int k=0;k<Jets.size();k++){                 // Loop through all jets in event
            num_trk=0;                                  // Reset num_trk at the beginning of each jet
            Score=0;                                    // Reset score at the beginning of each jet
            for(int j=E_start_idx;j<E_end_idx;j++){     // Loop through all tracks
                if(jet_ID[j]==Jets[k]){                 // Check if jet_ID matches current jet
                    num_trk++;                          // Count number of trks in each jet
                    if(trk_label[j]==0) Score++;       // If trk is HS add to the jet score
                }
            }
            Jet_Score.push_back(float(Score)/float(num_trk)); // After looping through all trks append score
            Jet_Score_h->Fill(float(Score)/float(num_trk));
        }

        // Iterate through jets and cut jets with low score in order to balance the dataset
        for (int k=0;k<Jets.size();k++){
            rng = ((float) rand() / (float) RAND_MAX);  // Generate a (psuedo)random number between 0 and 1
            if (Jet_Score[k]<0.05){
                if (rng<0.99) cut_labels.push_back(Jets[k]);
            }
            else if (Jet_Score[k]<0.1){
                if (rng<0.94) cut_labels.push_back(Jets[k]);
            }
            else if (Jet_Score[k]<0.15){
                if (rng<0.85) cut_labels.push_back(Jets[k]);
            }
            else if (Jet_Score[k]<0.2){
                if (rng<0.75) cut_labels.push_back(Jets[k]);
            }
            else if (Jet_Score[k]<0.3){
                if (rng<0.7) cut_labels.push_back(Jets[k]);
            }
            else if (Jet_Score[k]<0.35){
                if (rng<0.55) cut_labels.push_back(Jets[k]);
            }
            else if (Jet_Score[k]<0.45){
                if (rng<0.35) cut_labels.push_back(Jets[k]);
            }
            else continue;
        }

        // Loop through tracks in the event, check if jet_lable is an element of cut list, if so set isKept to zero.
        for(int j=E_start_idx;j<E_end_idx;j++){
            if (find(cut_labels.begin(), cut_labels.end(), jet_ID[j]) != cut_labels.end() ){
                isKept[j] = 0;      // Vertex label is found in cut list and trk is discarded
            };
        }
    }

    // Now that the dataset has been balanced, we split into testing and training datasets
    float split = 0.75;
    vector<int> Events;
    vector<int> TrainEvents;
    Events.assign(E.begin(),E.end());
    int split_idx = int(Events.size() * split);
    for(int i=0;i<split_idx;i++){
        TrainEvents.push_back(Events[i]);
    }
    for(int i=0;i<trk_label.size();i++){
        if (find(TrainEvents.begin(), TrainEvents.end(), Event_ID[i]) != TrainEvents.end() ){
            isTrain[i] = 1;
        }
        else isTest[i] = 1;
    }

    // Now Write Required information back into root file
    TFile output("output.root","recreate");
    TTree *t2 = new TTree("Data", "Balanced and Split Dataset");
    TTree *t3 = new TTree("train", "Balanced and Split Dataset for Training");
    TTree *t4 = new TTree("test", "Balanced and Split Dataset for Testing");

    int Keep;
    int Train;
    int Test;

    t2->Branch("features", &features);
    t2->Branch("jet_label", &j_label);
    t2->Branch("trk_label", &t_label);
    t2->Branch("trk_ID", &t_ID);
    t2->Branch("jet_ID", &j_ID);
    t2->Branch("Event_ID", &E_ID);
    t2->Branch("isKept", &Keep);
    t2->Branch("isTrain", &Train);
    t2->Branch("isTest", &Test);

    t3->Branch("features", &features);
    t3->Branch("jet_label", &j_label);
    t3->Branch("trk_label", &t_label);
    t3->Branch("trk_ID", &t_ID);
    t3->Branch("jet_ID", &j_ID);
    t3->Branch("Event_ID", &E_ID);

    t4->Branch("features", &features);
    t4->Branch("jet_label", &j_label);
    t4->Branch("trk_label", &t_label);
    t4->Branch("trk_ID", &t_ID);
    t4->Branch("jet_ID", &j_ID);
    t4->Branch("Event_ID", &E_ID);
    
    for(int i=0;i<trk_label.size();i++){
        tree->GetEntry(i);
        Keep=isKept[i];
        Train=isTrain[i];
        Test=isTest[i];
        t2->Fill();
        // Remap labels for test and train trees
        if(trk_label[i]==-999) continue;
        if(trk_label[i]==0) t_label=1;
        else t_label=0;
        if(jet_label[i]!=-1) j_label=0;
        else j_label=1;
        if(Train==1&&Keep==1) t3->Fill();
        if(Test==1&&Keep==1) t4->Fill();
    }

    t2->Write();
    t3->Write();
    t4->Write();

    Jet_Score_h->Write();

    ///////////////////////////////////////////
    /// The following is used for debugging ///
    ///////////////////////////////////////////

    // Recalulate Jet Score After Balancing
    E_start_idx=0;
    E_end_idx=0;
    for(int i=0;i<num_Events;i++){
        E_start_idx = E_end_idx;
        E_current_ID = Event_ID[E_start_idx];
        J.clear(); 
        Jets.clear();
        Jet_Score.clear();

        // Figure out the Event_ending_idx using a while loop
        while(E_current_ID == Event_ID[E_end_idx]){
            E_end_idx++;
        }

        // Figure out how mant jets per event
        for(int j=E_start_idx;j<E_end_idx;j++){
            J.insert(jet_ID[j]);
        }
        Jets.assign(J.begin(),J.end());
        for(int k=0;k<Jets.size();k++){
            num_trk=0;                
            Score=0;                 
            for(int j=E_start_idx;j<E_end_idx;j++){
                if(jet_ID[j]==Jets[k]&&isKept[j]==1){           
                    num_trk++;                   
                    if(jet_label[j]==-1) Score++;
                }
            }
            Jet_Score_AFTER_h->Fill(float(Score)/float(num_trk));
        }
    }

    int score_before_balance=0;
    int score_after_balance=0;
    int num_kept=0;
    for(int i=0;i<trk_label.size();i++){
        if(isKept[i]==1) num_kept++;
        if(trk_label[i]==0){
            score_before_balance++;
            if(isKept[i]==1){
                score_after_balance++;
            }
        }
    }
    cout << "Score Before Balance: " << (float) score_before_balance / (float)trk_label.size() << endl;
    cout << "Score After Balance: " << (float) score_after_balance / (float)num_kept << endl;

    cout << "Num trks Before Balance: " << trk_label.size() << endl;
    cout << "Num trks After Balance: " << (float)num_kept << endl;
    cout << "Reduction Factor: " << (float)num_kept / (float)trk_label.size()<< endl;

    Jet_Score_AFTER_h->Write();
    cout << "Train " << float(reduce(isTrain.begin(),isTrain.end())) / (float) isTrain.size() << endl;
    cout << "Test " << (float) reduce(isTest.begin(),isTest.end()) / (float) isTest.size() << endl;
    cout << "Kept " << (float) reduce(isKept.begin(),isKept.end()) / (float) isKept.size() << endl;
}

// TO-DO:
// -Convert to zero based indices
