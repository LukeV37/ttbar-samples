#include "get_data.cpp"
#include <vector>
#include <time.h>
#include <iostream>
#include <set>
#include <cstdlib>
#include <string>

using namespace std;

int main(){
  char FILENAME[] = "../../ZZ4nu_250k_output.txt";
  int nvars = 4;

  // get total number of lines in file
  int nlines = get_nlines(FILENAME);  
  //printf("Total number of lines: %d\n", nlines);

  //initlize data
  int maxjets = 100;
  int events = int(nlines / maxjets);

  //limit total number of events for debugging
  //events = 10;
  //nlines = events * maxjets;

  vector<vector<vector<string> > > data(events, vector<vector<string> >(maxjets, vector<string>(nvars)));
  data = get_data(FILENAME, nlines, maxjets, nvars);
  //cout << "Data Got" << endl;
  

  // balance dataset by cutting subgraphs
  // loop through each event and store vertex values
  // generate random number to see if vertex is cut
  srand( time( 0 ) );
  float threshold = 0.7;
  float rng = 0;
  set<string> s;
  vector<string> labels;
  string label;
  int unique_jet_ID = 1;


  for (int i = 0; i < events; i++){
    labels.clear();
    for (int j = 0; j < maxjets; j++){
      s.insert(data[i][j][3]);
      data[i][j].push_back(to_string(unique_jet_ID++));
    }
    labels.assign(s.begin(), s.end());
    //cout << "Label: " << labels[0] << endl;
    for (int k = 0; k < labels.size(); k++){
      label = labels[k];
      if(label=="-1") continue;
      rng = ((float) rand() / RAND_MAX);
      //cout << rng << endl;
      if (rng<threshold){
        for(int n=0;n<maxjets;n++){
          for(int m=0;m<4;m++){
            if(data[i][n][3]==label)
              data[i][n][m] = "0";
          }  
        }
      }
    }
  }

/*
  for (int i = 0; i < events; i++){
    for (int j = 0; j < maxjets; j++){
      cout << data[i][j][0] <<" "<< data[i][j][1] <<" "<< data[i][j][2] <<" "<< data[i][j][3] <<" "<< data[i][j][4] << endl;
    }
  }
*/
  

  //test train split
  float split = 0.75;
  int split_idx = int(events * split);

  ofstream train_file("train.txt");
  ofstream test_file("test.txt");

  for (int i=0;i<split_idx;i++){
    for (int j=0;j<maxjets;j++){
      train_file << data[i][j][0] <<" "<< data[i][j][1] <<" "<< data[i][j][2] <<" "<< data[i][j][3] <<" "<< data[i][j][4] << endl;
    }
  }      
  
  for (int i=split_idx;i<events;i++){
    for (int j=0;j<maxjets;j++){
      test_file << data[i][j][0] <<" "<< data[i][j][1] <<" "<< data[i][j][2] <<" "<< data[i][j][3] <<" "<< data[i][j][4] << endl;
    }
  }      

  cout << "Success: train.txt and test.txt have been written to disk." << endl;

  //edges will be implemented in separate script
}
