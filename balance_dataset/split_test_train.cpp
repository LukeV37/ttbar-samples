#include "get_data.cpp"
#include <vector>
#include <time.h>
#include <iostream>
#include <set>

using namespace std;

int main(){
  // get total number of lines in file
  char FILENAME[] = "/data/lvaughan/PU/ZZ4nu_250k_output.txt";
  int nlines = get_nlines(FILENAME);  
  //printf("Total number of lines: %d\n", nlines);

  //initlize data
  int maxjets = 100;
  int events = int(nlines / maxjets);
  events = 10;
  nlines = events * maxjets;

  vector<vector<vector<string> > > data(events, vector<vector<string> >(maxjets, vector<string>(4)));
  data = get_data(FILENAME, nlines, maxjets);
  //cout << "Data Got" << endl;
  


  // balance dataset by cutting subgraphs
  // loop through each event and store vertex values
  // generate random number to see if vertex is cut
  srand( (unsigned)time( NULL ) );
  float threshold = 0.7;
  float rng = 0;
  set<string> s;
  vector<string> labels;
  string label;


  for (int i = 0; i < events; i++){
    labels.clear();
    for (int j = 0; j < maxjets; j++){
      s.insert(data[i][j][3]);
      labels.assign(s.begin(), s.end());
      for (int k = 0; k < labels.size(); k++){
        label = labels[k];
        if(label=="-1") continue;
        rng = rand() / RAND_MAX;
        cout << rng << endl;
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
  }


  for (int i = 0; i < events; i++){
    for (int j = 0; j < maxjets; j++){
        cout << data[i][j][0] <<" "<< data[i][j][1] <<" "<< data[i][j][2] <<" "<< data[i][j][3] << endl;
    }
  }

  

  //test train split
  

  //edges will be implemented in separate script
}
