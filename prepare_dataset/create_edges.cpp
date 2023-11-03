#include "get_data.cpp"
#include <vector>
#include <iostream>
#include <string>

using namespace std;

int main(){
  string run_type = "train";
  char FILENAME[] = "./train.txt"; 
  int nvars = 5;

  // get total number of lines in file
  int nlines = get_nlines(FILENAME);
  //printf("Total number of lines: %d\n", nlines);

  //initlize data
  int maxjets = 100;
  int events = int(nlines / maxjets);

  vector<vector<vector<string> > > data(events, vector<vector<string> >(maxjets, vector<string>(nvars)));
  data = get_data(FILENAME, nlines, maxjets, nvars);
  //cout << "Data Got" << endl;

  //create edges and labels vectors
  struct Edge {
    string source;
	string neighbor;
  };

  struct Edge edge;
  vector<Edge> edges;
  vector<int> labels;

  int jet1_idx = 0, jet2_idx = 0;

  for (int i=0;i<events;i++){
    for (int jet1_idx=0;jet1_idx<maxjets;jet1_idx++){
      if (data[i][jet1_idx][3] == "-1") labels.push_back(int(1));
	  else labels.push_back(int(0));
      for (int jet2_idx=jet1_idx+1;jet2_idx<maxjets;jet2_idx++);{
	    //int p = data[i][jet1_idx][3].compare(data[i][jet2_idx][3]);
		//cout << data[i][jet1_idx][3] <<" "<< data[i][jet2_idx][3] <<" "<< p <<" "<< labels[labels.size() -1] << endl;
		if (data[i][jet1_idx][3].compare(data[i][jet2_idx][3])==0&&data[i][jet1_idx][3]!="0"){
          edge.source = data[i][jet1_idx][4];
		  edge.neighbor = data[i][jet2_idx][4];
		  edges.push_back(edge);
		}
	  }  
    }
  }
  //test train split
  ofstream data_file(run_type+"_data.txt");
  ofstream edges_file(run_type+"_edges.txt");

  for (int i=0;i<edges.size();i++){
    edges_file << edges[i].source <<" "<< edges[i].neighbor << endl;
  }
  /*
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
  */

  cout << "Success: train.txt and test.txt have been written to disk." << endl;
}
