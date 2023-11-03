#include "get_data.cpp"
#include <vector>
#include <iostream>
#include <string>
#include <numeric>

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

  for (int i=0;i<events;i++){
    for (int jet1_idx=0;jet1_idx<maxjets;jet1_idx++){
      if (data[i][jet1_idx][3] == "-1") labels.push_back(int(1));
	  else labels.push_back(int(0));
      for (int jet2_idx=jet1_idx+1;jet2_idx<maxjets;jet2_idx++){
	    //cout << "jet1_idx " << jet1_idx << " jet2_idx " << jet2_idx << endl;
		if (data[i][jet1_idx][3].compare(data[i][jet2_idx][3])==0&&data[i][jet1_idx][0]!="0"){
          edge.source = data[i][jet1_idx][4];
		  edge.neighbor = data[i][jet2_idx][4];
		  edges.push_back(edge);
		  edge.source = data[i][jet2_idx][4];
		  edge.neighbor = data[i][jet1_idx][4];
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

  int stored_jets=0;
  for (int i=0;i<events;i++){
    for (int j=0;j<maxjets;j++){
	  if (data[i][j][0]!="0"){
      data_file << data[i][j][0] <<" "<< data[i][j][1] <<" "<< data[i][j][2] <<" "<< to_string(labels[i*maxjets+j]) <<" "<< data[i][j][4] << endl;
	  stored_jets++;
	  }
    }
  }

  cout << "Success: "+run_type+"_data.txt and "+run_type+"_edges.txt have been written to disk." << endl;
  cout << "Mean of Labels: " << accumulate(labels.begin(), labels.end(), 0.0) / stored_jets << endl;
}
