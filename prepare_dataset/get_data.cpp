#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <string.h>

using namespace std;

int get_nlines(char filename[]){
  FILE* fp = fopen(filename, "r");
  if (fp == NULL)
      exit(EXIT_FAILURE);

  int nlines = 0;

  char* line = NULL;
  size_t len = 0;

  while ((getline(&line, &len, fp)) != -1) {
    // using printf() in all tests for consistency
    //printf("%s", line);
    nlines++;
  }
  fclose(fp);
  if (line)
    free(line);
  return nlines;
}

vector<string> get_line_data(string line_data){
  vector<string> nums;
  int N = line_data.length();
  int idx1 =0;
  int idx2 = 0;
  //int n=0;
  //cout << line_data << endl;

  for (int i=0; i<N; i++){
    if ((line_data[i]==' ')){
      idx2 = i;
      //cout << "N " << N << " idx1 " << idx1 << " idx2 " << idx2 << " n " << n << endl;
      nums.push_back(line_data.substr(idx1, idx2-idx1));
      idx1 = idx2+1;
      //n++;
    }
  }
  //cout << "N " << N << " idx1 " << idx1 << " idx2 " << idx2 << " n " << n << endl;
  nums.push_back(line_data.substr(idx1,N-idx1));
  return nums;
}

vector<vector<vector<string> > > get_data(char filename[], int nlines, int maxjets, int nvars){
  fstream myfile;
 
  // open file
  myfile.open(filename);
  //vector<string> g1(nlines);
  vector<string> g1;

  if (myfile.is_open()) {
    // checking if the file is open
    string str;
    int i=0;
 
    // read data from file object
    // and put it into string.
    //while (getline(myfile, str, ' ')) {
    while (getline(myfile, str)) {
      //g1[i++]=str;  
      g1.push_back(str);
    }
    // close the file object.
    myfile.close();
 }

  /* 
  cout << "\nVector elements are: " << endl;
  for (int i = 0; i < g1.size(); i++) {
        cout << g1[i] << endl;
   }
  */

  int events = int(nlines / maxjets);
  vector<vector<vector<string> > > data(events, vector<vector<string> >(maxjets, vector<string>(nvars)));

  vector<string> line_data(nvars);
  int idx;

  for (int i = 0; i < events; i++){
    for (int j = 0; j < maxjets; j++){
      line_data = get_line_data(g1[i*maxjets + j]);
      for (int k = 0; k < nvars; k++){
          data[i][j][k] = line_data[k];
          //data[i][j][k] = g1[k+j*4+i*maxjets*4];
        }
        //cout << line_data[0] <<" "<< line_data[1] <<" "<< line_data[2] <<" "<< line_data[3] << endl;
    }
  }

  return data;
}
