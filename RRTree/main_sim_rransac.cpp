#include<iostream>
#include <stdlib.h>
#include <time.h>
#include "Point.h"
#include "RRTree.h"
#include "ConsensusSet.h"
#include "Cluster.h"
using std::cout;
using std::endl;
using std::cin;

int main(int argc, char *argv[]) {
  srand(time(0));
  unsigned int tWindow;
  cout << "What is your time horizon (number of time steps)? ";
  cin >> tWindow;
  unsigned int numPtsPerTStep;
  cout << "How many data points per time step do you wish to enter? ";
  cin >> numPtsPerTStep;
  char chooseM;
  cout << "Do you wish to set the M (max nodes per level) value? (y/n) ";
  cin >> chooseM;
  if (chooseM == 'y') {
    int M;
    int m = 1;
    cout << "M = ";
    cin >> M;
    cout << endl << "Inserting " << tWindow*numPtsPerTStep << " points into R*-Tree" << endl << endl;
    RRTree* rrTree = new RRTree(M,m);
    ConsensusSet CS;
    Cluster ClustS;
    for (unsigned int k=0; k < numPtsPerTStep; k++) {
      for (unsigned int i=0; i < tWindow; i++) {
        Point dataPoint;
        float x = rand() % 50;
        float y = rand() % 50;
        float z = rand() % 50;
        int t = i;//rand() % 10; //Time horizon is 10
        vector<float> data{x,y,z};
        dataPoint.AddAllDimensionData(data);
        dataPoint.SetTStep(t);
        CS.AddPointToConsensusSet(dataPoint);
        ClustS.AddPointToCluster(dataPoint);
        rrTree->InsertData(dataPoint);
      }
    }
    rrTree->PrintTree();
    cout << endl;
    CS.PrintConsensusSet();
    cout << endl;
    ClustS.PrintCluster();
    delete rrTree;
  }/*
  else {
    cout << endl << "Inserting " << tWindow*numPtsPerTStep << " points into R*-Tree" << endl << endl;
    RRTree* rrTree = new RRTree();
    ConsensusSet CS;
    Cluster ClustS;
    for (unsigned int k=0; k < numPtsPerTStep; k++) {
      for (unsigned int i=0; i < tWindow; i++) {
        Point dataPoint;
        float x = rand() % 50;
        float y = rand() % 50;
        float z = rand() % 50;
        int t = i;//rand() % 10; //Time horizon is 10
        vector<float> data{x,y,z};
        dataPoint.AddAllDimensionData(data);
        dataPoint.SetTStep(t);
        CS.AddPointToConsensusSet(dataPoint);
        ClustS.AddPointToCluster(dataPoint);
        rrTree->InsertData(dataPoint);
      }
    }
    rrTree->PrintTree();
    cout << endl;
    CS.PrintConsensusSet();
    cout << endl;
    ClustS.PrintCluster();
    delete rrTree;
  }*/
  return 0;
}
