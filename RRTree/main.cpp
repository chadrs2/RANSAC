#include<iostream>
#include <stdlib.h>
#include "Point.h"
#include "RRTree.h"
using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
  unsigned int numPts = 5;
  if (argc > 1) {
    int M = *argv[1] - '0'; // *argv[0]=./run_file; *argv[1]=M; *argv[2]=m
    int m = *argv[2] - '0';
    //std::cout << M << " " << m << std::endl;
    RRTree* rrTree = new RRTree(M,m);
    for (unsigned int i=0; i < numPts; i++) {
      Point dataPoint;
      dataPoint.SetTStep((int)rand()%10);
      vector<float> data{(float)rand(),(float)rand(),(float)rand()};
      dataPoint.AddAllDimensionData(data);
      rrTree->InsertData(dataPoint);
      rrTree->printTree();
      /*cout << "User Point:" << endl;
      cout << "Cluster ID = " << dataPoint.GetClID() << endl;
      cout << "Time step = " << dataPoint.GetTStep() << endl;
      vector<float> newData = dataPoint.GetData();
      for (unsigned int i=0; i < newData.size(); i++) {
        cout << "Data in the " << i+1 << " dimension = " << newData.at(i) << endl;
      }*/
    }
    //rrTree->printTree();
    /*vector<float> transform{3,4.6,1.2};
    dataPoint.TransformData(transform);
    vector<float> transformedData = dataPoint.GetData();
    for (unsigned int i=0; i < transformedData.size(); i++) {
      cout << "Transformed Data in the " << i+1 << " dimension = " << transformedData.at(i) << endl;
    }*/
    delete rrTree;
  }
  else {
    RRTree* rrTree = new RRTree();
    for (unsigned int i=0; i < numPts; i++) {
      Point dataPoint;
      dataPoint.SetTStep((int)rand()%10);
      vector<float> data{(float)rand(),(float)rand(),(float)rand()};
      dataPoint.AddAllDimensionData(data);
      rrTree->InsertData(dataPoint);
      rrTree->printTree();
        /*cout << "User Point:" << endl;
        cout << "Cluster ID = " << dataPoint.GetClID() << endl;
        cout << "Time step = " << dataPoint.GetTStep() << endl;
        vector<float> newData = dataPoint.GetData();
        for (unsigned int i=0; i < newData.size(); i++) {
          cout << "Data in the " << i+1 << " dimension = " << newData.at(i) << endl;
        }*/
    }
    //rrTree->printTree();
    /*vector<float> transform{3,4.6,1.2};
    dataPoint.TransformData(transform);
    vector<float> transformedData = dataPoint.GetData();
    for (unsigned int i=0; i < transformedData.size(); i++) {
      cout << "Transformed Data in the " << i+1 << " dimension = " << transformedData.at(i) << endl;
    }*/
    delete rrTree;
  }
  return 0;
}
