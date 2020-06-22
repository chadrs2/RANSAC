#include<iostream>
#include <stdlib.h>
#include <time.h>
#include "Point.h"
#include "RRTree.h"
using std::cout;
using std::endl;
using std::cin;

int main(int argc, char *argv[]) {
  srand(time(0));
  unsigned int numPts = 5;//*argv[1] - '0';
  cout << endl << "Inserting " << numPts << " points into R*-Tree" << endl << endl;
  if (argc > 1) {
    int M = *argv[1] - '0'; // *argv[0]=./run_file; *argv[1]=M; *argv[2]=m
    int m = *argv[2] - '0';
    //std::cout << M << " " << m << std::endl;
    RRTree* rrTree = new RRTree(M,m);
    for (unsigned int i=0; i < numPts; i++) {
      Point dataPoint;
      dataPoint.SetTStep((int)rand()%10);
      float x = rand() % 50;
      float y = rand() % 50;
      float z = rand() % 50;
      vector<float> data{x,y,z};
      dataPoint.AddAllDimensionData(data);
      rrTree->InsertData(dataPoint);
    }
    rrTree->PrintTree();
    cout << endl;
    // Implementing query search
    vector<vector<float>> query{{1,15},{3,30},{44,50}};
    vector<Point*> ptsFoundInQuery = rrTree->EvaluateQuery(query);
    cout << "Found " << ptsFoundInQuery.size() << " points that fall within query" << endl;
    cout << "Query bounds are: {";
    for (unsigned int i=0; i < query.size(); i++) {
      if ((i+1) != query.size()) {
        cout << query.at(i).at(0) << "," << query.at(i).at(1) << "; ";
      }
      else {
        cout << query.at(i).at(0) << "," << query.at(i).at(1) << "}" << endl;
      }
    }
    for (unsigned int i=0; i < ptsFoundInQuery.size(); i++) {
      cout << "Point " << i+1 << ": {";
      vector<float> newData = ptsFoundInQuery.at(i)->GetData();
      for (unsigned int j=0; j < newData.size(); j++) {
        if ((j+1) != newData.size()) {
          cout << newData.at(j) << ",";
        }
        else {
          cout << newData.at(j);
        }
      }
      cout << "}" << endl;
    }
    cout << endl;
    //Remove point
    char userInput;
    cout << "Do you want to remove a data point? (y/n) ";
    cin >> userInput;
    while(userInput=='y') {
      cout << endl;
      Point dataPoint;
      for (unsigned int i=0; i < query.size(); i++) {
        int currInput;
        cout << "Type in " << i+1 << " input = ";
        cin >> currInput;
        dataPoint.AddData(currInput);
      }
      cout << endl;
      rrTree->RemovePoint(dataPoint);
      cout << endl;
      rrTree->PrintTree();
      cout << "Do you want to remove a data point? (y/n) ";
      cin >> userInput;
    }
    delete rrTree;
  }
  else {
    RRTree* rrTree = new RRTree();
    for (unsigned int i=0; i < numPts; i++) {
      Point dataPoint;
      dataPoint.SetTStep((int)rand()%10);
      float x = rand() % 50;
      float y = rand() % 50;
      float z = rand() % 50;
      vector<float> data{x,y,z};
      dataPoint.AddAllDimensionData(data);
      rrTree->InsertData(dataPoint);
    }
    rrTree->PrintTree();
    cout << endl;
    // Implementing query search
    vector<vector<float>> query{{1,15},{3,30},{44,50}};
    vector<Point*> ptsFoundInQuery = rrTree->EvaluateQuery(query);
    cout << "Found " << ptsFoundInQuery.size() << " points that fall within query" << endl;
    cout << "Query bounds are: {";
    for (unsigned int i=0; i < query.size(); i++) {
      if ((i+1) != query.size()) {
        cout << query.at(i).at(0) << "," << query.at(i).at(1) << "; ";
      }
      else {
        cout << query.at(i).at(0) << "," << query.at(i).at(1) << "}" << endl;
      }
    }
    for (unsigned int i=0; i < ptsFoundInQuery.size(); i++) {
      cout << "Point " << i+1 << ": {";
      vector<float> newData = ptsFoundInQuery.at(i)->GetData();
      for (unsigned int j=0; j < newData.size(); j++) {
        if ((j+1) != newData.size()) {
          cout << newData.at(j) << ",";
        }
        else {
          cout << newData.at(j);
        }
      }
      cout << "}" << endl;
    }
    cout << endl;
    //Remove point
    char userInput;
    cout << "Do you want to remove a data point? (y/n) ";
    cin >> userInput;
    while(userInput=='y') {
      cout << endl;
      Point dataPoint;
      for (unsigned int i=0; i < query.size(); i++) {
        int currInput;
        cout << "Type in " << i+1 << " input = ";
        cin >> currInput;
        dataPoint.AddData(currInput);
      }
      cout << endl;
      rrTree->RemovePoint(dataPoint);
      cout << endl;
      rrTree->PrintTree();
      cout << "Do you want to remove a data point? (y/n) ";
      cin >> userInput;
    }
    delete rrTree;
  }
  return 0;
}
