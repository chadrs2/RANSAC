#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include "Point.h"
#include "RRTree.h"
#include "ConsensusSet.h"
#include "Cluster.h"
using std::cout;
using std::endl;
using std::cin;

vector<vector<float>> FindQueryFromPoint(Point* currPoint, float Eps_t, float Eps_d) {
  vector<vector<float>> query;
  vector<float> ptData = currPoint->GetData();
  for (unsigned int i=0; i < ptData.size(); i++) {
    if (i+1 < ptData.size()) {
      vector<float> currDimData{ptData.at(i)-Eps_d,ptData.at(i)+Eps_d};
      query.push_back(currDimData);
    }
    else {
      vector<float> currDimData{ptData.at(i)-Eps_t,ptData.at(i)+Eps_t};
      query.push_back(currDimData);
    }
  }
  return query;
};
Point ConvertPointerPoint2StdPoint(Point* currPoint) {
    Point convertedPt;
    vector<float> ptData = currPoint->GetData();
    int ptTStep = currPoint->GetTStep();
    string ptClID = currPoint->GetClID();
    convertedPt.AddAllDimensionData(ptData);
    convertedPt.SetTStep(ptTStep);
    convertedPt.SetClID(ptClID);
    return convertedPt;
};

Cluster ExpandCluster(RRTree*& rrTree, Point* currPoint, string clusterID, float Eps_t, float Eps_d, unsigned int MinPts, int M) {
  Cluster tempCS(M,1);
  cout << "1" << endl;
  vector<Point*> seeds = rrTree->EvaluateQuery(FindQueryFromPoint(currPoint,Eps_t,Eps_d));
  if (seeds.size() < MinPts) {
    cout << "2" << endl;
    currPoint->SetClID("NOISE");
    return tempCS;
  }
  else {
    cout << "3" << endl;
    for (unsigned int i=0; i < seeds.size(); i++) {
      seeds.at(i)->SetClID(clusterID);
    }
    cout << "4" << endl;
    //seeds.erase(std::find(seeds.begin(),seeds.end(),currPoint));
    //vector<Point*>::iterator tempItr = std::find(SetOfPoints.begin(),SetOfPoints.end(),currPoint);
    //*tempItr = NULL;
    //Point tempPt = ConvertPointerPoint2StdPoint(currPoint);
    //tempCS.AddPointToCluster(tempPt);
    //rrTree->RemovePoint(tempPt);

    cout << "5" << endl;
    while(!seeds.empty()) {
      cout << "6" << endl;
      Point* currPt = seeds.front();
      cout << "seeds size = " << seeds.size() << endl;
      vector<Point*> result = rrTree->EvaluateQuery(FindQueryFromPoint(currPt,Eps_t,Eps_d));
      if (result.size() >= MinPts) {
        cout << "results size = " << result.size() << endl;
        for (unsigned int i=0; i < result.size(); i++) {
          cout << "7" << endl;
          if (result.at(i)->GetClID() == "UNCLASSIFIED" || result.at(i)->GetClID() == "NOISE") {
            if (result.at(i)->GetClID() == "UNCLASSIFIED") {
              seeds.push_back(result.at(i));
            }
            result.at(i)->SetClID(clusterID);
          }
        }
      }
      cout << "seeds size = " << seeds.size() << endl;
      cout << "8" << endl;
      //vector<Point*>::iterator tempItr = std::find(SetOfPoints.begin(),SetOfPoints.end(),seeds.front());
      //*tempItr = NULL;
      Point tempPt = ConvertPointerPoint2StdPoint(seeds.front());
      tempCS.AddPointToCluster(tempPt);
      cout << "9" << endl;
      seeds.erase(seeds.begin());
      rrTree->PrintTree();
      tempPt.PrintPoint();
      rrTree->RemovePoint(tempPt);
      cout << "Should've removed point" << endl;
      rrTree->PrintTree();
    }
    return tempCS;
  }
};

vector<Cluster> DBSCAN (RRTree*& rrTree, vector<Point*> SetOfPoints, float Eps_t, float Eps_d, int MinPts, int M) {
  vector<Cluster> SetOfClusters;
  string ClID = "1";
  for (unsigned int i=0; i < SetOfPoints.size(); i++) {
    Point* currPoint = SetOfPoints.at(i);
    currPoint->PrintPoint();
    if (currPoint->GetClID() == "UNCLASSIFIED") {
      cout << "0" << endl;
      Cluster currCluster = ExpandCluster(rrTree,currPoint,ClID,Eps_t,Eps_d,MinPts,M);
      currCluster.PrintCluster();
      cout << "10" << endl;
      if (!currCluster.IsClusterEmpty()) {
        cout << "20" << endl;
        SetOfClusters.push_back(currCluster);
        int currentClID = stoi(ClID);
        currentClID++;
        ClID = std::to_string(currentClID);
      }
    }
  }
  return SetOfClusters;
};

int main(int argc, char *argv[]) {
  // Initialization
  srand(time(0));
  unsigned int tWindow;
  cout << "What is your time horizon (number of time steps)? ";
  cin >> tWindow;
  unsigned int numPtsPerTStep;
  cout << "How many data points per time step do you wish to enter? ";
  cin >> numPtsPerTStep;
  int M;
  cout << "Set your M (max nodes per level) parameter (integer) = ";
  cin >> M;

  // Adding point to consensus set, cluster, or tree
  cout << endl << "Inserting " << tWindow*numPtsPerTStep << " points into R*-Tree" << endl << endl;
  RRTree* rrTree = new RRTree(M,1); //m = 1 (minimum number of data points per level)
  //ConsensusSet consensusSet; // Doesn't need parameters (M,1) because it is just a vector of vector Points
  //Cluster cluster(M,1); // Needs (M,1) parameter because it is a vector of RRTrees
  for (unsigned int k=0; k < numPtsPerTStep; k++) {
    for (unsigned int i=0; i < tWindow; i++) {
      Point dataPoint;
      float x = rand() % 50;
      float y = rand() % 50;
      float z = rand() % 50;
      float t = i;//rand() % 10; //Time horizon is 10
      vector<float> data{x,y,z,t}; //4th dimension is the time step
      dataPoint.AddAllDimensionData(data);
      dataPoint.SetTStep(t);
      /*
      if dataPoint is in model
        Add point to ConsensusSet       =>  consensusSet.AddPointToConsensusSet(dataPoint)
      else if dataPoint is in Cluster
        Add point to Cluster            =>  cluster.AddPointToCluster(dataPoint)
      else
        Add point to tree (see line 43) =>  rrTree->InsertData(dataPoint);
        Add point to SetOfPoints vector =>  setOfPoints.push_back(dataPoint);
      */
      rrTree->InsertData(dataPoint);
    }
  }
  vector<Point*> setOfPoints = rrTree->GetAllPointsInTree();
  cout << "Your current tree looks like this:" << endl << endl;
  rrTree->PrintTree();
  cout << endl << endl << endl;

  // Generating Clusters from RRTree
  cout << "--Generating Clusters from RRTree--" << endl << endl;
  float eps_t;
  float eps_d;
  int minPts;
  cout << "Enter your eps_t parameter = ";
  cin >> eps_t;
  cout << "Enter your eps_d parameter = ";
  cin >> eps_d;
  cout << "Enter your minPts parameter = ";
  cin >> minPts;
  vector<Cluster> SetOfClusters = DBSCAN(rrTree, setOfPoints, eps_t, eps_d, minPts, M);
  /*for (unsigned int i=0; i < setOfPoints.size(); i++) {
    delete setOfPoints.at(i);
  }
  setOfPoints.clear();*/
  cout << "Number of clusters = " << SetOfClusters.size() << endl;
  cout << "Printing each cluster:" << endl;
  for (unsigned int i=0; i < SetOfClusters.size(); i++) {
    SetOfClusters.at(i).PrintCluster();
  }
  cout << endl;
  cout << "Your new clustered tree looks like this:" << endl << endl;
  rrTree->PrintTree();
  cout << "100" << endl;
  delete rrTree;
  cout << "200" << endl;
  return 0;
};
