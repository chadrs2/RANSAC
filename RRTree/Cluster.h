#pragma once

#include "RRTree.h"

class Cluster {
public:
  Cluster(int userM) {
    M = userM;
    m = 1;
  };
  Cluster(int userM, int userm) {
    M = userM;
    m = userm;
  };
  ~Cluster() {
    for (unsigned int i=0; i < cluster.size(); i++) {
      cluster.at(i)->Clear();
    }
  };
  void AddPointToCluster(Point userPt) {
    unsigned int ptTStep = userPt.GetTStep();
    if (ptTStep >= cluster.size()) {
      RRTree* tempTree = new RRTree(M,m);
      tempTree->InsertData(userPt);
      cluster.push_back(tempTree);
    }
    else {
      cluster.at(ptTStep)->InsertData(userPt);
    }
  };
  bool IsClusterEmpty() {return cluster.empty();};
  vector<RRTree*> GetCluster() {return cluster;};
  void PrintCluster() {
    for (unsigned int i=0; i < cluster.size(); i++) {
      std::cout << "At Time " << i << ":" << std::endl;
      cluster.at(i)->PrintTree();
    }
  }
private:
  vector<RRTree*> cluster;
  int M;
  int m;
};
