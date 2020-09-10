#pragma once

#include "Point.h"

class ConsensusSet {
public:
  ConsensusSet() {};
  ~ConsensusSet() {};
  void AddPointToConsensusSet(Point userPt) {
    unsigned int ptTStep = userPt.GetTStep();
    if (ptTStep >= consensusSet.size()) {
      vector<Point> tempPtVect;
      tempPtVect.push_back(userPt);
      consensusSet.push_back(tempPtVect);
    }
    else {
      consensusSet.at(ptTStep).push_back(userPt);
    }
  };
  vector<vector<Point>> GetConsensusSet() {return consensusSet;};
  void PrintConsensusSet() {
    for (unsigned int i=0; i < consensusSet.size(); i++) {
      std::cout << "At Time " << i << ":" << std::endl;
      for (unsigned int j=0; j < consensusSet.at(i).size(); j++) {
        std::cout << "P" << j+1 << ": ";
        consensusSet.at(i).at(j).PrintPoint();
      }
    }
  }
private:
  vector<vector<Point>> consensusSet;
};
