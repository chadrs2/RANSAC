#pragma once

#include<iostream>
#include<vector>
#include<string>
using std::string;
using std::vector;

class Point {
public:
  Point() {
    ClID = "UNCLASSIFIED";
    tStep = 0;
  };
  Point(vector<float> userData, int userTStep) {
    ClID = "UNCLASSIFIED";
    tStep = userTStep;
    data = userData;
  };
  ~Point(){};
  void TransformData(vector<float> userTransform) {
    if (userTransform.size() == data.size()) {
      for (unsigned int i=0; i < data.size(); i++) {
        data.at(i) *= userTransform.at(i);
      }
    }
    else {
      std::cout << "Transform has different dimensions than data vector" << std::endl;
    }
  };
  // Setters/Adders
  void SetClID(string userClID) {ClID=userClID;};
  void SetTStep(int userTStep) {tStep=userTStep;};
  void AddData(float userData) {data.push_back(userData);};
  void AddAllDimensionData(vector<float> userData) {data = userData;};
  // Getters
  string GetClID() {return ClID;};
  int GetTStep() {return tStep;}
  vector<float> GetData() {return data;};
  int GetPointDimensions() {return data.size();}
private:
  string ClID;
  int tStep;
  vector<float> data;
};
