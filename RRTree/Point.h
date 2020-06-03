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
    data = NULL;
  };
  Point(vector<float> userData, int userTStep) {
    ClID = "UNCLASSIFIED";
    tStep = userTStep;
    data = userData;
  };
  ~Point(){};
  void TransformData(vector<float> userTransform) {
    if (userTransform.size() == data.size()) {
        data = userTransform * data;
    }
    else {
      std::cout << "Transform has different dimensions than data vector" << std::endl;
    }
  };
  // Setters/Adders
  void SetClID(string userClID) {ClID=userClID;};
  void SetTStep(int userTStep) {tStep=userTStep;};
  void AddData(float userData) {data.push_back(userData);};
  // Getters
  string GetClID() {return ClID;};
  int GetTStep() {return tStep;}
  vector<float> GetData() {return data;};
private:
  string ClID;
  int tStep;
  vector<float> data;
};
