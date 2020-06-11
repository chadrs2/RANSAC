#pragma once

#include "Point.h"

class RRTree {
public:
  struct Node {
    Node(){
      //nodeLevel=0;
      overflowTreatmentWasCalled=false;
      isLeafNode = true;
    };
    ~Node(){
      RemoveAllPointObjs();
    };
    void AddChildNode(Node* userNode) {childrenNodes.push_back(userNode);};
    /*void RemoveChildNode(Node* userChildNode2Remove) {
      childrenNodes.erase(std::remove(childrenNodes.begin(), childrenNodes.end(), userChildNode2Remove), childrenNodes.end());
    };*/
    void RemoveAllChildNodes() {
      childrenNodes.clear();
    }
    void RemoveAllPointObjs() {
      for (unsigned int i=0; i < pointObjs.size(); i++) {
        Point* tempPt = pointObjs.at(i);
        delete tempPt;
      }
      pointObjs.clear();
    }
    void AddNodePoint(vector<float> userData, int userTStep) {
      pointObjs.push_back(new Point(userData, userTStep));
      UpdateNodeBoundingBox();
    };
    void UpdateNodeBoundingBox() {
      vector<float> mins, maxs;
      if (pointObjs.size() > 0) {//If node is a leaf node
        //For each data point
        for (unsigned int ii=0; ii<pointObjs.size(); ii++) {
          //vector of dimensions where each dimension is a vector of min and max values
          vector<float> currPtData = pointObjs.at(ii)->GetData();
          //For each dimension
          for (unsigned int dim=0; dim<currPtData.size(); dim++) {
            if (ii != 0) {
              if (currPtData.at(dim) < mins.at(dim)) {
                mins.at(dim) = currPtData.at(dim);
              }
              if (currPtData.at(dim) > maxs.at(dim)) {
                maxs.at(dim) = currPtData.at(dim);
              }
            }
            else {
              mins.push_back(currPtData.at(dim));
              maxs.push_back(currPtData.at(dim));
            }
          }
        }
      }
      else {
        //For each childNode
        for (unsigned int ii=0; ii<childrenNodes.size(); ii++) {
          //vector of dimensions where each dimension is a vector of min and max values
          vector<vector<float>> currNodeBoundingBox = childrenNodes.at(ii)->GetNodeBounds();
          if (!currNodeBoundingBox.empty()) {
            //For each dimension
            for (unsigned int dim=0; dim<currNodeBoundingBox.size(); dim++) {
              if (ii != 0) {
                if (currNodeBoundingBox.at(dim).at(0) < mins.at(dim)) {
                  mins.at(dim) = currNodeBoundingBox.at(dim).at(0);
                }
                if (currNodeBoundingBox.at(dim).at(1) > maxs.at(dim)) {
                  maxs.at(dim) = currNodeBoundingBox.at(dim).at(1);
                }
              }
              else {
                mins.push_back(currNodeBoundingBox.at(dim).at(0));//0==min
                maxs.push_back(currNodeBoundingBox.at(dim).at(1));//1==max
              }
            }
          }
        }
      }
      vector<vector<float>> newNodeBoundingBox;
      for (unsigned int dim=0; dim<mins.size(); dim++) {
        vector<float> dimBounds;
        dimBounds.push_back(mins.at(dim));
        dimBounds.push_back(maxs.at(dim));
        newNodeBoundingBox.push_back(dimBounds);
      }
      nodeBoundingBox = newNodeBoundingBox;
    };
    void SetOverflowTreatment(bool userOTCalled) {overflowTreatmentWasCalled = userOTCalled;};
    // Getter functions
    //int GetNodeLevel() {return nodeLevel;};
    //void IncreaseNodeLevel() {nodeLevel++;};
    //void DecreaseNodeLevel() {nodeLevel--;};
    vector<Node*> GetNodeChildren() {return childrenNodes;};
    int GetNumberOfChildren() {return childrenNodes.size();};
    vector<vector<float>> GetNodeBounds() {return nodeBoundingBox;};
    vector<Point*> GetPointObjs() {return pointObjs;};
    int GetNumberOfPoinObjs() {return pointObjs.size();};
    bool WasOverflowTreatmentCalled() {return overflowTreatmentWasCalled;};
    bool IsNodeALeafNode() {
      if (pointObjs.size() > 0) {return true;}
      else if (childrenNodes.size() == 0) {return true;}
      else {return false;}
    };
    vector<float> GetCenterOfBoundingBox() {
      vector<float> boxOrigin;
      if (!nodeBoundingBox.empty()) {
        for (unsigned int dim=0; dim<nodeBoundingBox.size(); dim++) {
          boxOrigin.push_back((nodeBoundingBox.at(dim).at(0) + nodeBoundingBox.at(dim).at(1)) / 2);
        }
        return boxOrigin;
      }
      else {
        return boxOrigin;
      }
    };
    /*float nodePerimeter() {
      float perim = 0;
      for (unsigned int i = 0; i<nodeSize.size(); i++) {
      }
      return perim;
      perim += nodeSize.at(i).at(1)-nodeSize.at(i).at(0);
    };
    float nodeVolume() {
      float vol = 1;
      for (unsigned int i = 0; i<nodeSize.size(); i++) {
        vol *= nodeSize.at(i).at(1)-nodeSize.at(i).at(0);
      }
      return vol;
    };*/
    //Variables
    vector<Node*> childrenNodes;
    //int nodeLevel;
    vector<vector<float>> nodeBoundingBox; //1stVect: Dimension; 2ndVect: {min, max}
    vector<Point*> pointObjs;
    bool overflowTreatmentWasCalled;
    bool isLeafNode;
  };

  RRTree(); //Assumes M=2 and m=1
  RRTree(int userM, int userm);
  ~RRTree(); //Calls clear()
  //Empties RR*-tree
  void clear(); // Calls remove
  void removeNodes(Node *&tree); //Recursive function to delete all nodes in tree

  // Returns root node of tree
  //Node* GetRootNode() const;

  //True if data was able to be inserted -> Calls Insert
  void InsertData(Point newPoint);
  //True if data was able to be inserted.
  //Keeps nodes balanced and within M and m parameters -> Calls ChooseSubtree, OverflowTreatment
  void Insert(Node *&tree, Point newPoint);
  //Recursive function moving through tree until it reaches a leaf node -> Calls itself
  //Node* ChooseSubtree(Node *&tree, Point newPoint);
  void ChooseSubtree(Node *&tree, Point newPoint);
  //True if Split() was called.
  //Meant for dealing with a filled Node -> Calls ReInsert or Split depending on nodeLevel
  bool OverflowTreatment(Node *&tree);
  //Used to help avoid Split of nodes as often as possible....I'm a bit confused -> Calls Insert
  void ReInsert(Node *&tree);
  //Finds an even distribution of children Nodes and splits Node -> Calls ChooseSplitAxis, ChooseSplitIndex
  void Split(Node *&tree);
  //Finds axis to split on and returns axis index
  int ChooseSplitAxis(Node *&tree);
  //Finds index in which to divide children nodes at and returns index
  void ChooseSplitIndex(Node *&tree, int axsIDX);
  //Removes point from leaf node -> Calls ChooseSubtree
  void Remove(Node *& tree, Point pt2Remove);
  //Transforms point data
  void Transform(vector<float> userTransform, Point *&pt2Transform);
  //Builds and returns a temporary bounding box around list of Nodes
  vector<vector<float>> TempNodeBB(vector<Node*> userGroup);
  //Builds and returns a temporary bounding box around list of Points
  vector<vector<float>> TempPointBB(vector<Point*> userGroup);
  //Calculate Area-Value
  float CalculateAreaValue(vector<vector<float>> grp1BB, vector<vector<float>> grp2BB);
  //Calculate Margin-Value
  float CalculateMarginValue(vector<vector<float>> grp1BB, vector<vector<float>> grp2BB);
  //Calculate Overlap-Value
  float CalculateOverlapValue(vector<vector<float>> grp1BB, vector<vector<float>> grp2BB);
  //Starts at root and updates all nodes bounding boxes
  void updateNodeBoundingBoxes(Node *&tree);
  void printTree();
  void printEachLayer(vector<Node*> branches);
private:
  Node *root;
  int M; //Max number of nodes per level
  int m; //minimum number of nodes per level
};
