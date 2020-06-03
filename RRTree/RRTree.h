#pragma once

#include "Point.h"

class RRTree {
public:
  struct Node {
    Node(){
      nodeLevel=0;
      pointObj=NULL;
      overflowTreatmentWasCalled=false;
    };
    Node(Point newPoint) {
      nodeLevel=0;
      pointObj=newPoint;
      overflowTreatmentWasCalled=false;
    };
    ~Node(){};
    void AddChildNode(Node userNode) {childrenNodes.push_back(new Node(userNode));};
    void RemoveChildNode(Node userChildNode2Remove) {
      childrenNodes.erase(std::remove(childrenNodes.begin(), childrenNodes.end(), userChildNode2Remove), childrenNodes.end());
    };
    void AddNodePoint(vector<float> userData, int userTStep) {pointObj = new Point(userData, userTStep);};
    void UpdateNodeBoundingBox() {
      vector<float> mins, maxs;
      //For each childNode
      for (unsigned int ii=0; ii<childrenNodes.size(); ii++) {
        //vector of dimensions where each dimension is a vector of min and max values
        vector<vector<float>> currNodeBoundingBox = childrenNodes.at(ii).GetNodeBounds();
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
      vector<vector<float>> newNodeBoundingBox;
      for (unsigned int ii=0; ii<mins.size(); ii++) {
        vector<float> dimBounds;
        dimBounds.push_back(mins.at(ii));
        dimBounds.push_back(maxs.at(ii));
        newNodeBoundingBox.push_back(dimBounds);
      }
      nodeBoundingBox = newNodeBoundingBox;
    };
    // Getter functions
    int GetNodeLevel() {return nodeLevel;};
    void IncreaseNodeLevel() {nodeLevel++;};
    void DecreaseNodeLevel() {nodeLevel--;};
    vector<Node*> GetNodeChildren() {return childrenNodes;};
    vector<vector<float>> GetNodeBounds() {return nodeBoundingBox;};
    Point* GetNodePoint() {return pointObj;};
    bool WasOverflowTreatmentCalled() {return overflowTreatmentWasCalled;};


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
    int nodeLevel;
    vector<vector<float>> nodeBoundingBox;
    Point* pointObj;
    bool overflowTreatmentWasCalled;
  };

  RRTree(); //Assumes M=2 and m=1
  RRTree(int userM, int userm);
  ~RRTree(); //Calls clear()
  Node* getRootNode() const;

  //True if data was able to be inserted -> Calls Insert()
  bool InsertData(Point newPoint);

  //True if data was able to be inserted -> Calls ChooseSubtree, OverflowTreatment
  bool Insert(Node *&tree, Point newPoint);

  //True if Split() was called.
  //Meant for dealing with a full Node -> Calls ReInsert()/Split() depending on nodeLevel
  bool OverflowTreatment(Node *&tree, int currNodeLevel);

  //Used to help avoid Split() of nodes as often as possible.
  void ReInsert(Node *&tree);

  //
  void ChooseSubtree(Node *&tree, Point newPoint);

  //
  void CheckComp();

  //
  void Split();

  //
  void ChooseSplitAxis();

  //
  void ChooseSplitIndex();

  //
  void remove(Point pt2Remove);

  //
  void Transform(vector<double> userTransform, Point pt2Transform);

  //
  void clear();

  //
  void updateNodeBoundingBoxes(Node *&tree);
private:
  Node *root;
  int M; //Max number of nodes per level
  int m; //minimum number of nodes per level
}
