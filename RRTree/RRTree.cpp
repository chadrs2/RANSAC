#include "RRTree.h"
#include <cmath>
#include <map>
#include <set>
using std::set;
using std::map;

RRTree::RRTree() {
  M=2;
  m=1;
  root = NULL;
};
RRTree::RRTree(int userM, int userm) {
  M=userM;
  m=userm;
  root = NULL;
};
RRTree::~RRTree() {clear();}
Node* RRTree::GetRootNode() const {return root;};

//True if data was able to be inserted -> Calls Insert
bool RRTree::InsertData(Point newPoint) {
  return Insert(GetRootNode(), newPoint);
  return false;
};

//True if data was able to be inserted.
//Keeps nodes balanced and within M and m parameters -> Calls ChooseSubtree, OverflowTreatment
bool RRTree::Insert(Node *&tree, Point newPoint) {
  Node* node2InsertPt = ChooseSubtree(tree, newPoint);
  if (node2InsertPt->GetNumberOfChildren() < M) {
    node2InsertPt->AddChildNode(Node(newPoint));
    updateNodeBoundingBoxes(tree);
    return true;
  }
  else {
    node2InsertPt->AddChildNode(Node(newPoint));
    bool wasSplitPerformed = OverflowTreatment(node2InsertPt);
    if (wasSplitPerformed && (node2InsertPt->GetNodeLevel() == 0)) {
      //TODO: create new root
    }
    elseif (wasSplitPerformed) {
      //TODO: Propagate OverflowTreatment upwards??
    }
    updateNodeBoundingBoxes(tree);
    return true;
  }
  return false;
};

//Recursive function moving through tree until it reaches a leaf node -> Calls itself
Node* RRTree::ChooseSubtree(Node *&tree, Point newPoint) {
  if (tree->IsNodeALeafNode()) {
    return tree;
  }
  else {
    if (tree->GetNodeLevel() == 1) {//children are leaf nodes
      // TODO: Determine minimum overlap costs - ties: smallest area change
      ChooseSubtree(tree,newPoint);
    }
    else {
      // TODO: Determine minimum area costs - ties: smallest area
      ChooseSubtree(tree,newPoint);
    }
  }
};

//True if Split() was called.
//Meant for dealing with a filled Node -> Calls ReInsert or Split depending on nodeLevel
bool RRTree::OverflowTreatment(Node *&tree) {
  if (!(tree->WasOverflowTreatmentCalled()) && tree->GetNodeLevel() != 0) {
    ReInsert(tree);
  }
  else {
    Split(tree);
  }
};

//Used to help avoid Split of nodes as often as possible....I'm a bit confused -> Calls Insert
void RRTree::ReInsert(Node *&tree) {
  map<int, Node*, greater<int>> distanceBetweenNodeCenters;
  vector<Node*> childrenNodes = tree->GetNodeChildren();
  vector<float> treeCenterCoords = tree->GetCenterOfBoundingBox();
  for (unsigned int ii=0; ii<tree->GetNumberOfChildren(); ii++) {
    vector<float> currCenterCoords = childrenNodes.at(ii)->GetCenterOfBoundingBox();
    float dist=0.0;
    for (unsigned int jj=0; jj<currCenterCoords.size(); jj++) {
      dist += pow((treeCenterCoords.at(jj) - currCenterCoords.at(jj)),2);
    }
    distanceBetweenNodeCenters.insert(std::pair<int, Node*, greater<int>>(sqrt(dist),childrenNodes.at(ii)));
  }
  int p = floor(.3 * M);
  int i = 0;
  map<int, Node*, greater<int>>::iterator map_itr;
  for (map_itr = distanceBetweenNodeCenters.begin(); map_itr != distanceBetweenNodeCenters.end(); map_itr++) {
    i++;
    if (i <= p) {
      distanceBetweenNodeCenters.erase(map_itr);
      // TODO: Does this actually get rid of child node??
    }
  }
  tree->UpdateNodeBoundingBox();
  // TODO: R14 in paper ... call Insert with removed Nodes
};

//Finds an even distribution of children Nodes and splits Node -> Calls ChooseSplitAxis, ChooseSplitIndex
void RRTree::Split(Node *&tree) {
  int axisIDX = ChooseSplitAxis(tree);
  ChooseSplitIndex(tree, axisIDX);
};

//Finds axis to split on and returns axis index
int RRTree::ChooseSplitAxis(Node *&tree) {
  vector<Node*> currChildrenNodes = tree->GetNodeChildren();
  map<float, Node*, std::less<float>> sortedMapOfNodes;
  map<float, Node*, std::less<float>>::iterator itr;
  vector<Node*> sortedNodes;
  for (unsigned int axs=0; axs < tree->nodeBoundingBox.size(); axs++) {
    sortedMapOfNodes = {};
    sortedNodes = {};
    for (unsigned int entry; entry < currChildrenNodes.size(); entry++) {
      vector<vector<float>> currBB = currChildrenNodes.at(entry)->GetNodeBounds();
      sortedMapOfNodes.insert({currBB.at(axs).at(0),currChildrenNodes.at(entry)});
      // TODO: What happens if minimum value of bounding box is 0
    }
    for (itr=sortedMapOfNodes.begin(); itr != sortedMapOfNodes.end(); itr++) {
      sortedNodes.push_back(itr->second);
      sortedMapOfNodes.erase(itr);
    }
    vector<float> marginValues = {};
    for (unsigned int k=1; k <= (M-2*m+2); k++) {
      vector<Node*> grp1(sortedNodes.begin(),sortedNodes.begin()+k);
      vector<Node*> grp2(sortedNodes.begin()+k,sortedNodes.end());
      if (marginValues.empty()) {
          marginValues.push_back(CalculateMarginValue(grp1,grp2));
      }
      else {
        marginValues.at(axs) = marginValues.at(axs) + CalculateMarginValue(grp1,grp2);
      }
    }
  }
  // Return minimum margin value's index (b/c that idx correlates to the axis)
  return (std::min_element(marginValues.begin(), marginValues.end())-marginValues.begin());
};

//Finds index in which to divide children nodes at and adjusts tree pointer
void RRTree::ChooseSplitIndex(Node *&tree, int axsIDX) {
  vector<Node*> currChildrenNodes = tree->GetNodeChildren();
  map<float, Node*, std::less<float>> sortedMapOfNodes;
  map<float, Node*, std::less<float>>::iterator itr;
  vector<Node*> sortedNodes;
  for (unsigned int entry; entry < currChildrenNodes.size(); entry++) {
    vector<vector<float>> currBB = currChildrenNodes.at(entry)->GetNodeBounds();
    sortedMapOfNodes.insert({currBB.at(axsIDX).at(0),currChildrenNodes.at(entry)});
    // TODO: What happens if minimum value of bounding box is 0
  }
  for (itr=sortedMapOfNodes.begin(); itr != sortedMapOfNodes.end(); itr++) {
    sortedNodes.push_back(itr->second);
    sortedMapOfNodes.erase(itr);
  }
  vector<float> overlapValues;
  vector<float> areaValues;
  for (unsigned int k=1; k <= (M-2*m+2); k++) {
    vector<Node*> grp1(sortedNodes.begin(),sortedNodes.begin()+k);
    vector<Node*> grp2(sortedNodes.begin()+k,sortedNodes.end());
    vector<vector<float>> grp1BB = TempBB(grp1);
    vector<vector<float>> grp2BB = TempBB(grp2);
    overlapValues.push_back(CalculateOverlapValue(grp1BB,grp2BB));
    areaValues.push_back(CalculateAreaValue(grp1BB,grp2BB));
  }
  float currMin = *std::min_element(overlapValues.begin(), overlapValues.end());
  int currMinIDX = std::min_element(overlapValues.begin(), overlapValues.end())-overlapValues.begin();
  for (unsigned int i=0; i<overlapValues.size(); i++) {
    if ((i != currMinIDX) && (currMin == overlapValues.at(i))) {
      if (areaValues.at(currMinIDX) > areaValues.at(i)) {
        currMinIDX = i;
      }
    }
  }
  Node* grp1Node = new Node();
  vector<Node*> newGrp1(sortedNodes.begin(),sortedNodes.begin()+currMinIDX+1);
  for (unsigned int i=0; i < newGrp1.size(); i++) {
    grp1Node->AddChildNode(newGrp1.at(i));
  }
  Node* grp2Node = new Node();
  vector<Node*> newGrp2(sortedNodes.begin()+currMinIDX+1,sortedNodes.end());
  for (unsigned int i=0; i < newGrp2.size(); i++) {
    grp2Node->AddChildNode(newGrp2.at(i));
  }
  tree->RemoveAllChildNodes();
  grp1Node->UpdateNodeBoundingBox();
  grp2Node->UpdateNodeBoundingBox();
  tree->AddChildNode(grp1Node);
  tree->AddChildNode(grp2Node);
  // TODO: DO THIS WORK??
};

//Removes point from leaf node -> Calls ChooseSubtree
void RRTree::Remove(Point pt2Remove);

//Transforms point data
void RRTree::Transform(vector<double> userTransform, Point pt2Transform);

//Builds and returns a temporary bounding box around list of Nodes
vector<vector<float>> RRTree::TempBB(vector<Node*> userGroup) {
  vector<vector<float>> grpBB;
  vector<vector<float>> initBounds = userGroup.at(0)->GetNodeBounds();
  for (unsigned int dim=0; dim < initBounds.size(); dim++) {
    float min=0;
    float max=0;
    for (unsigned int i=0; i < userGroup.size(); i++) {
      vector<vector<float>> currBounds = userGroup.at(i)->GetNodeBounds();
      if (i==0) {
        min = *std::min_element(currBounds.at(dim).begin(),currBounds.at(dim).end());
        max = *std::max_element(currBounds.at(dim).begin(),currBounds.at(dim).end());
      }
      else {
        if(*std::min_element(currBounds.at(dim).begin(),currBounds.at(dim).end()) < min) {
          min = *std::min_element(currBounds.at(dim).begin(),currBounds.at(dim).end());
        }
        if (*std::max_element(currBounds.at(dim).begin(),currBounds.at(dim).end()) > max) {
          max = *std::max_element(currBounds.at(dim).begin(),currBounds.at(dim).end());
        }
      }
    }
    vector<float> boundsVect={};
    boundsVect.push_back(min);
    boundsVect.push_back(max);
    grpBB.push_back(boundsVect);
  }
  return grpBB;
};

//Calculate & Return Area-Value
float RRTree::CalculateAreaValue(vector<vector<float>> grp1BB, vector<vector<float>> grp2BB) {
  float area = 0.0;
  float tempArea = 1;
  for (unsigned int dim=0; dim<grp1BB.size(); dim++) {
    tempArea *= (grp1BB.at(dim).at(1) - grp1BB.at(dim).at(0));
  }
  area = tempArea;
  tempArea = 1;
  for (unsigned int dim=0; dim<grp2BB.size(); dim++) {
    tempArea *= (grp2BB.at(dim).at(1) - grp2BB.at(dim).at(0));
  }
  area += tempArea;
  return area;
};

//Calculate & Return Margin-Value
float RRTree::CalculateMarginValue(vector<Node*> userNodeGroup1, vector<Node*> userNodeGroup2) {
  float margin = 0.0;
  vector<vector<float>> grp1BB = TempBB(userNodeGroup1);
  vector<vector<float>> grp2BB = TempBB(userNodeGroup2);
  for (unsigned int dim=0; dim<grp1BB.size(); dim++) {
    margin += (grp1BB.at(dim).at(1) - grp1BB.at(dim).at(0));
  }
  for (unsigned int dim=0; dim<grp2BB.size(); dim++) {
    margin += (grp2BB.at(dim).at(1) - grp2BB.at(dim).at(0));
  }
  return margin;
};

//Calculate & Return Overlap-Value
float RRTree::CalculateOverlapValue(vector<vector<float>> grp1BB, vector<vector<float>> grp2BB) {
  float overlap = 1.0;
  bool isOverlap = true;
  for (unsigned int dim=0; dim < grp1BB.size(); dim++) {
    if (grp1BB.at(dim).at(1) < grp2BB.at(dim).at(0) || grp2BB.at(dim).at(1) < grp1BB.at(dim).at(0)) {
      isOverlap = false;
      break;
    }
  }
  if (isOverlap) {
    for (unsigned int dim=0; dim < grp1BB.size(); dim++) {
      if (grp1BB.at(dim).at(1) < grp2BB.at(dim).at(1)) {
        float maxPt = grp1BB.at(dim).at(1);
      }
      else {
        float maxPt = grp2BB.at(dim).at(1);
      }
      if (grp1BB.at(dim).at(0) > grp2BB.at(dim).at(0)) {
        float minPt = grp1BB.at(dim).at(0);
      }
      else {
        float minPt = grp2BB.at(dim).at(0);
      }
      overlap *= (maxPt-minPt);
    }
    return overlap;
  }
  else {
    return 0.0;
  }
};

//Empties RR*-tree
void RRTree::clear();

//Starts at root and updates all nodes bounding boxes
void RRTree::updateNodeBoundingBoxes(Node *&tree);
