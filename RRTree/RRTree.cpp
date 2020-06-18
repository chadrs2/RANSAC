#include "RRTree.h"
#include <cmath>
#include <map>
#include <set>
#include <algorithm>
using std::pow;
using std::sqrt;
using std::ceil;
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
RRTree::~RRTree() {Clear();}
//Empties RR*-tree
void RRTree::Clear() {
  if (root != NULL) {
    RemoveNodes(root);
  }
};
void RRTree::RemoveNodes(Node *&tree) {
  vector<Node*> currChildrenNodes = tree->GetNodeChildren();
  for(unsigned int i=0; i < currChildrenNodes.size(); i++) {
    RemoveNodes(currChildrenNodes.at(i));
  }
  currChildrenNodes.clear();
  delete tree;
};

//Finds all points inside query
vector<Point*> RRTree::EvaluateQuery(vector<vector<float>> userQuery) {
  vector<Node*> treeRoot{root};
  vector<Node*> foundNodesInQuery = SiblingNodesInQuery(treeRoot, userQuery); // Function in .h file
  vector<Point*> allPtsInFoundNodes = FindAllPointsInParentNodes(foundNodesInQuery);
  return allPtsInFoundNodes;
};
// Finds out if current node touches query bounds
bool RRTree::IsNodeInQuery(vector<vector<float>> userNode, vector<vector<float>> userQuery) {
  bool nodeIsInQuery=true;
  for (unsigned int dim=0; dim < userNode.size(); dim++) {
    if (!(((userNode.at(dim).at(0) > userQuery.at(dim).at(0)) && (userNode.at(dim).at(0) < userQuery.at(dim).at(1))) \
    || ((userNode.at(dim).at(1) < userQuery.at(dim).at(1)) && (userNode.at(dim).at(1) > userQuery.at(dim).at(0))))) {
      if (!(((userQuery.at(dim).at(0) > userNode.at(dim).at(0)) && (userQuery.at(dim).at(0) < userNode.at(dim).at(1))) \
      || ((userQuery.at(dim).at(1) < userNode.at(dim).at(1)) && (userQuery.at(dim).at(1) > userNode.at(dim).at(0))))) {
        nodeIsInQuery=false;
      }
    }
  }
  return nodeIsInQuery;
};
// Finds all points in passed in parent nodes
vector<Point*> RRTree::FindAllPointsInParentNodes(vector<Node*> parentNodes) {
  vector<Point*> allPts;
  for (unsigned int i=0; i < parentNodes.size(); i++) {
    if (parentNodes.at(i)->GetNumberOfChildren() > 0) {
      vector<Node*> currChildrenNodes = parentNodes.at(i)->GetNodeChildren();
      vector<Point*> tempPts = FindAllPointsInParentNodes(currChildrenNodes);
      for (unsigned int j=0; j < tempPts.size(); j++) {
        allPts.push_back(tempPts.at(j));
      }
    }
    else {
      vector<Point*> tempPts = parentNodes.at(i)->GetPointObjs();
      for (unsigned int j=0; j < tempPts.size(); j++) {
        allPts.push_back(tempPts.at(j));
      }
    }
  }
  return allPts;
};

//True if data was able to be inserted -> Calls Insert
void RRTree::InsertData(Point newPoint) {
  //std::cout << "InsertData is called" << std::endl;
  if (root == NULL) {
    root = new Node();
  }
  vector<Node*> treeRoot{root};
  Insert(treeRoot, 0, newPoint);
};

//Keeps nodes balanced and within M and m parameters -> Calls ChooseSubtree, OverflowTreatment
void RRTree::Insert(vector<Node*> &siblings, int childIDX, Point newPoint) {
//void RRTree::Insert(Node *&tree, Point newPoint) {
  if (siblings.at(childIDX) == NULL) {
    siblings.at(childIDX) = new Node();
  }
  // ChooseSubtree Algorithm
  if(!(siblings.at(childIDX)->IsNodeALeafNode())) {
    vector<Node*> currChildren = siblings.at(childIDX)->GetNodeChildren();
    if (currChildren.at(0)->IsNodeALeafNode()) {//children are leaf nodes
      //Determine minimum overlap costs - // TODO: ties: smallest area change
      Node* newNode = new Node();
      newNode->AddNodePoint(newPoint.GetData(),newPoint.GetTStep());
      vector<float> minOverlapCost;
      vector<float> minAreaCost;
      for (unsigned int i=0; i < currChildren.size(); i++) {
        vector<vector<float>> currBB = currChildren.at(i)->GetNodeBounds();
        vector<Node*> nodeAndPointNode; // TODO: Can I just reset vector to NULL or do I have to delete it first -------
        nodeAndPointNode.push_back(currChildren.at(i));
        nodeAndPointNode.push_back(newNode);
        vector<vector<float>> newBB = TempNodeBB(nodeAndPointNode);
        minOverlapCost.push_back(CalculateOverlapValue(currBB,newBB));
        minAreaCost.push_back(CalculateAreaValue(currBB,newBB));
        //Clean Up nodeAndPointNode vector
        nodeAndPointNode.clear();
      }
      int minOVIdx = std::min_element(minOverlapCost.begin(), minOverlapCost.end())-minOverlapCost.begin();
      delete newNode;
      //std::cout << "Enter " << minOVIdx << " branch (if)" << std::endl;
      Insert(currChildren, minOVIdx, newPoint);
      siblings.at(childIDX)->RemoveAllChildNodes();
      for (unsigned int j=0; j < currChildren.size(); j++) {
        siblings.at(childIDX)->AddChildNode(currChildren.at(j));
      }
      //std::cout << "CurrChildren of root: " << currChildren.size() << std::endl;
    }
    else {
      // Determine minimum area costs - // TODO: ties: smallest area
      Node* newNode = new Node();
      newNode->AddNodePoint(newPoint.GetData(),newPoint.GetTStep());
      vector<Node*> currChildren = siblings.at(childIDX)->GetNodeChildren();
      vector<float> minAreaCost;
      for (unsigned int i=0; i < currChildren.size(); i++) {
        vector<vector<float>> currBB = currChildren.at(i)->GetNodeBounds();
        vector<Node*> nodeAndPointNode; // TODO: Can I just reset vector to NULL or do I have to delete it first -------
        nodeAndPointNode.push_back(currChildren.at(i));
        nodeAndPointNode.push_back(newNode);
        vector<vector<float>> newBB = TempNodeBB(nodeAndPointNode);
        minAreaCost.push_back(CalculateAreaValue(currBB,newBB));
        //Clean Up nodeAndPointNode vector
        nodeAndPointNode.clear();
      }
      int minAreaIdx = std::min_element(minAreaCost.begin(), minAreaCost.end())-minAreaCost.begin();
      delete newNode;
      //std::cout << "Enter " << minAreaIdx << " branch (else)" << std::endl;
      Insert(currChildren, minAreaIdx, newPoint);
      siblings.at(childIDX)->RemoveAllChildNodes();
      for (unsigned int j=0; j < currChildren.size(); j++) {
        siblings.at(childIDX)->AddChildNode(currChildren.at(j));
      }
      //std::cout << "CurrChildren of root: " << currChildren.size() << std::endl;
    }
  }
  //Insert Algorithm
  //std::cout << "Number of points in current Node: " << siblings.at(childIDX)->GetNumberOfPoinObjs() << " when M is " << M << std::endl;
  //std::cout << "Number of children nodes in current node: " << siblings.at(childIDX)->GetNumberOfChildren() << std::endl;
  if ((siblings.at(childIDX)->GetNumberOfPoinObjs() < M) && (siblings.at(childIDX)->IsNodeALeafNode())) {
    //std::cout << "Add Point into current child node" << std::endl;
    siblings.at(childIDX)->AddNodePoint(newPoint.GetData(),newPoint.GetTStep());
  }
  else if ((siblings.at(childIDX)->GetNumberOfPoinObjs() == M) || (siblings.at(childIDX)->GetNumberOfChildren() > M)) {
    //std::cout << "Node is either full of nodes or points" << std::endl;
    if ((siblings.at(childIDX)->GetNumberOfPoinObjs() == M) && (siblings.at(childIDX)->IsNodeALeafNode())) {
      //std::cout << "Add point to this node" << std::endl;
      siblings.at(childIDX)->AddNodePoint(newPoint.GetData(),newPoint.GetTStep());
    }
    bool rootWasSplit = false;
    if (siblings.at(childIDX) == root) {
      rootWasSplit = true;
    }
    bool wasSplitPerformed = OverflowTreatment(siblings, childIDX);
    if (wasSplitPerformed && rootWasSplit) {
      //std::cout << "Need to create a new root" << std::endl;
      vector<Node*>::iterator itr = siblings.begin();
      siblings.erase(itr);
      //Node *tempRoot = root;
      for (unsigned int i=0; i < siblings.size(); i++) {
        root->AddChildNode(siblings.at(i));
      }
      root->UpdateNodeBoundingBox();
    }
  }
  updateNodeBoundingBoxes(siblings);
  //siblings.at(childIDX)->UpdateNodeBoundingBox();
};

//True if Split() was called.
//Meant for dealing with a filled Node -> Calls ReInsert or Split depending on nodeLevel
bool RRTree::OverflowTreatment(vector<Node*> &siblings, int childIDX) {
  //std::cout << "Overflow was called" << std::endl;
  if (!(siblings.at(childIDX)->WasOverflowTreatmentCalled()) && (siblings.at(childIDX) != root)) {
    siblings.at(childIDX)->SetOverflowTreatment(true);
    ReInsert(siblings, childIDX);
    return false; // Split wasn't called
  }
  else {
    Split(siblings, childIDX);
    for (unsigned int i=0; i < siblings.size(); i++) {
        siblings.at(i)->SetOverflowTreatment(true);
    }
    return true; // Split was called
  }
};

//Used to help avoid Split of nodes as often as possible....I'm a bit confused -> Calls Insert
void RRTree::ReInsert(vector<Node*> &siblings, int childIDX) {
  std::cout << "ReInsert was called" << std::endl;
  if (siblings.at(childIDX)->IsNodeALeafNode()) {
    map<int, Point*, std::greater<int>> distanceBetweenPointsAndCenters;
    vector<Point*> childrenPts = siblings.at(childIDX)->GetPointObjs();
    vector<float> treeCenterCoords = siblings.at(childIDX)->GetCenterOfBoundingBox();
    for (unsigned int ii=0; ii<childrenPts.size(); ii++) {
      vector<float> currCenterCoords = childrenPts.at(ii)->GetData();
      float dist=0.0;
      for (unsigned int jj=0; jj<currCenterCoords.size(); jj++) {
        dist += pow((treeCenterCoords.at(jj) - currCenterCoords.at(jj)),2);
      }
      distanceBetweenPointsAndCenters.insert({sqrt(dist),childrenPts.at(ii)});
    }
    int p = ceil(.3 * M);
    int i = 0;
    siblings.at(childIDX)->RemoveAllPointObjs(); //Clears pointers from tree to point objects...but doesn't get rid of children Nodes
    map<int, Point*, std::greater<int>>::iterator map_itr;
    for (map_itr = distanceBetweenPointsAndCenters.begin(); map_itr != distanceBetweenPointsAndCenters.end(); map_itr++) {
      i++;
      if (i <= p) {
        vector<Node*> tempSiblingNodes{root};
        Insert(tempSiblingNodes, 0, *(map_itr->second));
        distanceBetweenPointsAndCenters.erase(map_itr);
      }
      else {
        siblings.at(childIDX)->AddNodePoint(map_itr->second->GetData(),map_itr->second->GetTStep());
      }
    }
    siblings.at(childIDX)->UpdateNodeBoundingBox();
  }
  else {
    std::cout << "Shouldn't be here...should only reinsert points not nodes???" << std::endl;
    /*map<int, Node*, std::greater<int>> distanceBetweenNodeCenters;
    vector<Node*> childrenNodes = tree->GetNodeChildren();
    vector<float> treeCenterCoords = tree->GetCenterOfBoundingBox();
    for (unsigned int ii=0; ii<tree->GetNumberOfChildren(); ii++) {
      vector<float> currCenterCoords = childrenNodes.at(ii)->GetCenterOfBoundingBox();
      float dist=0.0;
      for (unsigned int jj=0; jj<currCenterCoords.size(); jj++) {
        dist += pow((treeCenterCoords.at(jj) - currCenterCoords.at(jj)),2);
      }
      distanceBetweenNodeCenters.insert(std::pair<int, Node*, std::greater<int>>(sqrt(dist),childrenNodes.at(ii)));
    }
    int p = ceil(.3 * M);
    int i = 0;
    siblings.at(childIDX)->RemoveAllChildNodes(); //Clears pointers from tree to children Nodes...but doesn't get rid of children Nodes
    map<int, Node*, std::greater<int>>::iterator map_itr;
    for (map_itr = distanceBetweenNodeCenters.begin(); map_itr != distanceBetweenNodeCenters.end(); map_itr++) {
      i++;
      if (i <= p) {
        vector<Node*> tempSiblingNodes{root};
        Insert(tempSiblingNodes, 0, *(map_itr->second));
        distanceBetweenNodeCenters.erase(map_itr);
      }
      else {
        siblings.at(childIDX)->AddChildNode(map_itr->second);
      }
    }
    siblings.at(childIDX)->UpdateNodeBoundingBox();*/
  }

};

//Finds an even distribution of children Nodes and splits Node -> Calls ChooseSplitAxis, ChooseSplitIndex
void RRTree::Split(vector<Node*> &siblings, int childIDX) {
  //std::cout << "Split was called" << std::endl;
  int axisIDX = ChooseSplitAxis(siblings.at(childIDX));
  //std::cout << "Split along the " << axisIDX+1 << " axis" << std::endl;
  ChooseSplitIndex(siblings, childIDX, axisIDX);
};

//Finds axis to split on and returns axis index
int RRTree::ChooseSplitAxis(Node *&tree) {
  //std::cout << "ChooseSplitAxis was called" << std::endl;
  if (!(tree->IsNodeALeafNode())) {
    //std::cout << "Current Node is not a leaf node" << std::endl;
    vector<Node*> currChildrenNodes = tree->GetNodeChildren();
    vector<Node*> sortedNodes;
    vector<float> marginValues;
    for (unsigned int axs=0; axs < tree->nodeBoundingBox.size(); axs++) {
      sortedNodes = SortNodesMin2Max(currChildrenNodes,axs);
      for (int k=1; k <= (M-2*m+2); k++) {
        vector<Node*> grp1(sortedNodes.begin(),sortedNodes.begin()+k);
        vector<Node*> grp2(sortedNodes.begin()+k,sortedNodes.end());
        if (marginValues.empty() || (marginValues.size()-1) < axs) {
          marginValues.push_back(CalculateMarginValue(TempNodeBB(grp1),TempNodeBB(grp2)));
        }
        else {
          marginValues.at(axs) = marginValues.at(axs) + CalculateMarginValue(TempNodeBB(grp1),TempNodeBB(grp2));
        }
      }
    }
    // Return minimum margin value's index (b/c that idx correlates to the axis)
    return (std::min_element(marginValues.begin(), marginValues.end())-marginValues.begin());
  }
  else {
    //std::cout << "Current Node is a leaf node" << std::endl;
    vector<Point*> currChildrenPoints = tree->GetPointObjs();
    vector<Point*> sortedPoints;
    vector<float> marginValues;
    for (unsigned int axs=0; axs < tree->nodeBoundingBox.size(); axs++) {
      sortedPoints = SortPointsMin2Max(currChildrenPoints,axs);
      for (int k=1; k <= (M-2*m+2); k++) {
        vector<Point*> grp1(sortedPoints.begin(),sortedPoints.begin()+k);
        vector<Point*> grp2(sortedPoints.begin()+k,sortedPoints.end());
        if (marginValues.empty() || (marginValues.size()-1) < axs) {
          marginValues.push_back(CalculateMarginValue(TempPointBB(grp1),TempPointBB(grp2)));
        }
        else {
          marginValues.at(axs) = marginValues.at(axs) + CalculateMarginValue(TempPointBB(grp1),TempPointBB(grp2));
        }
      }
    }
    // Return minimum margin value's index (b/c that idx correlates to the axis)
    return (std::min_element(marginValues.begin(), marginValues.end())-marginValues.begin());
  }
};

//Finds index in which to divide children nodes at and adjusts tree pointer
void RRTree::ChooseSplitIndex(vector<Node*> &siblings, int childIDX, int axsIDX) {
  //std::cout << "ChooseSplitIndex was called" << std::endl;
  if (!(siblings.at(childIDX)->IsNodeALeafNode())) {
    //std::cout << "Node is not a leaf node" << std::endl;
    vector<Node*> currChildrenNodes = siblings.at(childIDX)->GetNodeChildren();
    vector<Node*> sortedNodes;
    sortedNodes = SortNodesMin2Max(currChildrenNodes,axsIDX);
    vector<float> overlapValues;
    vector<float> areaValues;
    for (int k=1; k <= (M-2*m+2); k++) {
      vector<Node*> grp1(sortedNodes.begin(),sortedNodes.begin()+k);
      vector<Node*> grp2(sortedNodes.begin()+k,sortedNodes.end());
      overlapValues.push_back(CalculateOverlapValue(TempNodeBB(grp1),TempNodeBB(grp2)));
      areaValues.push_back(CalculateAreaValue(TempNodeBB(grp1),TempNodeBB(grp2)));
    }
    float currMin = *std::min_element(overlapValues.begin(), overlapValues.end());
    unsigned int currMinIDX = std::min_element(overlapValues.begin(), overlapValues.end())-overlapValues.begin();
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
    siblings.at(childIDX)->RemoveAllChildNodes();
    //Delete childIDX node here
    //std::cout << "About to delete node" << std::endl;
    vector<Node*>::iterator node2Remove_itr=siblings.begin();
    for (int i=0; i < childIDX; i++) {
      node2Remove_itr++;
    }
    if (siblings.at(childIDX) != root) {
      siblings.erase(node2Remove_itr);
    }
    //std::cout << "Deleted node" << std::endl;
    siblings.push_back(grp1Node);
    siblings.push_back(grp2Node);
    //tree->AddChildNode(grp1Node);
    //tree->AddChildNode(grp2Node);
  }
  else {//If children are Point Objects and not Nodes
    //std::cout << "Node is a leaf node" << std::endl;
    vector<Point*> currChildrenPoints = siblings.at(childIDX)->GetPointObjs();
    vector<Point*> sortedPoints;
    vector<float> marginValues;
    sortedPoints = SortPointsMin2Max(currChildrenPoints,axsIDX);
    vector<float> overlapValues;
    vector<float> areaValues;
    for (int k=1; k <= (M-2*m+2); k++) {
      vector<Point*> grp1(sortedPoints.begin(),sortedPoints.begin()+k);
      vector<Point*> grp2(sortedPoints.begin()+k,sortedPoints.end());
      overlapValues.push_back(CalculateOverlapValue(TempPointBB(grp1),TempPointBB(grp2)));
      areaValues.push_back(CalculateAreaValue(TempPointBB(grp1),TempPointBB(grp2)));
    }
    float currMin = *std::min_element(overlapValues.begin(), overlapValues.end());
    unsigned int currMinIDX = std::min_element(overlapValues.begin(), overlapValues.end())-overlapValues.begin();
    for (unsigned int i=0; i<overlapValues.size(); i++) {
      if ((i != currMinIDX) && (currMin == overlapValues.at(i))) {
        if (areaValues.at(currMinIDX) > areaValues.at(i)) {
          currMinIDX = i;
        }
      }
    }
    Node* grp1Node = new Node();
    vector<Point*> newGrp1(sortedPoints.begin(),sortedPoints.begin()+currMinIDX+1);
    for (unsigned int i=0; i < newGrp1.size(); i++) {
      grp1Node->AddNodePoint(newGrp1.at(i)->GetData(),newGrp1.at(i)->GetTStep());
    }
    Node* grp2Node = new Node();
    vector<Point*> newGrp2(sortedPoints.begin()+currMinIDX+1,sortedPoints.end());
    for (unsigned int i=0; i < newGrp2.size(); i++) {
      grp2Node->AddNodePoint(newGrp2.at(i)->GetData(),newGrp2.at(i)->GetTStep());
    }
    //std::cout << "About to delete node" << std::endl;
    siblings.at(childIDX)->RemoveAllPointObjs();
    //Delete childIDX node here
    vector<Node*>::iterator node2Remove_itr=siblings.begin();
    for (int i=0; i < childIDX; i++) {
      node2Remove_itr++;
    }
    if (siblings.at(childIDX) != root) {
      siblings.erase(node2Remove_itr);
    }
    //std::cout << "Deleted node" << std::endl;
    siblings.push_back(grp1Node);
    siblings.push_back(grp2Node);
    //std::cout << "# of siblings: " << siblings.size() << std::endl;
    //tree->AddChildNode(grp1Node);
    //tree->AddChildNode(grp2Node);
  }
};

//Removes point from leaf node -> Calls ChooseSubtree
void RRTree::Remove(Node *&tree, Point pt2Remove) {

};

//Transforms point data
void RRTree::Transform(vector<float> userTransform, Point *&pt2Transform) {
  //pt2Transform->TransformData(userTransform);
};

//Builds and returns a temporary bounding box around list of Nodes
vector<vector<float>> RRTree::TempNodeBB(vector<Node*> userGroup) {
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

//Builds and returns a temporary bounding box around list of Points
vector<vector<float>> RRTree::TempPointBB(vector<Point*> userGroup) {
  vector<vector<float>> grpBB;
  float min=0;
  float max=0;
  for (int dim=0; dim < userGroup.at(0)->GetPointDimensions(); dim++) {
    for (unsigned int i=0; i < userGroup.size(); i++) {
      vector<float> currPtData = userGroup.at(i)->GetData();
      if (i==0) {
        min = currPtData.at(dim);
        max = currPtData.at(dim);;
      }
      else {
        if(currPtData.at(dim) < min) {
          min = currPtData.at(dim);;
        }
        if (currPtData.at(dim) > max) {
          max = currPtData.at(dim);;
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
float RRTree::CalculateMarginValue(vector<vector<float>> grp1BB, vector<vector<float>> grp2BB) {
  float margin = 0.0;
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
      float maxPt=0.0;
      float minPt=0.0;
      if (grp1BB.at(dim).at(1) < grp2BB.at(dim).at(1)) {
        maxPt = grp1BB.at(dim).at(1);
      }
      else {
        maxPt = grp2BB.at(dim).at(1);
      }
      if (grp1BB.at(dim).at(0) > grp2BB.at(dim).at(0)) {
        minPt = grp1BB.at(dim).at(0);
      }
      else {
        minPt = grp2BB.at(dim).at(0);
      }
      overlap *= (maxPt-minPt);
    }
    return overlap;
  }
  else {
    return 0.0;
  }
  return overlap; //delete this line after uncommenting whole function
};

//Starts at root and updates all nodes bounding boxes
void RRTree::updateNodeBoundingBoxes(vector<Node*> &branches) {
  for (unsigned int i=0; i < branches.size(); i++) {
    branches.at(i)->UpdateNodeBoundingBox();
    vector<Node*> children = branches.at(i)->GetNodeChildren();
    if (children.size() > 0) {
      updateNodeBoundingBoxes(children);
    }
  }
};

//Print Tree by layer starting at the root Node
void RRTree::printTree() {
  vector<Node*> tempVectNodes{root};
  printEachLayer(tempVectNodes);
}

void RRTree::printEachLayer(vector<Node*> currLayer) {
  if (currLayer.size() > 1) {
    vector<Node*> nextLayer;
    for (unsigned int i=0; i < currLayer.size(); i++) {
      //Print out current Node
      std::cout << "Node " << i << " has bounds {";
      vector<vector<float>> treeBounds = currLayer.at(i)->GetNodeBounds();
      for (unsigned int dim=0; dim < treeBounds.size(); dim++) {
        std::cout << treeBounds.at(dim).at(0) << ",";
        if (dim+1 < treeBounds.size()) {
          std::cout << treeBounds.at(dim).at(1) << "; ";
        }
        else {
          std::cout << treeBounds.at(dim).at(1) << "}";
        }
      }
      std::cout << "\t";
      if (currLayer.at(i)->GetNumberOfChildren() > 0) {
        vector<Node*> currChildren = currLayer.at(i)->GetNodeChildren();
        for (unsigned int j=0; j < currChildren.size(); j++) {
          nextLayer.push_back(currChildren.at(j));
        }
      }
    }
    std::cout << std::endl;
    if (!nextLayer.empty()) {
      printEachLayer(nextLayer);
    }
    else {
      for (unsigned int i=0; i < currLayer.size(); i++) {
        if (currLayer.at(i)->GetNumberOfPoinObjs() > 0) {
          //Print out points in node
          vector<Point*> currPoints = currLayer.at(i)->GetPointObjs();
          for (unsigned int j=0; j < currPoints.size(); j++) {
            std::cout << "Point " << j << " belonging to " << i << " Node: Data {";
            vector<float> ptData = currPoints.at(j)->GetData();
            for (unsigned int k=0; k < ptData.size(); k++) {
              if (k+1 < ptData.size()) {
                std::cout << ptData.at(k) << ",";
              }
              else {
                std::cout << ptData.at(k) << "}; ";
              }
            }
            std::cout << "T_step = " << currPoints.at(j)->GetTStep();
            std::cout << "; Cluster ID = " << currPoints.at(j)->GetClID() << std::endl;
          }
        }
      }
    }
  }
  else if (currLayer.size() > 0) {
    //Print out current Node
    std::cout << "Root has bounds {";
    vector<vector<float>> treeBounds = currLayer.at(0)->GetNodeBounds();
    for (unsigned int dim=0; dim < treeBounds.size(); dim++) {
      std::cout << treeBounds.at(dim).at(0) << ",";
      if (dim+1 < treeBounds.size()) {
        std::cout << treeBounds.at(dim).at(1) << "; ";
      }
      else {
        std::cout << treeBounds.at(dim).at(1) << "}";
      }
    }
    std::cout << std::endl;
    if (currLayer.at(0)->GetNumberOfChildren() > 0) {
      printEachLayer(currLayer.at(0)->GetNodeChildren());
    }
    else if (currLayer.at(0)->GetNumberOfPoinObjs() > 0) {
      //Print out points in node
      vector<Point*> currPoints = currLayer.at(0)->GetPointObjs();
      for (unsigned int i=0; i < currPoints.size(); i++) {
        std::cout << "Point " << i << ": Data {";
        vector<float> ptData = currPoints.at(i)->GetData();
        for (unsigned int j=0; j < ptData.size(); j++) {
          if (j+1 < ptData.size()) {
            std::cout << ptData.at(j) << ",";
          }
          else {
            std::cout << ptData.at(j) << "}; ";
          }
        }
        std::cout << "T_step = " << currPoints.at(i)->GetTStep();
        std::cout << "; Cluster ID = " << currPoints.at(i)->GetClID() << std::endl;
      }
    }
    else {
      std::cout << "Root has no children" << std::endl;
    }
  }
  else {
    std::cout << "Shouldn't be called" << std::endl;
  }
};
