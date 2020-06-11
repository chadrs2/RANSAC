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
RRTree::~RRTree() {clear();}
//Empties RR*-tree
void RRTree::clear() {
  if (root != NULL) {
    removeNodes(root);
  }
};
void RRTree::removeNodes(Node *&tree) {
  vector<Node*> currChildrenNodes = tree->GetNodeChildren();
  for(unsigned int i=0; i < currChildrenNodes.size(); i++) {
    removeNodes(currChildrenNodes.at(i));
  }
  currChildrenNodes.clear();
  delete tree;
};
//Node* RRTree::GetRootNode() const {return root;};

//True if data was able to be inserted -> Calls Insert
void RRTree::InsertData(Point newPoint) {
  std::cout << "InsertData is called" << std::endl;
  Insert(root, newPoint);
};

//True if data was able to be inserted.
//Keeps nodes balanced and within M and m parameters -> Calls ChooseSubtree, OverflowTreatment
void RRTree::Insert(Node *&tree, Point newPoint) {
  if (tree == NULL) {
    tree = new Node();
  }
  // ChooseSubtree Algorithm
  if(!(tree->IsNodeALeafNode())) {
    vector<Node*> currChildren = tree->GetNodeChildren();
    if (tree->GetNodeChildren().at(0)->IsNodeALeafNode()) {//children are leaf nodes
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
      std::cout << "Enter " << minOVIdx << " branch (if)" << std::endl;
      Insert(currChildren.at(minOVIdx),newPoint);
    }
    else {
      // Determine minimum area costs - // TODO: ties: smallest area
      Node* newNode = new Node();
      newNode->AddNodePoint(newPoint.GetData(),newPoint.GetTStep());
      vector<Node*> children = tree->GetNodeChildren();
      vector<float> minAreaCost;
      for (unsigned int i=0; i < children.size(); i++) {
        vector<vector<float>> currBB = children.at(i)->GetNodeBounds();
        vector<Node*> nodeAndPointNode; // TODO: Can I just reset vector to NULL or do I have to delete it first -------
        nodeAndPointNode.push_back(children.at(i));
        nodeAndPointNode.push_back(newNode);
        vector<vector<float>> newBB = TempNodeBB(nodeAndPointNode);
        minAreaCost.push_back(CalculateAreaValue(currBB,newBB));
        //Clean Up nodeAndPointNode vector
        nodeAndPointNode.clear();
      }
      int minAreaIdx = std::min_element(minAreaCost.begin(), minAreaCost.end())-minAreaCost.begin();
      delete newNode;
      std::cout << "Enter " << minAreaIdx << " branch (else)" << std::endl;
      Insert(children.at(minAreaIdx),newPoint);
    }
  }
  //Insert Algorithm
  std::cout << "Number of points in current Node: " << tree->GetNumberOfPoinObjs() << " when M is " << M << std::endl;
  if (tree->GetNumberOfPoinObjs() < M) {
    if (tree->IsNodeALeafNode()) {
      tree->AddNodePoint(newPoint.GetData(),newPoint.GetTStep());
    }
  }
  else {
    tree->AddNodePoint(newPoint.GetData(),newPoint.GetTStep());
    bool wasSplitPerformed = OverflowTreatment(tree);
    //////////// NOT USED (BELOW THIS LINE) //////////////
    if (wasSplitPerformed && (tree->IsNodeALeafNode())) {
      std::cout << "Need to create a new root" << std::endl;
      //TODO: create new root & update each nodes level -----------------------------------------------------
    }
    else if (wasSplitPerformed) {
      std::cout << "Propagate OverflowTreatment upwards" << std::endl;
      //TODO: Propagate OverflowTreatment upwards?? ---------------------------------------------------------
    }
    //////////// NOT USED (ABOVE THIS LINE) ///////////////
  }
  //updateNodeBoundingBoxes(tree);
  tree->UpdateNodeBoundingBox();
};

//Recursive function moving through tree until it reaches a leaf node -> Calls itself
//Node* RRTree::ChooseSubtree(Node *&tree, Point newPoint) { //Keep getting erros with Node* returns
void RRTree::ChooseSubtree(Node *&tree, Point newPoint) {
  /*if (tree->IsNodeALeafNode()) {
    return tree;
  }
  else {
    if (tree->GetNodeLevel() == 1) {//children are leaf nodes
      //Determine minimum overlap costs - // TODO: ties: smallest area change
      Node* newNode = new Node(newPoint);
      vector<Node*> children = tree->GetNodeChildren();
      vector<float> minOverlapCost;
      vector<float> minAreaCost;
      for (unsigned int i=0; i < children.size(); i++) {
        vector<vector<float>> currBB = children.at(i)->GetNodeBounds();
        vector<Node*> nodeAndPointNode = NULL; // TODO: Can I just reset vector to NULL or do I have to delete it first -------
        nodeAndPointNode.push_back(children.at(i));
        nodeAndPointNode.push_back(newNode);
        vector<vector<float>> newBB = TempBB(nodeAndPointNode);
        minOverlapCost.push_back(CalculateOverlapValue(currBB,newBB));
        minAreaCost.push_back(CalculateAreaValue(currBB,newBB));
      }
      int minOVIdx = std::min_element(minOverlapCost.begin(), minOverlapCost.end())-minOverlapCost.begin();
      ChooseSubtree(children.at(minOVIdx),newPoint);
    }
    else {
      // Determine minimum area costs - // TODO: ties: smallest area
      Node* newNode = new Node(newPoint);
      vector<Node*> children = tree->GetNodeChildren();
      vector<float> minAreaCost;
      for (unsigned int i=0; i < children.size(); i++) {
        vector<vector<float>> currBB = children.at(i)->GetNodeBounds();
        vector<Node*> nodeAndPointNode = NULL; // TODO: Can I just reset vector to NULL or do I have to delete it first -------
        nodeAndPointNode.push_back(children.at(i));
        nodeAndPointNode.push_back(newNode);
        vector<vector<float>> newBB = TempBB(nodeAndPointNode);
        minAreaCost.push_back(CalculateAreaValue(currBB,newBB));
      }
      int minAreaIdx = std::min_element(minAreaCost.begin(), minAreaCost.end())-minAreaCost.begin();
      ChooseSubtree(children.at(minAreaIdx),newPoint);
    }
  }*/
};

//True if Split() was called.
//Meant for dealing with a filled Node -> Calls ReInsert or Split depending on nodeLevel
bool RRTree::OverflowTreatment(Node *&tree) {
  std::cout << "Overflow was called" << std::endl;
  if (!(tree->WasOverflowTreatmentCalled()) && (tree != root)) {
    tree->SetOverflowTreatment(true);
    ReInsert(tree);
    return false; // Split wasn't called
  }
  else {
    tree->SetOverflowTreatment(true);
    Split(tree);
    return true; // Split was called
  }
};

//Used to help avoid Split of nodes as often as possible....I'm a bit confused -> Calls Insert
void RRTree::ReInsert(Node *&tree) {
  std::cout << "ReInsert was called" << std::endl;
  if (tree->IsNodeALeafNode()) {
    map<int, Point*, std::greater<int>> distanceBetweenNodeCenters;
    vector<Point*> childrenPts = tree->GetPointObjs();
    vector<float> treeCenterCoords = tree->GetCenterOfBoundingBox();
    for (unsigned int ii=0; ii<childrenPts.size(); ii++) {
      vector<float> currCenterCoords = childrenPts.at(ii)->GetData();
      float dist=0.0;
      for (unsigned int jj=0; jj<currCenterCoords.size(); jj++) {
        dist += pow((treeCenterCoords.at(jj) - currCenterCoords.at(jj)),2);
      }
      distanceBetweenNodeCenters.insert({sqrt(dist),childrenPts.at(ii)});
    }
    int p = ceil(.3 * M);
    int i = 0;
    tree->RemoveAllChildNodes(); //Clears pointers from tree to children Nodes...but doesn't get rid of children Nodes
    map<int, Point*, std::greater<int>>::iterator map_itr;
    for (map_itr = distanceBetweenNodeCenters.begin(); map_itr != distanceBetweenNodeCenters.end(); map_itr++) {
      i++;
      if (i <= p) {
        Insert(root, *(map_itr->second));
        distanceBetweenNodeCenters.erase(map_itr);
      }
      else {
        tree->AddNodePoint(map_itr->second->GetData(),map_itr->second->GetTStep());
      }
    }
    tree->UpdateNodeBoundingBox();
  }
  else {
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
    map<int, Node*, std::greater<int>>::iterator map_itr;
    for (map_itr = distanceBetweenNodeCenters.begin(); map_itr != distanceBetweenNodeCenters.end(); map_itr++) {
      i++;
      if (i <= p) {
        distanceBetweenNodeCenters.erase(map_itr);
        // TODO: Does this actually get rid of child node?? ---------------------------------------------------------
      }
    }
    tree->UpdateNodeBoundingBox();
    // TODO: R14 in paper ... call Insert with removed Nodes --------------------------------------------------------
    */
  }

};

//Finds an even distribution of children Nodes and splits Node -> Calls ChooseSplitAxis, ChooseSplitIndex
void RRTree::Split(Node *&tree) {
  std::cout << "Split was called" << std::endl;
  int axisIDX = ChooseSplitAxis(tree);
  std::cout << "Split along the " << axisIDX+1 << " axis" << std::endl;
  ChooseSplitIndex(tree, axisIDX);
};

//Finds axis to split on and returns axis index
int RRTree::ChooseSplitAxis(Node *&tree) {
  std::cout << "ChooseSplitAxis was called" << std::endl;
  if (tree != root) {
    std::cout << "Current Node is not the root" << std::endl;
    /*vector<Node*> currChildrenNodes = tree->GetNodeChildren();
    map<float, Node*, std::less<float>> sortedMapOfNodes;
    map<float, Node*, std::less<float>>::iterator itr;
    vector<Node*> sortedNodes;
    vector<float> marginValues;
    for (unsigned int axs=0; axs < tree->nodeBoundingBox.size(); axs++) {
      sortedMapOfNodes = {};
      sortedNodes = {};
      for (unsigned int entry=0; entry < currChildrenNodes.size(); entry++) {
        vector<vector<float>> currBB = currChildrenNodes.at(entry)->GetNodeBounds();
        sortedMapOfNodes.insert({currBB.at(axs).at(0),currChildrenNodes.at(entry)});
        // TODO: What happens if minimum value of bounding box is 0 ////////////////////////////////
      }
      for (itr=sortedMapOfNodes.begin(); itr != sortedMapOfNodes.end(); itr++) {
        sortedNodes.push_back(itr->second);
        sortedMapOfNodes.erase(itr);
      }
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
    return (std::min_element(marginValues.begin(), marginValues.end())-marginValues.begin());*/
    return 1;
  }
  else {
    std::cout << "Current Node is the root" << std::endl;
    vector<Point*> currChildrenPoints = tree->GetPointObjs();
    map<float, Point*, std::less<float>> sortedMapOfPoints;
    map<float, Point*, std::less<float>>::iterator pt_itr;
    vector<Point*> sortedPoints;
    vector<float> marginValues;
    for (unsigned int axs=0; axs < tree->nodeBoundingBox.size(); axs++) {
      sortedMapOfPoints = {};
      sortedPoints = {};
      for (unsigned int entry=0; entry < currChildrenPoints.size(); entry++) {
        vector<float> currPtsData = currChildrenPoints.at(entry)->GetData();
        sortedMapOfPoints.insert({currPtsData.at(axs),currChildrenPoints.at(entry)});
      }
      for (pt_itr=sortedMapOfPoints.begin(); pt_itr != sortedMapOfPoints.end(); pt_itr++) {
        sortedPoints.push_back(pt_itr->second);
        sortedMapOfPoints.erase(pt_itr);
      }
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
void RRTree::ChooseSplitIndex(Node *&tree, int axsIDX) {
  std::cout << "ChooseSplitIndex was called" << std::endl;
  if (tree != root) {
    /*vector<Node*> currChildrenNodes = tree->GetNodeChildren();
    map<float, Node*, std::less<float>> sortedMapOfNodes;
    map<float, Node*, std::less<float>>::iterator itr;
    vector<Node*> sortedNodes;
    for (unsigned int entry; entry < currChildrenNodes.size(); entry++) {
      vector<vector<float>> currBB = currChildrenNodes.at(entry)->GetNodeBounds();
      sortedMapOfNodes.insert({currBB.at(axsIDX).at(0),currChildrenNodes.at(entry)});
      // TODO: What happens if minimum value of bounding box is 0 //////////////////////////////////
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
    // TODO: DO THIS WORK?? ----------------------------------------------------------------------------------
    */
  }
  else {//If children are Point Objects and not Nodes
    vector<Point*> currChildrenPoints = tree->GetPointObjs();
    map<float, Point*, std::less<float>> sortedMapOfPoints;
    map<float, Point*, std::less<float>>::iterator pt_itr;
    vector<Point*> sortedPoints;
    vector<float> marginValues;
    for (unsigned int entry=0; entry < currChildrenPoints.size(); entry++) {
      vector<float> currPtsData = currChildrenPoints.at(entry)->GetData();
      sortedMapOfPoints.insert({currPtsData.at(axsIDX),currChildrenPoints.at(entry)});
    }
    for (pt_itr=sortedMapOfPoints.begin(); pt_itr != sortedMapOfPoints.end(); pt_itr++) {
      sortedPoints.push_back(pt_itr->second);
      sortedMapOfPoints.erase(pt_itr);
    }
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
    tree->RemoveAllPointObjs();
    tree->AddChildNode(grp1Node);
    tree->AddChildNode(grp2Node);
    tree->UpdateNodeBoundingBox();
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
void RRTree::updateNodeBoundingBoxes(Node *&tree) {
  /*tree->UpdateNodeBoundingBox();
  vector<Node*> children = tree->GetNodeChildren();
  for (unsigned int i=0; i < children.size(); i++) {
    updateNodeBoundingBoxes(children.at(i));
  }*/
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
        std::cout << treeBounds.at(dim).at(0) << ", ";
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
                std::cout << ptData.at(k) << ", ";
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
    std::cout << "Node has bounds {";
    vector<vector<float>> treeBounds = currLayer.at(0)->GetNodeBounds();
    for (unsigned int dim=0; dim < treeBounds.size(); dim++) {
      std::cout << treeBounds.at(dim).at(0) << ", ";
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
            std::cout << ptData.at(j) << ", ";
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
