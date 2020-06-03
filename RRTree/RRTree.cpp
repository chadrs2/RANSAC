#include "RRTree.h"
#include <cmath>
#include <map>
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
Node* RRTree::getRootNode() const {return root;};


bool InsertData(Point newPoint) {
  return Insert(root, newPoint);
  return false;
};

bool RRTree::Insert(Node *&tree, Point newPoint) {
  Node* insertNode = ChooseSubtree(tree, newPoint);
  if (insertNode.getChildren().size() < M) {
    N = new Node(newPoint);
    insertNode->childrenNodes.push_back(N);
    insertNode->nodeLevel++;
  }
  if (insertNode.getChildren().size() == M) {

    ///////////////////
  }
  updateNodeBoundingBoxes(root);
  return true;
};

void RRTree::OverflowTreatment();

void RRTree::ReInsert();

Node* RRTree::ChooseSubtree(Node *&tree, Point newPoint) {
  /*vector<int> COV;
  vector<Node*> currChildrenVect = tree.getChildren();
  map<float,Node*> nodePerims;
  map<float,Node*> nodeVols;

  for (unsigned int i=0; i<currChildrenVect.size(); i++) {
//    nodePerims.insert(std::pair<float,Node*>(currChildrenVect.at(i).nodePerimeter(),currChildrenVect.at(i)));
//    nodeVols.insert(std::pair<float,Node*>(currChildrenVect.at(i).nodeVolume(),currChildrenVect.at(i)));
    vector<vector<float>> nodeSize = currChildrenVect.at(i).getNodeSize();
    vector<float> pointData = newPoint.getData();
    bool pointFitsInSubtree;
    for (unsigned int dim=0; dim<nodeSize.size(); dim++) {
      // First value is the minimum in that dimension and second is the maximum
      if (nodeSize.at(dim).at(0) > pointData.at(dim) && nodeSize.at(dim).at(1) < pointData.at(dim)) {
        pointFitsInSubtree = false;
        break;
      }
      else {
        pointFitsInSubtree = true;
      }
    }
    if (pointFitsInSubtree) {
      nodePerims.insert(std::pair<float,Node*>(currChildrenVect.at(i).nodePerimeter(),currChildrenVect.at(i)));
      nodeVols.insert(std::pair<float,Node*>(currChildrenVect.at(i).nodeVolume(),currChildrenVect.at(i)));
      COV.push_back(i);
    }
  }
  if (COV.size() != 0) {
    if (nodeVols.find(0) != nodeVols.end()) { //there exists a node with zero volume
      return nodePerims.begin()->second;
    }
    else {
      return nodeVols.begin()->second;
    }
  }*/


  //Node* N = tree;
  /*if (N->nodeLevel == 0) {
    return N;
  }
  else {
    if (N.getChildren().at(0).nodeLevel == 0) {

    }
    else {

    }
  }*/
};

void RRTree::CheckComp();

void RRTree::Split();

void RRTree::ChooseSplitAxis();

void RRTree::ChooseSplitIndex();

void RRTree::remove(Point pt2Remove);

void RRTree::Transform(vector<double> userTransform, Point pt2Transform)

void RRTree::clear();

/*void updateNodeBoundingBoxes(Node *&tree) {
  vector<Node*> currChildrenVect = tree->childrenNodes;
  vector<float> mins;
  vector<float> maxs;
  for (unsigned int ii=0; ii<currChildrenVect.size(); ii++) {
    Node currNode = currChildrenVect.at(ii);
    vector<float> currNodeSize = currNode.getNodeSize();
    if (currNodeSize != NULL) {
      for (unsigned int jj=0; jj<currNodeSize.size(); jj+=2) {
        if (ii != 0) {
          if (currNodeSize.at(jj) < mins.at(jj)) {
            mins.at(jj) = currNodeSize.at(jj);
          }
          if (currNodeSize.at(jj+1) > maxs.at(jj+1)) {
            maxs.at(jj+1) = currNodeSize.at(jj+1);
          }
        }
        else {
          mins.push_back(currNodeSize.at(jj));
          maxs.push_back(currNodeSize.at(jj+1));
        }
      }
    }
  }
  vector<float> currNodeSize;
  for (unsigned int ii=0; ii<mins.size(); ii++) {
    currNodeSize.push_back(mins.at(ii));
    currNodeSize.push_back(maxs.at(ii));
  }
  tree->nodeSize = currNodeSize;
  for (unsigned int ii=0; ii<currChildrenVect.size(); ii++) {
    updateNodeBoundingBoxes(currChildrenVect.at(ii));
  }
};*/
