#include<iostream>
#include "Point.h"
#include "RRTree.h"
using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
  if (argc > 1) {
    int M = *argv[1] - '0';
    int m = *argv[2] - '0';
    //std::cout << M << " " << m << std::endl;
    RRTree rrTree(M,m);
  }
  else {
    RRTree rrTree();
  }
  return 0;
}
