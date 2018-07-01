#ifndef CONTOUREXTRACTION_H
#define CONTOUREXTRACTION_H

#include <utility>
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include "Matrix.h"
#include "Utility.h"

using namespace std;

vector<IntPairVector> InitializeClkNeighbors();
inline IntPairVector operator+(const IntPairVector p1, const IntPairVector p2);
inline IntPairVector operator-(const IntPairVector p1, const IntPairVector p2);
int findIndex(vector<IntPairVector> clkNeighbors, IntPairVector p);
bool OutOfBound(IntPairVector P, unsigned int xdim, unsigned int ydim);
void findStartPoint(Matrix<int> p, IntPairVector& StartPnt, IntPairVector& PrevNeigh);
vector <IntPairVector> getNeighbors(IntPairVector P, IntPairVector PrevNeigh, unsigned int xdim, unsigned int ydim, vector<IntPairVector> clkNeighbors);
Matrix<int> MooreNeighborTracing(Matrix<int> p);
vector<IntPairVector> MooreNeighborTracingIndx(Matrix<int> p);

#endif // CONTOUREXTRACTION_H
