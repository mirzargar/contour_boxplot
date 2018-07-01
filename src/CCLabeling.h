#ifndef CCLABELING_H
#define CCLABELING_H

#include <utility>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include "Matrix.h"

using namespace std;

typedef pair<int, int> IntPairVector;

bool sort_pred(IntPairVector a, IntPairVector b);
int getBCCLabelIntArray(Matrix<int> p, int labelNum);

bool PairComp(IntPairVector a, IntPairVector b);
Matrix<int> CCLabeling(Matrix<int> p);
Matrix<int> getBCC(Matrix<int> Connected);

#endif // CCLABELING_H
