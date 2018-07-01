#ifndef UTILITY_H
#define UTILITY_H

#include <sstream>
#include "Matrix.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkConnectedComponentImageFilter.h"

#include "itksys/SystemTools.hxx"
#include <sstream>

#define _USE_MATH_DEFINES

#define StartPath 4
#define FinishPath 3
#define OnPath 2
#define OnBorder 1

typedef pair<int, int> IntPairVector;
typedef pair<double, double> DoublePairVect;
typedef pair<double, int> my_pair;


//////////////////
///// map projection
/////////////////

DoublePairVect LambertConformalConicSphere1SP(double R, double lon, double lon0, double lat, double lat0, double SP1);
vector<DoublePairVect> getMapProjSphere(vector<DoublePairVect> pset, double cntre_lat, double cntre_lon, double SP);
DoublePairVect getMapProjSphere(double lon, double lat, double cntre_lat, double cntre_lon, double SP);

DoublePairVect LambertConformalConicEllipse1SP(double a, double e, double lon, double lon0, double lat, double lat0, double SP1);
vector<DoublePairVect> getMapProjEllipsoid(vector<DoublePairVect> pset);

/////////////////
// Combinatorics
/////////////////
unsigned long choose(unsigned long n, unsigned long k);
template <typename Iterator>
bool next_combination(const Iterator first, Iterator k, const Iterator last);
Matrix<int> getCombinationMatrix(unsigned int n, unsigned int k);

Matrix<int> zeroCrossing(Matrix<double> M);

Matrix<int> getLevelSetMask(Matrix<double> m, double levelSet);
Matrix<int> getBinaryUnion(Matrix<int> m, Matrix<int> n);
Matrix<int> getBinaryIntersect(Matrix<int> m, Matrix<int> n);
double getBinaryIntersectSum(const Matrix<int>& m, const Matrix<int>& n);
double getPercentBinarySubset(Matrix<int> m, Matrix<int> n);
bool IsOnBorder(IntPairVector p, IntPairVector border);
unsigned int getState(IntPairVector p1, IntPairVector p2, IntPairVector border);
Matrix<int> getCCLabeling(Matrix<int> p, bool recip);
bool CCWithHole(Matrix<int> p);
bool ContourHitsEdge(vector<IntPairVector> contour, unsigned int Xdim, unsigned int Ydim);

string StringCat(string s1, string s2);

#endif // UTILITY_H
