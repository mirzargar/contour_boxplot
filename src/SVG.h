#ifndef SVG_H
#define SVG_H

#include "Utility.h"
#include "ContourExtraction.h"
#include "CCLabeling.h"
#include "Matrix.h"

void InitializeSVG(string SVGfname, double xmin, double ymin, double xmax, double ymax);

void AddSVGBG(char* SVGfname, char* BGfname, int x, int y, int width, int height);

void DrawSVGIntPath(string SVGfname, vector<IntPairVector> Pset, char* strokeColor, double strokeWidth, char* StrokeType, char* fill, bool smooth, Matrix<double> lon, Matrix<double> lat, double xmin, double ymin);
void DrawSVGIntPath(string SVGfname, vector<IntPairVector> Pset, char* strokeColor, double strokeWidth, char* StrokeType, char* fill, bool smooth, Matrix<double> lon, Matrix<double> lat, double xmin, double ymin, unsigned int Xdim, unsigned int Ydim);

vector<IntPairVector> SmoothPath(vector<IntPairVector> v);
vector<DoublePairVect> SmoothPath(vector<DoublePairVect> v);

void FillPolygon(string SVGfname, vector<IntPairVector> Pset, char* fillColor, char* StrokeColor, double StrokeWidth, double opacity, bool smooth);

void DrawConnectedComp(string SVGfname, Matrix<int> m, char* StrokeColor, double StrokeWidth, double opacity, Matrix<double> lon, Matrix<double> lat, double xmin, double ymin, bool smooth);

void DrawCCSilhouette(string SVGfname, Matrix<int> m, char* StrokeColor, double StrokeWidth, char *StrokeType, char *fill, Matrix<double> lon, Matrix<double> lat, double xmin, double ymin, bool smooth, unsigned int Xdim, unsigned int Ydim);

void CloseSVGFile(string SVGfname);

#endif // SVG_H
