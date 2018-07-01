#ifndef CBD_H
#define CBD_H

#include "Utility.h"
#include "EnsembleHandler.h"
#include "CCLabeling.h"
#include "ContourExtraction.h"
#include "SVG.h"

using namespace std;

class CBD{
private:
    vector<double> Probs;
    double PThreshold;
    Matrix<double> UPercents;
    Matrix<double> IPercents;
    unsigned int numJ;
    Ensemble E;
    //unsigned int MRun;
    string Dir;
    Matrix<double> lon;
    Matrix<double> lat;
    double shift_x;
    double shift_y;

public:

    CBD(string src_dir, unsigned int j, double levelSet, string des_dir, unsigned int ESize, unsigned int Xdim, unsigned int Ydim);

    // getter & setters
    Matrix<int> getUBand(const Ensemble e, const vector<int>& indx);
    Matrix<int> getIBand(const Ensemble e, const vector<int>& indx);

    void getPercentsFast();
    void getPercentProbabilities();
    friend bool sort_pred(const my_pair& left, const my_pair& right);
    friend vector<int> getSortedIndex(const vector<double>& v);
    void WriteSortedProbabilities(string fileName);

    void getTimePercents(unsigned int ESize);
    void genSVG(unsigned int outlierNum, string des_dir, double xmin, double ymin, double xmax, double ymax);
    //void genAnimation(string des_dir, double xmin, double ymin, double xmax, double ymax);
    Matrix<int> getBandMask(const Ensemble e, vector<int> v);
    double getThreshold(int outlierNum);

    Matrix<int> getTrimmedMean(Ensemble e, vector<int> v);
    void Load_lonLat(string fileName);

    friend void WriteLog();
};

#endif // CBD_H
