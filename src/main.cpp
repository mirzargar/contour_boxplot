#include <iostream>
#include <time.h>
#include <cstdlib>
#include "EnsembleHandler.h"
#include "CCLabeling.h"
#include "ContourExtraction.h"
#include "CBD.h"

//#include "QuickView.h"

using namespace std;

//./contourboxplot test 21 258 test 185 129
int main (int argc, char *argv[]){
    try{

        if(argc < 6){
            cout << "Corret usage: ./contourboxplot <src_dir> <ens_size> <level_set> <des_dir> <Xdim> <Ydim>" << endl;
            cout << "NOTE: files are considered to be linearized - if change needed consult function LoadFieldFileWithoutCoord in EnsembleHandler" << endl;
            cout << "      the files should be named as mem_c_0 to mem_c_ens_size - if change needed consult LoadFieldEnsemble in EnsembleHandler" << endl;
            cout << "      the value of j has been hard coded in the main file (j=2) - it also has been hard coded partially in getTimePercents function in CBD" << endl;
            cout << "      the outlier num is specified as the number of members with zero depth value. If every member has a non-zero depth value then it has been set to three (see CBD constructor)" << endl;
            cout << "      The SVG generator has been tuned for SREF ensemble - it WIL NOT work for other dataset most probably" << endl;
            exit(-1);
        }

        //Default settings
        unsigned int j = 2;


        string src_dir = argv[1], des_dir = argv[4];
        double levelset = atof(argv[3]);
        unsigned int ESize = atoi(argv[2]);
        unsigned int Xdim = atoi(argv[5]), Ydim = atoi(argv[6]);

        cout << endl
             << "    source directory: " << src_dir <<  endl <<
                "    dimensions: [" << Xdim << ", " << Ydim << "]" << endl <<
                "    ensemble size: " << ESize << endl << endl;

        CBD(src_dir, j, levelset, des_dir, ESize, Xdim, Ydim);

        exit(1);


    }catch(int i){
        if( i == 1)
            cout<<"***Dimension mismatch***"<<endl;
        system ("PAUSE");
        return 0;
    }
    exit(1);
}
