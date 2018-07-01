#include "EnsembleHandler.h"

using namespace std;

Ensemble::Ensemble(){
    EnsembleSize = -1;
    Xdim = -1;
    Ydim = -1;
    LevelSet = -1.0;
}

Ensemble::Ensemble(const Ensemble& e)
:EnsembleSize(e.EnsembleSize){
    EnsembleSize = e.EnsembleSize;
    Xdim = e.Xdim;
    Ydim = e.Ydim;
    FieldAvg = e.FieldAvg;
    LevelSet = e.LevelSet;
    //gridRes = e.gridRes;

    Members.clear();
    Orders.clear();
    MemFileName.clear();
    for (unsigned int i=0; i<EnsembleSize; i++){
        Members.push_back(e.Members[i]);
        Orders.push_back(e.Orders[i]);
        MemFileName.push_back(e.MemFileName[i]);
    }
}

Ensemble::~Ensemble(){
    Members.clear();
}

void Ensemble::InitializeEnsemble(unsigned int ESize, const double levelSet, unsigned int xdim, unsigned int ydim){
    EnsembleSize = ESize;
    Xdim = xdim;
    Ydim = ydim;

    LevelSet = levelSet;
    Members.resize(ESize, Matrix<int>(Xdim, Ydim, 0));
    MemFileName.resize(ESize, "");
    FieldAvg = Matrix<double>(Xdim, Ydim, 0.0);
    for(unsigned int i=0; i<ESize; i++)
        Orders.push_back(i);
}

string getFileName(string fname, unsigned int fileNumber){
    stringstream ss;
    string s, fileName;
    s = fname;
    ss.str("");
    ss << s << "_" << fileNumber << ".txt";
    fileName = ss.str();
    return fileName;
}

Matrix<double> Ensemble::LoadFieldFileWithCoord(Matrix<double> &toReturn, string fileName){
    //TODO: Not stable
    ifstream data(fileName.c_str());
cout << "+" << fileName << endl;

    if(!data.is_open()){
    printf("Bad Vector Filename: %s\n", fileName.c_str());
    exit(-1);
    }

    double d1, d2, d3;
    char separator;
    int counter = 0, i =0, j = 0;
    while(counter < Xdim*Ydim){
        data >> d1 >> separator>> d2 >> separator >> d3;
        i = counter%Xdim; j = counter/Xdim;
        toReturn(i,Ydim-1-j) = d3;
        counter++;

    }
    data.close();
    return toReturn;
}

Matrix<double> Ensemble::LoadFieldFileWithoutCoord(Matrix<double>& toReturn, string fileName){

    ifstream data(fileName.c_str());

    if(!data.is_open()){
    printf("Bad Vector Filename: %s\n", fileName.c_str());
    exit(-1);
    }

    double d;
    int counter = 0, i =0, j = 0;
    while(counter < Xdim*Ydim){
        data >> d;
        i = counter%Xdim; j = counter/Xdim;
        toReturn(i,Ydim-1-j) = d;
        counter++;
    }
    data.close();
    return toReturn;
}

void Ensemble::LoadFieldEnsemble(string dir){

    Matrix<double> mem = Matrix<double>(Xdim, Ydim);
    Matrix<double> Fieldsum = Matrix<double>(Xdim, Ydim, 0.0);

    stringstream ss;
    int counter = 0;

    for (int e = 0; e < EnsembleSize; ++e){
        ss.str(string());
        ss.clear();
        ss << dir << "/mem_c_" << e << ".txt";
        MemFileName[counter] = ss.str();

        LoadFieldFileWithoutCoord(mem, ss.str());
        Fieldsum += mem;
        if(counter > EnsembleSize){
            cerr << "wrong size of ensemble: "<< counter << ", " << EnsembleSize << endl;
            exit(-1);
        }

        Members[counter++] = getLevelSetMask(mem, LevelSet);
    }

    if(counter != EnsembleSize){
        cerr << "size mismatch: " << counter << ", " << EnsembleSize <<endl;
        exit(-1);
    }

    FieldAvg = Matrix<double>(Xdim, Ydim);
    for(unsigned int i=0; i<Xdim; ++i)
        for(unsigned int j=0; j<Ydim; ++j)
            FieldAvg(i,j) = Fieldsum(i,j)/EnsembleSize;

}

int Ensemble::getXdim(){
    return Xdim;
}

void Ensemble::setXdim(int Xdim){
    Xdim = Xdim;
}

void Ensemble::setYdim(int Ydim){
    Ydim = Ydim;
}

int Ensemble::getYdim(){
    return Ydim;
}

Matrix<int> Ensemble::getMember(unsigned int index){
    if(index < 0 || index > EnsembleSize){
        throw 1;
        return Matrix<int>();
    }
    return Members[index];
}

void Ensemble::printSize(){
    cout<<"[" << EnsembleSize << "]" <<endl;
}

unsigned int Ensemble::getSize() const{
    return EnsembleSize;
}


ostream& operator<<(ostream& out, const Ensemble& e){
    //TODO: maybe expand
    out<<endl;
    out<<"Ensemble Size: "<<e.EnsembleSize<<endl;
    out<<"Ensemble member dimension: ["<<e.Xdim<<", "<<e.Ydim<<"]"<<endl;
    out<<"Field Average: "<<e.FieldAvg<<endl;
    return out;
}

void WriteEnsembleLog(const Ensemble& e, string logName){

    ofstream out(logName.c_str());
    out<<"Ensemble Size: "<<e.EnsembleSize<<endl;
    out<<"Level Set: " << e.LevelSet << endl;
    out<<"Ensemble member dimension: ["<<e.Xdim<<", "<<e.Ydim<<"]"<<endl;

    out<<"Ensemble Field Average: "<<(!(e.FieldAvg.IsEmpty()) ? "Set" : "Not Set")<<endl;

    for(int i=0; i<e.EnsembleSize; i++)
        out << e.MemFileName[i] << ", " << e.Orders[i] << endl;
}


