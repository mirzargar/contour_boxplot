#include "ContourExtraction.h"

vector<IntPairVector> InitializeClkNeighbors(){
    vector<IntPairVector> clkNeighbors;
    clkNeighbors.push_back(make_pair(-1,-1)); clkNeighbors.push_back(make_pair(-1,0));
    clkNeighbors.push_back(make_pair(-1,1)); clkNeighbors.push_back(make_pair(0,1));
    clkNeighbors.push_back(make_pair(1,1)); clkNeighbors.push_back(make_pair(1,0));
    clkNeighbors.push_back(make_pair(1,-1)); clkNeighbors.push_back(make_pair(0,-1));
    return clkNeighbors;
}


inline IntPairVector operator+(const IntPairVector p1, const IntPairVector p2){
    IntPairVector result;
    result.first = p1.first + p2.first;
    result.second = p1.second + p2.second;
    return result;
}

inline IntPairVector operator-(const IntPairVector p1, const IntPairVector p2){
    IntPairVector result;
    result.first = p1.first - p2.first;
    result.second = p1.second - p2.second;
    return result;
}

int findIndex(vector<IntPairVector> clkNeighbors, IntPairVector p){
    int index = 100;

    for(unsigned int i=0; i<clkNeighbors.size(); i++)
        if((clkNeighbors[i].first == p.first) && (clkNeighbors[i].second == p.second))
            index = i;
        else
            continue;
    if(index == 100){
        cerr << "Bad vector: " << p.first << ", " << p.second << endl;
        exit(-1);
    }
    return index;
}

bool OutOfBound(IntPairVector P, unsigned int xdim, unsigned int ydim){
    if((P.first < 0) || (P.second < 0) || (P.first >= xdim) || (P.second >= ydim))
        return true;
    else
        return false;
}

// works with binary images assuming 0 is background and 1 is foreground
void findStartPoint(Matrix<int> p, IntPairVector& StartPnt, IntPairVector& PrevNeigh){
    bool found = false;
    PrevNeigh = make_pair(0,0);
    for(unsigned int i=0; i<p.getXdim(); i++){
        for(unsigned int j=0; j<p.getYdim(); j++){
            if(p(i,j) == 0){ //background
                PrevNeigh = make_pair(i,j);
                continue;
            }else if(p(i,j) == 1){ //foreground
                found = true;
                StartPnt.first = i; StartPnt.second = j;
                break;
            }else{
                cerr << "unexpected value: " << p(i,j) << " at position: (" << i << ", " << j << ")" << endl;
                exit(-1);
            }
        }
        if(found)
            break;
    }

    //taking care of boundary case
    if((StartPnt.second == 0))
        if(StartPnt.first != 0){
            PrevNeigh.first = StartPnt.first-1;
            PrevNeigh.second = 0;
        }
}

vector <IntPairVector> getNeighbors(IntPairVector P, IntPairVector PrevNeigh, unsigned int xdim, unsigned int ydim, vector<IntPairVector> clkNeighbors){
    vector<IntPairVector> Neighbors = clkNeighbors;

    int index = findIndex(clkNeighbors, PrevNeigh - P);
    rotate(Neighbors.begin(), Neighbors.begin()+index, Neighbors.end());

    return Neighbors;
}

//works with binary images assuming 0 is background and 1 is foreground
Matrix<int> MooreNeighborTracing(Matrix<int> p){

    //initializations
    vector<IntPairVector> clkNeighbors = InitializeClkNeighbors();
    IntPairVector StartPnt, CurrentPnt, PrevNeigh, P;
    Matrix<int> labels(p.getXdim(), p.getYdim(), 0);
    vector<IntPairVector> ContourIndx;
    vector<IntPairVector> Neighbors = clkNeighbors;

    Matrix<int> temp(p.getXdim()+2, p.getYdim()+2, 0);

    //padd the image to ensure no special boundary checking is required
    for(unsigned int i=0; i<p.getXdim(); i++)
        for(unsigned int j=0; j<p.getYdim(); j++)
            temp(i+1,j+1) = p(i,j);

    findStartPoint(temp, StartPnt, PrevNeigh);
    ContourIndx.push_back((StartPnt));
    int itr = 1;

    Neighbors = getNeighbors(StartPnt, PrevNeigh, p.getXdim(), p.getYdim(), clkNeighbors);

    CurrentPnt = Neighbors[itr]+StartPnt;
    P.first = StartPnt.first; P.second = StartPnt.second;

    itr++;

    while(!((CurrentPnt.first == StartPnt.first) && (CurrentPnt.second == StartPnt.second))){//TODO: apply Jacob's stopping criteria
        if(temp(CurrentPnt.first, CurrentPnt.second) == 1){
            P.first = CurrentPnt.first; P.second = CurrentPnt.second;
            ContourIndx.push_back((P));

            Neighbors = getNeighbors(P, PrevNeigh, p.getXdim(), p.getYdim(), clkNeighbors);

            itr = 1;
            CurrentPnt = Neighbors[itr] + CurrentPnt;
        }else{
            PrevNeigh.first = CurrentPnt.first; PrevNeigh.second = CurrentPnt.second;
            CurrentPnt = P + Neighbors[itr];
            itr++;
        }
    }

    for(unsigned int k=0; k<ContourIndx.size(); k++)
        labels(ContourIndx[k].first-1, ContourIndx[k].second-1) = 1;

    return labels;
}

//works with binary images assuming 0 is background and 1 is foreground
vector<IntPairVector> MooreNeighborTracingIndx(Matrix<int> p){

    //initializations
    vector<IntPairVector> clkNeighbors = InitializeClkNeighbors();
    IntPairVector StartPnt, CurrentPnt, PrevNeigh, P;
    vector<IntPairVector> ContourIndx, toReturn;
    vector<IntPairVector> Neighbors = clkNeighbors;

    Matrix<int> temp(p.getXdim()+2, p.getYdim()+2, 0);

    //padd the image to ensure no special boundary checking is required
    for(unsigned int i=0; i<p.getXdim(); i++)
        for(unsigned int j=0; j<p.getYdim(); j++)
            temp(i+1,j+1) = p(i,j);

    findStartPoint(temp, StartPnt, PrevNeigh);
    ContourIndx.push_back((StartPnt));
    int itr = 1;

    Neighbors = getNeighbors(StartPnt, PrevNeigh, p.getXdim(), p.getYdim(), clkNeighbors);

    CurrentPnt = Neighbors[itr]+StartPnt;
    P.first = StartPnt.first; P.second = StartPnt.second;

    itr++;

    while(!((CurrentPnt.first == StartPnt.first) && (CurrentPnt.second == StartPnt.second))){//TODO: apply Jacob's stopping criteria
        if(itr > Neighbors.size())
            break;
        if(temp(CurrentPnt.first, CurrentPnt.second) == 1){
            P.first = CurrentPnt.first; P.second = CurrentPnt.second;
            ContourIndx.push_back((P));

            Neighbors = getNeighbors(P, PrevNeigh, p.getXdim(), p.getYdim(), clkNeighbors);

            itr = 1;
            //PrevNeigh.first = CurrentPnt.first; PrevNeigh.second = CurrentPnt.second;
            CurrentPnt = Neighbors[itr] + CurrentPnt;
            //PrevNeigh.first = Neighbors[0].first + CurrentPnt.first; PrevNeigh.second = Neighbors[0].second + CurrentPnt.second;
        }else{
            PrevNeigh.first = CurrentPnt.first; PrevNeigh.second = CurrentPnt.second;
            CurrentPnt = P + Neighbors[itr];
            itr++;
        }
    }

    for(unsigned int k=0; k<ContourIndx.size(); k++)
        toReturn.push_back(make_pair(ContourIndx[k].first-1, ContourIndx[k].second-1));
        //labels(ContourIndx[k].first-1, ContourIndx[k].second-1) = 1;

    return toReturn;
}
