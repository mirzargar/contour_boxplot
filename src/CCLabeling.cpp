#include "CCLabeling.h"


bool sort_pred(IntPairVector a, IntPairVector b){
    return (a.first < b.first);
}

int getBCCLabelIntArray(Matrix<int> p, int labelNum){
    vector<IntPairVector> labelCount(labelNum,make_pair(0,0));
    for(int i=0; i<labelNum; i++)
        labelCount[i].second = i;

    for(int i=0; i<p.getXdim(); i++)
        for(int j=0; j<p.getYdim(); j++)
            labelCount[p(i,j)].first++;
    sort(labelCount.begin(), labelCount.end(), sort_pred);

    //not to get the background
    if(labelCount.back().second != 0)
        return labelCount.back().second;
    else
        return labelCount[labelCount.size()-2].second;
}

bool PairComp(IntPairVector a, IntPairVector b){
    if (a.first == b.first)
        return (a.second<b.second);
    else
        return (a.first<b.first);
}

//Works on binary images
Matrix<int> CCLabeling(Matrix<int> p){
    int Value = 1, Increment = 1,
            Mark = Value, Difference = Increment,
            M = p.getXdim()+2, N = p.getYdim()+2, NumOfObj = 0;
    vector<IntPairVector> Offset, Neighbors, Index;

    Matrix<int> connected(p.getXdim(), p.getYdim()), image(M, N);

    Offset.push_back(make_pair(-1,-1)); Offset.push_back(make_pair(0,-1)); Offset.push_back(make_pair(1,-1));
    Offset.push_back(make_pair(-1,0)); Offset.push_back(make_pair(1,0));
    Offset.push_back(make_pair(-1,1)); Offset.push_back(make_pair(0,1)); Offset.push_back(make_pair(1,1));

    for(int i=1; i<M-1; i++)
        for(int j=1; j<N-1; j++)
            image(i,j) = p(i-1, j-1);

    for(int i=1; i<M-1; i++)
        for(int j=1; j<N-1; j++){
            if(image(i,j) == 1){
                NumOfObj++;
                Index.clear();
                Index.push_back(make_pair(i,j));
                connected(i-1, j-1) = Mark;
                int bkt = 0;
                while(!Index.empty()){
                    bkt++;
                    for(int k=0; k<Index.size(); k++)
                        image(Index[k].first, Index[k].second) = 0;
                    Neighbors.clear();
                    for(int k=0; k<Offset.size(); k++)
                        for(int l=0; l<Index.size(); l++)
                        Neighbors.push_back(make_pair(Index[l].first+Offset[k].first, Index[l].second+Offset[k].second));

                    sort(Neighbors.begin(), Neighbors.end(), PairComp);
                    vector<IntPairVector>::iterator it;
                    it = unique(Neighbors.begin(), Neighbors.end());
                    Neighbors.resize(distance(Neighbors.begin(), it));

                    Index.clear();
                    for(int k=0; k<Neighbors.size(); k++)
                        if(image(Neighbors[k].first, Neighbors[k].second)){
                            Index.push_back(Neighbors[k]);
                            connected(Neighbors[k].first-1, Neighbors[k].second-1) = Mark;
                        }
                }
                Mark += Difference;
            }
        }
    return connected;
}

//Works on binary images
Matrix<int> getBCC(Matrix<int> Connected){

    Matrix<int> Connected1(Connected.getXdim(), Connected.getYdim(), 0),
            Mask1(Connected.getXdim(), Connected.getYdim(), 1),
            Mask2(Connected.getXdim(), Connected.getYdim(), 1);
    int labelNum, maxBCCVal;

    labelNum = Connected.getMax();
    maxBCCVal = getBCCLabelIntArray(Connected, labelNum+1);

    for(int i=0; i<Connected.getXdim(); i++)
        for(int j=0; j<Connected.getYdim(); j++)
            if(Connected(i,j) == maxBCCVal)
                Mask1(i,j) = 0;

    Connected1 = CCLabeling(Mask1);
    labelNum = Connected1.getMax();
    maxBCCVal = getBCCLabelIntArray(Connected1, labelNum+1);

    for(int i=0; i<Connected1.getXdim(); i++)
        for(int j=0; j<Connected1.getYdim(); j++)
            if(Connected1(i,j) == maxBCCVal)
                Mask2(i,j) = 0;
    return Mask2;
}


