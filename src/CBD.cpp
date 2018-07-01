#include "CBD.h"


CBD::CBD(string src_dir, unsigned int j, double levelSet, string des_dir, unsigned int ESize, unsigned int Xdim, unsigned int Ydim){

    E.InitializeEnsemble(ESize, levelSet, Xdim, Ydim);

    cout << "*** about to load src_dir/mem_c_0.txt ... src_dir/mem_c_ens_size.txt" << endl;

    E.LoadFieldEnsemble(src_dir);

    cout << "*** ensmeble loaded successfully" << endl;

    // NOTE: This the coordinate file for SREF ensemble - only required if you are working with SREF
    Load_lonLat("lonLat_212.txt");

    double xmin = lon(0,0), ymax = lat(0,0);
    double ymin = lat(0,lat.getYdim()-1);
    double xmax = lon(lon.getXdim()-1, 0);

    cout << "*** about to compute the depth values" << endl;

    numJ = j;
    Dir = src_dir;
    getTimePercents(ESize);
    // computing the threshold based on the argument made in original "contour boxplot" paper
    // in this case number 3 specifies the original guess for the number of outliers
    // alternatively (and probably better to use Tukey's rule)
    PThreshold = getThreshold(3);

    getPercentProbabilities();

    cout << "*** computing depth values are done!" << endl;
    cout << "*** writing the depth values and orders: analysis_CBD.txt" << endl;

    // NOTE: the number of outliers has been specified as the number of members with zero depth value
    // if every member has a non-zero depth value then it has been set to three - one can alternatively use
    // Tukey's rule
    unsigned int outlierNum = 0;
    for (auto bd : Probs)
        if (bd <= 1e-12)
            outlierNum++;

    outlierNum = outlierNum == 0 ? 3 : outlierNum;
    cout << "    outlier number: " << outlierNum << endl;

    WriteSortedProbabilities("analysis_CBD.txt");

    cout << "*** generating the SVG file: des_dir/CBP.svg" << endl;

    // NOTE: this code is not handling holes properly
    genSVG(outlierNum, des_dir, xmin, ymin, xmax, ymax);

    cout << "*** SVG file generated successfully!!!" << endl;
}

bool sort_predD(const my_pair& left, const my_pair& right){
    return left.first < right.first;
}

void CBD::Load_lonLat(string fileName){
    //TODO: Not stable
    ifstream data(fileName.c_str());
    DoublePairVect prj;

    if(!data.is_open()){
    printf("Bad Vector Filename: %s\n", fileName.c_str());
    exit(-1);
    }
    lon = Matrix<double>(E.getXdim(), E.getYdim(), 0.0);
    lat = Matrix<double>(E.getXdim(), E.getYdim(), 0.0);

    double d1, d2;
    char separator;
    int counter = 0, i =0, j = 0;
    double cntre_lat, cntre_lon, SP;

    while(counter < E.getXdim()*E.getYdim()){
        data >> d1 >> separator>> d2;
        i = counter%E.getXdim(); j = counter/E.getXdim();
            cntre_lat = 35.0;
            cntre_lon = -95.0;
            SP = 25.0;
            //NEEDS TO BE CHANGE IF MAP CHANGE
            shift_x = -123.721;
            shift_y = -111.456;
        prj = getMapProjSphere(-(95.0+265.0-d1), d2, cntre_lat, cntre_lon, SP);
        //the shift required to have the center to be center of the map (the center for the map is 95W, 35N)
        lon(i,E.getYdim()-1-j) = prj.first;//-(95.0+265.0-d1);
        lat(i,E.getYdim()-1-j) = prj.second;//2;
        counter++;
    }

    data.close();
}


vector<int> getSortedIndex(const vector<double>& v){
    vector<my_pair> VIndx;
    for(unsigned int i=0; i<v.size(); i++)
        VIndx.push_back(my_pair(v[i],i));
    sort(VIndx.begin(), VIndx.end(), sort_predD);
    vector<int> indx;
    for(unsigned int i=0; i<VIndx.size(); i++)
        indx.push_back(VIndx[i].second);
    return indx;
}

Matrix<int> CBD::getUBand(const Ensemble e, const vector<int>& indx){
    Matrix<int> temp(e.Xdim, e.Ydim);
    temp = e.Members[indx[0]];
    for(unsigned int i=1;i<indx.size();i++)
        temp = getBinaryUnion(temp, e.Members[indx[i]]);
    return temp;
}

Matrix<int> CBD::getIBand(const Ensemble e, const vector<int>& indx){
    Matrix<int> temp(e.Xdim, e.Ydim);
    temp = e.Members[indx[0]];
    for(unsigned int i=1;i<indx.size();i++)
        temp = getBinaryIntersect(temp, e.Members[indx[i]]);
    return temp;
}


void CBD::getPercentProbabilities(){
    if(Probs.size() != 0)
        cout<<"The probabilities are about to be overwritten"<<endl;
    Probs.clear();
    double temp;
    PThreshold = PThreshold;
    for(unsigned int i=0; i<UPercents.getXdim(); ++i){
        temp = 0.0;
        for(unsigned int j=0; j<UPercents.getYdim(); ++j)
            temp += (((UPercents(i,j) >= PThreshold) && (IPercents(i,j) >= PThreshold))? 1.0 : 0.0);
        Probs.push_back(temp/double(UPercents.getYdim()));
    }
}

void CBD::getTimePercents(unsigned int ESize){
    if((E.getSize())%ESize != 0){
        cerr << "size mismatch: " << E.getSize() << ", " << ESize << endl;
        exit(-1);
    }

    Matrix<int> Combinations = getCombinationMatrix(ESize, numJ);
    Matrix<int> UBD = Matrix<int>(E.Xdim, E.Ydim);
    Matrix<int> IBD = Matrix<int>(E.Xdim, E.Ydim);

    Matrix<double> UEnumerator = Matrix<double>(ESize, Combinations.getXdim(), 0.0), UDenumerator = Matrix<double>(ESize, Combinations.getXdim(), 0.0);
    Matrix<double> IDenumerator = Matrix<double>(ESize, Combinations.getXdim(), 0.0), IEnumerator = Matrix<double>(ESize, Combinations.getXdim(), 0.0);
    UPercents = Matrix<double>(ESize, Combinations.getXdim(), 0.0);
    IPercents = Matrix<double>(ESize, Combinations.getXdim(), 0.0);

    //#pragma omp parallel for schedule(dynamic,1)
        for(unsigned int j=0; j<Combinations.getXdim(); ++j){
            // NOTE: hard coded j = 2
                vector<int> JSet; JSet.push_back(Combinations(j,0)-1); JSet.push_back(Combinations(j,1)-1);
                UBD = getUBand(E, JSet);
                IBD = getIBand(E, JSet);
                for(unsigned int k=0; k<ESize; ++k){
                    // NOTE: hard coded j = 2
                        if((JSet[0] == k) || (JSet[1] == k))
                            continue;
                        Matrix<int> current = E.Members[k];
                        UEnumerator(k,j) += double(getBinaryIntersect(current, UBD).getSum());
                        UDenumerator(k,j) += double(current.getSum());

                        IEnumerator(k,j) += double(getBinaryIntersect(IBD, current).getSum());
                        IDenumerator(k,j) += double(IBD.getSum());
                }
        }
    //}
    for(int i=0; i<ESize; i++)
        for(int j=0; j<Combinations.getXdim(); j++){
            if(UDenumerator(i,j) == 0.0 || IDenumerator(i,j) == 0){
                UPercents(i,j) = 0.0;
                IPercents(i,j) = 0.0;
            }else{
                UPercents(i,j) = double(UEnumerator(i,j))/double(UDenumerator(i,j));
                IPercents(i,j) = double(IEnumerator(i,j))/double(IDenumerator(i,j));
            }
        }
}

double CBD::getThreshold(int outlierNum){
    double toReturn = 0.0;
    vector<double> UV(UPercents.getXdim(), 0.0), IV(IPercents.getXdim(), 0.0);

    for(int i=0; i<UPercents.getXdim(); i++)
        for(int j=0; j<UPercents.getYdim(); j++){
            if(UPercents(i,j) > UV[i])
                UV[i] = UPercents(i,j);
            if(IPercents(i,j) > IV[i])
                IV[i] = IPercents(i,j);
    }
    sort(UV.begin(), UV.end());
    sort(IV.begin(), IV.end());

    if(outlierNum > UV.size()){
        cerr << "Outlier num is bigger than the size of ensemble: " << outlierNum << ", "  << UV.size() << endl;
        exit(-1);
    }

    toReturn = min(UV[outlierNum-1], IV[outlierNum-1]);
    return toReturn;
}

void CBD::WriteSortedProbabilities(string fileName){
    ofstream data(fileName.c_str());
    int counter = Probs.size();
    vector<int> SortedIndices = getSortedIndex(Probs);
    vector<double> SortedProbs = Probs;
    sort(SortedProbs.begin(), SortedProbs.end());
    data << "depth value, " << "member_indx" << endl;
    while(--counter>= 0)
        data << SortedProbs[counter] << " " << SortedIndices[counter] << endl;;
}

Matrix<int> getMask(Matrix<int> m){
    Matrix<int> toReturn = Matrix<int>(m.getXdim(), m.getYdim());
    for(int i=0; i<toReturn.getXdim(); i++)
        for(int j=0; j<toReturn.getYdim(); j++)
            toReturn(i,j) = int((m(i,j)-1)*.5);
    return toReturn;
}

Matrix<int> CBD::getBandMask(const Ensemble e, vector<int> v){
    Matrix<int> U = getUBand(e, v);
    Matrix<int> I = getIBand(e, v);
    Matrix<int> toReturn = Matrix<int>(U.getXdim(), U.getYdim(), 0);

    for(int i=0; i<toReturn.getXdim(); i++)
        for(int j=0; j<toReturn.getYdim(); j++){
            toReturn(i,j) = int((U(i,j) - I(i,j)));
        }
    return toReturn;
}

Matrix<int> CBD::getTrimmedMean(Ensemble e, vector<int> v){
    Matrix<double> TMean = Matrix<double>(e.getXdim(), e.getYdim(), 0.0);
    Matrix<int> toReturn = Matrix<int>(e.getXdim(), e.getYdim(), 0);
    Matrix<int> current = Matrix<int>(e.getXdim(), e.getYdim(), 0);
    for(int k=0; k<v.size(); k++){
        current = e.Members[v[k]];
        for(int i=0; i<e.getXdim(); i++)
            for(int j=0; j<e.getYdim(); j++)
                TMean(i,j) = TMean(i,j) + double(current(i,j))/double(v.size());
    }

    for(int i=0; i<e.getXdim(); i++)
        for(int j=0; j<e.getYdim(); j++){
            if(TMean(i,j) >= .5)
                toReturn(i,j) = 1;
        }
    return toReturn;
}

void CBD::genSVG(unsigned int outlierNum, string des_dir, double xmin, double ymin, double xmax, double ymax){
    Ensemble currentEnsemble;
    stringstream ss;
    string filename = "";
    unsigned int num = E.EnsembleSize;

    currentEnsemble.InitializeEnsemble(num, E.LevelSet, E.Xdim, E.Ydim);
    currentEnsemble.LoadFieldEnsemble(Dir);
    ss.str(string());
    ss.clear();
    ss << des_dir << "/CBP.svg";
    filename = ss.str();
    InitializeSVG(filename, 0, 0, xmax-xmin, ymax-ymin);

        vector<int> SortedIndices = getSortedIndex(Probs);
        vector<double> SortedProbs = Probs;
        sort(SortedProbs.begin(), SortedProbs.end());

        Matrix<int> mask = Matrix<int>(E.getXdim(), E.getYdim());

        cout << "   about to draw" << endl;

        //*******************
        //***Drawing outliers
        //*******************

        for(int i=0; i<outlierNum; i++){
            mask = currentEnsemble.Members[SortedIndices[i]];
            DrawCCSilhouette(filename, mask, "red", 1, "dash", "none", lon, lat, xmin, ymin, true, E.getXdim(), E.getYdim());

        }
        cout<<"      done drawing outliers"<<endl;

        //********************
        //***Draw envelop
        //********************
        vector<int> envelop (SortedIndices.begin()+outlierNum, SortedIndices.end());
        mask = getBandMask(currentEnsemble, envelop);
        DrawConnectedComp(filename, mask, "#bbbbff", 1, .5, lon, lat, xmin, ymin, true);
        cout<<"      done drawing envelop"<<endl;


        //********************
        //***Drawing band
        //********************
        vector<int> BWidth (SortedIndices.end()-floor(SortedIndices.size()/2), SortedIndices.end());
        mask = getBandMask(currentEnsemble, BWidth);//getUBand(BWidth);
        DrawConnectedComp(filename, mask, "#8888bb", 1, .5, lon, lat, xmin, ymin, true);
        cout<<"      done drawing band"<<endl;


        //********************
        //***Drawing Trimmead mean
        //********************
        vector<int> vv (SortedIndices.begin(), SortedIndices.end());
        mask = getTrimmedMean(currentEnsemble, vv);
        DrawCCSilhouette(filename, mask, "#cc00ff", 1, "fill", "none", lon, lat, xmin, ymin, true, E.getXdim(), E.getYdim());
        cout<<"      done drawing the trimmed mean"<<endl;


        //******************
        //***Drawing median
        //******************
        mask = currentEnsemble.Members[SortedIndices.back()];
        DrawCCSilhouette(filename, mask, "yellow", 1, "fill", "none", lon, lat, xmin, ymin, true, E.getXdim(), E.getYdim());


        CloseSVGFile(filename);
}

