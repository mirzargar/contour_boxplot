#include "Utility.h"

//////////////////
// Combinatorics
//////////////////
unsigned long choose(unsigned long n, unsigned long k) {
    if (k > n) {
        return 0;
    }
    unsigned long r = 1;
    for (unsigned long d = 1; d <= k; ++d) {
        r *= n--;
        r /= d;
    }
    return r;
}

template <typename Iterator>
bool next_combination(const Iterator first, Iterator k, const Iterator last){
   /* Credits: Thomas Draper */
   if ((first == last) || (first == k) || (last == k))
      return false;
   Iterator itr1 = first;
   Iterator itr2 = last;
   ++itr1;
   if (last == itr1)
      return false;
   itr1 = last;
   --itr1;
   itr1 = k;
   --itr2;
   while (first != itr1){
      if (*--itr1 < *itr2){
         Iterator j = k;
         while (!(*itr1 < *j)) ++j;
         std::iter_swap(itr1,j);
         ++itr1;
         ++j;
         itr2 = k;
         std::rotate(itr1,j,last);
         while (last != j){
            ++j;
            ++itr2;
         }
         std::rotate(k,itr2,last);
         return true;
      }
   }
   std::rotate(first,k,last);
   return false;
}

Matrix<int> getCombinationMatrix(unsigned int n, unsigned int k){
    std::vector<int> ints;
    for (unsigned int i = 1; i <= n; ints.push_back(i++));
    unsigned long chooseCount = choose(n,k);

    Matrix<int> combinations(chooseCount, k);

    int counter = 0;
    do{
       for (unsigned int i = 0; i < k; ++i)
            combinations(counter,i) = ints[i];
       counter++;
    }while(next_combination(ints.begin(),ints.begin() + k,ints.end()));
    return combinations;
}

Matrix<int> zeroCrossing(Matrix<double> M){
    Matrix<int> toReturn(M.getXdim(), M.getYdim());
    for(unsigned int i=0; i<M.getXdim()-1;++i)
        for(unsigned int j=0; j<M.getYdim()-1;++j){
            double nx = M(i+1, j);
            double ny = M(i, j+1);
            double current = M(i,j);
            if((nx>0 && current<0) || (nx<0 && current>0) || (nx==0 && current>0) || (nx>0 && current==0)){
                if(fabs(nx)<fabs(current))
                    toReturn(i+1,j) = 1;
                else
                    toReturn(i,j) = 1;
            }
            if((ny>0 && current<0) || (ny<0 && current>0) || (ny==0 && current>0) || (ny>0 && current ==0)){
                if(fabs(ny)<fabs(current))
                    toReturn(i, j+1) = 1;
                else
                    toReturn(i, j) = 1;
            }
        }
    return toReturn;
}

Matrix<int> getLevelSetMask(Matrix<double> m, double levelSet){
    Matrix<int> mask = Matrix<int>(m.getXdim(), m.getYdim(), 0);
    for(unsigned int i=0; i<m.getXdim(); ++i)
        for(unsigned int j=0; j<m.getYdim(); ++j){
            mask(i,j) = ((m(i,j) <= levelSet) ? 0 : 1);
        }
        return mask;
}

Matrix<int> getBinaryUnion(Matrix<int> m, Matrix<int> n){
    Matrix<int> temp(m.getXdim(), m.getYdim());
    temp = MaxMatrix(m, n);
    return temp;
}

double getBinaryIntersectSum(const Matrix<int>& m, const Matrix<int>& n){
    if((m.getXdim() != n.getXdim()) || (m.getYdim() != n.getYdim())){
        throw 1;
        return -1;
    }

    double toReturn = 0.0;
    for(unsigned int i=0; i<m.getXdim();++i)
        for(unsigned int j=0; j<m.getYdim();++j){
                if(m(i,j) == 1 && n(i,j) == 1)
                    toReturn += 1;
        }
    return toReturn;
}

Matrix<int> getBinaryIntersect(Matrix<int> m, Matrix<int> n){
    if((m.getXdim() != n.getXdim()) || (m.getYdim() != n.getYdim())){
        throw 1;
        return Matrix<int>();
    }

    Matrix<int> temp(m.getXdim(), m.getYdim());
    for(unsigned int i=0; i<temp.getXdim();++i)
        for(unsigned int j=0; j<temp.getYdim();++j){
            if(m(i,j) == 1 && n(i,j) == 1)
                temp(i,j) = 1;
        }
    return temp;
}

double getPercentBinarySubset(Matrix<int> m, Matrix<int> n){
    if((m.getXdim() != n.getXdim()) || (m.getYdim() != n.getYdim())){
        throw 1;
        return -1.0;
    }

    Matrix<int> diffMatrix(m.getXdim(), m.getYdim());
    diffMatrix = getBinaryIntersect(m, n);
    double temp = 0.0;
    temp = double(diffMatrix.getSum())/double(m.getSum());
    return fabs(temp);
}


string StringCat(string s1, string s2){
    stringstream ss;
    string s3;
    ss.str("");
    ss << s1 << s2;
    s3 = ss.str();
    return s3;
}


/***************************** Map Projectioin Utilities *****************************/
DoublePairVect LambertConformalConicSphere1SP(double R, double lon, double lon0, double lat, double lat0, double SP1){
    //NOT TESTED
    DoublePairVect toReturn;

    double n = sin(SP1);
    double theta = n*(lon-lon0);
    double F =  (cos(SP1)*pow(tan(M_PI/4.0+SP1/2.0), n))/n;
    double rho = (R*F)/(pow(tan(M_PI/4.0+lat/2.0), n));
    double rho0 = (R*F)/(pow(tan(M_PI/4.0+lat0/2.0), n));

    toReturn.first = rho*sin(theta);
    toReturn.second = rho0 - rho*cos(theta);

    return toReturn;
}

DoublePairVect getMapProjSphere(double lon, double lat, double cntre_lat, double cntre_lon, double SP){
    //hard coded Lambeert confomal conic projection for NCEP sphere
    //NOT TESTED

    DoublePairVect toReturn;
    //DoublePairVect proj;

    // NCEP specifications
    double R = 637.1200; //radius of earth
    double lon0 = cntre_lon*M_PI/180; //centre longitude
    double lat0 = cntre_lat*M_PI/180; // centre latitude
    double SP1 = SP*M_PI/180;  // standard parallel

        toReturn = LambertConformalConicSphere1SP(R, M_PI*(lon)/180, lon0, M_PI*(lat)/180, lat0, SP1);
    return toReturn;
}

vector<DoublePairVect> getMapProjSphere(vector<DoublePairVect> pset, double cntre_lat, double cntre_lon, double SP){
    //hard coded Lambeert confomal conic projection for NCEP sphere
    //NOT TESTED

    vector<DoublePairVect> toReturn;
    DoublePairVect proj;


    // NCEP specifications
    double R = 637.1200; //radius of earth
    double lon0 = cntre_lon*M_PI/180; //centre longitude
    double lat0 = cntre_lat*M_PI/180; // centre latitude
    double SP1 = SP*M_PI/180;  // standard parallel

    for(int i=0; i<pset.size(); i++){
        proj = LambertConformalConicSphere1SP(R, M_PI*(pset[i].first)/180, lon0, M_PI*(pset[i].second)/180, lat0, SP1);
        toReturn.push_back(proj);
    }
    return toReturn;
}

DoublePairVect LambertConformalConicEllipse1SP(double a, double e, double lon, double lon0, double lat, double lat0, double SP1){
    //Lambert Conformal Conic map projection for ellipsoide using 1 standard parallel

    DoublePairVect toReturn;

    //map constants
    double m1 = cos(SP1)/sqrt(1-(e*sin(SP1))*(e*sin(SP1)));

    double t0 = sqrt(((1-sin(lat0))/(1+sin(lat0)))*pow((1+e*sin(lat0))/(1-e*sin(lat0)), e));
    double t1 = sqrt(((1-sin(SP1))/(1+sin(SP1)))*pow((1+e*sin(SP1))/(1-e*sin(SP1)), e));

    double n = sin(SP1);

    double F = m1/(n*(pow(t1, n)));

    double rho0 = a*F*(pow(t0, n));

    // point dependent spec.
    double t = sqrt(((1-sin(lat))/(1+sin(lat)))*pow((1+e*sin(lat))/(1-e*sin(lat)), e));

    double rho = a*F*(pow(t, n));

    double theta = n*(lon-lon0);


    toReturn.first = rho*sin(theta);
    toReturn.second = rho0 - rho*cos(theta);

    return toReturn;
}

vector<DoublePairVect> getMapProjEllipsoid(vector<DoublePairVect> pset){
    //hard coded Lambeert confomal conic projection for ellipe having 1 standard parallel
    //the map projection specifications are matching the NCEP specifications
    // using WGS 84 datum  specification coming with natural earth shapefile

    vector<DoublePairVect> toReturn;
    DoublePairVect proj;

    //WGS 84 datum information
    double a = 637.8137; //in kilometers
    double f = double(1/298.257223563);
    double e = sqrt(2*f-f*f);

    // NCEP specifications
    double lon0 = -95*M_PI/180; //centre longitude
    double lat0 = 35*M_PI/180; // centre latitude
    double SP1 = 25*M_PI/180;  // standard parallel

    for(int i=0; i<pset.size(); i++){
        proj = LambertConformalConicEllipse1SP(a, e, M_PI*(pset[i].first)/180, lon0, M_PI*(pset[i].second)/180, lat0, SP1);
        toReturn.push_back(proj);
    }
    return toReturn;
}

bool IsOnBorder(IntPairVector p, IntPairVector border){
    if(p.first == 0 || p.first == border.first)
        return true;
    if(p.second == 0 || p.second == border.second)
        return true;
    return false;
}

unsigned int getState(IntPairVector p1, IntPairVector p2, IntPairVector border){
    if(IsOnBorder(p1, border) && IsOnBorder(p2, border))
        //going along the edge
        return OnBorder;
    if(!IsOnBorder(p1, border) && !IsOnBorder(p2, border)){
        //going along the path
        return OnPath;
    }
    if(!IsOnBorder(p1, border) && IsOnBorder(p2, border))
        //finished a path
        return FinishPath;
    if(IsOnBorder(p1, border) && !IsOnBorder(p2, border))
        //starting a path
        return StartPath;
}


Matrix<int> getCCLabeling(Matrix<int> p, bool recip){

    Matrix<int> m = Matrix<int>(p.getXdim(), p.getYdim(), 1);
    if(recip){
        for(int i=0; i<p.getXdim(); i++)
            for(int j=0; j<p.getYdim(); j++)
                if(p(i,j) == 1)
                    m(i,j) = 0;
    }else{
        for(int i=0; i<p.getXdim(); i++)
            for(int j=0; j<p.getYdim(); j++)
                m(i,j) = p(i,j);
    }


    const unsigned int Dimension = 2;
    typedef unsigned char                       PixelType;
    typedef itk::RGBPixel<unsigned char>         RGBPixelType;
    typedef itk::Image<PixelType, 2>     ImageType;
    typedef itk::Image<RGBPixelType, 2>  RGBImageType;

    //fill ITK image type
    ImageType::Pointer image = ImageType::New();
    typename ImageType::IndexType start;
    start[0] = 0;
    start[1] = 0;

    typename ImageType::SizeType size;
    unsigned int NumRows = m.getXdim();
    unsigned int NumCols = m.getYdim();
    size[0] = NumRows;
    size[1] = NumCols;

    ImageType::RegionType region;// region(start, size);
    region.SetSize(size);
    region.SetIndex(start);

    image->SetRegions(region);
    image->Allocate();


    for(typename ImageType::IndexValueType r = 0; r < m.getXdim(); r++)
      {
        for(typename ImageType::IndexValueType c = 0; c < m.getYdim(); c++)
        {
        typename ImageType::IndexType pixelIndex = {{r,c}};
          if(m(r,c) == 1)
              image->SetPixel(pixelIndex, 255);
          else
              image->SetPixel(pixelIndex, 0);
        //image->SetPixel(pixelIndex, 255);

        }
      }

    typedef itk::Image< unsigned short, Dimension > OutputImageType;

      typedef itk::ConnectedComponentImageFilter <ImageType, OutputImageType >
        ConnectedComponentImageFilterType;

      ConnectedComponentImageFilterType::Pointer connected =
        ConnectedComponentImageFilterType::New ();
      connected->SetInput(image);
      connected->Update();


      Matrix<int> out = Matrix<int>(m.getXdim(), m.getYdim(), 0);


      for(typename ImageType::IndexValueType r = 0; r < m.getXdim(); r++)
        {
          for(typename ImageType::IndexValueType c = 0; c < m.getYdim(); c++)
          {
          typename ImageType::IndexType pixelIndex = {{r,c}};
              out(r,c) =  connected->GetOutput()->GetPixel(pixelIndex);
          }
        }
    return out;
}

bool CCWithHole(Matrix<int> p){
    Matrix<int> cc = getCCLabeling(p, true);
    int maxVal = cc.getMax();
    if(maxVal >1 )
        return true;
    return false;
}

bool ContourHitsEdge(vector<IntPairVector> contour, unsigned int Xdim, unsigned int Ydim){
    for(int i=0; i<contour.size(); i++){
        if(contour[i].first == Xdim-1 || contour[i].first == 0)
            return true;
        if(contour[i].second == Ydim-1 || contour[i].second == 0)
            return true;
    }
    return false;
}
