#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <string.h>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>

using namespace std;

template <class T>
class Matrix{

private:
    unsigned int Xdim, Ydim;
    T **p;
    void allocArrays(){
        p = new T*[Xdim];
        for (unsigned int i=0; i<Xdim; i++)
            p[i] = new T[Ydim];
    }

public:
    friend class Ensemble;

    //constructors
    Matrix();
    Matrix<T>(unsigned int sizeX, unsigned int sizeY);
    Matrix<T>(unsigned int sizeX, unsigned int sizeY, T val);
    ~Matrix<T>();
    Matrix<T>(const Matrix<T>& m);

    //operators
    Matrix<T>& operator=(const Matrix<T>& rhs);
    Matrix<T> operator+(const Matrix<T>& m);
    Matrix<T> operator-(const Matrix<T>& m);
    Matrix<T>& operator+=(const Matrix<T>& m);
    Matrix<T>& operator-=(const Matrix<T>& m);
    T& operator()(int x, int y)const;

    //getters
    void getSize() const;
    unsigned int getXdim() const;
    unsigned int getYdim() const;
    T getSum() const;
    bool IsEmpty() const;
    T getMax();
    T getMin();
    T get(unsigned int i, unsigned int j) const;
    void set(T val);

    // operations
    Matrix<T> Transpose();
    Matrix<T> SubCols(const vector<int>& c);
    Matrix<T> SubRows(const vector<int>& r);
    void loadMatrixFile(unsigned int xdim, unsigned int ydim, string fileName);

    // friends
    friend ostream& operator<<(ostream& out, const Matrix<T>& m){
        out<<endl;
        for(unsigned int i=0;i<m.Xdim;++i){
            for(unsigned int j=0;j<m.Ydim;++j)
                out<<m.p[i][j]<<" ";
            out<<endl;
        }
        return out;
    }

    // NOTE: all of them have one to one relation
    friend Matrix<T> operator*(const Matrix<T>& m1, const Matrix<T>& m2){
        Matrix<T> prod(m1.Xdim, m2.Ydim);
        for(unsigned int i=0;i<prod.Xdim;++i)
            for(unsigned int j=0;j<prod.Ydim;++j)
                for(unsigned int k=0;k<m1.Ydim;++k)
                    prod.p[i][j] += m1.p[i][k]*m2.p[k][j];
        return prod;
    }

    friend Matrix<T> operator*(T c, const Matrix<T>& m2){
        Matrix<T> prod(m2);
        for(unsigned int i=0;i<m2.Xdim;++i)
            for(unsigned int j=0;j<m2.Ydim;++j)
                prod.p[i][j] = c*m2.p[i][j];
        return prod;
    }

    friend Matrix<T> MaxMatrix(const Matrix<T>& m, const Matrix<T>& n){
        if((m.Xdim != n.Xdim) || (m.Ydim != n.Ydim)){
            throw 1;
            return Matrix<T>();
        }

        Matrix<T> temp(m.Xdim, m.Ydim);
        for(unsigned int i=0; i<temp.Xdim;++i)
            for(unsigned int j=0; j<temp.Ydim;++j)
                temp.p[i][j] = max(m.p[i][j], n.p[i][j]);
        return temp;
    }

    friend Matrix<T> MinMatrix(const Matrix<T>& m, const Matrix<T>& n){
        if((m.Xdim != n.Xdim) || (m.Ydim != n.Ydim)){
            throw 1;
            return Matrix<T>();
        }

        Matrix<T> temp(m.Xdim, m.Ydim);
        for(unsigned int i=0; i<temp.Xdim;++i)
            for(unsigned int j=0; j<temp.Ydim;++j)
                temp.p[i][j] = min(m.p[i][j], n.p[i][j]);
        return temp;
    }

    friend bool IsEqual(const Matrix& m, const Matrix& n){
        if((m.Xdim != n.Xdim) || (m.Ydim != n.Ydim)){
            throw 1;
            return false;
        }

        bool flag = true;
        for(unsigned int i=0; i<m.Xdim;++i)
            for(unsigned int j=0; j<m.Ydim;++j)
                if(m.p[i][j] != n.p[i][j]){
                    flag = false;
                    return flag;
                }
        return flag;
    }

    friend Matrix<T> AddConstant(const Matrix<T>& m, T d){
        //Matrix temp(*this); // copy constructor
        Matrix<T> temp(m.Xdim, m.Ydim);
        for (unsigned int i=0;i<temp.Xdim;++i)
            for (unsigned int j=0;j<temp.Ydim;j++)
                temp.p[i][j] = m.p[i][j] + d;
        return (temp);
    }

    friend void WriteMatrixFile(const Matrix<T>& m, const string fileName){
        ofstream data(fileName.c_str());
        unsigned int counter = 0;
        while(counter < m.Xdim*m.Ydim){
            int i = counter/m.Ydim, j = counter%m.Ydim;
            if(counter%m.Ydim == m.Ydim-1)
                data << m.p[i][j] << endl;
            else
                data << m.p[i][j] << " ";
            counter++;
        }
        data.close();
    }

    friend void WriteMatrixSREF(const Matrix<T>& m, const string fileName){
        ofstream data(fileName.c_str());
        for(int j=m.Ydim-1; j>=0; j--)
        for(int i=0; i<m.Xdim; i++)
                data << m.p[i][j] << endl;
        data.close();
    }
};

////////////////////
// constructors
///////////////////
template <class T>
Matrix<T>::Matrix(){
    Xdim=0;
    Ydim=0;
    allocArrays();
}

template <class T>
Matrix<T>::Matrix(unsigned int sizeX, unsigned int sizeY):
Xdim(sizeX),Ydim(sizeY){
    allocArrays();
    for (unsigned int i=0;i<Xdim;i++)
        for(unsigned int j=0;j<Ydim;j++)
            p[i][j] = 0;
}

template <class T>
Matrix<T>::Matrix(unsigned int sizeX, unsigned int sizeY, T val):
Xdim(sizeX),Ydim(sizeY){
    allocArrays();
    for (unsigned int i=0;i<Xdim;i++)
        for(unsigned int j=0;j<Ydim;j++)
            p[i][j] = val;
}

template <class T>
Matrix<T>::Matrix(const Matrix& m)
:Xdim(m.Xdim),Ydim(m.Ydim){
    allocArrays();
    for(unsigned int i=0;i<Xdim;++i)
        for(unsigned int j=0;j<Ydim;++j)
            p[i][j] = m.p[i][j];
}

template <class T>
Matrix<T>::~Matrix(){
    for(unsigned int i=0;i<Xdim;i++)
        delete[] p[i];
    delete [] p;
}

////////////////
// getters
////////////////

template <class T>
void Matrix<T>::getSize() const{
    cout<<"[" << Xdim << "," << Ydim << "]" <<endl;
}

template <class T>
unsigned int Matrix<T>::getXdim() const{
    return Xdim;
}

template <class T>
unsigned int Matrix<T>::getYdim() const{
    return Ydim;
}

template <class T>
T Matrix<T>::getSum() const{
    T temp = 0;
    for(unsigned int i=0; i<Xdim;++i)
        for(unsigned int j=0; j<Ydim;++j)
            temp += p[i][j];
    return temp;
}

template <class T>
bool Matrix<T>::IsEmpty() const{
    if(Xdim == 0 && Ydim == 0)
        return true;
    else
        return false;
}

template <class T>
T Matrix<T>::getMax(){
    T maxVal = -1000;//TODO: put NAN
    for(int i=0; i<Xdim; i++)
        for(int j=0; j<Ydim; j++)
            if(p[i][j] > maxVal)
                maxVal = p[i][j];
    return maxVal;
}

template <class T>
T Matrix<T>::getMin(){
    T minVal = 1000;//TODO: put NAN
    for(int i=0; i<Xdim; i++)
        for(int j=0; j<Ydim; j++)
            if(p[i][j] < minVal)
                minVal = p[i][j];
    return minVal;
}

template<class T>
T Matrix<T>::get(unsigned int i, unsigned int j)const{
    if(i<0 || i> Xdim ||
       j<0 || j> Ydim){
        cerr << "out of bound indices: " << i << ", " << j << endl;
        exit(-1);
    }
    return p[i][j];
}

template<class T>
void Matrix<T>::set(T val){
    for(unsigned int i=0; i<Xdim; i++)
        for(unsigned int j=0; j<Ydim; j++)
                p[i][j] = val;
}

/////////////////
// operators
/////////////////

template <class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m){
    if(this == &m){
        return *this;
    }else{
        if(Xdim!=m.Xdim || Ydim!=m.Ydim){
            if(Xdim != 0 && Ydim != 0)
                this->~Matrix();
            Xdim = m.Xdim; Ydim = m.Ydim;
            allocArrays();
        }
        for(unsigned int i=0;i<Xdim;i++)
            for(unsigned int j=0;j<Ydim;j++)
                p[i][j] = m.p[i][j];
    }
    return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& m){
    for(unsigned int i=0;i<Xdim;++i)
        for(unsigned int j=0;j<Ydim;++j)
            p[i][j] += m.p[i][j];
    return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& m){
    Matrix<T> temp(*this); // copy constructor
    return (temp += m);
}

template <class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& m){
    for(unsigned int i=0;i<Xdim;++i)
        for(unsigned int j=0;j<Ydim;++j)
            p[i][j] -= m.p[i][j];
    return *this;
}

template <class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& m){
    Matrix<T> temp(*this); // copy constructor
    return (temp -= m);
}

template <class T>
T& Matrix<T>::operator()(int i, int j)const{
    return p[i][j];
}

//////////////////
// operations
//////////////////
template <class T>
Matrix<T> Matrix<T>::Transpose(){
    Matrix<T> tm(Ydim, Xdim);
    for(unsigned int i=0;i<tm.Xdim;++i)
        for(unsigned int j=0;j<tm.Ydim;++j)
            tm.p[i][j] = p[j][i];
    return tm;
}

template <class T>
Matrix<T> Matrix<T>::SubCols(const vector<int>& c){
    Matrix<T> m(Xdim, c.size());
    for(std::vector<int>::size_type j = 0; j != c.size(); ++j)
        for(unsigned int i=0;i<m.Xdim;++i)
            m.p[i][j] = p[i][c[j]-1];
    return m;
}

template <class T>
Matrix<T> Matrix<T>::SubRows(const vector<int>& r){
    Matrix<T> m(r.size(), Ydim);
    for(std::vector<int>::size_type i=0; i !=r.size(); ++i)
        for(unsigned int j=0;j<m.Ydim;++j)
            m.p[i][j] = p[r[i]-1][j];
    return m;
}

template <class T>
void Matrix<T>::loadMatrixFile(unsigned int xdim, unsigned int ydim, string fileName){

    ifstream data(fileName.c_str());

    if(!data.is_open()){
    printf("Bad Vector Filename: %s\n", fileName.c_str());
    exit(-1);
    }

    int d;
    int counter = 0;
    while(counter < xdim*ydim){
        data >> d;
        int i = counter/ydim, j = counter%ydim;
        //m(i,j)=d;
        p[i][j] = d;
        counter++;
    }
    data.close();
}
#endif // MATRIX_H
