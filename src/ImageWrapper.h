#ifndef __IMAGE_H_INCLUDED__
#define __IMAGE_H_INCLUDED__

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkCastImageFilter.h"

//#include "QuickView.h"
#include <string.h>
#include <iostream>
//#include <tr1/memory>

using namespace std;


///////////////////
/// ITK typdefs
////////////////////
typedef itk::Image< double, 2 >         ImageType;
typedef itk::ImageFileReader<ImageType> ReaderType;

typedef ImageType::SizeType ImageSizeType;
typedef ImageType::Pointer ITKImage;
typedef  itk::ImageFileWriter< ImageType  > WriterType;

typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> filterType;
typedef itk::ImageDuplicator<ImageType> DuplicatorType;
typedef itk::CastImageFilter<ImageType, ImageType> CastType;


class ImageWrapper{
private:
    ITKImage image_;
    unsigned int rows_, cols_;

public:

    //TODO: maybe make ensemble a friend
    ImageWrapper(): image_(ImageType::New()), rows_(0), cols_(0){}
    ImageWrapper(string FileName);
    ImageWrapper(const ImageWrapper& img)
       : image_(img.image_), rows_(img.rows_), cols_(img.cols_){}
    ImageWrapper(const ITKImage& img, int rows, int cols)
        : image_(img), rows_(rows), cols_(cols){}
    ImageWrapper(int rows, int cols, double val);
    ~ImageWrapper(){}

    // operators
    inline double& operator() (int row, int col);
    inline double operator() (int row, int col) const;
    inline ImageWrapper& operator= (const ImageWrapper& img);

    void LoadImage(string FileName);
    void Clear();

    // getter & setters
    void set(int row, int col, double val);
    void set(double val);
    int GetRowSize() const;
    int GetColSize() const;
    double getMaxVal() const;
    ImageSizeType GetSize() const;
    ITKImage& GetImage();
    void getSubImage(std::shared_ptr< ImageWrapper >& subImage, int width, int height, int middleIndxX, int middleIndxY) const; //it's only half width/height
    void getSubImage(std::shared_ptr< ImageWrapper >& subImage, int radius, int middleIndxX, int middleIndxY) const;

    bool IsZeroSubImage(int radius, int middleIndxX, int middleIndxY) const;
    void WriteImage(const string FileName);

    // friends
    friend ostream& operator<< (ostream& OutputStream, const ImageWrapper& Image);
    //friend void WriteImage(const string FileName);
    friend void VisImage(ImageWrapper& Image);
};

inline
double& ImageWrapper::operator() (int row, int col){
    if (row >= rows_ || col >= cols_ || col < 0 || row < 0){
        cout << "ImageWrapper::operator() -> Image subscript out of bound" << endl;
        throw out_of_range("Image subscript out of bound");
    }

    ImageType::IndexType Indx;
    Indx[0] = row;
    Indx[1] = col;

    return image_->GetPixel(Indx);
}

inline
double ImageWrapper::operator() (int row, int col) const{
    if (row >= rows_ || col >= cols_ || col < 0 || row < 0){
        cout << "ImageWrapper::operator() -> mage subscript out of bound" << endl;
        throw out_of_range("Image subscript out of bound");
    }

    ImageType::IndexType Indx;
    Indx[0] = row;
    Indx[1] = col;

    return image_->GetPixel(Indx);
}

inline
ImageWrapper& ImageWrapper::operator= (const ImageWrapper& img){
    if (this != &img){
        image_ = img.image_;
        rows_ = img.rows_;
        cols_ = img.cols_;
    }
    return *this;
}

#endif
