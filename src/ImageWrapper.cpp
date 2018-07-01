#include "ImageWrapper.h"

//////////////////
/// constructors
//////////////////
///
ImageWrapper::ImageWrapper(string FileName){

    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(FileName.c_str());

    reader->Update();
    image_ = reader->GetOutput();

    ImageType::RegionType region = image_->GetLargestPossibleRegion();
    ImageType::SizeType size = region.GetSize();
    rows_ = size[0];
    cols_ = size[1];
}

void ImageWrapper::WriteImage(const string FileName){
    ImageType::RegionType region;
    ImageType::IndexType start;
    start[0] = 0;
    start[1] = 0;

    ImageType::SizeType size;
    size[0] = this->GetRowSize();
    size[1] = this->GetColSize();

    //typedef unsigned char     PixelType;
    //  const     unsigned int    Dimension = 2;
    //  typedef itk::Image< PixelType, Dimension >  ImageType;

      region.SetSize(size);
      region.SetIndex(start);

      ImageType::Pointer image = ImageType::New();
      image->SetRegions(region);
      image->Allocate();

      ImageType::IndexType Indx;
      for (int i=0; i<size[0]; i++){
          for (int j=0; j<size[1]; j++){
              Indx[0] = i;
              Indx[1] = j;
              image->SetPixel(Indx, image_->GetPixel(Indx));
          }
      }

      typedef  itk::ImageFileWriter< ImageType  > WriterType;
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName(FileName);
      writer->SetInput(image);
      writer->Update();
}

ImageWrapper::ImageWrapper(int rows, int cols, double val) : image_(ImageType::New()){
    //NOTE: operator= was not overloaded so, I couldn't initialized image_ by image_ = ImageType::New();
    // so the only solution was initialization list
    ImageType::RegionType region;
    ImageType::IndexType start;
    start[0] = 0;
    start[1] = 0;

    if (rows <0 || cols <0)
        throw out_of_range("image dimensions can not be negative");

    rows_ = rows;
    cols_ = cols;

    ImageType::SizeType size;
    size[0] = rows_;
    size[1] = cols_;

    typedef unsigned char     PixelType;
      const     unsigned int    Dimension = 2;
      typedef itk::Image< PixelType, Dimension >  ImageType;

      region.SetSize(size);
      region.SetIndex(start);

      image_->SetRegions(region);
      image_->Allocate();

      ImageType::IndexType Indx;
      for (int i=0; i<size[0]; i++){
          for (int j=0; j<size[1]; j++){
              Indx[0] = i;
              Indx[1] = j;
              image_->SetPixel(Indx, val);
          }
      }
}

//////////////////////
/// getter & setter
//////////////////////

ImageSizeType ImageWrapper::GetSize() const{
    ImageSizeType Indx;
    Indx[0] = rows_;
    Indx[1] = cols_;
    return Indx;
}

int ImageWrapper::GetRowSize() const{
    return rows_;
}

int ImageWrapper::GetColSize() const{
    return cols_;
}

void ImageWrapper::set(int row, int col, double val){
    if (row >= rows_ || col >= cols_ || row < 0 || col < 0){
        cout << "ImageWrapper::set -> image subscript out of bound " << row << ", "
             << col << endl;
        throw out_of_range("image subscript out of bound");
    }
    ImageType::IndexType Indx;
    Indx[0] = row;
    Indx[1] = col;
    image_->SetPixel(Indx, val);
}

double ImageWrapper::getMaxVal() const{
    double maxVal = -numeric_limits<double>::max();

     ImageType::IndexType Indx;
     double tmp = -1.0;
    for (int row=0; row<rows_; row++)
        for (int col=0; col<cols_; col++){
            Indx[0] = row;
            Indx[1] = col;
            tmp = image_->GetPixel(Indx);
            if (tmp > maxVal)
                maxVal = tmp;
        }
    return maxVal;
}

void ImageWrapper::set(double val){
    for (int row=0; row<rows_; row++)
        for(int col=0; col<cols_; col++){
            ImageType::IndexType Indx;
            Indx[0] = row;
            Indx[1] = col;
            image_->SetPixel(Indx, val);
        }
}

ITKImage& ImageWrapper::GetImage(){return image_;}

void ImageWrapper::getSubImage(std::shared_ptr< ImageWrapper >& subImage, int width, int height, int middleIndxX, int middleIndxY) const{
    // this function gets the middle of the subImage (middleIndxX, middleIndxY) and *half* of the width and heigth
    // and return the subimage - in case, the subimage goes beyond the boudary of the original image this function
    // only fills in the portion inside the boundary of the image
    subImage = std::shared_ptr<ImageWrapper>(new ImageWrapper(2*width+1, 2*height+1, 0.0));

    ImageType::IndexType Indx;

    for (int row=0; row<2*width+1; row++)
        for(int col=0; col<2*height+1; col++){
            Indx[0] = (middleIndxX) + row - width;
            Indx[1] = (middleIndxY) + col - height;
            if (Indx[0] >= 0 && Indx[0] < rows_ && Indx[1] >=0 && Indx[1] < cols_)
                (*subImage)(row,col) = image_->GetPixel(Indx);
        }
}

void ImageWrapper::getSubImage(std::shared_ptr< ImageWrapper >& subImage, int radius, int middleIndxX, int middleIndxY) const{
    // this function gets the middle of the subImage (middleIndxX, middleIndxY) and the radius of the subImage circle
    // and return the subimage - in case, the subimage goes beyond the boudary of the original image this function
    // only fills in the portion inside the boundary of the image
    subImage = std::shared_ptr<ImageWrapper>(new ImageWrapper(2*radius+1, 2*radius+1, 0.0));

    ImageType::IndexType Indx;

    for (int row=0; row<2*radius+1; row++)
        for(int col=0; col<2*radius+1; col++){
            if (pow(col-radius, 2.0) + pow(row-radius, 2.0) <= pow(radius, 2.0)){
                Indx[0] = (middleIndxX) + row - radius;
                Indx[1] = (middleIndxY) + col - radius;
                if (Indx[0] >= 0 && Indx[0] < rows_ && Indx[1] >=0 && Indx[1] < cols_)
                    (*subImage)(row,col) = image_->GetPixel(Indx);
            }
        }
}

bool ImageWrapper::IsZeroSubImage(int radius, int middleIndxX, int middleIndxY) const{
    // this function check whether all the values in the current is zero

    bool toReturn = true;
    ImageType::IndexType Indx;

    for (int row=0; row<2*radius+1; row++)
        for(int col=0; col<2*radius+1; col++){
            if (pow(col-radius, 2.0) + pow(row-radius, 2.0) <= pow(radius, 2.0)){
                Indx[0] = (middleIndxX) + row - radius;
                Indx[1] = (middleIndxY) + col - radius;
                if (Indx[0] >= 0 && Indx[0] < rows_ && Indx[1] >=0 && Indx[1] < cols_){
                    //(*subImage)(row,col) = image_->GetPixel(Indx);
                    if (image_->GetPixel(Indx) > 0)
                        toReturn = false;
                }
            }
        }
    return toReturn;
}

/////////////////////////
/// other mem funcstions
/////////////////////////

void ImageWrapper::LoadImage(string FileName){//TODO: catch bad file exception
    // read the image file
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(FileName.c_str());

    reader->Update();
    image_ = reader->GetOutput();

    // set image size
    ImageType::RegionType region = image_->GetLargestPossibleRegion();
    ImageType::SizeType size = region.GetSize();
    rows_ = size[0];
    cols_ = size[1];
}

void ImageWrapper::Clear(){
    for (int i=0; i<rows_; i++)
        for (int j=0; j<cols_; j++)
            this->set(i,j, 0.0);
}
