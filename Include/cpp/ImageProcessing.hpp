
#ifndef IMAGE_H
#define IMAGE_H
 
#include "stdio.h"
#include "memory.h"
#include "image_filtering.h"
#include <iostream>
#include <fstream>
#include <typeinfo>
#include "Vector.h"
#include "minmax.h"
 
using namespace std;
 
enum collapse_type{collapse_average,collapse_max,collapse_min};
enum color_type{RGB,BGR,DATA,GRAY};
 
// template class for image
template <class T>
class Image
{
public:
    T* pData;
    int imWidth,imHeight,nChannels;
    int nPixels,nElements;
    bool IsDerivativeImage;
    color_type colorType;
public:
    Image(void);
    Image(int width,int height,int nchannels=1);
    Image(const T& value,int _width,int _height,int _nchannels=1);
    Image(const Image<T>& other);
    ~Image(void);
 
    virtual Image<T>& operator=(const Image<T>& other);
 
    virtual inline void computeDimension(){nPixels=imWidth*imHeight;nElements=nPixels*nChannels;};
 
    virtual void allocate(int width,int height,int nchannels=1);
 
    template <class T1>
    void allocate(const Image<T1>& other);
 
    virtual void clear();
    virtual void reset();
    virtual void copy_data(const Image<T>& other);
 
    template <class T1>
    void copy(const Image<T1>& other);
 
    inline T*& data(){return pData;};
    inline const T*& data() const{return (const T*&)pData;};
    inline int width() const {return imWidth;};
    inline int height() const {return imHeight;};
    inline int nchannels() const {return nChannels;};
    inline int total_pixels() const {return nPixels;};
    inline bool isDerivativeImage() const {return IsDerivativeImage;};
 
    inline color_type colortype() const{return colorType;};
 
    void setColorType(int colorVal) {
    switch (colorVal) {
      case 1: colorType = GRAY; break;
      default: colorType = RGB;
    }
    return;
    }
 
    bool IsFloat () const;
 
    template <class T1>
    bool matchDimension  (const Image<T1>& image) const;
 
    bool matchDimension (int width,int height,int nchannels) const;
 
    inline void setDerivative(bool isDerivativeImage=true){IsDerivativeImage=isDerivativeImage;};
 
    // function of basic image operations
    virtual bool imresize(double ratio);
 
    template <class T1>
    void imresize(Image<T1>& result,double ratio) const;
 
    void imresize(int dstWidth,int dstHeight);
 
    // image IO's
    virtual bool saveImage(const char* filename) const;
    virtual bool loadImage(const char* filename);
    virtual bool saveImage(ofstream& myfile) const;
    virtual bool loadImage(ifstream& myfile);
 
    template <class T1>
    void calc_horizontal_derivative(Image<T1>& image,bool IsAdvancedFilter=false) const;
 
    template <class T1>
    void calc_vertical_derivative(Image<T1>& image,bool IsAdvancedFilter=false) const;
 
    template <class T1>
    void GaussianSmoothing(Image<T1>& image,double sigma,int fsize) const;
 
    template <class T1>
    void imfilter_h(Image<T1>& image,double* filter,int fsize) const;
 
    template <class T1>
    void imfilter_v(Image<T1>& image,double* filter,int fsize) const;
 
    template <class T1>
    void horizontal_vertical_filtering(Image<T1>& image,const double* hfilter,int hfsize,const double* vfilter,int vfsize) const;
 
    // function to desaturating
    template <class T1>
    void desaturate(Image<T1>& image) const;
 
    template <class T1>
    void collapse(Image<T1>& image,collapse_type type = collapse_average) const;
 
    template <class T1,class T2,class T3>
    void mult(const Image<T1>& image1,const Image<T2>& image2,const Image<T3>& image3);
 
    void scale(double value);
 
    template <class T1,class T2>
    void add(const Image<T1>& image1,const Image<T2>& image2);
 
    template <class T1>
    void add_double_value(const Image<T1>& image1,const double value);
 
    template <class T1>
    void add(const Image<T1>& image1);
 
    template <class T1,class T2>
    void sub(const Image<T1>& image1,const Image<T2>& image2);
 
    // function to threshold an image
    void threshold();
 
    // function for bicubic image interpolation
    template <class T1>
    inline void BicubicCoeff(double a[][4],const T* pIm,const T1* pImDx,const T1* pImDy,const T1* pImDxDy,const int offsets[][2]) const;
 
    template <class T1,class T2>
    void warpImageBicubicRef(const Image<T>& ref,Image<T>& output,const Image<T1>& imdx,const Image<T1>& imdy, const Image<T1>& imdxdy,const Image<T2>& vx,const Image<T2>& vy) const;
 
    template <class T1>
    void warpImageBicubicRef(const Image<T>& ref,Image<T>& output,const Image<T1>& vx,const Image<T1>& vy) const;
};
 
 
typedef Image<double> DImage;
 
//------------------------------------------------------------------------------------------
// constructor
//------------------------------------------------------------------------------------------
template <class T>
Image<T>::Image()
{
    pData=NULL;
    imWidth=imHeight=nChannels=nPixels=nElements=0;
    colorType=RGB;
    IsDerivativeImage=false;
}
 
//------------------------------------------------------------------------------------------
// constructor with specified dimensions
//------------------------------------------------------------------------------------------
template <class T>
Image<T>::Image(int width,int height,int nchannels)
{
    imWidth=width;
    imHeight=height;
    nChannels=nchannels;
    computeDimension();
    pData=NULL;
    pData=new T[nElements];
    if(nElements>0)
        memset(pData,0,sizeof(T)*nElements);
    IsDerivativeImage=false;
}
 
template <class T>
Image<T>::Image(const T& value,int _width,int _height,int _nchannels)
{
    pData=NULL;
    allocate(_width,_height,_nchannels);
    setValue(value);
}
 
template <class T>
void Image<T>::allocate(int width,int height,int nchannels)
{
    clear();
    imWidth=width;
    imHeight=height;
    nChannels=nchannels;
    computeDimension();
    pData=NULL;
 
    if(nElements>0)
    {
        pData=new T[nElements];
        memset(pData,0,sizeof(T)*nElements);
    }
}
 
template <class T>
template <class T1>
void Image<T>::allocate(const Image<T1> &other)
{
    allocate(other.width(),other.height(),other.nchannels());
    IsDerivativeImage = other.isDerivativeImage();
    colorType = other.colortype();
}
 
//------------------------------------------------------------------------------------------
// copy constructor
//------------------------------------------------------------------------------------------
template <class T>
Image<T>::Image(const Image<T>& other)
{
    imWidth=imHeight=nChannels=nElements=0;
    pData=NULL;
    copyData(other);
}
 
//------------------------------------------------------------------------------------------
// destructor
//------------------------------------------------------------------------------------------
template <class T>
Image<T>::~Image()
{
    if(pData!=NULL)
        delete []pData;
}
 
//------------------------------------------------------------------------------------------
// clear the image
//------------------------------------------------------------------------------------------
template <class T>
void Image<T>::clear()
{
    if(pData!=NULL)
        delete []pData;
    pData=NULL;
    imWidth=imHeight=nChannels=nPixels=nElements=0;
}
 
//------------------------------------------------------------------------------------------
// reset the image (reset the buffer to zero)
//------------------------------------------------------------------------------------------
template <class T>
void Image<T>::reset()
{
    if(pData!=NULL)
        memset(pData,0,sizeof(T)*nElements);
}
 
//------------------------------------------------------------------------------------------
// copy from other image
//------------------------------------------------------------------------------------------
template <class T>
void Image<T>::copy_data(const Image<T>& other)
{
    imWidth=other.imWidth;
    imHeight=other.imHeight;
    nChannels=other.nChannels;
    nPixels=other.nPixels;
    IsDerivativeImage=other.IsDerivativeImage;
    colorType = other.colorType;
 
    if(nElements!=other.nElements)
    {
        nElements=other.nElements;
        if(pData!=NULL)
            delete []pData;
        pData=NULL;
        pData=new T[nElements];
    }
    if(nElements>0)
        memcpy(pData,other.pData,sizeof(T)*nElements);
}
 
template <class T>
template <class T1>
void Image<T>::copy(const Image<T1>& other)
{
    clear();
 
    imWidth=other.width();
    imHeight=other.height();
    nChannels=other.nchannels();
    computeDimension();
 
    IsDerivativeImage=other.isDerivativeImage();
    colorType = other.colortype();
 
    pData=NULL;
    pData=new T[nElements];
    const T1*& srcData=other.data();
    for(int i=0;i<nElements;i++)
        pData[i]=srcData[i];
}
 
//------------------------------------------------------------------------------------------
// override equal operator
//------------------------------------------------------------------------------------------
template <class T>
Image<T>& Image<T>::operator=(const Image<T>& other)
{
    copy_data(other);
    return *this;
}
 
template <class T>
bool Image<T>::IsFloat() const
{
    if(typeid(T)==typeid(float) || typeid(T)==typeid(double) || typeid(T)==typeid(long double))
        return true;
    else
        return false;
}
 
template <class T>
template <class T1>
bool Image<T>::matchDimension(const Image<T1>& image) const
{
    if(imWidth==image.width() && imHeight==image.height() && nChannels==image.nchannels())
        return true;
    else
        return false;
}
 
template <class T>
bool Image<T>::matchDimension(int width, int height, int nchannels) const
{
    if(imWidth==width && imHeight==height && nChannels==nchannels)
        return true;
    else
        return false;
}
 
//------------------------------------------------------------------------------------------
// resize the image
//------------------------------------------------------------------------------------------
template <class T>
bool Image<T>::imresize(double ratio)
{
    if(pData==NULL)
        return false;
 
    T* pDstData;
    int DstWidth,DstHeight;
    DstWidth=(double)imWidth*ratio;
    DstHeight=(double)imHeight*ratio;
    pDstData=new T[DstWidth*DstHeight*nChannels];
 
    image_filtering::img_resizing(pData,pDstData,imWidth,imHeight,nChannels,ratio);
 
    delete []pData;
    pData=pDstData;
    imWidth=DstWidth;
    imHeight=DstHeight;
    computeDimension();
    return true;
}
 
template <class T>
template <class T1>
void Image<T>::imresize(Image<T1>& result,double ratio) const
{
    int DstWidth,DstHeight;
    DstWidth=(double)imWidth*ratio;
    DstHeight=(double)imHeight*ratio;
    if(result.width()!=DstWidth || result.height()!=DstHeight || result.nchannels()!=nChannels)
        result.allocate(DstWidth,DstHeight,nChannels);
    else
        result.reset();
    image_filtering::img_resizing(pData,result.data(),imWidth,imHeight,nChannels,ratio);
}
 
template <class T>
void Image<T>::imresize(int dstWidth,int dstHeight)
{
    DImage foo(dstWidth,dstHeight,nChannels);
    image_filtering::img_resizing(pData,foo.data(),imWidth,imHeight,nChannels,dstWidth,dstHeight);
    copy_data(foo);
}
 
//------------------------------------------------------------------------------------------
// function of reading or writing images (uncompressed)
//------------------------------------------------------------------------------------------
template <class T>
bool Image<T>::saveImage(const char *filename) const
{
    ofstream myfile(filename,ios::out | ios::binary);
    if(myfile.is_open())
    {
        bool foo = saveImage(myfile);
        myfile.close();
        return foo;
    }
    else
        return false;
}
 
template <class T>
bool Image<T>::saveImage(ofstream& myfile) const
{
    char type[16];
    sprintf(type,"%s",typeid(T).name());
    myfile.write(type,16);
    myfile.write((char *)&imWidth,sizeof(int));
    myfile.write((char *)&imHeight,sizeof(int));
    myfile.write((char *)&nChannels,sizeof(int));
    myfile.write((char *)&IsDerivativeImage,sizeof(bool));
    myfile.write((char *)pData,sizeof(T)*nElements);
    return true;
}
 
template <class T>
bool Image<T>::loadImage(const char *filename)
{
    ifstream myfile(filename, ios::in | ios::binary);
    if(myfile.is_open())
    {
        bool foo = loadImage(myfile);
        myfile.close();
        return foo;
    }
    else
        return false;
}
 
template <class T>
bool Image<T>::loadImage(ifstream& myfile)
{
    char type[16];
    myfile.read(type,16);
 
    if(strcasecmp(type,"uint16")==0)
        sprintf(type,"unsigned short");
    if(strcasecmp(type,"uint32")==0)
        sprintf(type,"unsigned int");
    if(strcasecmp(type,typeid(T).name())!=0)
    {
        cout<<"The type of the image is different from the type of the object!"<<endl;
        return false;
    }
 
    int width,height,nchannels;
    myfile.read((char *)&width,sizeof(int));
    myfile.read((char *)&height,sizeof(int));
    myfile.read((char *)&nchannels,sizeof(int));
 
    if(!matchDimension(width,height,nchannels))
    {
        allocate(width,height,nchannels);
    }
    
    myfile.read((char *)&IsDerivativeImage,sizeof(bool));
    myfile.read((char *)pData,sizeof(T)*nElements);
 
    return true;
}
 
template <class T>
template <class T1>
void Image<T>::calc_horizontal_derivative(Image<T1>& result,bool IsAdvancedFilter) const
{
    if(matchDimension(result)==false)
        result.allocate(imWidth,imHeight,nChannels);
    result.reset();
    result.setDerivative();
    T1*& data=result.data();
    int i,j,k,offset;
    if(IsAdvancedFilter==false)
        for(i=0;i<imHeight;i++)
            for(j=0;j<imWidth-1;j++)
            {
                offset=i*imWidth+j;
                for(k=0;k<nChannels;k++)
                    data[offset*nChannels+k]=(T1)pData[(offset+1)*nChannels+k]-pData[offset*nChannels+k];
            }
    else
    {
        double xFilter[5]={1,-8,0,8,-1};
        for(i=0;i<5;i++)
            xFilter[i]/=12;
        image_filtering::horizontal_filtering(pData,data,imWidth,imHeight,nChannels,xFilter,2);
    }
}
 
//------------------------------------------------------------------------------------------
// function to get y-derivative of the image
//------------------------------------------------------------------------------------------
template <class T>
template <class T1>
void Image<T>::calc_vertical_derivative(Image<T1>& result,bool IsAdvancedFilter) const
{
    if(matchDimension(result)==false)
        result.allocate(imWidth,imHeight,nChannels);
    result.setDerivative();
    T1*& data=result.data();
    int i,j,k,offset;
    if(IsAdvancedFilter==false)
        for(i=0;i<imHeight-1;i++)
            for(j=0;j<imWidth;j++)
            {
                offset=i*imWidth+j;
                for(k=0;k<nChannels;k++)
                    data[offset*nChannels+k]=(T1)pData[(offset+imWidth)*nChannels+k]-pData[offset*nChannels+k];
            }
    else
    {
        double yFilter[5]={1,-8,0,8,-1};
        for(i=0;i<5;i++)
            yFilter[i]/=12;
        image_filtering::vertical_filtering(pData,data,imWidth,imHeight,nChannels,yFilter,2);
    }
}
 
template <class T>
template <class T1>
void Image<T>::GaussianSmoothing(Image<T1>& image,double sigma,int fsize) const
{
    Image<T1> foo;
    // constructing the 1D gaussian filter
    double* gFilter;
    gFilter=new double[fsize*2+1];
    double sum=0;
    sigma=sigma*sigma*2;
    for(int i=-fsize;i<=fsize;i++)
    {
        gFilter[i+fsize]=exp(-(double)(i*i)/sigma);
        sum+=gFilter[i+fsize];
    }
    for(int i=0;i<2*fsize+1;i++)
        gFilter[i]/=sum;
 
    // apply filtering
    horizontal_vertical_filtering(image,gFilter,fsize,gFilter,fsize);
 
    delete[] gFilter;
}
 
template <class T>
template <class T1>
void Image<T>::imfilter_h(Image<T1>& image,double* filter,int fsize) const
{
    if(matchDimension(image)==false)
        image.allocate(imWidth,imHeight,nChannels);
    image_filtering::horizontal_filtering(pData,image.data(),imWidth,imHeight,nChannels,filter,fsize);
}
 
template <class T>
template <class T1>
void Image<T>::imfilter_v(Image<T1>& image,double* filter,int fsize) const
{
    if(matchDimension(image)==false)
        image.allocate(imWidth,imHeight,nChannels);
    image_filtering::vertical_filtering(pData,image.data(),imWidth,imHeight,nChannels,filter,fsize);
}
 
 
template <class T>
template <class T1>
void Image<T>::horizontal_vertical_filtering(Image<T1> &image, const double *hfilter, int hfsize, const double *vfilter, int vfsize) const
{
    if(matchDimension(image)==false)
        image.allocate(imWidth,imHeight,nChannels);
    T1* pTempBuffer;
    pTempBuffer=new T1[nElements];
    image_filtering::horizontal_filtering(pData,pTempBuffer,imWidth,imHeight,nChannels,hfilter,hfsize);
    image_filtering::vertical_filtering(pTempBuffer,image.data(),imWidth,imHeight,nChannels,vfilter,vfsize);
    delete[] pTempBuffer;
}
 
//------------------------------------------------------------------------------------------
//   function for desaturation
//------------------------------------------------------------------------------------------
template <class T>
template <class T1>
void Image<T>::desaturate(Image<T1> &image) const
{
    if(nChannels!=3)
    {
        collapse(image);
        return;
    }
    if(!(image.width()==imWidth && image.height()==imHeight && image.nChannels==1))
        image.allocate(imWidth,imHeight,1);
    T1* data=image.data();
    int offset;
    for(int i=0;i<nPixels;i++)
    {
        offset=i*3;
        if(colorType == RGB)
            data[i]=(double)pData[offset]*.299+pData[offset+1]*.587+pData[offset+2]*.114;
        else
            data[i]=(double)pData[offset]*.114+pData[offset+1]*.587+pData[offset+2]*.299;
    }
}
 
template <class T>
template <class T1>
void Image<T>::collapse(Image<T1> &image,collapse_type type) const
{
    if(!(image.width()==imWidth && image.height()==imHeight && image.nChannels==1))
        image.allocate(imWidth,imHeight,1);
    image.IsDerivativeImage = IsDerivativeImage;
    if(nChannels == 1)
    {
        image.copy(*this);
        return;
    }
    T1* data=image.data();
    int offset;
    double temp;
    for(int i=0;i<nPixels;i++)
    {
        offset=i*nChannels;
        switch(type){
            case collapse_average:
                temp=0;
                for(int j=0;j<nChannels;j++)
                    temp+=pData[offset+j];
                data[i]=temp/nChannels;
                break;
            case collapse_max:
                data[i] = pData[offset];
                for(int j=1;j<nChannels;j++)
                    data[i] = my_max(data[i],pData[offset+j]);
                break;
            case collapse_min:
                data[i] = pData[offset];
                for(int j = 1;j<nChannels;j++)
                    data[i]=my_min(data[i],pData[offset+j]);
                break;
        }
    }
}
 
//------------------------------------------------------------------------------------------
// function to multiply image1, image2 and image3 to the current image
//------------------------------------------------------------------------------------------
template <class T>
template <class T1,class T2,class T3>
void Image<T>::mult(const Image<T1>& image1,const Image<T2>& image2,const Image<T3>& image3)
{
    if(image1.matchDimension(image2)==false || image2.matchDimension(image3)==false)
    {
        cout<<"Error in image dimensions--function Image<T>::Multiply()!"<<endl;
        return;
    }
    if(matchDimension(image1)==false)
        allocate(image1);
 
    const T1*& pData1=image1.data();
    const T2*& pData2=image2.data();
    const T3*& pData3=image3.data();
 
    for(int i=0;i<nElements;i++)
        pData[i]=pData1[i]*pData2[i]*pData3[i];
}
 
template <class T>
void Image<T>::scale(double value)
{
    for(int i=0;i<nElements;i++)
        pData[i]*=value;
}
 
//------------------------------------------------------------------------------------------
// function to add image2 to image1 to the current image
//------------------------------------------------------------------------------------------
template <class T>
template <class T1,class T2>
void Image<T>::add(const Image<T1>& image1,const Image<T2>& image2)
{
    if(image1.matchDimension(image2)==false)
    {
        cout<<"Error in image dimensions--function Image<T>::Add()!"<<endl;
        return;
    }
    if(matchDimension(image1)==false)
        allocate(image1);
 
    const T1*& pData1=image1.data();
    const T2*& pData2=image2.data();
    for(int i=0;i<nElements;i++)
        pData[i]=pData1[i]+pData2[i];
}
 
template <class T>
template <class T1>
void Image<T>::add_double_value(const Image<T1>& image1,const double ratio)
{
    if(matchDimension(image1)==false)
    {
        cout<<"Error in image dimensions--function Image<T>::Add()!"<<endl;
        return;
    }
    const T1*& pData1=image1.data();
    for(int i=0;i<nElements;i++)
        pData[i]+=pData1[i]*ratio;
}
 
template <class T>
template <class T1>
void Image<T>::add(const Image<T1>& image1)
{
    if(matchDimension(image1)==false)
    {
        cout<<"Error in image dimensions--function Image<T>::Add()!"<<endl;
        return;
    }
    const T1*& pData1=image1.data();
    for(int i=0;i<nElements;i++)
        pData[i]+=pData1[i];
}
 
//------------------------------------------------------------------------------------------
// function to subtract image2 from image1
//------------------------------------------------------------------------------------------
template <class T>
template <class T1,class T2>
void Image<T>::sub(const Image<T1> &image1, const Image<T2> &image2)
{
    if(image1.matchDimension(image2)==false)
    {
        cout<<"Error in image dimensions--function Image<T>::Subtract()!"<<endl;
        return;
    }
    if(matchDimension(image1)==false)
        allocate(image1);
 
    const T1*& pData1=image1.data();
    const T2*& pData2=image2.data();
    for(int i=0;i<nElements;i++)
        pData[i]=(T)pData1[i]-pData2[i];
}
 
template <class T>
void Image<T>::threshold()
{
    T ImgMax;
    if(IsFloat())
        ImgMax = 1;
    else
        ImgMax = 255;
    for(int i = 0;i<nPixels*nChannels;i++)
        pData[i] = my_min(my_max(pData[i],0),ImgMax);
}
 
template <class T>
template <class T1>
void Image<T>::BicubicCoeff(double a[][4],const T* pIm,const T1* pImDx,const T1* pImDy,const T1* pImDxDy,const int offsets[][2]) const
{
        a[0][0] = pIm[offsets[0][0]];
        a[1][0] = pImDx[offsets[0][0]];
        a[2][0] = -3*pIm[offsets[0][0]] + 3*pIm[offsets[1][0]] -2*pImDx[offsets[0][0]] - pImDx[offsets[1][0]];
        a[3][0] =   2*pIm[offsets[0][0]] -  2*pIm[offsets[1][0]] +   pImDx[offsets[0][0]] +pImDx[offsets[1][0]];
 
        a[0][1] = pImDy[offsets[0][0]];
        a[1][1] = pImDxDy[offsets[0][0]];
        a[2][1] = -3*pImDy[offsets[0][0]] + 3*pImDy[offsets[1][0]] - 2*pImDxDy[offsets[0][0]] - pImDxDy[offsets[1][0]];
        a[3][1] = 2*pImDy[offsets[0][0]] - 2*pImDy[offsets[1][0]] + pImDxDy[offsets[0][0]] + pImDxDy[offsets[1][0]];
 
        a[0][2] =      -3*pIm[offsets[0][0]]      + 3*pIm[offsets[0][1]]       -2*pImDy[offsets[0][0]]        - pImDy[offsets[0][1]];
        a[1][2] = -3*pImDx[offsets[0][0]] + 3*pImDx[offsets[0][1]] -2*pImDxDy[offsets[0][0]] - pImDxDy[offsets[0][1]];
        a[2][2] =            9*pIm[offsets[0][0]]      -        9*pIm[offsets[1][0]]     -        9*pIm[offsets[0][1]]     +    9*pIm[offsets[1][1]] +
                                6*pImDx[offsets[0][0]]   +    3*pImDx[offsets[1][0]]   -     6*pImDx[offsets[0][1]] -    3*pImDx[offsets[1][1]] +
                                6*pImDy[offsets[0][0]]   -     6*pImDy[offsets[1][0]] +      3*pImDy[offsets[0][1]] -    3*pImDy[offsets[1][1]] +
                            4*pImDxDy[offsets[0][0]] + 2*pImDxDy[offsets[1][0]] + 2*pImDxDy[offsets[0][1]] + pImDxDy[offsets[1][1]];
        a[3][2] =           -6*pIm[offsets[0][0]]      +      6*pIm[offsets[1][0]]     +       6*pIm[offsets[0][1]]     -     6*pIm[offsets[1][1]] +
                            (-3)*pImDx[offsets[0][0]]   -     3*pImDx[offsets[1][0]]   +    3*pImDx[offsets[0][1]] +   3*pImDx[offsets[1][1]] +
                            (-4)*pImDy[offsets[0][0]]   +    4*pImDy[offsets[1][0]]    -    2*pImDy[offsets[0][1]] +   2*pImDy[offsets[1][1]] +
                        (-2)*pImDxDy[offsets[0][0]]  - 2*pImDxDy[offsets[1][0]]   -    pImDxDy[offsets[0][1]]   -  pImDxDy[offsets[1][1]];
 
        a[0][3] =      2*pIm[offsets[0][0]]        - 2*pIm[offsets[0][1]]       + pImDy[offsets[0][0]]        + pImDy[offsets[0][1]];
        a[1][3] = 2*pImDx[offsets[0][0]]  - 2*pImDx[offsets[0][1]]  + pImDxDy[offsets[0][0]] + pImDxDy[offsets[0][1]];
        a[2][3] =           -6*pIm[offsets[0][0]]      +      6*pIm[offsets[1][0]]     +       6*pIm[offsets[0][1]]     -     6*pIm[offsets[1][1]] +
                            (-4)*pImDx[offsets[0][0]]   -     2*pImDx[offsets[1][0]]   +    4*pImDx[offsets[0][1]] +   2*pImDx[offsets[1][1]] +
                            (-3)*pImDy[offsets[0][0]]   +    3*pImDy[offsets[1][0]]    -    3*pImDy[offsets[0][1]] +   3*pImDy[offsets[1][1]] +
                        (-2)*pImDxDy[offsets[0][0]]  -     pImDxDy[offsets[1][0]] -  2*pImDxDy[offsets[0][1]]   -  pImDxDy[offsets[1][1]];
        a[3][3] =            4*pIm[offsets[0][0]]      -        4*pIm[offsets[1][0]]     -        4*pIm[offsets[0][1]]     +    4*pIm[offsets[1][1]] +
                                2*pImDx[offsets[0][0]]   +    2*pImDx[offsets[1][0]]   -     2*pImDx[offsets[0][1]] -    2*pImDx[offsets[1][1]] +
                                2*pImDy[offsets[0][0]]   -     2*pImDy[offsets[1][0]] +      2*pImDy[offsets[0][1]] -    2*pImDy[offsets[1][1]] +
                                pImDxDy[offsets[0][0]] +     pImDxDy[offsets[1][0]] +      pImDxDy[offsets[0][1]] + pImDxDy[offsets[1][1]];
}
 
template <class T>
template <class T1>
void Image<T>::warpImageBicubicRef(const Image<T>& ref,Image<T>& output,const Image<T1>& vx,const Image<T1>& vy) const
{
    double dfilter[3] = {-0.5,0,0.5};
    DImage imdx,imdy,imdxdy;
    imfilter_h(imdx,dfilter,1);
    imfilter_v(imdy,dfilter,1);
    imdx.imfilter_v(imdxdy,dfilter,1);
    warpImageBicubicRef(ref,output,imdx,imdy,imdxdy,vx,vy);
}
 
template <class T>
template <class T1,class T2>
void Image<T>::warpImageBicubicRef(const Image<T>& ref,Image<T>& output,const Image<T1>& imdx,const Image<T1>& imdy,const Image<T1>& imdxdy,
                                                                        const Image<T2>& vx,const Image<T2>& vy) const
{
    T* pIm = pData;
    const T1* pImDx = imdx.data();
    const T1* pImDy = imdy.data();
    const T1* pImDxDy = imdxdy.data();
    int width = vx.width();
    int height = vx.height();
    if(!output.matchDimension(width,height,nChannels))
        output.allocate(width,height,nChannels);
    double a[4][4];
    int offsets[2][2];
 
    T ImgMax;
    if(IsFloat())
        ImgMax = 1;
    else
        ImgMax = 255;
 
    for(int i  = 0; i<height; i++)
        for(int j = 0;j<width;j++)
        {
            int offset = i*width+j;
            double x = j + vx.pData[offset];
            double y = i + vy.pData[offset];
            if(x<0 || x>imWidth-1 || y<0 || y>imHeight-1)
            {
                for(int k = 0; k<nChannels;k++)
                    output.pData[offset*nChannels+k] = ref.pData[offset*nChannels+k];
                continue;
            }
            int x0 = x;
            int y0 = y;
            int x1 = x0+1;
            int y1 = y0+1;
            x0 = my_min(my_max(x0,0),imWidth-1);
            x1 = my_min(my_max(x1,0),imWidth-1);
            y0 = my_min(my_max(y0,0),imHeight-1);
            y1 = my_min(my_max(y1,0),imHeight-1);
 
            double dx = x - x0;
            double dy = y- y0;
            double dx2 = dx*dx;
            double dy2 = dy*dy;
            double dx3 = dx*dx2;
            double dy3 = dy*dy2;
 
 
            for(int k = 0;k<nChannels;k++)
            {
                offsets[0][0] = (y0*imWidth+x0)*nChannels + k;
                offsets[1][0] = (y0*imWidth+x1)*nChannels + k;
                offsets[0][1] = (y1*imWidth+x0)*nChannels + k;
                offsets[1][1] = (y1*imWidth+x1)*nChannels + k;
 
                // set the sampling coefficients
                BicubicCoeff(a,pIm,pImDx,pImDy,pImDxDy,offsets);
 
                // now use the coefficients for interpolation
                output.pData[offset*nChannels+k] = a[0][0] +          a[0][1]*dy +          a[0][2]*dy2 +           a[0][3]*dy3+
                                                                                        a[1][0]*dx +   a[1][1]*dx*dy   + a[1][2]*dx*dy2   + a[1][3]*dx*dy3 +
                                                                                        a[2][0]*dx2 + a[2][1]*dx2*dy + a[2][2]*dx2*dy2 + a[2][3]*dx2*dy3+
                                                                                        a[3][0]*dx3 + a[3][1]*dx3*dy + a[3][2]*dx3*dy2 + a[3][3]*dx3*dy3;
            }
        }
}
 
#endif
