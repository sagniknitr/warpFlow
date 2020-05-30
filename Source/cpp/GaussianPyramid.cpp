#include <math.h>
#include "GaussianPyramid.hpp"
 
gaussian_pyramid::gaussian_pyramid(void)
{
    image_pyramid=NULL;
}
 
gaussian_pyramid::~gaussian_pyramid(void)
{
    /*if(image_pyramid!=NULL)
        delete []image_pyramid;*/
 
    delete []image_pyramid;
}
 
void gaussian_pyramid::form_image_pyramid(const DImage &image, double ratio, int minWidth)
{
    /* Ensure that the ratio is within the allowed range */
    if(ratio > 0.98 || ratio < 0.4)
    {
        ratio=0.75;
    }
 
    /* Currently, the levels are chosen based on the image width and the parameter ratio */
    pyramid_levels=log((double)minWidth/image.width())/log(ratio);
    if(image_pyramid != NULL)
    {
        delete []image_pyramid;
    }
 
    image_pyramid = new DImage[pyramid_levels];
    image_pyramid[0].copy_data(image);
 
    double sigma = (1/ratio-1);
    const double gain = -1.3862943611198906;
    int n = gain / log(ratio);
    double sigma_count = sigma * n;
 
    for(int i = 1; i < pyramid_levels; ++i)
    {
        DImage temp;
        
        if(i <= n)
        {
            double s = sigma*i;
            image.GaussianSmoothing(temp, s, s*3);
            temp.imresize(image_pyramid[i], pow(ratio, i));
        }
        else
        {
            image_pyramid[i-n].GaussianSmoothing(temp, sigma_count, sigma_count*3);
            double rate=(double)pow(ratio, i) * image.width() / temp.width();
            temp.imresize(image_pyramid[i], rate);
        }
    }
}
