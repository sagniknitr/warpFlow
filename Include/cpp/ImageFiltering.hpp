#ifndef IMAGE_FILTERING_H
#define IMAGE_FILTERING_H
 
#include <stdio.h>
#include <math.h>
#include "minmax.h"
 
class image_filtering
{
public:
    image_filtering(void);
    ~image_filtering(void);
    
    template <class T>
    static inline T put_within_range(const T& x,const int& max_val) 
    {
        return my_min(my_max(x,0),max_val-1);
    };
 
    template <class T1,class T2>
    static inline void bilinear_interpolation(const T1* p_input,
                                                        int width,
                                                        int height,
                                                        int channels,
                                                        double x,
                                                        double y,
                                                        T2* p_output);
 
    template <class T1,class T2>
    static void img_resizing(const T1* p_input_img,
                                        T2* p_output_img,
                                        int input_width,
                                        int input_height,
                                        int channels,
                                        double ratio);
 
    template <class T1,class T2>
    static void img_resizing(const T1* pSrcImage,
                                    T2* pDstImage,
                                    int SrcWidth,
                                    int SrcHeight,
                                    int nChannels,
                                    int DstWidth,
                                    int DstHeight);
 
    template <class T1,class T2>
    static void horizontal_filtering(const T1* p_inp_img,
                                                    T2* p_outp_img,
                                                    int width,
                                                    int height,
                                                    int channels,
                                                    const double* p_filter,
                                                    int filter_order);
 
    template <class T1,class T2>
    static void vertical_filtering(const T1* pSrcImage,
                                            T2* pDstImage,
                                            int width,
                                            int height,
                                            int nChannels,
                                            const double* pfilter1D,
                                            int fsize);
};
 
template <class T1,class T2>
inline void image_filtering::bilinear_interpolation(const T1* p_input,
                                                              int width,
                                                              int height,
                                                              int channels,
                                                              double x,
                                                              double y,
                                                              T2* p_output)
{
    int x_int = (int)x;
    int y_int = (int)y;
 
    double horizontal_derivative = my_max(my_min(x - x_int, 1), 0);
    double vertical_derivative = my_max(my_min(y - y_int, 1), 0);
 
    int horizontal_index = my_min(my_max(x_int, 0), width - 1);
    int vertical_index_1 = my_min(my_max(y_int, 0), height - 1);
    int index = (vertical_index_1*width + horizontal_index) * channels;
 
    double slope = fabs(1-horizontal_derivative) * fabs(1-vertical_derivative);
 
    for(int k = 0; k < channels; ++k)
    {
        p_output[k] += p_input[index+k]*slope;
    }
            
    int vertical_index_2 = my_min(my_max(y_int + 1, 0), height - 1);
 
    index = (vertical_index_2*width+horizontal_index)*channels;
    slope = fabs(1-horizontal_derivative) * fabs(-vertical_derivative);
 
    for(int k = 0; k < channels; ++k)
    {
        p_output[k]+=p_input[index+k]*slope;
    }
 
    horizontal_index = my_min(my_max(x_int + 1, 0), width - 1);
    
    index = (vertical_index_1*width+horizontal_index)*channels;
 
    slope = fabs(-horizontal_derivative) * fabs(1-vertical_derivative);
 
    for(int k = 0; k < channels; ++k)
    {
        p_output[k] += p_input[index+k]*slope;
    }
 
    index = (vertical_index_2*width+horizontal_index)*channels;
    slope = fabs(-horizontal_derivative) * fabs(-vertical_derivative);
 
    for(int k = 0; k < channels; ++k)
    {
        p_output[k]+=p_input[index+k]*slope;
    }
}
 
template <class T1,class T2>
void image_filtering::img_resizing(const T1* p_input_img,
                                        T2* p_output_img,
                                        int input_width,
                                        int input_height,
                                        int channels,
                                        double ratio)
{
    double owidth  =  (double)input_width * ratio;
    double oheight =  (double)input_height * ratio;
 
    int output_width  =  (int)floor(owidth);
    int output_height =  (int)floor(oheight);
    
    memset(p_output_img, 0 , sizeof(T2)*output_width*output_height*channels);
 
    for(int i=0;i<output_height;i++)
    {
        for(int j=0;j<output_width;j++)
        {
            double horizontal_index=(double)(j+1)/ratio-1;
            double vertical_index=(double)(i+1)/ratio-1;
 
            bilinear_interpolation(p_input_img,
                                   input_width,
                                   input_height,
                                   channels,
                                   horizontal_index,
                                   vertical_index,
                                   p_output_img+(i*output_width+j)*channels);
        }
    }
}
 
template <class T1,class T2>
void image_filtering::img_resizing(const T1 *p_input_img, 
                                        T2 *p_output_img, 
                                        int input_width, 
                                        int input_height, 
                                        int channels, 
                                        int output_width, 
                                        int output_height)
{
    double horizontal_ratio=(double)output_width/input_width;
    double vertical_ratio=(double)output_height/input_height;
    memset(p_output_img,sizeof(T2)*output_width*output_height*channels,0);
 
    for(int i = 0;i < output_height; ++i)
    {
        for(int j = 0;j < output_width; ++j)
        {
            double horizontal_indx = (double)(j + 1) / horizontal_ratio - 1;
            double vertical_indx = (double)(i + 1) / vertical_ratio - 1;
 
            bilinear_interpolation(p_input_img,
                                    input_width,
                                    input_height,
                                    channels,
                                    horizontal_indx,
                                    vertical_indx,
                                    p_output_img + (i*output_width + j) * channels);
        }
    }
}
 
template <class T1,class T2>
void image_filtering::horizontal_filtering(const T1* p_inp_img,
                                                    T2* p_outp_img,
                                                    int width,
                                                    int height,
                                                    int channels,
                                                    const double* p_filter,
                                                    int filter_order)
{
    memset(p_outp_img, 0, sizeof(T2)*width*height*channels);
 
    for(int i = 0;i < height;++i)
    {
        for(int j = 0; j < width; ++j)
        {
            int index = i * width * channels;
            
            T2* p_output = p_outp_img+index+j*channels;
 
            for(int l = -filter_order;l<=filter_order;l++)
            {
                double filter_coef = p_filter[l+filter_order];
                
                int filter_offset_indx = my_min(my_max(j+l, 0), width - 1);
 
                for(int k = 0; k < channels; ++k)
                {
                    p_output[k] += p_inp_img[index + filter_offset_indx*channels + k] * filter_coef;
                }
            }
        }
    }
}
 
template <class T1,class T2>
void image_filtering::vertical_filtering(const T1* p_inp_img,
                                                 T2* p_outp_img,
                                                 int width,
                                                 int height,
                                                 int channels,
                                                 const double* p_filter,
                                                 int filter_order)
{
    memset(p_outp_img, 0, sizeof(T2)*width*height*channels);
    
    for(int i = 0;i < height; ++i)
    {
        for(int j=0; j < width; ++j)
        {
            T2* p_output = p_outp_img+(i*width+j)*channels;
 
            for(int l=-filter_order; l <= filter_order; ++l)
            {
                double filter_coef = p_filter[l+filter_order];
                
                int filter_offset_indx = put_within_range(i+l,height);
 
                for(int k = 0; k < channels; ++k)
                {
                    p_output[k] += p_inp_img[(filter_offset_indx*width + j) * channels + k] * filter_coef;
                }
            }
        }
    }
}
 
#endif
