#include "CalcOpticalFlow.hpp"
#include "ImageFiltering.hpp"
#include "GaussianPyramid.hpp"
#include <cstdlib>
#include <iostream>
#include <math.h>
 
using namespace std;
 
Vector<double> calc_optical_flow::lplacian_parameter;
 
calc_optical_flow::calc_optical_flow(void)
{
}
 
calc_optical_flow::~calc_optical_flow(void)
{
}
 
void calc_optical_flow::calculate_gradients(DImage &horizontal_derivative_image, 
                                                    DImage &vertical_derivative_image, 
                                                    DImage &temporal_gradient_image, 
                                                    const DImage &prev_frame, 
                                                    const DImage &curr_frame)
{
    const double smoothing_filter_coefs[5] = 
    {
        0.02,
        0.11,
        0.74,
        0.11,
        0.02
    };
 
    DImage prev_frame_smoothed, curr_frame_smoothed, temp;
 
    prev_frame.horizontal_vertical_filtering(prev_frame_smoothed,smoothing_filter_coefs,2,smoothing_filter_coefs,2);
    curr_frame.horizontal_vertical_filtering(curr_frame_smoothed,smoothing_filter_coefs,2,smoothing_filter_coefs,2);
    temp.copy_data(prev_frame_smoothed);
    temp.scale(0.4);
    temp.add_double_value(curr_frame_smoothed, 0.6);
 
    temp.calc_horizontal_derivative(horizontal_derivative_image, true);
    temp.calc_vertical_derivative(vertical_derivative_image, true);
    temporal_gradient_image.sub(curr_frame_smoothed, prev_frame_smoothed);
}
 
void calc_optical_flow::generate_pixel_mask_in_boundary(DImage &mask, 
                                                                      const DImage &horizontal_flow, 
                                                                      const DImage &vertical_flow,
                                                                      int interval)
{
    int width,height;
    width = horizontal_flow.width();
    height = horizontal_flow.height();
    if(mask.matchDimension(horizontal_flow)==false)
        mask.allocate(width,height);
 
    const double *p_horizontal_flow,*p_vertical_flow;
    
    p_horizontal_flow = horizontal_flow.data();
    p_vertical_flow   = vertical_flow.data();
 
    mask.reset();
    double *p_mask=mask.data();
 
    for(int i = 0; i < height; ++i)
    {
        for(int j = 0; j < width; ++j)
        {
            int offset               = i*width + j;
            double horizontal_offset = j + p_vertical_flow[offset];
            double vertical_offset   = i + p_horizontal_flow[offset];
            
 
            if((horizontal_offset < interval)          || 
                (horizontal_offset > width-1-interval) || 
                (vertical_offset < interval)           || 
                (vertical_offset>height-1-interval))
                {
                    continue;
                }
 
            p_mask[offset]=1;
        }
    }
}
 
void calc_optical_flow::compute_optical_flow(const DImage &prev_frame, 
                                                      const DImage &curr_frame, 
                                                      DImage &warped_curr_frame, 
                                                      DImage &horizontal_flow, 
                                                      DImage &vertical_flow,
                                                      double alpha, 
                                                      int outer_iterations_count, 
                                                      int inner_iterations_count, 
                                                      int sor_iterations_count)
{
    DImage mask, horz_grad_img, ver_grad_img, temp_grad_img;
    
    int width    = prev_frame.width();
    int height   = prev_frame.height();
    int channels = prev_frame.nchannels();
    int total_pixels   = width * height;
 
    for(int outer_iter_indx = 0; outer_iter_indx < outer_iterations_count; ++outer_iter_indx)
    {
        DImage hor_flow_field_deriv(width,height);
        DImage ver_flow_field_deriv(width,height);
        
        calculate_gradients(horz_grad_img, 
                            ver_grad_img,temp_grad_img, 
                            prev_frame,warped_curr_frame);
        
        generate_pixel_mask_in_boundary(mask,horizontal_flow,vertical_flow); 
        
        hor_flow_field_deriv.reset();
        ver_flow_field_deriv.reset();
 
        for(int inner_iter_indx = 0; inner_iter_indx < inner_iterations_count; ++inner_iter_indx)
        {
            DImage hor_deriv_curr_flow_field(width,height);
            DImage ver_deriv_curr_flow_field(width,height);
 
            DImage hor_sec_deriv(width,height);
            DImage ver_sec_deriv(width,height);
            
            DImage hor_deriv_frm_ver_deriv(width,height);
            DImage ver_deriv_frm_hor_deriv(width,height);
 
            DImage phi_img(width,height);
            DImage psi_img(width,height,channels);
            
            DImage derived_img1;
            DImage derived_img2;
            DImage derived_img3;
            DImage derived_img4;
            DImage derived_img5;
            DImage derived_img6; 
            DImage derived_img7;
            DImage derived_img8;
            DImage derived_img9;
            DImage derived_img10;
            
            if(inner_iter_indx==0)
            {
                hor_deriv_curr_flow_field.copy_data(horizontal_flow);
                ver_deriv_curr_flow_field.copy_data(vertical_flow);
            }
            else
            {
                hor_deriv_curr_flow_field.add(horizontal_flow, hor_flow_field_deriv);
                ver_deriv_curr_flow_field.add(vertical_flow, ver_flow_field_deriv);
            }
            
            hor_deriv_curr_flow_field.calc_horizontal_derivative(hor_sec_deriv);
            hor_deriv_curr_flow_field.calc_vertical_derivative(ver_sec_deriv);
            ver_deriv_curr_flow_field.calc_horizontal_derivative(hor_deriv_frm_ver_deriv);
            ver_deriv_curr_flow_field.calc_vertical_derivative(ver_deriv_frm_hor_deriv);
 
            phi_img.reset();
            
            double* phi_img_pixels = phi_img.data();
            double temp;
 
            const double *p_horz_deriv_pixels;
            const double *p_ver_deriv_pixels;
            const double *p_sec_hor_deriv_pixels;
            const double *p_sec_ver_deriv_pixels;
 
            DImage lap_flt_hor_out;
            DImage lap_flt_ver_out;
 
            const double epsilon     = 1e-06;
            const double power_alpha = 0.5;
            const double noise_threhold = 1e-20;
 
            p_horz_deriv_pixels=hor_sec_deriv.data();
            p_ver_deriv_pixels=ver_sec_deriv.data();
            
            p_sec_hor_deriv_pixels=hor_deriv_frm_ver_deriv.data();
            p_sec_ver_deriv_pixels=ver_deriv_frm_hor_deriv.data();
 
            for(int i = 0; i < total_pixels; ++i)
            {
                temp = p_horz_deriv_pixels[i] * p_horz_deriv_pixels[i] + 
                       p_ver_deriv_pixels[i] * p_ver_deriv_pixels[i] + 
                       p_sec_hor_deriv_pixels[i] * p_sec_hor_deriv_pixels[i] + 
                       p_sec_ver_deriv_pixels[i] * p_sec_ver_deriv_pixels[i];
                
                phi_img_pixels[i] = power_alpha / sqrt(temp + epsilon);
            }
 
            psi_img.reset();
 
            double* psi_img_pixels = psi_img.data();
            
            const double *p_tmp_grad_pixels;
            const double *p_hor_flow_field_deriv_pixels;
            const double *p_ver_flow_field_deriv_pixels;
 
            p_horz_deriv_pixels = horz_grad_img.data();
            p_ver_deriv_pixels  = ver_grad_img.data();
            p_tmp_grad_pixels   = temp_grad_img.data();
 
            p_hor_flow_field_deriv_pixels = hor_flow_field_deriv.data();
            p_ver_flow_field_deriv_pixels = ver_flow_field_deriv.data();
 
            if(channels==1)
            {
                /* Grayscale */
                for(int i = 0;i < total_pixels;++i)
                {
                    temp = p_tmp_grad_pixels[i] + 
                           p_horz_deriv_pixels[i] * p_hor_flow_field_deriv_pixels[i] + 
                           p_ver_deriv_pixels[i] * p_ver_flow_field_deriv_pixels[i];
 
                    temp *= temp;
 
                    if(lplacian_parameter[0] >= noise_threhold)
                    {
                        psi_img_pixels[i] = 1.0 / (2.0 * sqrt(temp + epsilon));
                    }
                }
            }
            else
            {
                /* Colour Image */
                for(int i = 0;i < total_pixels; ++i)
                {
                    for(int j = 0;j < channels; ++j)
                    {
                        int index = i*channels + j;
                        
                        temp = p_tmp_grad_pixels[index] + 
                             p_horz_deriv_pixels[index] * p_hor_flow_field_deriv_pixels[i] + 
                             p_ver_deriv_pixels[index] * p_ver_flow_field_deriv_pixels[i];
                        
                        temp *= temp;
 
                        if(lplacian_parameter[j] >= noise_threhold)
                        {
                            psi_img_pixels[index]=1.0 / (2.0 * sqrt(temp + epsilon));
                        }
                    }
                }
            }
 
            derived_img1.mult(psi_img, horz_grad_img, ver_grad_img);
            derived_img2.mult(psi_img, horz_grad_img, horz_grad_img);
            derived_img3.mult(psi_img, ver_grad_img, ver_grad_img);
            derived_img4.mult(psi_img, horz_grad_img, temp_grad_img);
            derived_img5.mult(psi_img, ver_grad_img, temp_grad_img);
 
            if(channels > 1)
            {
                /* Colour Image */
                derived_img1.collapse(derived_img6);
                derived_img2.collapse(derived_img7);
                derived_img3.collapse(derived_img8);
                derived_img4.collapse(derived_img9);
                derived_img5.collapse(derived_img10);
            }
            else
            {
                /* Gray Scale Image */
                derived_img6.copy_data(derived_img1);
                derived_img7.copy_data(derived_img2);
                derived_img8.copy_data(derived_img3);
                derived_img9.copy_data(derived_img4);
                derived_img10.copy_data(derived_img5);
            }
            
            /* laplacian filtering of the current flow field */
            laplacian_filter(lap_flt_hor_out, horizontal_flow, phi_img);
            laplacian_filter(lap_flt_ver_out, vertical_flow, phi_img);
 
            for(int i=0;i<total_pixels;i++)
            {
                derived_img9.data()[i] = -derived_img9.data()[i] - alpha*lap_flt_hor_out.data()[i];
                derived_img10.data()[i] = -derived_img10.data()[i] - alpha*lap_flt_ver_out.data()[i];
            }
 
            /* SOR iterations */
            const double omega = 1.8;
 
            hor_flow_field_deriv.reset();
            ver_flow_field_deriv.reset();
 
            for(int i = 0; i < sor_iterations_count; ++i)
            {
                for(int j = 0; j < height; ++j)
                {
                    for(int k = 0; k < width; ++k)
                    {
                        int index = j*width + k;
                        double sigma = 0;
                        double tau = 0;
                        double scale = 0;
                        double offset = 0;
 
                        if(k > 0)
                        {
                            offset = phi_img_pixels[index-1];
                            sigma += offset*hor_flow_field_deriv.data()[index-1];
                            tau += offset*ver_flow_field_deriv.data()[index-1];
                            scale += offset;
                        }
                        
                        if(k < width-1)
                        {
                            offset = phi_img_pixels[index];
                            sigma += offset*hor_flow_field_deriv.data()[index+1];
                            tau += offset*ver_flow_field_deriv.data()[index+1];
                            scale += offset;
                        }
 
                        if(j > 0)
                        {
                            offset = phi_img_pixels[index-width];
                            sigma += offset*hor_flow_field_deriv.data()[index-width];
                            tau += offset*ver_flow_field_deriv.data()[index-width];
                            scale += offset;
                        }
 
                        if(j < height-1)
                        {
                            offset = phi_img_pixels[index];
                            sigma += offset*hor_flow_field_deriv.data()[index+width];
                            tau += offset*ver_flow_field_deriv.data()[index+width];
                            scale += offset;
                        }
 
                        sigma *= -alpha;
                        tau *= -alpha;
                        scale *= alpha;
 
                        sigma += derived_img6.data()[index] * 
                            ver_flow_field_deriv.data()[index];
 
                        hor_flow_field_deriv.data()[index] = (1.0 - omega) * 
                            hor_flow_field_deriv.data()[index] + 
                            omega / (derived_img7.data()[index] + alpha*0.05 + scale) * 
                            (derived_img9.data()[index] - sigma);
 
                        tau += derived_img6.data()[index] * hor_flow_field_deriv.data()[index];
 
                        ver_flow_field_deriv.data()[index] = (1.0 - omega) * 
                            ver_flow_field_deriv.data()[index] + 
                            omega / (derived_img8.data()[index] + alpha*0.05 + scale) * 
                            (derived_img10.data()[index] - tau);
                    }
                }
            }
        }
 
        horizontal_flow.add(hor_flow_field_deriv);
        vertical_flow.add(ver_flow_field_deriv);
        
        curr_frame.warpImageBicubicRef(prev_frame,warped_curr_frame,horizontal_flow,vertical_flow);
        warped_curr_frame.threshold();
 
        estimate_laplacian_noise(prev_frame,warped_curr_frame, lplacian_parameter);
    }
 
}
 
void calc_optical_flow::estimate_laplacian_noise(const DImage& prev_frame,const DImage& curr_frame, Vector<double>& estimated_noise)
{
    int channels = prev_frame.nchannels();
 
    if(estimated_noise.dim()!=channels)
    {
        estimated_noise.init(channels);
    }
    else
    {
        estimated_noise.reset();
    }
    
    Vector<double> histogram(channels);
    
    for(int i = 0;i < channels; ++i)
    {
        histogram[i] = 0.0;
    }
 
    for(int i = 0; i < prev_frame.total_pixels(); ++i)
    {
        for(int j = 0; j < channels; ++j)
        {
            int index = i*channels + j;
            double residue = fabs(prev_frame.data()[index] - 
                                  curr_frame.data()[index]);
            
            const double min_residue = 0.0;
            const double max_residue = 1000000.0;
            
            if(residue > min_residue && residue < max_residue)
            {
                estimated_noise[j] += residue;
                histogram[j]++;
            }
        }
    }
 
    for(int i = 0; i < channels; ++i)
    {
        if(histogram[i] == 0)
        {
            cout<<"Histogram value cannot be zero!..."<<endl;
            cout<<"Results are wrong!..."<<endl;
            
            return;
        }
        else
        {
            estimated_noise[i] /= histogram[i];
        }
    }
 
    return;
}
 
void calc_optical_flow::laplacian_filter(DImage &output, 
                                                const DImage &input, 
                                                const DImage& weighting_img)
{
    if(output.matchDimension(input) == false)
    {
        output.allocate(input);
    }
    
    output.reset();
 
    if(input.matchDimension(weighting_img) == false)
    {
        cout<<"Image dimension does not match!...\n"<<endl;
        return;
    }
 
    const double *p_in_pixels     = input.data();
    const double *p_weight_pixels = weighting_img.data();
    
    int width  = input.width();
    int height = input.height();
 
    DImage weighted_residual_image(width,height);
    double *p_weighted_residual_pixels = weighted_residual_image.data();
    double *p_out_pixels = output.data();
 
    for(int i = 0;i < height; ++i)
    {
        for(int j = 0; j < width-1; ++j)
        {
            int index = i*width + j;
            p_weighted_residual_pixels[index] = (p_in_pixels[index + 1] - p_in_pixels[index]) * 
                                                 p_weight_pixels[index];
        }
    }
        
    for(int i = 0; i < height; ++i)
    {
        for(int j = 0; j < width; ++j)
        {
            int index = i*width + j;
            
            if(j < (width-1))
            {
                p_out_pixels[index] -= p_weighted_residual_pixels[index];
            }
 
            if(j > 0)
            {
                p_out_pixels[index] += p_weighted_residual_pixels[index - 1];
            }
        }
    }
        
    weighted_residual_image.reset();
 
    
    for(int i = 0; i < height-1; ++i)
    {
        for(int j = 0; j < width; ++j)
        {
            int index = i*width + j;
            p_weighted_residual_pixels[index] = (p_in_pixels[index + width] - p_in_pixels[index]) * 
                                                  p_weight_pixels[index];
        }
    }
 
    for(int i = 0; i < height; ++i)
    {
        for(int j = 0; j < width; ++j)
        {
            int index = i*width + j;
            
            if(i < height-1)
            {
                p_out_pixels[index] -= p_weighted_residual_pixels[index];
            }
            
            if(i > 0)
            {
                p_out_pixels[index] += p_weighted_residual_pixels[index-width];
            }
        }
    }
}
 
void calc_optical_flow::generate_features(DImage &generated_features, const DImage &input)
{
    int width     = input.width();
    int height    = input.height();
    int channels = input.nchannels();
    
    if(channels==1)
    {
        /* Gray Scale */
        generated_features.allocate(input.width(),input.height(), 3);
        DImage horizontal_derivative, vertical_derivative;
 
        input.calc_horizontal_derivative(horizontal_derivative, true);
        input.calc_vertical_derivative(vertical_derivative, true);
        
        double* p_features = generated_features.data();
 
        for(int i = 0; i < height; ++i)
        {
            for(int j = 0; j < width; ++j)
            {
                int index = i*width + j;
                
                p_features[index*3]     = input.data()[index];
                p_features[index*3 + 1] = horizontal_derivative.data()[index];
                p_features[index*3 + 2] = vertical_derivative.data()[index];
            }
        }
    }
    else if(channels == 3)
    {
        DImage desaturated_input;
        input.desaturate(desaturated_input);
 
        generated_features.allocate(input.width(), input.height(), 5);
 
        DImage horizontal_derivative;
        DImage vertical_derivative;
 
        desaturated_input.calc_horizontal_derivative(horizontal_derivative,true);
        desaturated_input.calc_vertical_derivative(vertical_derivative,true);
 
        double* p_features=generated_features.data();
 
        for(int i = 0;i < height; ++i)
        {
            for(int j = 0;j < width; ++j)
            {
                int index=i*width + j;
 
                p_features[index*5]     = desaturated_input.data()[index];
                p_features[index*5 + 1] = horizontal_derivative.data()[index];
                p_features[index*5 + 2] = vertical_derivative.data()[index];
                p_features[index*5 + 3] = input.data()[index*3 + 1] - input.data()[index*3];
                p_features[index*5 + 4] = input.data()[index*3 + 1] - input.data()[index*3 + 2];
            }
        }
    }
    else
    {
        cout<<"Unsupported Channel Format!..."<<endl;
    }
}
 
 
void calc_optical_flow::pyramidal_warping_optical_flow(DImage &horizontal_flow_vector, 
                                                                     DImage &vertical_flow_vector, 
                                                                     DImage &warped_img,
                                                                     const DImage &prev_frame, 
                                                                     const DImage &curr_frame, 
                                                                     double alpha, 
                                                                     double ratio, 
                                                                     int min_width,
                                                                     int outer_iteartions_count, 
                                                                     int inner_iterations_count, 
                                                                     int sor_iterations_factor)
{
    gaussian_pyramid pyramid_prev_frame;
    gaussian_pyramid pyramid_curr_frame;
 
    pyramid_prev_frame.form_image_pyramid(prev_frame,ratio,min_width);
    pyramid_curr_frame.form_image_pyramid(curr_frame,ratio,min_width);
 
    DImage prev_frame_features;
    DImage curr_frame_features;
    DImage warped_curr_frame;
 
    lplacian_parameter.init(prev_frame.nchannels()+2);
 
    for(int i = 0; i < lplacian_parameter.dim(); ++i)
    {
        lplacian_parameter[i] = 0.02;
    }
 
    int pyramid_levels = pyramid_prev_frame.get_pyramid_levels() - 1;
 
    for(int i = pyramid_levels; i >= 0; i--)
    {
        int width=pyramid_prev_frame.Image(i).width();
        int height=pyramid_prev_frame.Image(i).height();
        
        generate_features(prev_frame_features, pyramid_prev_frame.Image(i));
        generate_features(curr_frame_features, pyramid_curr_frame.Image(i));
 
        if(i == pyramid_levels)
        {
            /* Highest Pyramid Level */
            horizontal_flow_vector.allocate(width,height);
            vertical_flow_vector.allocate(width,height);
            
            warped_curr_frame.copy_data(curr_frame_features);
        }
        else
        {
 
            horizontal_flow_vector.imresize(width,height);
            horizontal_flow_vector.scale(1/ratio);
            vertical_flow_vector.imresize(width,height);
            vertical_flow_vector.scale(1/ratio);
 
            curr_frame_features.warpImageBicubicRef(prev_frame_features, 
                                                    warped_curr_frame, 
                                                    horizontal_flow_vector, 
                                                    vertical_flow_vector);  
        }
 
        compute_optical_flow(prev_frame_features,
                             curr_frame_features,
                             warped_curr_frame,
                             horizontal_flow_vector,
                             vertical_flow_vector,
                             alpha,
                             outer_iteartions_count + i,
                             inner_iterations_count,
                             sor_iterations_factor + i*3);
 
    }
 
    curr_frame.warpImageBicubicRef(prev_frame,warped_img,horizontal_flow_vector,vertical_flow_vector);
}
 
 

