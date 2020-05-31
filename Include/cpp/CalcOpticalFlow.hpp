#ifndef CALC_OPTICAL_FLOW_H
#define CALC_OPTICAL_FLOW_H
 
#include "Image.hpp"
 
class calc_optical_flow
{
public:
    calc_optical_flow(void);
    ~calc_optical_flow(void);
 
    static Vector <double> lplacian_parameter;
 
    static void calculate_gradients(DImage &horizontal_derivative_image, 
                                            DImage &vertical_derivative_image, 
                                            DImage &temporal_gradient_image, 
                                            const DImage &prev_frame, 
                                            const DImage &curr_frame);
    
    static void generate_pixel_mask_in_boundary(DImage &mask, 
                                                              const DImage &horizontal_flow, 
                                                              const DImage &vertical_flow,
                                                              int interval = 0);
 
    static void compute_optical_flow(const DImage &prev_frame, 
                                      const DImage &curr_frame, 
                                      DImage &warped_curr_frame, 
                                      DImage &horizontal_flow, 
                                      DImage &vertical_flow,
                                      double alpha, 
                                      int outer_iterations_count, 
                                      int inner_iterations_count, 
                                      int SOR_iterations_count);
 
    static void estimate_laplacian_noise(const DImage& prev_frame,
                                                    const DImage& curr_frame,
                                                    Vector<double>& estimated_noise);
    
    static void laplacian_filter(DImage& output,
                                        const DImage& input,
                                        const DImage& weighting_img);
 
    static void pyramidal_warping_optical_flow(DImage &horizontal_flow_vector, 
                                                             DImage &vertical_flow_vector, 
                                                             DImage &warped_img,
                                                             const DImage &prev_frame, 
                                                             const DImage &curr_frame, 
                                                             double alpha, 
                                                             double ratio, 
                                                             int min_width,
                                                             int outer_iteartions_count, 
                                                             int inner_iterations_count, 
                                                             int sor_iterations_factor);
 
    static void generate_features(DImage& imfeature,const DImage& im);
};
 
#endif //CALC_OPTICAL_FLOW_H
