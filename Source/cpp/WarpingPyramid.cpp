#include "WarpingPyramid.hpp"
#include "Image.hpp"
#include "CalcOpticalFlow.hpp"
using namespace std;
 
void warping_pyramid(double * p_horizontal_flow, 
                            double * p_vertical_flow, 
                            double * p_warped_image,
                            const double * p_prev_frame, 
                            const double * p_curr_frame,
                            double alpha, 
                            double ratio, 
                            int min_width,
                            int outer_iteartions_count, 
                            int inner_iterations_count,
                            int sor_iterations_count,
                            int colour_format,
                            int height, 
                            int width, 
                            int channels) {
  DImage prev_frame_flattened;
  DImage curr_frame_flattened;
  DImage horizontal_flow;
  DImage vertical_flow;
  DImage warped_image;
 
  prev_frame_flattened.allocate(width, height, channels);
  curr_frame_flattened.allocate(width, height, channels);
 
  memcpy(prev_frame_flattened.pData, p_prev_frame, height * width * channels * sizeof(double));
  memcpy(curr_frame_flattened.pData, p_curr_frame, height * width * channels * sizeof(double));
 
  prev_frame_flattened.setColorType(colour_format);
  curr_frame_flattened.setColorType(colour_format);
 
  calc_optical_flow::pyramidal_warping_optical_flow(horizontal_flow,
                                                    vertical_flow,
                                                    warped_image,
                                                    prev_frame_flattened,
                                                    curr_frame_flattened,
                                                    alpha,
                                                    ratio, 
                                                    min_width,
                                                    outer_iteartions_count, 
                                                    inner_iterations_count,
                                                    sor_iterations_count);
 
 
  memcpy(p_horizontal_flow, horizontal_flow.pData, height * width * sizeof(double));
  memcpy(p_vertical_flow, vertical_flow.pData, height * width * sizeof(double));
  memcpy(p_warped_image, warped_image.pData, height * width * channels * sizeof(double));
 
  prev_frame_flattened.clear();
  curr_frame_flattened.clear();
  horizontal_flow.clear();
  vertical_flow.clear();
  warped_image.clear();
 
  return;
}
 




 
 
 
 

