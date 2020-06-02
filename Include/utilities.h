#ifndef __UTILITIES_H__
#define __UTILITIES_H__

#include "common_datatypes.h"

#ifdef __cplusplus
extern "C" {
#endif

/*Fast logarithm for Image data using LUT */
typedef struct __log_LUT {
  bool isLutpreset = false;
  uint8_t log_image_lut[256];
 }log_LUT_u8;

/*generic log implementation*/
 extern float32_t log(float32_t);

typedef enum _interpolation_type {
	 BILIENAR = 0,
	 BICUBIC, 
	 CUBICSPLINE
 } interpolation_type;

typedef enum _noisemodel {
	LAPLACIAN = 0
} noisemodel;



#ifdef __cplusplus
}
#endif




#endif /*__UTILITIES_H__*/