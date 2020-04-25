#include "warpFlow.h"

inline void bilinear_interpolation(const uint8_t pu8_image[], 
	                               const int32_t s32_width,
                                   int32_t s32_height,
	                               int32_t nChannels,
                                   float32_t  x, float32_t  y,
                                   uint8_t result[]) {

  int32_t xx, yy, m, n, u, v, l, offset;
  xx = x;
  yy = y;
  float32_t  dx, dy, s;
  dx = __max(__min(x - xx, 1), 0);
  dy = __max(__min(y - yy, 1), 0);

  for (m = 0; m <= 1; m++)
    for (n = 0; n <= 1; n++) {
      u = EnforceRange(xx + m, s32_width);
      v = EnforceRange(yy + n, s32_height);
      offset = (v * s32_width + u) * nChannels;
      s = fabs(1 - m - dx) * fabs(1 - n - dy);
      for (l = 0; l < nChannels; l++) result[l] += pu8_image[offset + l] * s;
    }
}



/*-------------------------------------------------------------------------------------------------*/
// function to interplate multi-channel image plane for (x,y)
/* --------------------------------------------------------------------------------------------------*/

inline void bilinear_interpolation_transpose(
    const uint8_t pInput[],
	int32_t width,
	int32_t height,
	int32_t nChannels,
	float32_t  x, 
	float32_t  y,
    uint8_t pDstImage[]) {

  int32_t xx, yy, m, n, u, v, l, offset;
  xx = x;
  yy = y;
  float32_t  dx, dy, s;
  dx = __max(__min(x - xx, 1), 0);
  dy = __max(__min(y - yy, 1), 0);

  for (m = 0; m <= 1; m++)
    for (n = 0; n <= 1; n++) {
      u = EnforceRange(xx + m, width);
      v = EnforceRange(yy + n, height);
      offset = (v * width + u) * nChannels;
      s = fabs(1 - m - dx) * fabs(1 - n - dy);
      for (l = 0; l < nChannels; l++) pDstImage[offset + l] += pInput[l] * s;
    }
}

void resize_Image(const uint8_t pSrcImage[],
	uint8_t pDstImage[],
    int32_t SrcWidth,
	int32_t SrcHeight, 
	int32_t nChannels,
    float32_t  Ratio) {

  int32_t DstWidth, DstHeight;
  DstWidth = (float32_t )SrcWidth * Ratio;
  DstHeight = (float32_t )SrcHeight * Ratio;
  memset(pDstImage, 0, sizeof(uint8_t) * DstWidth * DstHeight * nChannels);

  float32_t  x, y;

  for (int32_t i = 0; i < DstHeight; i++)
    for (int32_t j = 0; j < DstWidth; j++) {
      x = (float32_t )(j + 1) / Ratio - 1;
      y = (float32_t )(i + 1) / Ratio - 1;

      // bilinear interpolation
      bilinear_interpolation(pSrcImage, SrcWidth, SrcHeight, nChannels, x, y,
                          pDstImage + (i * DstWidth + j) * nChannels);
    }
}


void resize_image_adaptive(const uint8_t pSrcImage,
	                        uint8_t pDstImage,
                            int32_t SrcWidth,
	                         int32_t SrcHeight,
	                        int32_t nChannels,
                                  int32_t DstWidth,
	int32_t DstHeight) {
  float32_t  xRatio = (float32_t )DstWidth / SrcWidth;
  float32_t  yRatio = (float32_t )DstHeight / SrcHeight;
  memset(pDstImage, sizeof(uint8_t) * DstWidth * DstHeight * nChannels, 0);

  float32_t  x, y;

  for (int32_t i = 0; i < DstHeight; i++)
    for (int32_t j = 0; j < DstWidth; j++) {
      x = (float32_t )(j + 1) / xRatio - 1;
      y = (float32_t )(i + 1) / yRatio - 1;

      // bilinear interpolation
      bilinear_interpolation(pSrcImage, SrcWidth, SrcHeight, nChannels, x, y,
                          pDstImage + (i * DstWidth + j) * nChannels);
    }
}

//------------------------------------------------------------------------------------------------------------
//  horizontal direction filtering
//------------------------------------------------------------------------------------------------------------
void horizontal_filtering(const uint8_t pSrcImage[], uint8_t pDstImage[], int32_t width,
                                 int32_t height, int32_t nChannels,
                                 const float32_t * pfilter1D, int32_t fsize) {
  memset(pDstImage, 0, sizeof(uint8_t) * width * height * nChannels);
  uint8_t * pBuffer;
  float32_t  w;
  int32_t i, j, l, k, offset, jj;
  for (i = 0; i < height; i++)
    for (j = 0; j < width; j++) {
      offset = i * width * nChannels;
      pBuffer = pDstImage + offset + j * nChannels;
      for (l = -fsize; l <= fsize; l++) {
        w = pfilter1D[l + fsize];
        jj = EnforceRange(j + l, width);
        for (k = 0; k < nChannels; k++)
          pBuffer[k] += pSrcImage[offset + jj * nChannels + k] * w;
      }
    }
}

//------------------------------------------------------------------------------------------------------------
//  horizontal direction filtering transpose
//------------------------------------------------------------------------------------------------------------
void horizontal_filtering_transpose(const uint8_t pSrcImage[], uint8_t pDstImage[],
                                           int32_t width, int32_t height, int32_t nChannels,
                                           const float32_t * pfilter1D, int32_t fsize) {
  memset(pDstImage, 0, sizeof(uint8_t) * width * height * nChannels);
  const uint8_t* pBuffer;
  float32_t  w;
  int32_t i, j, l, k, offset, jj;
  for (i = 0; i < height; i++)
    for (j = 0; j < width; j++) {
      int32_t offset0 = i * width * nChannels;
      pBuffer = pSrcImage + (i * width + j) * nChannels;
      for (l = -fsize; l <= fsize; l++) {
        w = pfilter1D[l + fsize];
        jj = EnforceRange(j + l, width);
        offset = offset0 + jj * nChannels;
        for (k = 0; k < nChannels; k++) pDstImage[offset + k] += pBuffer[k] * w;
      }
    }
}
//------------------------------------------------------------------------------------------------------------
// fast filtering algorithm for laplacian
//------------------------------------------------------------------------------------------------------------

void laplacian_filtering(const uint8_t pSrcImage[], uint8_t pDstImage[], int32_t width,
                                int32_t height, int32_t nChannels) {
  int32_t LineWidth = width * nChannels;
  int32_t nElements = width * height * nChannels;
  // first treat the corners
  for (int32_t k = 0; k < nChannels; k++) {
    pDstImage[k] =
        pSrcImage[k] * 2 - pSrcImage[nChannels + k] - pSrcImage[LineWidth + k];
    pDstImage[LineWidth - nChannels + k] =
        pSrcImage[LineWidth - nChannels + k] * 2 -
        pSrcImage[LineWidth - 2 * nChannels + k] -
        pSrcImage[2 * LineWidth - nChannels + k];
    pDstImage[nElements - LineWidth + k] =
        pSrcImage[nElements - LineWidth + k] * 2 -
        pSrcImage[nElements - LineWidth + nChannels + k] -
        pSrcImage[nElements - 2 * LineWidth + k];
    pDstImage[nElements - nChannels + k] =
        pSrcImage[nElements - nChannels + k] * 2 -
        pSrcImage[nElements - 2 * nChannels + k] -
        pSrcImage[nElements - LineWidth - nChannels + k];
  }
  // then treat the borders
  for (int32_t i = 1; i < width - 1; i++)
    for (int32_t k = 0; k < nChannels; k++) {
      pDstImage[i * nChannels + k] = pSrcImage[i * nChannels + k] * 3 -
                                     pSrcImage[(i - 1) * nChannels + k] -
                                     pSrcImage[(i + 1) * nChannels + k] -
                                     pSrcImage[i * nChannels + LineWidth + k];
      pDstImage[nElements - LineWidth + i * nChannels + k] =
          pSrcImage[nElements - LineWidth + i * nChannels + k] * 3 -
          pSrcImage[nElements - LineWidth + (i - 1) * nChannels + k] -
          pSrcImage[nElements - LineWidth + (i + 1) * nChannels + k] -
          pSrcImage[nElements - 2 * LineWidth + i * nChannels + k];
    }
  for (int32_t i = 1; i < height - 1; i++)
    for (int32_t k = 0; k < nChannels; k++) {
      pDstImage[i * LineWidth + k] = pSrcImage[i * LineWidth + k] * 3 -
                                     pSrcImage[i * LineWidth + nChannels + k] -
                                     pSrcImage[(i - 1) * LineWidth + k] -
                                     pSrcImage[(i + 1) * LineWidth + k];
      pDstImage[(i + 1) * LineWidth - nChannels + k] =
          pSrcImage[(i + 1) * LineWidth - nChannels + k] * 3 -
          pSrcImage[(i + 1) * LineWidth - 2 * nChannels + k] -
          pSrcImage[i * LineWidth - nChannels + k] -
          pSrcImage[(i + 2) * LineWidth - nChannels + k];
    }
  // now the interior
  for (int32_t i = 1; i < height - 1; i++)
    for (int32_t j = 1; j < width - 1; j++) {
      int32_t offset = (i * width + j) * nChannels;
      for (int32_t k = 0; k < nChannels; k++)
        pDstImage[offset + k] = pSrcImage[offset + k] * 4 -
                                pSrcImage[offset + nChannels + k] -
                                pSrcImage[offset - nChannels + k] -
                                pSrcImage[offset - LineWidth + k] -
                                pSrcImage[offset + LineWidth + k];
    }
}

//------------------------------------------------------------------------------------------------------------
// vertical direction filtering
//------------------------------------------------------------------------------------------------------------
void vertical_filtering(const uint8_t pSrcImage[], uint8_t pDstImage[], int32_t width,
                                 int32_t height, int32_t nChannels,
                                 const float32_t * pfilter1D, int32_t fsize) {
  memset(pDstImage, 0, sizeof(uint8_t) * width * height * nChannels);
  uint8_t* pBuffer;
  float32_t  w;
  int32_t i, j, l, k, offset, ii;
  for (i = 0; i < height; i++)
    for (j = 0; j < width; j++) {
      pBuffer = pDstImage + (i * width + j) * nChannels;
      for (l = -fsize; l <= fsize; l++) {
        w = pfilter1D[l + fsize];
        ii = EnforceRange(i + l, height);
        for (k = 0; k < nChannels; k++)
          pBuffer[k] += pSrcImage[(ii * width + j) * nChannels + k] * w;
      }
    }
}

//------------------------------------------------------------------------------------------------------------
// vertical direction filtering transpose
//------------------------------------------------------------------------------------------------------------
void vertical_filtering_transpose(const uint8_t pSrcImage[], uint8_t pDstImage[],
                                           int32_t width, int32_t height, int32_t nChannels,
                                           const float32_t * pfilter1D, int32_t fsize) {
  memset(pDstImage, 0, sizeof(uint8_t) * width * height * nChannels);
  const uint8_t* pBuffer;
  float32_t  w;
  int32_t i, j, l, k, offset, ii;
  for (i = 0; i < height; i++)
    for (j = 0; j < width; j++) {
      pBuffer = pSrcImage + (i * width + j) * nChannels;
      for (l = -fsize; l <= fsize; l++) {
        w = pfilter1D[l + fsize];
        ii = EnforceRange(i + l, height);
        offset = (ii * width + j) * nChannels;
        for (k = 0; k < nChannels; k++)
          // pBuffer[k]+=pSrcImage[(ii*width+j)*nChannels+k]*w;
          pDstImage[offset + k] += pBuffer[k] * w;
      }
    }
}

//------------------------------------------------------------------------------------------------------------
// 2d filtering
//------------------------------------------------------------------------------------------------------------
void filtering_2d(const uint8_t pSrcImage[], uint8_t pDstImage[], int32_t width,
                                int32_t height, int32_t nChannels,
                                const float32_t * pfilter2D, int32_t fsize) {
  float32_t  w;
  int32_t i, j, u, v, k, ii, jj, wsize, offset;
  wsize = fsize * 2 + 1;
  float32_t * pBuffer = malloc(sizeof(float32_t) * nChannels);
  for (i = 0; i < height; i++)
    for (j = 0; j < width; j++) {
      for (k = 0; k < nChannels; k++) pBuffer[k] = 0;
      for (u = -fsize; u <= fsize; u++)
        for (v = -fsize; v <= fsize; v++) {
          w = pfilter2D[(u + fsize) * wsize + v + fsize];
          ii = EnforceRange(i + u, height);
          jj = EnforceRange(j + v, width);
          offset = (ii * width + jj) * nChannels;
          for (k = 0; k < nChannels; k++)
            pBuffer[k] += pSrcImage[offset + k] * w;
        }
      offset = (i * width + j) * nChannels;
      for (k = 0; k < nChannels; k++) pDstImage[offset + k] = pBuffer[k];
    }
   free(pBuffer);
}

//------------------------------------------------------------------------------------------------------------
// 2d filtering transpose
//------------------------------------------------------------------------------------------------------------
void filtering_2d_transpose(const uint8_t pSrcImage[], uint8_t pDstImage[],
                                          int32_t width, int32_t height, int32_t nChannels,
                                          const float32_t * pfilter2D, int32_t fsize) {
  float32_t  w;
  int32_t i, j, u, v, k, ii, jj, wsize, offset;
  wsize = fsize * 2 + 1;
  memset(pDstImage, 0, sizeof(uint8_t) * width * height * nChannels);
  for (i = 0; i < height; i++)
    for (j = 0; j < width; j++) {
      int32_t offset0 = (i * width + j) * nChannels;
      for (u = -fsize; u <= fsize; u++)
        for (v = -fsize; v <= fsize; v++) {
          w = pfilter2D[(u + fsize) * wsize + v + fsize];
          ii = EnforceRange(i + u, height);
          jj = EnforceRange(j + v, width);
          int32_t offset = (ii * width + jj) * nChannels;
          for (k = 0; k < nChannels; k++)
            pDstImage[offset + k] += pSrcImage[offset0 + k] * w;
        }
    }
}

//------------------------------------------------------------------------------------------------------------
// function to sample a patch from the source image
//------------------------------------------------------------------------------------------------------------
void get_patch(const uint8_t pSrcImage[], uint8_t pPatch[], int32_t width,
                               int32_t height, int32_t nChannels, float32_t  x0, float32_t  y0,
                               int32_t wsize) {
  // suppose pPatch has been allocated and cleared before calling the function
  int32_t wlength = wsize * 2 + 1;
  float32_t  x, y;
  for (int32_t i = -wsize; i <= wsize; i++)
    for (int32_t j = -wsize; j <= wsize; j++) {
      y = y0 + i;
      x = x0 + j;
      if (x < 0 || x > width - 1 || y < 0 || y > height - 1) continue;
      bilinear_interpolation(
          pSrcImage, width, height, nChannels, x, y,
          pPatch + ((i + wsize) * wlength + j + wsize) * nChannels);
    }
}

//------------------------------------------------------------------------------------------------------------
// function to warp an image with respect to flow field
// pWarpIm2 has to be allocated before hands
//------------------------------------------------------------------------------------------------------------
void warp_image(uint8_t pWarpIm2[], const uint8_t pIm1[], const uint8_t pIm2[],
                                const uint8_t pVx[], const uint8_t pVy[], int32_t width,
                                int32_t height, int32_t nChannels) {
  memset(pWarpIm2, 0, sizeof(uint8_t) * width * height * nChannels);
  for (int32_t i = 0; i < height; i++)
    for (int32_t j = 0; j < width; j++) {
      int32_t offset = i * width + j;
      float32_t  x, y;
      y = i + pVy[offset];
      x = j + pVx[offset];
      offset *= nChannels;
      if (x < 0 || x > width - 1 || y < 0 || y > height - 1) {
        for (int32_t k = 0; k < nChannels; k++)
          pWarpIm2[offset + k] = pIm1[offset + k];
        continue;
      }
      BilinearInterpolate(pIm2, width, height, nChannels, x, y,
                          pWarpIm2 + offset);
    }
}

void warp_image_flow(uint8_t pWarpIm2[], const uint8_t pIm1[],
                                    const uint8_t pIm2, const float32_t pFlow[], int32_t width,
                                    int32_t height, int32_t nChannels) {
  memset(pWarpIm2, 0, sizeof(uint8_t) * width * height * nChannels);
  for (int32_t i = 0; i < height; i++)
    for (int32_t j = 0; j < width; j++) {
      int32_t offset = i * width + j;
      float32_t  x, y;
      y = i + pFlow[offset * 2 + 1];
      x = j + pFlow[offset * 2];
      offset *= nChannels;
      if (x < 0 || x > width - 1 || y < 0 || y > height - 1) {
        for (int32_t k = 0; k < nChannels; k++)
          pWarpIm2[offset + k] = pIm1[offset + k];
        continue;
      }
      BilinearInterpolate(pIm2, width, height, nChannels, x, y,
                          pWarpIm2 + offset);
    }
}

void warpImage(uint8_t pWarpIm2, const uint8_t pIm2, const T2* pVx,
                                const T2* pVy, int32_t width, int32_t height,
                                int32_t nChannels) {
  memset(pWarpIm2, 0, sizeof(uint8_t) * width * height * nChannels);
  for (int32_t i = 0; i < height; i++)
    for (int32_t j = 0; j < width; j++) {
      int32_t offset = i * width + j;
      float32_t  x, y;
      y = i + pVy[offset];
      x = j + pVx[offset];
      offset *= nChannels;
      if (x < 0 || x > width - 1 || y < 0 || y > height - 1) continue;
      BilinearInterpolate(pIm2, width, height, nChannels, x, y,
                          pWarpIm2 + offset);
    }
}

void warpImage_transpose(uint8_t pWarpIm2, const uint8_t pIm2,
                                          const T2* pVx, const T2* pVy,
                                          int32_t width, int32_t height,
                                          int32_t nChannels) {
  memset(pWarpIm2, 0, sizeof(uint8_t) * width * height * nChannels);
  for (int32_t i = 0; i < height; i++)
    for (int32_t j = 0; j < width; j++) {
      int32_t offset = i * width + j;
      float32_t  x, y;
      y = i + pVy[offset];
      x = j + pVx[offset];
      offset *= nChannels;
      if (x < 0 || x > width - 1 || y < 0 || y > height - 1) continue;
      // BilinearInterpolate(pIm2,width,height,nChannels,x,y,pWarpIm2+offset);
      BilinearInterpolate_transpose(pIm2 + offset, width, height, nChannels, x,
                                    y, pWarpIm2);
    }
}

//////////////////////////////////////////////////////////////////////////////////////
// different format
//////////////////////////////////////////////////////////////////////////////////////
void warpImage(uint8_t pWarpIm2, const uint8_t pIm2, const T2* flow,
                                int32_t width, int32_t height, int32_t nChannels) {
  memset(pWarpIm2, 0, sizeof(uint8_t) * width * height * nChannels);
  for (int32_t i = 0; i < height; i++)
    for (int32_t j = 0; j < width; j++) {
      int32_t offset = i * width + j;
      float32_t  x, y;
      y = i + flow[offset * 2 + 1];
      x = j + flow[offset * 2];
      offset *= nChannels;
      if (x < 0 || x > width - 1 || y < 0 || y > height - 1) continue;
      BilinearInterpolate(pIm2, width, height, nChannels, x, y,
                          pWarpIm2 + offset);
    }
}

void warpImage_transpose(uint8_t pWarpIm2, const uint8_t pIm2,
                                          const T2* flow, int32_t width, int32_t height,
                                          int32_t nChannels) {
  memset(pWarpIm2, 0, sizeof(uint8_t) * width * height * nChannels);
  for (int32_t i = 0; i < height; i++)
    for (int32_t j = 0; j < width; j++) {
      int32_t offset = i * width + j;
      float32_t  x, y;
      y = i + flow[offset * 2 + 1];
      x = j + flow[offset * 2];
      offset *= nChannels;
      if (x < 0 || x > width - 1 || y < 0 || y > height - 1) continue;
      // BilinearInterpolate(pIm2,width,height,nChannels,x,y,pWarpIm2+offset);
      BilinearInterpolate_transpose(pIm2 + offset, width, height, nChannels, x,
                                    y, pWarpIm2);
    }
}

void warpImage(uint8_t pWarpIm2, T3* pMask, const uint8_t pIm1,
                                const uint8_t pIm2, const T2* pVx, const T2* pVy,
                                int32_t width, int32_t height, int32_t nChannels) {
  memset(pWarpIm2, 0, sizeof(uint8_t) * width * height * nChannels);
  for (int32_t i = 0; i < height; i++)
    for (int32_t j = 0; j < width; j++) {
      int32_t offset = i * width + j;
      float32_t  x, y;
      y = i + pVy[offset];
      x = j + pVx[offset];
      offset *= nChannels;
      if (x < 0 || x > width - 1 || y < 0 || y > height - 1) {
        for (int32_t k = 0; k < nChannels; k++)
          pWarpIm2[offset + k] = pIm1[offset + k];
        pMask[i * width + j] = 0;
        continue;
      }
      pMask[i * width + j] = 1;
      BilinearInterpolate(pIm2, width, height, nChannels, x, y,
                          pWarpIm2 + offset);
    }
}

//------------------------------------------------------------------------------------------------------------
// function to crop an image from the source
// assume that pDstImage has been allocated
// also Left and Top must be valid, DstWidth and DstHeight should ensure that
// the image lies inside the image boundary
//------------------------------------------------------------------------------------------------------------
void crop_image(const uint8_t pSrcImage[], int32_t SrcWidth,
                                int32_t SrcHeight, int32_t nChannels, uint8_t pDstImage[],
                                int32_t Left, int32_t Top, int32_t DstWidth,
                                int32_t DstHeight) {
  if (typeid(uint8_t) == typeid(T2)) {
    for (int32_t i = 0; i < DstHeight; i++)
      memcpy(pDstImage + i * DstWidth * nChannels,
             pSrcImage + ((i + Top) * SrcWidth + Left) * nChannels,
             sizeof(uint8_t) * DstWidth * nChannels);
    return;
  }
  int32_t offsetSrc, offsetDst;
  for (int32_t i = 0; i < DstHeight; i++)
    for (int32_t j = 0; j < DstWidth; j++) {
      offsetSrc = ((i + Top) * SrcWidth + Left + j) * nChannels;
      offsetDst = (i * DstWidth + j) * nChannels;
      for (int32_t k = 0; k < nChannels; k++)
        pDstImage[offsetDst + k] = pSrcImage[offsetSrc + k];
    }
}

//------------------------------------------------------------------------------------------------------------
// function to generate a 2D Gaussian image
// pImage must be allocated before calling the function
//------------------------------------------------------------------------------------------------------------
template <class T>
void ImageProcessing::generate2DGaussian(T*& pImage, int32_t wsize, float32_t  sigma) {
  if (sigma == -1) sigma = wsize / 2;
  float32_t  alpha = 1 / (2 * sigma * sigma);
  int32_t winlength = wsize * 2 + 1;
  if (pImage == NULL) pImage = new T[winlength * winlength];
  float32_t  total = 0;
  for (int32_t i = -wsize; i <= wsize; i++)
    for (int32_t j = -wsize; j <= wsize; j++) {
      pImage[(i + wsize) * winlength + j + wsize] =
          exp(-(float32_t )(i * i + j * j) * alpha);
      total += pImage[(i + wsize) * winlength + j + wsize];
    }
  for (int32_t i = 0; i < winlength * winlength; i++) pImage[i] /= total;
}

//------------------------------------------------------------------------------------------------------------
// function to generate a 1D Gaussian image
// pImage must be allocated before calling the function
//------------------------------------------------------------------------------------------------------------
template <class T>
void ImageProcessing::generate1DGaussian(T*& pImage, int32_t wsize, float32_t  sigma) {
  if (sigma == -1) sigma = wsize / 2;
  float32_t  alpha = 1 / (2 * sigma * sigma);
  int32_t winlength = wsize * 2 + 1;
  if (pImage == NULL) pImage = new T[winlength];
  float32_t  total = 0;
  for (int32_t i = -wsize; i <= wsize; i++) {
    pImage[i + wsize] = exp(-(float32_t )(i * i) * alpha);
    total += pImage[i + wsize];
  }
  for (int32_t i = 0; i < winlength; i++) pImage[i] /= total;
}