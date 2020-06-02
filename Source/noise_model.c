#ifdef __cplusplus
extern "C" {
#endif

#include "common_datatypes.h"

#define PI 3.1415926535897932384626433832

typedef struct _gaussian_mixture {
  int32_t s32_num_of_channels;
  float32_t* f32_alpha;
  float32_t* f32_sigma;
  float32_t* f32_beta;
  float32_t* f32_sigma_square;
  float32_t* f32_beta_square;
} gaussian_mixture;

s8_status init_gaussian_mixture(int32_t s32_channels) {
  gaussian_mixture* gmm = nullptr;
  gmm->s32_num_of_channels = s32_channels;
  gaussian_mixture_allocate(gmm, s32_channels);
  for (int i = 0; i < s32_channels; i++) {
    gmm->f32_alpha[i] = 0.95;
    gmm->f32_sigma[i] = 0.05;
    gmm->f32_beta[i] = 0.5;
  }
  gaussian_mixture_square(gmm, s32_channels);
}

s8_status gaussian_mixture_square(gaussian_mixture* gmm, int32_t s32_channel) {
  for (uint32_t u32_c = 0; u32_c < (uint32_t)s32_channel; u32_c++) {
    gmm->f32_sigma_square[u32_c] =
        gmm->f32_sigma[u32_c] * gmm->f32_sigma[u32_c];
    gmm->f32_beta_square[u32_c] = gmm->f32_beta[u32_c] * gmm->f32_beta[u32_c];
  }
}

int32_t gaussian_mixture_allocate(gaussian_mixture* gmm, int32_t s32_channels) {
  gmm->f32_alpha = malloc(s32_channels * sizeof(float32_t));
  gmm->f32_sigma = malloc(s32_channels * sizeof(float32_t));
  gmm->f32_beta = malloc(s32_channels * sizeof(float32_t));
  gmm->f32_sigma_square = malloc(s32_channels * sizeof(float32_t));
  gmm->f32_beta_square = malloc(s32_channels * sizeof(float32_t));
}

uint8_t gaussian_mixture_free(gaussian_mixture* gmm) {
  if (!gmm->f32_alpha) free(gmm->f32_alpha);
  if (!gmm->f32_sigma) free(gmm->f32_sigma);
  if (!gmm->f32_beta) free(gmm->f32_beta);
  if (!gmm->f32_sigma_square) free(gmm->f32_sigma_square);
  if (!gmm->f32_beta_square) free(gmm->f32_beta_square);
  gmm->f32_alpha = gmm->f32_sigma = gmm->f32_beta = gmm->f32_sigma_square =
      gmm->f32_beta_square = nullptr;
}

uint8_t get_gaussian(gaussian_mixture* gmm, float32_t f32_x, int32_t s32_i,
                     int32_t s32_k) {
  if (s32_i == 0)
    return exp(-f32_x / (2 * gmm->f32_sigma_square[s32_k])) /
           (2 * PI * gmm->f32_sigma[s32_k]);
  else
    return exp(-f32_x / (2 * gmm->f32_beta_square[s32_k])) /
           (2 * PI * gmm->f32_beta[s32_k]);
}

gaussian_mixture* shrink_gaussian_mixture(gaussian_mixture* gmm,
                                          int32_t s32_channels) {
  gaussian_mixture* GM = nullptr;
  gaussian_mixture_allocate(GM, s32_channels);
  for (int32_t s32_i = 0; s32_i < s32_channels; s32_i++) {
    GM->f32_alpha[s32_i] = gmm->f32_alpha[s32_i];
    GM->f32_sigma[s32_i] = gmm->f32_sigma[s32_i];
    GM->f32_beta[s32_i] = gmm->f32_beta[s32_i];
  }
  gaussian_mixture_square(GM, s32_channels);
  return GM;
}

float32_t get_exponent(float32_t f32_input) {
  float32_t f32_y = 0.0f;
  float32_t f32_b = 0.5f;
  const int32_t precision = 24;

  const float32_t LOG2 = 0.693147180;

  while (f32_input < 1.0) {
    f32_input /= 2;
    f32_y -= 1;
  }

  for (int i = 0; i < precision; i++) {
    f32_input = f32_input * f32_input;

    if (f32_input >= 2.0) {
      f32_input /= 2.0;
      f32_y += f32_b;
    }
    f32_b /= 2.0;
  }
}
gaussian_mixture* shrink_gaussian_mixture(gaussian_mixture* gmm, int32_t s32_channels) {

	gaussian_mixture* GM = nullptr;
    gaussian_mixture_allocate(GM, s32_channels);
    for (int32_t s32_i = 0; s32_i < s32_channels; s32_i++) {
      GM->f32_alpha[s32_i] = gmm->f32_alpha[s32_i];
      GM->f32_sigma[s32_i] = gmm->f32_sigma[s32_i];
      GM->f32_beta[s32_i] = gmm->f32_beta[s32_i];
    }
    gaussian_mixture_square(GM, s32_channels);
    return GM;
  
}


#ifdef __cplusplus
}
#endif
