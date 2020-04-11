#include "warpFlow.h"

#define PI 3.1415926535897932384626433832

typedef struct _gaussian_mixture {
  int32_t  s32_num_of_channels;
  float32_t* f32_alpha;
  float32_t* f32_sigma;
  float32_t* f32_beta;
  float32_t* f32_sigma_square;
  float32_t* f32_beta_square;
} gaussian_mixture;


u8_status init_gaussian_mixture(int32_t s32_channels) {
  
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

u8_status  gaussian_mixture_square(gaussian_mixture *gmm, int32_t s32_channel) {
  for (uint32_t u32_c = 0; u32_c < (uint32_t)s32_channel; u32_c++) {
    gmm->f32_sigma_square[u32_c] = gmm->f32_sigma[u32_c] * gmm->f32_sigma[u32_c];
    gmm->f32_beta_square[u32_c] = gmm->f32_beta[u32_c] * gmm->f32_beta[u32_c];
  }
}


int32_t gaussian_mixture_allocate(gaussian_mixture *gmm, int32_t s32_channels) {
  gmm->f32_alpha = malloc(s32_channels * sizeof(float32_t));
  gmm->f32_sigma = malloc(s32_channels * sizeof(float32_t));
  gmm->f32_beta = malloc(s32_channels * sizeof(float32_t));
  gmm->f32_sigma_square =malloc(s32_channels * sizeof(float32_t));
  gmm->f32_beta_square = malloc(s32_channels * sizeof(float32_t));
}



class GaussianMixture {
 public:
  int nChannels;
  double* alpha;
  double* sigma;
  double* beta;
  double* sigma_square;
  double* beta_square;

 public:
  GaussianMixture() {
    nChannels = 0;
    alpha = sigma = beta = sigma_square = beta_square = NULL;
  }
  GaussianMixture(int _nChannels) {
    nChannels = _nChannels;
    allocate();
    for (int i = 0; i < nChannels; i++) {
      alpha[i] = 0.95;
      sigma[i] = 0.05;
      beta[i] = 0.5;
    }
    square();
  }
  GaussianMixture(const GaussianMixture& GM) {
    clear();
    copy(GM);
  }
  void copy(const GaussianMixture& GM) {
    nChannels = GM.nChannels;
    allocate();
    for (int i = 0; i < nChannels; i++) {
      alpha[i] = GM.alpha[i];
      sigma[i] = GM.sigma[i];
      beta[i] = GM.beta[i];
    }
    square();
  }
  void operator=(const GaussianMixture& GM) {
    clear();
    copy(GM);
  }
  GaussianMixture shrink(int N) {
    GaussianMixture GM(N);
    for (int i = 0; i < N; i++) {
      GM.alpha[i] = alpha[i];
      GM.sigma[i] = sigma[i];
      GM.beta[i] = beta[i];
    }
    GM.square();
    return GM;
  }
  void allocate() {
    alpha = new double[nChannels];
    sigma = new double[nChannels];
    beta = new double[nChannels];
    sigma_square = new double[nChannels];
    beta_square = new double[nChannels];
  }
  void clear() {
    if (!alpha) delete[] alpha;
    if (!sigma) delete[] sigma;
    if (!beta) delete[] beta;
    if (!sigma_square) delete[] sigma_square;
    if (!beta_square) delete[] beta_square;
    alpha = sigma = beta = sigma_square = beta_square = NULL;
  }
  void reset() {
    // for(int i = 0;i<nChannels;i++)
    //	alpha[i] = sigma[i] = beta[i] = sigma_square[i] = beta_square[i] = 0;
    for (int i = 0; i < nChannels; i++) {
      alpha[i] = 0.95;
      sigma[i] = 0.05;
      beta[i] = 0.5;
    }
    square();
  }
  void reset(int _nChannels) {
    clear();
    nChannels = _nChannels;
    allocate();
    reset();
  }
  double Gaussian(double x, int i, int k) const {
    if (i == 0)
      return exp(-x / (2 * sigma_square[k])) / (2 * PI * sigma[k]);
    else
      return exp(-x / (2 * beta_square[k])) / (2 * PI * beta[k]);
  }
  ~GaussianMixture() { clear(); }
  void square() {
    for (int i = 0; i < nChannels; i++) {
      sigma_square[i] = sigma[i] * sigma[i];
      beta_square[i] = beta[i] * beta[i];
    }
  }
  void display() {
    for (int i = 0; i < nChannels; i++)
      cout << "alpha: " << alpha[i] << " sigma: " << sigma[i]
           << " beta: " << beta[i] << " sigma^2: " << sigma_square[i]
           << " beta^2: " << beta_square[i] << endl;
  }
  bool write(const char* filename) {
    ofstream myfile(filename, ios::out | ios::binary);
    if (myfile.is_open()) {
      bool foo = write(myfile);
      myfile.close();
      return foo;
    }
    return false;
  }
  bool write(ofstream& myfile) {
    myfile.write((char*)&nChannels, sizeof(int));
    myfile.write((char*)alpha, sizeof(double) * nChannels);
    myfile.write((char*)sigma, sizeof(double) * nChannels);
    myfile.write((char*)beta, sizeof(double) * nChannels);
    return true;
  }
  bool read(const char* filename) {
    ifstream myfile(filename, ios::in | ios::binary);
    if (myfile.is_open()) {
      bool foo = read(myfile);
      myfile.close();
      square();
      return foo;
    }
    return false;
  }
  bool read(ifstream& myfile) {
    myfile.read((char*)&nChannels, sizeof(int));
    allocate();
    myfile.read((char*)alpha, sizeof(double) * nChannels);
    myfile.read((char*)sigma, sizeof(double) * nChannels);
    myfile.read((char*)beta, sizeof(double) * nChannels);
    square();
    return true;
  }
};

// class Laplacian
//{
// public:
//	int nChannels;
//	Vector<double> scale;
// public:
//	Laplacian()
//	{
//	}
//	Laplacian(int _nChannels)
//	{
//		nChannels = _nChannels;
//		scale.allocate(nChannels);
//	}
//	Laplacian(const Laplacian
//
//};