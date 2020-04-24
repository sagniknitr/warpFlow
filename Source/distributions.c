#include "math.h"
#include "memory.h"
#include "stdlib.h"
#include "warpFlow.h"

static inline uint8_t _abs(int8_t x) { return (x >= 0) ? x : -x; }

#ifndef PI
#define PI 3.1415927
#endif

enum SortType { SortAscending, SortDescending };

void get_mean_var(uint8_t signal[], int32_t length, float32_t* mean,
                  float32_t* var) {
  float32_t m_mean = 0, m_var = 0;

  int32_t i;
  for (i = 0; i < length; i++) m_mean += signal[i];
  m_mean /= length;
  for (i = 0; i < length; i++)
    m_var += (signal[i] - m_mean) * (signal[i] - m_mean);
  m_var /= length - 1;
  *mean = m_mean;
  *var = m_var;
}


uint32_t sum(int32_t num_data, uint8_t pData[]) {
  uint32_t sum = 0;
  int32_t i;
  for (i = 0; i < num_data; i++) sum += pData[i];
  return sum;
}

void normalize(int32_t num_data, uint8_t pData[]) {
  int32_t i;
  uint32_t Sum;
  Sum = sum(num_data, pData);
  for (i = 0; i < num_data; i++) pData[i] /= Sum;
}



////////////////////////////////////////////////////////////
// sort data in descending order
template <class T>
void CStochastic::sort(int32_t Number, T* pData, int32_t* pIndex, SortType m_SortType) {
  int32_t i, j, offset_extreme, *flag;
  float32_t extreme;
  flag = new int32_t[Number];
  memset(flag, 0, sizeof(int32_t) * Number);
  for (i = 0; i < Number; i++) {
    if (m_SortType == SortDescending)
      extreme = -1E100;
    else
      extreme = 1E100;
    offset_extreme = 0;
    for (j = 0; j < Number; j++) {
      if (flag[j] == 1) continue;
      if ((m_SortType == SortDescending && extreme < pData[j]) ||
          (m_SortType == SortAscending && extreme > pData[j])) {
        extreme = pData[j];
        offset_extreme = j;
      }
    }
    pIndex[i] = offset_extreme;
    flag[offset_extreme] = 1;
  }
  delete flag;
}

template <class T>
T CStochastic::Min(int32_t NumData, T* pData) {
  int32_t i;
  T result = pData[0];
  for (i = 1; i < NumData; i++) result = __min(result, pData[i]);
  return result;
}

template <class T>
T CStochastic::Min(int32_t NumData, T* pData1, T* pData2) {
  int32_t i;
  T result = pData1[0] + pData2[0];
  for (i = 1; i < NumData; i++) result = __min(result, pData1[i] + pData2[i]);
  return result;
}

template <class T>
T CStochastic::Max(int32_t NumData, T* pData) {
  int32_t i;
  T result = pData[0];
  for (i = 1; i < NumData; i++) result = __max(result, pData[i]);
  return result;
}

template <class T>
int32_t CStochastic::FindMax(int32_t NumData, T* pData) {
  int32_t i, index;
  T result = pData[0];
  index = 0;
  for (i = 1; i < NumData; i++)
    if (pData[i] > result) {
      index = i;
      result = pData[i];
    }
  return index;
}

template <class T1, class T2>
void CStochastic::ComputeMeanCovariance(int32_t Dim, int32_t NumData, T1* pData,
                                        T2* pMean, T2* pCovariance,
                                        float32_t* pWeight) {
  int32_t i, j, k;
  memset(pMean, 0, sizeof(T2) * Dim);
  memset(pCovariance, 0, sizeof(T2) * Dim * Dim);

  bool IsWeightLoaded = false;
  float32_t Sum;
  if (pWeight != NULL) IsWeightLoaded = true;

  // compute mean first
  Sum = 0;
  if (IsWeightLoaded)
    for (i = 0; i < NumData; i++) {
      if (pWeight[i] == 0) continue;
      for (j = 0; j < Dim; j++) pMean[j] += pData[i * Dim + j] * pWeight[i];
      Sum += pWeight[i];
    }
  else {
    for (i = 0; i < NumData; i++)
      for (j = 0; j < Dim; j++) pMean[j] += pData[i * Dim + j];
    Sum = NumData;
  }
  for (j = 0; j < Dim; j++) pMean[j] /= Sum;

  // compute covariance;
  T2* pTempVector;
  pTempVector = new T2[Dim];

  for (i = 0; i < NumData; i++) {
    for (j = 0; j < Dim; j++) pTempVector[j] = pData[i * Dim + j] - pMean[j];
    if (IsWeightLoaded) {
      if (pWeight[i] == 0) continue;
      for (j = 0; j < Dim; j++)
        for (k = 0; k <= j; k++)
          pCovariance[j * Dim + k] +=
              pTempVector[j] * pTempVector[k] * pWeight[i];
    } else
      for (j = 0; j < Dim; j++)
        for (k = 0; k <= j; k++)
          pCovariance[j * Dim + k] += pTempVector[j] * pTempVector[k];
  }
  for (j = 0; j < Dim; j++)
    for (k = j + 1; k < Dim; k++)
      pCovariance[j * Dim + k] = pCovariance[k * Dim + j];

  for (j = 0; j < Dim * Dim; j++) pCovariance[j] /= Sum;

  delete[] pTempVector;
}

template <class T1, class T2>
void CStochastic::ComputeVectorMean(int32_t Dim, int32_t NumData, T1* pData, T2* pMean,
                                    float32_t* pWeight) {
  int32_t i, j;
  memset(pMean, 0, sizeof(T2) * Dim);
  bool IsWeightLoaded;
  float32_t Sum;
  if (pWeight = NULL)
    IsWeightLoaded = false;
  else
    IsWeightLoaded = true;

  Sum = 0;
  if (IsWeightLoaded)
    for (i = 0; i < NumData; i++) {
      if (pWeight[i] == 0) continue;
      for (j = 0; j < Dim; j++) pMean[j] += pData[i * Dim + j] * pWeight[i];
      Sum += pWeight[i];
    }
  else {
    for (i = 0; i < NumData; i++)
      for (j = 0; j < Dim; j++) pMean[j] += pData[i * Dim + j];
    Sum = NumData;
  }
  for (j = 0; j < Dim; j++) pMean[j] /= Sum;
}

template <class T1, class T2>
float32_t CStochastic::VectorSquareDistance(int32_t Dim, T1* pVector1, T2* pVector2) {
  float32_t result = 0, temp;
  int32_t i;
  for (i = 0; i < Dim; i++) {
    temp = pVector1[i] - pVector2[i];
    result += temp * temp;
  }
  return result;
}

template <class T1>
void CStochastic::KMeanClustering(int32_t Dim, int32_t NumData, int32_t NumClusters,
                                  T1* pData, int32_t* pPartition,
                                  float32_t** pClusterMean, int32_t MaxIterationNum,
                                  int32_t MinClusterSampleNumber) {
  int32_t i, j, k, l, Index, ClusterSampleNumber;
  float32_t MinDistance, Distance;
  float32_t** pCenters;
  pCenters = new float32_t*[NumClusters];
  for (i = 0; i < NumClusters; i++) pCenters[i] = new float32_t[Dim];

  // generate randome guess of the partition
_CStochastic_KMeanClustering_InitializePartition:
  for (i = 0; i < NumClusters; i++) {
    Index = UniformSampling(NumData);
    for (j = 0; j < Dim; j++) pCenters[i][j] = pData[Index * Dim + j];
  }

  for (k = 0; k < MaxIterationNum; k++) {
    // step 1. do partition
    for (i = 0; i < NumData; i++) {
      MinDistance = 1E100;
      for (j = 0; j < NumClusters; j++) {
        Distance = VectorSquareDistance(Dim, pData + i * Dim, pCenters[j]);
        if (Distance < MinDistance) {
          MinDistance = Distance;
          Index = j;
        }
      }
      pPartition[i] = Index;
    }
    // step 2. compute mean
    for (i = 0; i < NumClusters; i++) {
      memset(pCenters[i], 0, sizeof(float32_t) * Dim);
      ClusterSampleNumber = 0;
      for (j = 0; j < NumData; j++)
        if (pPartition[j] == i) {
          for (l = 0; l < Dim; l++) pCenters[i][l] += pData[j * Dim + l];
          ClusterSampleNumber++;
        }
      // maybe the initial partition is bad
      // if so just do initial partition again
      if (ClusterSampleNumber < MinClusterSampleNumber)
        goto _CStochastic_KMeanClustering_InitializePartition;
      for (l = 0; l < Dim; l++) pCenters[i][l] /= ClusterSampleNumber;
    }
  }
  // output the final partition if necessary
  if (pClusterMean != NULL)
    for (i = 0; i < NumClusters; i++)
      for (l = 0; l < Dim; l++) pClusterMean[i][l] = pCenters[i][l];
  // free buffer
  for (i = 0; i < NumClusters; i++) delete pCenters[i];
  delete[] pCenters;
}


float32_t norm(uint8_t X[], int32_t Dim) {
  float32_t result = 0;
  int32_t i;
  for (i = 0; i < Dim; i++) result += X[i] * X[i];
  result = sqrt(result);
  return result;
}

template <class T1, class T2>
int32_t CStochastic::FindClosestPoint(T1* pPointSet, int32_t NumPoints, int32_t nDim,
                                  T2* QueryPoint) {
  int32_t i, j, Index = 0, offset;
  T1 MinDistance, Distance, x;
  MinDistance = 0;
  for (j = 0; j < nDim; j++) MinDistance += _abs(pPointSet[j] - QueryPoint[j]);
  for (i = 1; i < NumPoints; i++) {
    Distance = 0;
    offset = i * nDim;
    for (j = 0; j < nDim; j++) {
      x = pPointSet[offset + j] - QueryPoint[j];
      Distance += _abs(x);
    }
    if (Distance < MinDistance) {
      MinDistance = Distance;
      Index = i;
    }
  }
  return Index;
}

template <class T1, class T2>
void CStochastic::GaussianFiltering(T1* pSrcArray, T2* pDstArray, int32_t NumPoints,
                                    int32_t nChannels, int32_t size, float32_t sigma) {
  int32_t i, j, u, l;
  float32_t *pGaussian, temp;
  pGaussian = new float32_t[2 * size + 1];
  Generate1DGaussian(pGaussian, size, sigma);
  for (i = 0; i < NumPoints; i++)
    for (l = 0; l < nChannels; l++) {
      temp = 0;
      for (j = -size; j <= size; j++) {
        u = i + j;
        u = __max(__min(u, NumPoints - 1), 0);
        temp += pSrcArray[u * nChannels + l] * pGaussian[j + size];
      }
      pDstArray[i * nChannels + l] = temp;
    }
  delete pGaussian;
}
