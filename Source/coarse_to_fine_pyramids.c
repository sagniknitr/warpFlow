#include "warpFlow.h"
#include "utilities.h"


s8_status free_pyramid(void) {

}

void construct_pyramid(const uint8_t image[],uint32_t u32_width,  float32_t ratio, uint32_t minWidth)
{
	// the ratio cannot be arbitrary numbers
	if(ratio>0.98 || ratio<0.4)
		ratio=0.75;
	// first decide how many levels
	uint32_t nLevels=log((float32_t)minWidth/u32_width)/log(ratio);
	if(ImPyramid!=NULL)
		delete []ImPyramid;
	ImPyramid=new DImage[nLevels];
	ImPyramid[0].copyData(image);
	double baseSigma=(1/ratio-1);
	int n=log(0.25)/log(ratio);
	double nSigma=baseSigma*n;
	for(int i=1;i<nLevels;i++)
	{
		DImage foo;
		if(i<=n)
		{
			double sigma=baseSigma*i;
			image.GaussianSmoothing(foo,sigma,sigma*3);
			foo.imresize(ImPyramid[i],pow(ratio,i));
		}
		else
		{
			ImPyramid[i-n].GaussianSmoothing(foo,nSigma,nSigma*3);
			double rate=(double)pow(ratio,i)*image.width()/foo.width();
			foo.imresize(ImPyramid[i],rate);
		}
	}
}
