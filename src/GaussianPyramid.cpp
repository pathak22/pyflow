// Author: Ce Liu (c) Dec, 2009; celiu@mit.edu
// Modified By: Deepak Pathak (c) 2016; pathak@berkeley.edu

#include "GaussianPyramid.h"
#include "math.h"

GaussianPyramid::GaussianPyramid(void)
{
	ImPyramid=NULL;
}

GaussianPyramid::~GaussianPyramid(void)
{
	if(ImPyramid!=NULL)
		delete []ImPyramid;
}

//---------------------------------------------------------------------------------------
// function to construct the pyramid
// this is the slow way
//---------------------------------------------------------------------------------------
/*void GaussianPyramid::ConstructPyramid(const DImage &image, double ratio, int minWidth)
{
	// the ratio cannot be arbitrary numbers
	if(ratio>0.98 || ratio<0.4)
		ratio=0.75;
	// first decide how many levels
	nLevels=log((double)minWidth/image.width())/log(ratio);
	if(ImPyramid!=NULL)
		delete []ImPyramid;
	ImPyramid=new DImage[nLevels];
	ImPyramid[0].copyData(image);
	double baseSigma=(1/ratio-1);
	for(int i=1;i<nLevels;i++)
	{
		DImage foo;
		double sigma=baseSigma*i;
		image.GaussianSmoothing(foo,sigma,sigma*2.5);
		foo.imresize(ImPyramid[i],pow(ratio,i));
	}
}//*/

//---------------------------------------------------------------------------------------
// function to construct the pyramid
// this is the fast way
//---------------------------------------------------------------------------------------
void GaussianPyramid::ConstructPyramid(const DImage &image, double ratio, int minWidth)
{
	// the ratio cannot be arbitrary numbers
	if(ratio>0.98 || ratio<0.4)
		ratio=0.75;
	// first decide how many levels
	nLevels=log((double)minWidth/image.width())/log(ratio);
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

void GaussianPyramid::ConstructPyramidLevels(const DImage &image, double ratio, int _nLevels)
{
	// the ratio cannot be arbitrary numbers
	if(ratio>0.98 || ratio<0.4)
		ratio=0.75;
	nLevels = _nLevels;
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

void GaussianPyramid::displayTop(const char *filename)
{
	// ImPyramid[nLevels-1].imwrite(filename);
}
