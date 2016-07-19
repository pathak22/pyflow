// Author: Ce Liu (c) Dec, 2009; celiu@mit.edu
// Modified By: Deepak Pathak (c) 2016; pathak@berkeley.edu

#include "Stochastic.h"
#include "time.h"
#include "stdlib.h"
#include "stdio.h"

CStochastic::CStochastic(void)
{
}

CStochastic::~CStochastic(void)
{
}

void CStochastic::ConvertInt2String(int x,char* string,int BitNumber)
{
	int i,Base=1;
	for(i=1;i<BitNumber;i++)
		Base*=10;
	for(i=0;i<BitNumber;i++)
	{
		string[i]=x/Base+'0';
		x%=Base;
		Base/=10;
	}
	string[i]='\0';
}

double CStochastic::UniformSampling()
{
	return (double)rand()/((double)RAND_MAX+(double)1);
}

int CStochastic::UniformSampling(int R)
{
	int Index=(double)UniformSampling()*R;
	if(Index>R-1)
		Index=R-1;
	return Index;
}

double CStochastic::GaussianSampling()
{
	int i;
	double result=0;
	for (i=0;i<12;i++)
		result+=UniformSampling();
	result-=6;
	return result;
}


double CStochastic::GetMean(double* signal,int length)
{
	double mean=0;
	int i;
	for(i=0;i<length;i++)
		mean+=signal[i];
	mean/=length;
	return mean;
}

int CStochastic::Sampling(double* Density,int NumSamples)
{
	double RandNumber=UniformSampling();
	int i;
	double sum=0;
	for(i=0;i<NumSamples;i++)
	{
		sum+=Density[i];
		if(sum>=RandNumber)
			return i;
	}
	return NumSamples-1;
}

void CStochastic::Generate1DGaussian(double* pGaussian,int size,double sigma)
{
	int i;
	if(sigma==0)
		sigma=size/2;
	for(i=-size;i<=size;i++)
		pGaussian[i+size]=exp(-(double)i*i/(2*sigma));
}

void CStochastic::Generate2DGaussian(double* pGaussian,int WinSize,double sigma)
{
	int i,j,WinLength=WinSize*2+1;
	double Sigma;
	if(sigma==0)
		Sigma=WinSize;
	else
		Sigma=sigma;
	Sigma*=Sigma;
	for (i=-WinSize;i<=WinSize;i++)
		for(j=-WinSize;j<=WinSize;j++)
			pGaussian[(i+WinSize)*WinLength+j+WinSize]=exp(-(double)(i*i+j*j)/(2*Sigma));
	Normalize(WinLength*WinLength,pGaussian);
}

double CStochastic::entropy(double* pDensity,int n)
{
	double result=0;
	int i;
	for(i=0;i<n;i++)
		result-=log(__max(pDensity[i],1E-6))*pDensity[i];
	return result;
}
