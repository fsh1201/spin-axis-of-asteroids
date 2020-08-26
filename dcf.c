#include "dcf.h"

double DCF(double** arr1, double** arr2, int na1, int na2, double dt, double dtau)
{
	double d = 0;	//dcf
	double* ar1, * ar2;
	ar1 = (double*)malloc(na1 * sizeof(double));
	ar2 = (double*)malloc(na2 * sizeof(double));
	for (int i = 0; i < na1; i++)
	{
		ar1[i] = arr1[i][1];
	}
	for (int i = 0; i < na2; i++)
	{
		ar2[i] = arr2[i][1];
	}
	double s1 = stdd(ar1, na1);
	double s2 = stdd(ar2, na2);
	double m1 = mean(ar1, na1);
	double m2 = mean(ar2, na2);
	int M = 0;
	for (int i = 0; i < na1; i++)
	{
		for (int j = 0; j < na2; j++)
		{
			if (arr2[j][0] - arr1[i][0] >= dt - dtau && arr2[j][0] - arr1[i][0] <= dt + dtau)
			{
				d = d + (arr1[i][1] - m1) * (arr2[j][1] - m2) / (s1 * s2);
				M++;
			}
		}
	}
	d = d / M;
	free(ar1);
	free(ar2);
	return d;
}

double timedelay_DCF(double** arr1, double** arr2, int n1, int n2, double tsyn)
{
	int ss = (int)(((arr2[n2 - 1][0] - arr1[0][0]) - (arr2[0][0] - arr1[n1 - 1][0])) / (0.01 * tsyn) + 1);
	double** s = (double**)malloc( ss * sizeof(double));
	for (int i = 0; i < ss; i++)
	{
		s[i] = (double*)malloc(2 * sizeof(double));
	}
	int sss = 0;
	for (double t = arr2[0][0] - arr1[n1 - 1][0]; t <= arr2[n2 - 1][0] - arr1[0][0]; t = t + 0.01 * tsyn)
	{
		s[sss][0] = t;
		s[sss][1] = DCF(arr1, arr2, n1, n2, t, 0.01 * tsyn);
		sss++;
	}
	double** sf = avef(s, ss);
	double T = 0;
	double temp = -1e26;
	for (int i = 0; i < ss; i++)
	{
		if (temp < sf[i][1])
		{
			temp = sf[i][1];
			T = sf[i][0];
		}
	}
	for (int i = 0; i < ss; i++)
	{
		free(s[i]);
		free(sf[i]);
	}
	free(s);
	free(sf);
	if (temp > 0.95 && (min(arr1[n1 - 1][0], arr2[n2 - 1][0] - T) - max(arr1[0][0], arr2[0][0] - T)) * 24 > 1 && (Abs(mod(T, tsyn) / tsyn - 1) < 0.2 || Abs(mod(T, tsyn) / tsyn - 1) > 0.8))
	{
		return T;
	}
	else
	{
		return 0;
	}
}