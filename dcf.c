#include "dcf.h"
#include "juvf.h"
#include <stdlib.h>

double dcf(double** arr1, double** arr2, int na1, int na2, double tsyn)
{
	double temp = -1e50;
	double T = 0;	// ±º‰—”≥Ÿ
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
	for (double tau = arr2[0][0] - arr1[na1 - 1][0]; tau < arr2[na2 - 1][0] - arr1[0][0]; tau = tau + 0.01 * tsyn)
	{
		d = 0;
		M = 0;
		for (int i = 0; i < na1; i++)
		{
			for (int j = 0; j < na2; j++)
			{
				if (arr2[j][0] - arr1[i][0] >= tau - 0.01 * tsyn / 2 && arr2[j][0] - arr1[i][0] <= tau + 0.01 * tsyn / 2)
				{
					d = d + (arr1[i][1] - m1) * (arr2[j][1] - m2) / (s1 * s2);
					M++;
				}
			}
		}
		d = d / M;
		if (temp < d && M > 30)
		{
			temp = d;
			T = tau;
		}
	}
	free(ar1);
	free(ar2);
	if (temp > 0.98 && (min(arr1[na1 - 1][0], arr2[na2 - 1][0] - T) - max(arr1[0][0], arr2[0][0] - T)) * 24 > 2.074 && T <= 90 && (abs(mod(T, tsyn) / tsyn - 1) < 0.2 || abs(mod(T, tsyn) / tsyn - 1) > 0.8))
	{
		return T;
	}
	else
	{
		return 0;
	}
}