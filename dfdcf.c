#include "dfdcf.h"
#define dw 50
#define ws 50
#define we 10000

double** DFT(double** arr, int n)
{
	double** dft;	//功率谱
	dft = (double**)malloc((we/dw) * sizeof(double));
	for (int i = 0; i < we/dw; i++)
	{
		dft[i] = (double*)malloc(2 * sizeof(double));
	}
	double r_cos, r_sin;
	int m = 0;
	for (double w = ws; w <= we; w = w + dw)
	{
		r_cos = 0;
		r_sin = 0;
		for (int i = 0; i < n - 1; i++)
		{
			r_cos = r_cos + pow(arr[i][1] * cos(w * arr[i][0]) * (arr[i + 1][0] - arr[i][0]), 2);
			r_sin = r_sin + pow(arr[i][1] * sin(w * arr[i][0]) * (arr[i + 1][0] - arr[i][0]), 2);
		}
		dft[m][0] = w;	//圆频率
		dft[m][1] = r_cos + r_sin;	//功率
		m++;
	}
	return dft;
}

double DFDCF(double** arr1, double** arr2, int na1, int na2, double dt, double tsyn)
{
	double** ar1, ** ar2;
	ar1 = (double**)malloc(na1 * sizeof(double));
	ar2 = (double**)malloc(na2 * sizeof(double));
	for (int i = 0; i < na1; i++)
	{
		ar1[i] = (double*)malloc(2 * sizeof(double));
	}
	for (int i = 0; i < na2; i++)
	{
		ar2[i] = (double*)malloc(2 * sizeof(double));
	}
	double ts, te;
	ts = max(arr1[0][0], arr2[0][0] - dt);
	te = min(arr1[na1 - 1][0], arr2[na2 - 1][0] - dt);
	int m1 = 0;
	int m2 = 0;
	for (int i = 0; i < na1; i++)
	{
		if (arr1[i][0] >= ts && arr1[i][0] <= te)
		{
			ar1[m1][0] = arr1[i][0];
			ar1[m1][1] = arr1[i][1];
			m1++;
		}
	}
	for (int i = 0; i < na2; i++)
	{
		if (arr2[i][0] - dt >= ts && arr2[i][0] - dt <= te)
		{
			ar2[m2][0] = arr2[i][0] - dt;
			ar2[m2][1] = arr2[i][1];
			m2++;
		}
	}
	double** dft1, ** dft2;
	dft1 = DFT(ar1, m1);
	dft2 = DFT(ar2, m2);
	double dcf_r;
	dcf_r = DCF(dft1, dft2, we/dw, we/dw, 0, dw/2);
	for (int i = 0; i < na1; i++)
	{
		free(ar1[i]);
	}
	free(ar1);
	for (int i = 0; i < na2; i++)
	{
		free(ar2[i]);
	}
	free(ar2);
	for (int i = 0; i < we / dw; i++)
	{
		free(dft1[i]);
		free(dft2[i]);
	}
	free(dft1);
	free(dft2);

	return dcf_r;
}

double timedelay_DFDCF(double** arr1, double** arr2, int n1, int n2, double tsyn)
{
	double temp = -1e26;
	double dcf_r;
	double T = 0;
	for (double t = arr2[0][0] - arr1[n1 - 1][0] - 0.04; t <= arr2[n2 - 1][0] - arr1[0][0] + 0.04; t = t + 0.01 * tsyn)
	{
		dcf_r = DFDCF(arr1, arr2, n1, n2, t, 0.01 * tsyn);
		if (temp < dcf_r)
		{
			temp = dcf_r;
			T = t;
		}
	}
	if (temp > 0.98 && (min(arr1[n1 - 1][0], arr2[n2 - 1][0] - T) - max(arr1[0][0], arr2[0][0] - T)) * 24 > 2.074 && T <= 90 && (Abs(mod(T, tsyn) / tsyn - 1) < 0.2 || Abs(mod(T, tsyn) / tsyn - 1) > 0.8))
	{
		return T;
	}
	else
	{
		return 0;
	}
}