#include "juvf.h"
#include <math.h>
#include <stdlib.h>

/*行列式*/
double hhlx(double** arr, int na)
{
	if (na == 1)
	{
		return arr[0][0];
	}
	else
	{
		double s = 0;
		for (int i = 0; i < na; i++)
		{
			double** arr1;
			arr1 = (double**)malloc((na - 1) * sizeof(double));
			for (int i = 0; i < na - 1; i++)
			{
				arr1[i] = (double*)malloc((na - 1) * sizeof(double));
			}
			for (int j = 1; j < na; j++)
			{
				for (int k = 0; k < na - 1; k++)
				{
					if (k >= i)
					{
						arr1[j - 1][k] = arr[j][k + 1];
					}
					else
					{
						arr1[j - 1][k] = arr[j][k];
					}
				}
			}
			s = s + hhlx(arr1, na - 1) * pow(-1, i) * arr[0][i];
			for (int i = 0; i < na - 1; i++)
			{
				free(arr1[i]);
			}
			free(arr1);
		}
		return s;
	}
}

/*逆矩阵*/
double** inv(double** a, int n)
{
	double det = hhlx(a, n);
	double** as;
	as = (double**)malloc(n * sizeof(double));
	for (int i = 0; i < n; i++)
	{
		as[i] = (double*)malloc(n * sizeof(double));
	}
	for (int is = 0; is < n; is++)
	{
		for (int js = 0; js < n; js++)
		{
			double** ab;
			ab = (double**)malloc((n - 1) * sizeof(double));
			for (int i = 0; i < n - 1; i++)
			{
				ab[i] = (double*)malloc((n - 1) * sizeof(double));
			}
			for (int i = 0; i < n - 1; i++)
			{
				for (int j = 0; j < n - 1; j++)
				{
					if (i >= is)
					{
						if (j >= js)
							ab[i][j] = a[i + 1][j + 1];
						else
							ab[i][j] = a[i + 1][j];
					}
					else
					{
						if (j >= js)
							ab[i][j] = a[i][j + 1];
						else
							ab[i][j] = a[i][j];
					}
				}
			}
			as[js][is] = pow(-1, (double)is + (double)js) * hhlx(ab, n - 1);

			for (int i = 0; i < n - 1; i++)
			{
				free(ab[i]);
			}
			free(ab);
		}
	}
	double** ai;
	ai = (double**)malloc(n * sizeof(double));
	for (int i = 0; i < n; i++)
	{
		ai[i] = (double*)malloc(n * sizeof(double));
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			ai[i][j] = as[i][j] / det;
		}
	}
	return ai;
}

/*矩阵相乘*/
double** AB(double** a, int ma, int na, double** b, int mb, int nb)
{
	if (na != mb)
	{
		printf("计算错误！");
		return NULL;
	}
	else
	{
		double** ab;
		ab = (double**)malloc(ma * sizeof(double));
		for (int i = 0; i < ma; i++)
		{
			ab[i] = (double*)malloc(nb * sizeof(double));
		}
		for (int i = 0; i < ma; i++)
		{
			for (int j = 0; j < nb; j++)
			{
				ab[i][j] = 0;
				for (int k = 0; k < na; k++)
				{
					ab[i][j] = ab[i][j] + a[i][k] * b[k][j];
				}
			}
		}
		return ab;
	}
}

/*矩阵转置*/
double** TA(double** a, int ma, int na)
{
	double** ta;
	ta = (double**)malloc(na * sizeof(double));
	for (int i = 0; i < na; i++)
	{
		ta[i] = (double*)malloc(ma * sizeof(double));
	}
	for (int i = 0; i < na; i++)
	{
		for (int j = 0; j < ma; j++)
		{
			ta[i][j] = a[j][i];
		}
	}
	return ta;
}

/*取余*/
double mod(double a, double b)
{
	while (a<0)
	{
		a = a + b;
	}
	return a - (double)(((int)(a / b)) * b);
}

/*求数组平均值*/
double mean(double* xulx, int na)
{
	double s = 0;
	for (int i = 0; i < na; i++)
	{
		s = s + xulx[i];
	}
	return s / na;
}

/*求数组标准差*/
double stdd(double* xulx, int na)
{
	double s = 0;
	double m = mean(xulx, na);
	for (int i = 0; i < na; i++)
	{
		s = s + pow((xulx[i] - m), 2);
	}
	s = sqrt(s / na);
	return s;
}

/*数组最大值*/
double amax(double* xulx, int na)
{
	double aa = -1e26;
	for (int i = 0; i < na; i++)
	{
		if (aa < xulx[i])
		{
			aa = xulx[i];
		}
	}
	return aa;
}

/*数组最小值*/
double amin(double* xulx, int na)
{
	double aa = 1e26;
	for (int i = 0; i < na; i++)
	{
		if (aa > xulx[i])
		{
			aa = xulx[i];
		}
	}
	return aa;
}

/*点乘*/
double dotpro(double* xulx1, double* xulx2, int na)
{
	double s = 0;
	for (int i = 0; i < na; i++)
	{
		s = s + xulx1[i] * xulx2[i];
	}
	return s;
}

/*叉乘*/
double* crosspro(double* xulx1, double* xulx2)
{
	double* s;
	s = (double*)malloc(3 * sizeof(double));
	s[0] = xulx1[1] * xulx2[2] - xulx1[2] * xulx2[1];
	s[1] = xulx1[2] * xulx2[0] - xulx1[0] * xulx2[2];
	s[2] = xulx1[0] * xulx2[1] - xulx1[1] * xulx2[0];
	return s;
}

/*向量归一化*/
double* nor(double* xulx, int na)
{
	double* s;
	s = (double*)malloc(na * sizeof(double));
	double a;
	a = dotpro(xulx, xulx, na);
	for (int i = 0; i < na; i++)
	{
		s[i] = xulx[i] / sqrt(a);
	}
	return s;
}

double r2d(double rad)
{
	return rad * 180 / pi;
}