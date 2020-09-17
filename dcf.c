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
	double** co = (double**)malloc((n1 + n2) * sizeof(double));
	for (int i = 0; i < n1 + n2; i++)
	{
		co[i] = (double*)malloc(2 * sizeof(double));
	}
	double t2 = min(arr1[n1 - 1][0], arr2[n2 - 1][0] - T);
	double t1 = max(arr1[0][0], arr2[0][0] - T);
	double Vm2 = 0;
	int nn = 0;
	if ((t2 - t1) * 24 > 1)
	{
		//printf("%f %f\n", t2, t1);
		
		for (int i = 0; i < n1; i++)
		{
			if (arr1[i][0] > t1 && arr1[i][0] < t2)
			{
				co[nn][0] = arr1[i][0];
				co[nn][1] = arr1[i][1];
				//printf("%f ", co[nn][1]);
				nn++;
			}
			
		}
		for (int i = 0; i < n2; i++)
		{
			if (arr2[i][0] - T > t1 && arr2[i][0] - T < t2)
			{
				co[nn][0] = arr2[i][0] - T;
				co[nn][1] = arr2[i][1];
				//printf("%f ", co[nn][1]);
				nn++;
			}
		}
		double** br = (double**)malloc(3 * sizeof(double));
		for (int i = 0; i < 3; i++)
		{
			br[i] = (double*)malloc((n1 + n2) * sizeof(double));
		}
		/*for (int i = 0; i < nn; i++)
		{
			printf("%f ", co[i][0]-t1);
		}*/
		int n[3] = { 0 };
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < nn; j++)
			{
				//printf("%f %f\n", co[i][0] - t1, (t2 - t1) / 3);
				if ((int)((co[j][0] - t1) / ((t2 - t1) / 3)) == i)
				{
					br[i][n[i]] = co[j][1];
					//printf("%f ", br[i][n[i]]);
					n[i]++;
				}
			}
		}
		//printf("\n");
		double X_[3] = { 0 };
		double Vl2[3] = { 0 };
		for (int i = 0; i < 3; i++)
		{
			//printf("%d ", n[i]);
			//X_[i] = mean(br[i], n[i]);
			//printf("%f ", X_[i]);
			for (int j = 0; j < n[i]; j++)
			{
				Vl2[i] += pow(br[i][j], 2);
			}
			Vl2[i] = Vl2[i] - n[i] * pow(X_[i], 2);
			Vm2 += Vl2[i];
		}
		//printf("\n");
		for (int i = 0; i < 3; i++)
		{
			free(br[i]);
		}
		free(br);
	}
	for (int i = 0; i < n1 + n2; i++)
	{
		free(co[i]);
	}
	free(co);

	if (temp > 0.98 && (min(arr1[n1 - 1][0], arr2[n2 - 1][0] - T) - max(arr1[0][0], arr2[0][0] - T)) * 24 > 4 && Vm2/nn < 1.5)
	{
		//printf("%f\n", Vm2/nn);
		return T;
	}
	else
	{
		return 0;
	}
}

double DCF_jur(double** arr1, double** arr2, int na1, int na2, double dt, double dtau)
{
	double temp = -1e50;
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
	if (M < 5)
	{
		return -13;
	}
	return d;
}