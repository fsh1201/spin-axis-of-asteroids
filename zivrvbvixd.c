#pragma warning(disable:4996)
#include <stdio.h>
#include "juvf.h"
#include "dcf.h"
#include "dfdcf.h"
#include "lagrange.h"
#include <math.h>

#define tmin 5.9	//周期起始值
#define tmax 6	//周期结束值
#define tstep 0.0001	//周期步长
#define dp 1e-5
#define dlp 1e-3
#define dbp 1e-3
#define stopc 0	//0：停止条件为误差，1：停止条件为迭代次数
#define stopin 200	//最大迭代次数
#define eqN 10
#define Tsyn 0.3031695
#define T_asteroid 1378.9	//小行星公转周期
//#define newton	//牛顿迭代
#define tra	//遍历
#define fun	//矩阵函数
#define man	//

/*jurkevich方法寻找周期*/
double jur(double** arr, int na)
{
	int m = 8;	//jurkevich方法分组数量
	double** jur;
	jur = (double**)malloc(m * sizeof(double));
	for (int i = 0; i < m; i++)
	{
		jur[i] = (double*)malloc(na * sizeof(double));
	}
	int* mn;	//每段里的数据点数量
	mn = (int*)malloc(m * sizeof(int));
	double* Xl;	//每段平均值
	Xl = (double*)malloc(m * sizeof(double));
	double* Vl2;	//每段离差平方和
	Vl2 = (double*)malloc(m * sizeof(double));
	double Vm2 = 0;
	double temp = 1e26;
	double T = 0;
	for (double t = tmin; t < tmax; t = t + tstep)
	{
		Vm2 = 0;
		for (int i = 0; i < m; i++)
		{
			mn[i] = 0;
			Xl[i] = 0;
			Vl2[i] = 0;
		}
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < na; j++)
			{
				if ((int)(mod(arr[j][0], t/24) / (t/24 / m)) == i)
				{
					jur[i][mn[i]] = arr[j][1];
					mn[i]++;
				}
			}
		}
		for (int i = 0; i < m; i++)
		{
			Xl[i] = mean(jur[i], mn[i]);
			for (int j = 0; j < mn[i]; j++)
			{
				Vl2[i] = Vl2[i] + pow(jur[i][j], 2);
			}
			Vl2[i] = Vl2[i] - mn[i] * pow(Xl[i], 2);
			Vm2 = Vm2 + Vl2[i];
		}
		if (temp >= Vm2)
		{
			temp = Vm2;
			T = t;
		}
	}
	for (int i = 0; i < m; i++)
	{
		free(jur[i]);
	}
	free(jur);
	free(Xl);
	free(Vl2);
	free(mn);
	return T;
}

double dN(double ti, double tj, double tsyn)
{
	double a;
	a = (double)((int)((ti - tj) / tsyn - 0.5));
	return a;
}

/*二分线赤道坐标*/
double* rc(double* re, double* rs, double lp, double bp)
{
	double* ren, * rsn;	//归一化后坐标
	ren = nor(re, 3);
	//rsn = nor(rs, 3);
	rsn = (double*)malloc(3 * sizeof(double));
	double R = sqrt(pow(rs[0] - re[0], 2) + pow(rs[1] - re[1], 2) + pow(rs[2] - re[2], 2));
	double es[3];
	for (int i = 0; i < 3; i++)
	{
		es[i] = (rs[i] - re[i]) / R;
	}
	double cosls = es[0] / sqrt(es[0] * es[0] + es[1] * es[1]);
	double sinls = es[1] / sqrt(es[0] * es[0] + es[1] * es[1]);
	double Rh = sqrt(rs[0] * rs[0] + rs[1] * rs[1] + rs[2] * rs[2]);
	double Rg = sqrt(re[0] * re[0] + re[1] * re[1] + re[2] * re[2]);
	rsn[0] = R / Rh * cosls + Rg / Rh * re[0];
	rsn[1] = R / Rh * sinls + Rg / Rh * re[1];
	rsn[2] = Rg / Rh * re[2];
	double* rse;	//角平分线坐标
	rse = (double*)malloc(3 * sizeof(double));
	for (int i = 0; i < 3; i++)
	{
		rse[i] = ren[i] + rsn[i];
	}
	double* rh;	//角平分线黄道坐标
	rh = (double*)malloc(3 * sizeof(double));
	double rr;
	rr = sqrt(dotpro(rse, rse, 3));
	for (int i = 0; i < 3; i++)
	{
		rh[i] = rse[i] / rr;
	}
	double* rc;
	rc = (double*)malloc(3 * sizeof(double));
	rc[0] = rh[0] * sin(bp) * cos(lp) + rh[1] * sin(bp) * sin(lp) - rh[2] * cos(bp);
	rc[1] = -rh[0] * sin(lp) + rh[1] * cos(lp);
	rc[2] = rh[0] * cos(bp) * cos(lp) + rh[1] * cos(bp) * sin(lp) + rh[2] * sin(bp);

	free(ren);
	free(rsn);
	free(rse);
	free(rh);

	return rc;
}

/*二分线赤经*/
double LL(double* r)
{
	double l;
	l = atan(r[1] / r[0]);
	if (r[0] > 0 && r[1] < 0)
	{
		l = l + 2 * pi;
	}
	if (r[0] < 0)
	{
		l = l + pi;
	}
	return l;
}

/*二分线赤纬*/
double BB(double* r)
{
	double b;
	b = asin(r[2] / sqrt(dotpro(r, r, 3)));
	return b;
}

double f(double a, double Ti, double Tj, double Psid, double *rs1, double *re1,double *rs2,double*re2, double lp, double bp,double Psyn,double NN,double nn)
{
	double* r1, * r2;
	r1 = rc(re1, rs1, lp, bp);
	r2 = rc(re2, rs2, lp, bp);
	double L1, L2;
	L1 = LL(r1);
	L2 = LL(r2);
	free(r1);
	free(r2);
	double f;
	f = NN + a * ((L1 - L2) / (2 * pi) + nn) - (Ti - Tj) / Psid;
	return f;
}

double dLdlp(double* r, double lp, double bp)
{
	return (r[2] * sin(bp) - sin(BB(r))) * cos(LL(r)) / (cos(bp) * cos(BB(r))) - sin(LL(r)) * sin(LL(r)) * sin(bp);
}

double dLdbp(double *r, double lp, double bp)
{
	return -sin(LL(r)) * tan(BB(r));
}

double dfdPsid(double Ti, double Tj, double Psid)
{
	return (Ti - Tj) / pow(Psid, 2);
}

double dfdlp(double a, double* rs1, double *re1, double *rs2,double *re2,double lp, double bp)
{
	double* ri, * rj;
	ri = rc(re1, rs1, lp, bp);
	rj = rc(re2, rs2, lp, bp);
	double ddd= a / (2 * pi) * (dLdlp(ri, lp, bp) - dLdlp(rj, lp, bp));

	free(ri);
	free(rj);
	
	return ddd;
}

double dfdbp(double a, double* rs1, double* re1, double* rs2, double* re2,double lp,double bp)
{
	double* ri, * rj;
	ri = rc(re1, rs1, lp, bp);
	rj = rc(re2, rs2, lp, bp);
	double ddd = a / (2 * pi) * (dLdbp(ri, lp, bp) - dLdbp(rj, lp, bp));

	free(ri);
	free(rj);

	return ddd;
}



int main()
{
	double data[12][11] = {
		{2434569.7359, -3.237879954E-01, 9.636850742E-01, -1.526119477E-06, 5.152208989E-01, -1.155420138E+00, 2.743755859E-01, 288.1, 1.9, 0, 0},
	{2436599.7916, 6.403060498E-01, -7.487962495E-01, -7.145046607E-07, -1.492194507E+00, 8.200193030E-01, 5.930770902E-02, 151.20948, 1.9949335, 6696, 1},
	{2436616.7652, 8.354539750E-01, -5.277044849E-01, -4.602186595E-06, -1.412829949E+00, 9.246156086E-01, 1.025083114E-01, 146.79757, 3.4741631, 6752, 1},
	{2438552.1604, 2.661605233E-01, 9.792867301E-01, -3.917984512E-06, -4.952022877E-01, -1.463244005E+00, 6.161114354E-01, 251.30276, 21.743968, 13136, 3},
	{2439624.7835, 6.144010303E-01, 8.024949683E-01, 1.457588037E-06, -2.163927676E+00, 1.566890287E-01, 6.117868233E-01, 175.85847, 15.747511, 16674, 3},
	{2440018.8962, 1.719112504E-01, 1.000813546E+00, -3.419637746E-06, 3.048981565E-01, -1.402875449E+00, 4.529973251E-01, 282.26186, 17.512564, 17974, 4},
	{2440019.8052, 1.567415984E-01, 1.003399120E+00, -3.472548405E-06, 2.991846935E-01, -1.398028289E+00, 4.508948236E-01, 282.07936, 17.504279, 17977, 4},
	{2440083.7619, -8.003170343E-01, 6.204654438E-01, 6.240321544E-06, -8.060535542E-03, -1.521685988E+00, 2.863758558E-01, 269.69650, 10.658051, 18188, 4},
	{2441086.6687, 6.027677336E-01, 8.114798672E-01, -4.846188216E-06, -1.884530755E+00, -5.116899464E-01, 6.951379321E-01, 0, 0, 21496, 5},
	{2444550.0866, -7.140871391E-01, -6.871650167E-01, 7.453019583E-06, 5.446338612E-01, 7.615896451E-01, -5.064734592E-01, 0, 0, 32920, 7},
	{2445495.6473, 2.066596420E-01, 9.940093353E-01, -1.818511293E-06, -8.321825356E-02, -1.470939411E+00, 5.369818470E-01, 0, 0, 36039, 8}

	};

	double*** rs, *** re;	//相同特征时间点的太阳坐标与地球坐标
	double** tsyn;	//会合时间
	tsyn = (double**)malloc(eqN  * sizeof(double));
	rs = (double***)malloc(eqN * sizeof(double));
	re = (double***)malloc(eqN* sizeof(double));
	for (int i = 0; i < eqN; i++)
	{
		tsyn[i] = (double*)malloc(2 * sizeof(double));
		rs[i] = (double**)malloc(2 * sizeof(double));
		re[i] = (double**)malloc(2 * sizeof(double));
		for (int j = 0; j < 2; j++)
		{
			rs[i][j] = (double*)malloc(3 * sizeof(double));
			re[i][j] = (double*)malloc(3 * sizeof(double));
		}
	}

	double* nn;	//时间延迟
	nn = (double*)malloc(eqN * sizeof(double));
	double* NN;	//时间延迟
	NN = (double*)malloc(eqN * sizeof(double));
	int nxy = eqN;	//方程个数

	for (int i = 0; i < nxy; i++)
	{
		tsyn[i][0] = data[i][0];
		tsyn[i][1] = data[i + 1][0];
		for (int j = 0; j < 3; j++)
		{
			rs[i][0][j] = data[i][j + 1] - data[i][j + 4];
			rs[i][1][j] = data[i+1][j + 1] - data[i+1][j + 4];
			re[i][0][j] = -data[i][j + 1];
			re[i][1][j] = -data[i+1][j + 1];
		}
		nn[i] = data[i][10] - data[i + 1][10];
		NN[i] = data[i][9] - data[i + 1][9];
	}

	for (int i = 0; i < nxy; i++)
	{
		printf("%f ", tsyn[i][0]);
		for (int j = 0; j < 3; j++)
		{
			printf("%f ", rs[i][0][j]);
		}
		for (int j = 0; j < 3; j++)
		{
			printf("%f ", re[i][0][j]);
		}
		printf("%f ", tsyn[i][1]);
		for (int j = 0; j < 3; j++)
		{
			printf("%f ", rs[i][1][j]);
		}
		for (int j = 0; j < 3; j++)
		{
			printf("%f ", re[i][1][j]);
		}
		printf("%d ", NN[i]);
		printf("%d\n", nn[i]);
	}

#ifdef newton
	/*自转轴指向初值*/
	double chi2 = 0;
	double temp = 1e26;

	int fsh = 0;

	double la, be;
	double dx[3] = { 0 };
	int in = 0;
	double psid = Tsyn;
	double X[3] = { 0 };
	for (double a = -1.0; a < 2.0; a = a + 2.0)
	{
		psid = Tsyn;
		temp = 1e26;
		for (double l = 0; l <= 2 * pi; l = l + 0.1)
		{
			for (double b = -pi / 2; b < pi / 2; b = b + 0.1)
			{
				chi2 = 0;
				for (int i = 0; i < nxy; i++)
				{
					chi2 = chi2 + pow(f(a, tsyn[i][0], tsyn[i][1], psid, rs[i][0], re[i][0], rs[i][1], re[i][1], l, b, psid,(double)NN[i],(double)nn[i]), 2);
				}
				//printf("%f\n", chi2);
				if (temp > chi2)
				{
					temp = chi2;
					la = l;
					be = b;
				}
			}
		}
		printf("%f %f\n", la, be);

		in = 0;
		temp = 1e26;
		chi2 = 0;
		for (int i = 0; i < 3; i++)
		{
			dx[i] = 0;
		}
		do
		{
			in++;
			double** F, ** J, ** Jt, ** JtJ, ** JtF, ** JtJin, ** epsi;
			F = (double**)malloc(nxy * sizeof(double));
			J = (double**)malloc(nxy * sizeof(double));
			for (int i = 0; i < nxy; i++)
			{
				F[i] = (double*)malloc(sizeof(double));
				J[i] = (double*)malloc(3 * sizeof(double));
			}
			psid = psid - dx[0];
			la = la - dx[1];
			be = be - dx[2];
			chi2 = 0;
			for (int i = 0; i < nxy; i++)
			{
				F[i][0] = f(a, tsyn[i][0], tsyn[i][1], psid, rs[i][0], re[i][0], rs[i][1], re[i][1], la, be, Tsyn, (double)NN[i], (double)nn[i]);
				chi2 = chi2 + pow(F[i][0], 2);
				J[i][0] = dfdPsid(tsyn[i][0], tsyn[i][1], psid);
				J[i][1] = dfdlp(a, rs[i][0], re[i][0], rs[i][1], re[i][1], la, be);
				J[i][2] = dfdbp(a, rs[i][0], re[i][0], rs[i][1], re[i][1], la, be);
			}


#ifdef fun
			Jt = TA(J, nxy, 3);
			JtJ = AB(Jt, 3, nxy, J, nxy, 3);
			JtJin = inv(JtJ, 3);
			JtF = AB(Jt, 3, nxy, F, nxy, 1);
			epsi = AB(JtJin, 3, 3, JtF, 3, 1);
#endif	//fun


			for (int i = 0; i < 3; i++)
			{
				dx[i] = epsi[i][0];
			}
			if (temp > chi2)
			{
				temp = chi2;
				X[0] = psid;
				X[1] = la;
				X[2] = be;
			}
			for (int i = 0; i < nxy; i++)
			{
				free(F[i]);
				free(J[i]);
			}
			free(F);
			free(J);
			for (int i = 0; i < 3; i++)
			{
				free(Jt[i]);
				free(JtJ[i]);
				free(JtF[i]);
				free(epsi[i]);
				free(JtJin[i]);
			}
			free(Jt);
			free(JtJ);
			free(JtF);
			free(epsi);
			free(JtJin);
			printf("%15.10f %15.10f %15.10f\n", dx[0], dx[1], dx[2]);
			printf("%f %f %f\n", psid, r2d(la), r2d(be));
			fsh = 1;
			if (stopc == 0)
			{
				if (Abs(dx[0]) < dp && Abs(dx[1]) < dlp && Abs(dx[2]) < dbp)
				{
					fsh = 0;
				}
			}
			if (stopc == 1)
			{
				if (in > stopin)
					fsh = 0;
			}
			if (in > 2000)
				break;
		} while (fsh);
		printf("x残差小：%f %f %f\n", X[0], r2d(mod(X[1], 2 * pi)), r2d(asin(sin(X[2]))));
		printf("x迭代后：%f %f %f\n\n", psid, r2d(mod(la, 2 * pi)), r2d(asin(sin(be))));
		printf("%d\n", in);
	}
#endif //newton


#ifdef tra
	printf("%.6f %.6f\n", T_asteroid / (T_asteroid / (Tsyn)+atan(1 / 3) / (2 * pi) + 1), T_asteroid / (T_asteroid / (Tsyn)-atan(1 / 3) / (2 * pi) - 1));
	double X[3];
	double temp;
	double chi2;
	for (double a = -1.0; a < 2.0; a = a + 2.0)
	{
		temp = 1e26;
		for (int i = 0; i < 3; i++)
		{
			X[i] = 0;
		}
		for (double psid = T_asteroid / (T_asteroid / (Tsyn) + atan(1 / 3) / (2 * pi) + 1); psid < T_asteroid / (T_asteroid / (Tsyn) - atan(1 / 3) / (2 * pi) - 1); psid = psid + 1e-7)
		{
			for (double la = 0; la < 2 * pi; la = la + 0.017)
			{
				for (double be = -pi / 2; be < pi / 2; be = be + 0.017)
				{
					chi2 = 0;
					for (int i = 0; i < nxy; i++)
					{
						chi2 += pow(f(a, tsyn[i][0], tsyn[i][1], psid, rs[i][0], re[i][0], rs[i][1], re[i][1], la, be, Tsyn,NN[i],nn[i]), 2);
					}
					if (temp > chi2)
					{
						temp = chi2;
						X[0] = psid;
						X[1] = r2d(la);
						X[2] = r2d(be);
					}
				}
			}
		}
		printf("%f %.6f %f %f\n", a, X[0], X[1], X[2]);
	}
#endif //tra


	return 0;
}