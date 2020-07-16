#pragma warning(disable:4996)
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "juvf.h"
#include "dcf.h"
#define pi 3.1415926535897932	//圆周率
#define e 2.71828182845904523	//自然常数
#define tmin 4.1	//周期起始值
#define tmax 4.2	//周期结束值
#define tstep 0.0001	//周期步长

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
	a = (double)((int)((ti - tj) / tsyn));
	return a;
}

/*二分线赤道坐标*/
double* rc(double* re, double* rs, double lp, double bp)
{
	double* ren, * rsn;	//归一化后坐标
	ren = nor(re, 3);
	rsn = nor(rs, 3);
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

double f(double a, double Ti, double Tj, double Psid, double *rs1, double *re1,double *rs2,double*re2, double lp, double bp,double Psyn)
{
	double dn = dN(Ti, Tj, Psyn);
	double* r1, * r2;
	r1 = rc(re1, rs1, lp, bp);
	r2 = rc(re2, rs2, lp, bp);
	double L1, L2;
	L1 = LL(r1);
	L2 = LL(r2);
	double f;
	f = dn + a * (L1 - L2) / (2 * pi) - (Ti - Tj) / Psid;

	free(r1);
	free(r2);

	return f;
}

double dLdlp(double* r, double lp, double bp)
{
	return (r[2] * sin(bp) - sin(BB(r))) * cos(LL(r)) / (cos(bp) * cos(BB(r))) - sin(LL(r)) * sin(LL(r)) * sin(bp);
	//return ((-r[0] * cos(lp) - r[1] * sin(lp)) * (r[0] * sin(bp) * cos(lp) + r[1] * sin(bp) * sin(lp) - r[2] * cos(bp)) - (-r[0] * sin(lp) + r[1] * cos(lp)) * (-r[0] * sin(bp) * sin(lp) + r[1] * sin(bp) * cos(lp)))
		// (pow(r[0] * sin(bp) * cos(lp) + r[1] * sin(bp) * sin(lp) - r[2] * cos(bp), 2) + pow(-r[0] * sin(lp) + r[1] * cos(lp), 2));
}

double dLdbp(double *r, double lp, double bp)
{
	return -sin(LL(r)) * tan(BB(r));
	//return -(-r[0] * sin(lp) + r[1] * cos(lp)) * (r[0] * cos(bp) * cos(lp) + r[1] * cos(bp) * sin(lp) + r[2] * sin(bp))
		// (pow(r[0] * sin(bp) * cos(lp) + r[1] * sin(bp) * sin(lp) - r[2] * cos(bp), 2) + pow(-r[0] * sin(lp) + r[1] * cos(lp), 2));
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
	char *olcname;
	olcname = (char*)malloc(100 * sizeof(char));
	printf("输入光变曲线文件路径：");
	scanf("%s", olcname);
	FILE* olc;
	olc = fopen(olcname, "r");

	int np = 0;
	int nlc;	//光变曲线数量
	(void)fscanf(olc, "%d", &nlc);
	int** nlp;	//数据点数量
	nlp = (int**)malloc(nlc * sizeof(int));
	for (int i = 0; i < nlc; i++)
	{
		nlp[i] = (int*)malloc(2 * sizeof(int));
	}
	double*** lp;	//数据点
	lp = (double***)malloc(nlc * sizeof(double));

	double* lmax, * lmin;
	lmax = (double*)malloc(nlc * sizeof(double));
	lmin = (double*)malloc(nlc * sizeof(double));
	for (int i = 0; i < nlc; i++)
	{
		lmax[i] = -1e9;
		lmin[i] = 1e9;
		(void)fscanf(olc, "%d %d", &nlp[i][0], &nlp[i][1]);
		lp[i] = (double**)malloc(nlp[i][0] * sizeof(double));
		for (int j = 0; j < nlp[i][0]; j++)
		{
			lp[i][j] = (double*)malloc(8 * sizeof(double));
			for (int k = 0; k < 8; k++)
			{
				(void)fscanf(olc, "%lf", &lp[i][j][k]);
			}
			if (lmax[i] < lp[i][j][1])
			{
				lmax[i] = lp[i][j][1];
			}
			if (lmin[i] > lp[i][j][1])
			{
				lmin[i] = lp[i][j][1];
			}
			np++;
		}
	}

	int npc = 0;
	double** lpc;
	lpc = (double**)malloc(np * sizeof(double));
	for (int i = 0; i < np; i++)
	{
		lpc[i] = (double*)malloc(8 * sizeof(double));
	}
	for (int i = 0; i < nlc; i++)
	{
		for (int j = 0; j < nlp[i][0]; j++)
		{
			for (int k = 0; k < 8; k++)
			{
				lpc[npc][k] = lp[i][j][k];
			}
			lpc[npc][1] = (lpc[npc][1] - lmin[i]) / (lmax[i] - lmin[i]);
			npc++;
		}
	}

	/*测试jurkevich*/
	/*double** xy;
	xy = (double**)malloc(1000 * sizeof(double));
	for (int i = 0; i < 1000; i++)
	{
		xy[i] = (double*)malloc(2 * sizeof(double));
	}
	for (int i = 0; i < 1000; i++)
	{
		xy[i][0] = 0.1 * i;
		xy[i][1] = sin(0.5 * xy[i][0]);
	}
	printf("%f\n", jur(xy, 1000));*/

	double Tsyn = jur(lpc, npc);	//会合周期
	printf("%f\n", Tsyn);

	double*** rs, *** re;	//相同特征时间点的太阳坐标与地球坐标
	double** tsyn;	//会合时间
	tsyn = (double**)malloc(nlc * (nlc - 1) / 2 * sizeof(double));
	rs = (double***)malloc(nlc * (nlc - 1) / 2 * sizeof(double));
	re = (double***)malloc(nlc * (nlc - 1) / 2 * sizeof(double));
	for (int i = 0; i < nlc * (nlc - 1) / 2; i++)
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

	int** xy;	//会合光变曲线
	xy = (int**)malloc(nlc * (nlc - 1) / 2 * sizeof(int));
	for (int i = 0; i < nlc * (nlc - 1) / 2; i++)
	{
		xy[i] = (int*)malloc(2 * sizeof(int));
	}
	double* timedelay;	//时间延迟
	timedelay = (double*)malloc(nlc * (nlc - 1) / 2 * sizeof(double));
	int nxy = 0;	//方程个数
	int flag = 0;	//检验是否线性相关
	for (int i = 0; i < nlc - 1; i++)
	{
		if (lp[i][nlp[i][0] - 1][0] - lp[i][0][0] < 1)	//一段光变曲线数据在一天之内
		{
			for (int j = i + 1; j < nlc; j++)
			{
				if (lp[j][nlp[j][0] - 1][0] - lp[j][0][0] < 1)	//一段光变曲线数据在一天之内
				{
					double ddcf = dcf(lp[i], lp[j], nlp[i][0], nlp[j][0], Tsyn / 24);
					if (ddcf != 0)
					{
						flag = 0;
						for (int xyi = 0; xyi < nxy; xyi++)
						{
							for (int xyj = 0; xyj < nxy; xyj++)
							{
								if (xy[xyi][0] == xy[xyj][0] && xy[xyi][1] == i && xy[xyj][1] == j)	//线性相关
								{
									flag = 1;
								}
							}
						}
						if (!flag)	//线性无关
						{
							xy[nxy][0] = i;
							xy[nxy][1] = j;
							timedelay[nxy] = ddcf;
							tsyn[nxy][0] = (min(lp[i][nlp[i][0] - 1][0], lp[j][nlp[j][0] - 1][0] - ddcf) - max(lp[i][0][0], lp[j][0][0] - ddcf)) / 2 + max(lp[i][0][0], lp[j][0][0] - ddcf);
							tsyn[nxy][1] = tsyn[nxy][0] + ddcf;
							for (int k = 0; k < nlp[i][0] - 1; k++)
							{
								if (lp[i][k][0] <= tsyn[nxy][0] && lp[i][k + 1][0] >= tsyn[nxy][0])
								{
									for (int ki = 0; ki < 3; ki++)
									{
										rs[nxy][0][ki] = (lp[i][k + 1][ki + 2] - lp[i][k][ki + 2]) / (lp[i][k + 1][0] - lp[i][k][0]) * (tsyn[nxy][0] - lp[i][k][0]) + lp[i][k][ki + 2];
										re[nxy][0][ki] = (lp[i][k + 1][ki + 5] - lp[i][k][ki + 5]) / (lp[i][k + 1][0] - lp[i][k][0]) * (tsyn[nxy][0] - lp[i][k][0]) + lp[i][k][ki + 5];
									}
									break;
								}
							}
							for (int k = 0; k < nlp[j][0] - 1; k++)
							{
								if (lp[j][k][0] <= tsyn[nxy][1] && lp[j][k + 1][0] >= tsyn[nxy][1])
								{
									for (int ki = 0; ki < 3; ki++)
									{
										rs[nxy][1][ki] = (lp[j][k + 1][ki + 2] - lp[j][k][ki + 2]) / (lp[j][k + 1][0] - lp[j][k][0]) * (tsyn[nxy][1] - lp[j][k][0]) + lp[j][k][ki + 2];
										re[nxy][1][ki] = (lp[j][k + 1][ki + 5] - lp[j][k][ki + 5]) / (lp[j][k + 1][0] - lp[j][k][0]) * (tsyn[nxy][1] - lp[j][k][0]) + lp[j][k][ki + 5];
									}
									break;
								}
							}
							nxy++;
						}
					}
				}
			}
		}
	}
	printf("%d\n", nxy);
	for (int i = 0; i < nxy; i++)
	{
		printf("%d %d %f %f %f\n", xy[i][0], xy[i][1], timedelay[i], tsyn[i][0], tsyn[i][1]);
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				printf("%f ", rs[i][j][k]);
			}
			for (int k = 0; k < 3; k++)
			{
				printf("%f ", re[i][j][k]);
			}
			printf("\n");
		}
	}

	/*自转轴指向初值*/
	double lai, bei;
	double chi2 = 0;
	double temp = 1e26;

	double la, be;
	double dx[3] = { 0 };	//迭代步长
	int in = 0;	//迭代次数
	double psid = Tsyn / 24;
	double X[3] = { 0 };	//解
	for (double a = -1.0; a < 2.0; a = a+2.0)
	{
		/*遍历初值*/
		temp = 1e26;
		for (double lp = 0; lp < 2 * pi; lp = lp + 0.05)
		{
			for (double bp = 0; bp < pi / 2; bp = bp + 0.05)
			{
				for (int i = 0; i < nxy; i++)
				{
					chi2 = chi2 + pow(f(a, tsyn[i][0], tsyn[i][1], Tsyn/24, rs[i][0], re[i][0], rs[i][1], re[i][1], lp, bp, Tsyn/24), 2);
				}
				if (temp > chi2)
				{
					temp = chi2;
					lai = lp;
					bei = bp;
				}
			}
		}
		in = 0;
		la = lai;
		be = bei;

		for (int i = 0; i < 3; i++)
		{
			dx[i] = 0;
		}
		temp = 1e26;
		do
		{
			double** J;	//雅可比矩阵
			double** F;
			F = (double**)malloc(nxy * sizeof(double));
			J = (double**)malloc(nxy * sizeof(double));
			for (int i = 0; i < nxy; i++)
			{
				J[i] = (double*)malloc(3 * sizeof(double));
				F[i] = (double*)malloc(sizeof(double));
			}
			double** Jt;
			double** JtJ;
			double** JtF;
			double** JtJin;
			double** epsi;
			epsi = (double**)malloc(3 * sizeof(double));
			for (int i = 0; i < 3; i++)
			{
				epsi[i] = (double*)malloc(sizeof(double));
			}
			for (int i = 0; i < 3; i++)
			{
				epsi[i][0] = 0;
			}

			chi2 = 0;
			la = la - dx[1];
			be = be - dx[2];
			psid = psid - dx[0];
			printf("x迭代中：%f %f %f\t", psid, la, be);
			for (int i = 0; i < nxy; i++)
			{
				J[i][0] = dfdPsid(tsyn[i][0], tsyn[i][1], psid);
				J[i][1] = dfdlp(a, rs[i][0], re[i][0], rs[i][1], re[i][1], la, be);
				J[i][2] = dfdbp(a, rs[i][0], re[i][0], rs[i][1], re[i][1], la, be);
				F[i][0] = f(a, tsyn[i][0], tsyn[i][1], psid, rs[i][0], re[i][0], rs[i][1], re[i][1], la, be, Tsyn / 24);
				chi2 = chi2 + pow(f(a, tsyn[i][0], tsyn[i][1], psid, rs[i][0], re[i][0], rs[i][1], re[i][1], la, be, Tsyn / 24), 2);
			}
			if (temp > chi2)
			{
				temp = chi2;
				X[0] = psid;
				X[1] = la;
				X[2] = be;
			}
			in++;
			Jt = TA(J, nxy, 3);
			JtJ = AB(Jt, 3, nxy, J, nxy, 3);
			JtJin = inv(JtJ, 3);
			JtF = AB(Jt, 3, nxy, F, nxy, 1);
			epsi = AB(JtJin, 3, 3, JtF, 3, 1);
			for (int i = 0; i < 3; i++)
			{
				dx[i] = epsi[i][0];
			}
			for (int i = 0; i < 3; i++)
			{
				free(epsi[i]);
				free(JtJin[i]);
				free(JtJ[i]);
				free(JtF[i]);
				free(Jt[i]);
			}
			free(epsi);
			free(JtJin);
			free(JtJ);
			free(JtF);
			free(Jt);
			for (int i = 0; i < nxy; i++)
			{
				free(J[i]);
				free(F[i]);
			}
			free(J);
			free(F);
			printf("dx:%f %f %f\n", dx[0],dx[1],dx[2]);
		} while (abs(dx[0]) > 0.00001 || abs(dx[1]) > 0.1 || abs(dx[2]) > 0.1);
		printf("x迭代后:%f %f %f\n", r2d(mod(X[1], 2 * pi)), r2d(asin(sin(X[2]))), X[0] * 24);
		printf("x残差最小：%f %f %f\n", psid*24, r2d(mod(la, 2 * pi))), r2d(asin(sin(be)));
	}


	return 0;
}