#pragma warning(disable:4996)
#include <stdio.h>
#include <math.h>
#include "juvf.h"
#include "dcf.h"
#include "lagrange.h"

#define tmin 0.349	//周期起始值
#define tmax 0.351	//周期结束值
#define tstep 0.000001	//周期步长
#define dp 1e-5
#define dlp 1e-3
#define dbp 1e-3
#define alpha 1	//学习率
#define stopc 0	//0：停止条件为误差，1：停止条件为迭代次数
#define stopin 200	//停止条件迭代次数
#define MAX_ITER 300	//最大迭代次数
//#define La_inter	//拉格朗日插值
#define Li_inter	//线性插值
//#define eq_out	//输出方程
//#define print_EQ	//打印方程
#define dcf_test	//输出相关的曲线

//#define Newton_iter	//牛顿法迭代求解
//#define Tra	//遍历求解
//#define Tra_int	//取整方法遍历
#define siga	//求恒星周期方差法

#define T_asteroid 1825	//小行星公转周期
//#define Fi	//先遍历i
#define Fj	//先遍历曲线，判断i
//#define T_test	//周期判定
#define T_JUR	//Jurkevich
//#define T_JUR_DCF	//Jurkevich_DCF
//#define lcf_out	//光变曲线输出
//#define lcf_out_jur_dcf	//光变曲线输出
#define spin_axis	//求自转轴
#define MAX_T 3600

/*jurkevich方法寻找周期*/
double jur(double** arr, int na)
{
	int m = 8;	//jurkevich方法分组数量
	double** jur;
	jur = (double**)malloc(m * sizeof(double));
	if (jur == NULL)
		exit(20);
	for (int i = 0; i < m; i++)
	{
		jur[i] = (double*)malloc(na * sizeof(double));
		if (jur[i] == NULL)
			exit(20);
	}
	int* mn;	//每段里的数据点数量
	mn = (int*)malloc(m * sizeof(int));
	if (mn == NULL)
		exit(20);
	double* Xl;	//每段平均值
	Xl = (double*)malloc(m * sizeof(double));
	if (Xl == NULL)
		exit(20);
	double* Vl2;	//每段离差平方和
	Vl2 = (double*)malloc(m * sizeof(double));
	if (Vl2 == NULL)
		exit(20);
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


#ifdef Fi
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < na; j++)
			{
				if ((int)(mod(arr[j][0], t) / (t / m)) == i)
				{
					jur[i][mn[i]] = arr[j][1];
					mn[i]++;
				}
			}
		}
#endif	//Fi


#ifdef Fj
		for (int j = 0; j < na; j++)
		{
			int i = (int)(mod(arr[j][0], t) / (t / m));
			switch (i)
			{
			case 0:
			{
				jur[0][mn[0]] = arr[j][1];
				mn[0]++;
				break;
			}
			case 1:
			{
				jur[1][mn[1]] = arr[j][1];
				mn[1]++;
				break;
			}
			case 2:
			{
				jur[2][mn[2]] = arr[j][1];
				mn[2]++;
				break;
			}
			case 3:
			{
				jur[3][mn[3]] = arr[j][1];
				mn[3]++;
				break;
			}
			case 4:
			{
				jur[4][mn[4]] = arr[j][1];
				mn[4]++;
				break;
			}
			case 5:
			{
				jur[5][mn[5]] = arr[j][1];
				mn[5]++;
				break;
			}
			case 6:
			{
				jur[6][mn[6]] = arr[j][1];
				mn[6]++;
				break;
			}
			case 7:
			{
				jur[7][mn[7]] = arr[j][1];
				mn[7]++;
				break;
			}
			default:
				break;
			}
		}
#endif //Fj


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


#ifdef T_test
	double jp[20] = { 0 };
	int k = 0;
	for (int i = 0; i < 20; i++)
	{
		k = 0;
		for (int j = 0; j < na; j++)
		{
			if ((int)(mod(arr[j][0], T) / (T / 20)) == i)
			{
				jp[i] += arr[j][1];
				k++;
			}
		}
		jp[i] = jp[i] / k;
	}
	int change = 0;
	for (int i = 0; i < 20 - 2; i++)
	{
		if ((jp[i] - jp[i + 1] < 0 && jp[i + 1] - jp[i + 2]>0) || (jp[i] - jp[i + 1] > 0 && jp[i + 1] - jp[i + 2] < 0))
		{
			change++;
		}
	}
	printf("%d\n", change);
	if (change > 5)
	{
		T = T / 2;
	}
	if (change < 3)
	{
		T = T * 2;
	}
	if (change > 8)
	{
		T = T / 3;
	}
#endif	//T_test


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

double Jur_DCF(double*** lp, double*** lpc, int nlc, int** nlp)
{
	int M = 0;
	double dcf_r = 0;
	double d = 0;
	double temp = -1e26;
	double T = 0;
	for (double t = tmin; t < tmax; t = t + tstep)
	{
		dcf_r = 0;
		M = 0;
		d = 0;
		for (int i = 0; i < nlc; i++)
		{
			for (int j = 0; j < nlp[i][0]; j++)
			{
				lpc[i][j][0] = lp[i][j][0] - (double)((int)(lp[i][j][0] / t) * t);
			}
		}
		for (int i = 0; i < nlc - 1; i++)
		{
			for (int j = 1; j < nlc; j++)
			{
				if (min(lpc[j][nlp[j][0] - 1][0], lpc[i][nlp[i][0] - 1][0]) - max(lpc[j][0][0], lpc[i][0][0] > 0.01))
				{
					d = DCF_jur(lpc[i], lpc[j], nlp[i][0], nlp[j][0], 0, t * 0.01);
					if (d != -13)
					{
						dcf_r += d;
						M++;
					}
				}
			}
		}
		dcf_r = dcf_r / M;
		if (temp < dcf_r)
		{
			temp = dcf_r;
			T = t;
		}
	}
	return T;
}

double dN(double ti, double tj, double tsyn)
{
	double a;
	a = (double)((int)((ti - tj) / tsyn-0.5));
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
	if (rse == NULL)
		exit(20);
	for (int i = 0; i < 3; i++)
	{
		rse[i] = ren[i] + rsn[i];
	}
	double* rh;	//角平分线黄道坐标
	rh = (double*)malloc(3 * sizeof(double));
	if (rh == NULL)
		exit(20);
	double rr;
	rr = sqrt(dotpro(rse, rse, 3));
	for (int i = 0; i < 3; i++)
	{
		rh[i] = rse[i] / rr;
	}
	double* rc;
	rc = (double*)malloc(3 * sizeof(double));
	if (rc == NULL)
		exit(20);
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

double f(double a, double Ti, double Tj, double Psid, double* rs1, double* re1, double* rs2, double* re2, double lp, double bp, double Psyn)
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

	if (L1 - L2 < -pi)
	{
		return f + 1;
	}
	else
	{
		if (L1 - L2 > pi)
		{
			return f - 1;
		}
		else
		{
			return f;
		}
	}
}

double dLdlp(double* r, double lp, double bp)
{
	return (r[2] * sin(bp) - sin(BB(r))) * cos(LL(r)) / (cos(bp) * cos(BB(r))) - sin(LL(r)) * sin(LL(r)) * sin(bp);
}

double dLdbp(double* r, double lp, double bp)
{
	return -sin(LL(r)) * tan(BB(r));
}

double dfdPsid(double Ti, double Tj, double Psid)
{
	return (Ti - Tj) / pow(Psid, 2);
}

double dfdlp(double a, double* rs1, double* re1, double* rs2, double* re2, double lp, double bp)
{
	double* ri, * rj;
	ri = rc(re1, rs1, lp, bp);
	rj = rc(re2, rs2, lp, bp);
	double ddd = a / (2 * pi) * (dLdlp(ri, lp, bp) - dLdlp(rj, lp, bp));

	free(ri);
	free(rj);

	return ddd;
}

double dfdbp(double a, double* rs1, double* re1, double* rs2, double* re2, double lp, double bp)
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
	/*char* olcname;
	olcname = (char*)malloc(100 * sizeof(char));
	if (olcname == NULL)
		exit(20);
	printf("输入光变曲线文件路径：");
	(void)scanf("%s", olcname);
	FILE* olc;
	olc = fopen(olcname, "r");*/

	FILE* olc;
	olc = fopen("D:\\wuli\\小行星反演\\version_0.2.1_副本\\monilc\\monilc.txt", "r");

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
	double*** lpc_j;
	lpc_j = (double***)malloc(nlc * sizeof(double));

	double* lmax, * lmin;	//亮度最大值,最小值
	lmax = (double*)malloc(nlc * sizeof(double));
	lmin = (double*)malloc(nlc * sizeof(double));

	double* delt;	//最大时间间隔
	delt = (double*)malloc(nlc * sizeof(double));

	int* qs = (int*)malloc(nlc * sizeof(int));	//去重的数量

	double jdmax = -1e26;	//光变曲线时间最大值
	double jdmin = 1e26;	//光变曲线时间最小值

	for (int i = 0; i < nlc; i++)
	{
		lmax[i] = -1e9;
		lmin[i] = 1e9;
		delt[i] = -1e26;
		qs[i] = 0;
		(void)fscanf(olc, "%d %d", &nlp[i][0], &nlp[i][1]);
		lp[i] = (double**)malloc(nlp[i][0] * sizeof(double));
		lpc_j[i] = (double**)malloc(nlp[i][0] * sizeof(double));
		for (int j = 0; j < nlp[i][0]; j++)
		{
			lp[i][j] = (double*)malloc(8 * sizeof(double));
			lpc_j[i][j] = (double*)malloc(8 * sizeof(double));
			for (int k = 0; k < 8; k++)
			{
				(void)fscanf(olc, "%lf", &lp[i][j][k]);
				lpc_j[i][j][k] = lp[i][j][k];
			}
			if (lmax[i] < lp[i][j][1])
			{
				lmax[i] = lp[i][j][1];
			}
			if (lmin[i] > lp[i][j][1])
			{
				lmin[i] = lp[i][j][1];
			}
			if (jdmax < lp[i][j][0])
			{
				jdmax = lp[i][j][0];
			}
			if (jdmin > lp[i][j][0])
			{
				jdmin = lp[i][j][0];
			}
			np++;
			if (j > 0)
			{
				if (delt[i] < lp[i][j][0] - lp[i][j - 1][0])
				{
					delt[i] = lp[i][j][0] - lp[i][j - 1][0];
				}
				/*if (lp[i][j][0] - lp[i][j - 1][0] < 0.0005)	//去重
				{
					np--;
					j--;
					nlp[i][0]--;
					qs[i]++;
				}*/
			}
		}
	}

	fclose(olc);

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


#ifdef T_JUR
	double Tsyn = jur(lpc, npc);	//会合周期
	printf("会合周期：%f天\n", Tsyn);
#endif	//T_JUR


#ifdef T_JUR_DCF
	double Tsyn_J_D = Jur_DCF(lp, lpc_j, nlc, nlp);
	printf("%f\n", Tsyn_J_D);
#endif	//T_JUR_DCF


#ifdef lcf_out
	FILE* lcf;
	lcf = fopen("E:\\lcf.txt", "w");
	for (int i = 0; i < nlc; i++)
	{
		for (int j = 0; j < nlp[i][0]; j++)
		{
			fprintf(lcf, "%f %f\n", mod(lp[i][j][0], Tsyn), lp[i][j][1]);
		}
	}
#endif // lcf_out


#ifdef lcf_out_jur_dcf
	FILE* lcf;
	lcf = fopen("E:\\lcf.txt", "w");
	for (int i = 0; i < nlc; i++)
	{
		for (int j = 0; j < nlp[i][0]; j++)
		{
			fprintf(lcf, "%f %f\n", lp[i][j][0] - (lp[i][0][0] - mod(lp[i][0][0], Tsyn_JUR)), lp[i][j][1]);
		}
	}
	FILE* lcf1;
	lcf1 = fopen("E:\\lcf_JD.txt", "w");
	for (int i = 0; i < nlc; i++)
	{
		for (int j = 0; j < nlp[i][0]; j++)
		{
			fprintf(lcf1, "%f %f\n", lp[i][j][0] - (lp[i][0][0] - mod(lp[i][0][0], Tsyn_J_D)), lp[i][j][1]);
		}
	}
#endif // lcf_out_jur_dcf



#ifdef spin_axis
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
	double* lag;	//时间延迟
	lag = (double*)malloc(nlc * (nlc - 1) / 2 * sizeof(double));
	int nxy = 0;	//方程个数
	int flag = 0;	//检验是否线性相关
	double* par;	//拉格朗日插值结果
	for (int i = 0; i < nlc - 1; i++)
	{
		if (delt[i] < 0.02 && lp[i][nlp[i][0] - 1][0] - lp[i][0][0] > 0.1)	//观测间隔在15分钟内
		{
			for (int j = i + 1; j < nlc; j++)
			{
				flag = 0;
				for (int xyi = 0; xyi < nxy; xyi++)
				{
					for (int xyj = 0; xyj < nxy; xyj++)
					{
						if ((xy[xyi][0] == xy[xyj][0] && xy[xyi][1] == i && xy[xyj][1] == j) || (xy[xyi][1] == xy[xyj][1] && xy[xyi][0] == i && xy[xyj][0] == j))	//线性相关
						{
							flag = 1;
						}
					}
				}
				if (lp[j][0][0] - lp[i][0][0] > MAX_T)	//时间延迟大于1000天
					break;
				if (!flag && delt[j] < 0.02 && lp[j][nlp[j][0] - 1][0] - lp[j][0][0] > 0.1 && lp[i][0][0] < lp[j][0][0])	//观测间隔在15分钟内且时间延迟最大值为180天
				{
					double ddcf = timedelay_strf(lp[i], lp[j], nlp[i][0], nlp[j][0], Tsyn);	//时间延迟
					if (ddcf != 0)
					{
						xy[nxy][0] = i;
						xy[nxy][1] = j;
						lag[nxy] = ddcf;
						tsyn[nxy][0] = (min(lp[i][nlp[i][0] - 1][0], lp[j][nlp[j][0] - 1][0] - ddcf) - max(lp[i][0][0], lp[j][0][0] - ddcf)) / 2 + max(lp[i][0][0], lp[j][0][0] - ddcf);
						tsyn[nxy][1] = tsyn[nxy][0] + ddcf;


#ifdef La_inter
						for (int k = 0; k < nlp[i][0] - 1; k++)
						{
							if (lp[i][k][0] <= tsyn[nxy][0] && lp[i][k + 1][0] >= tsyn[nxy][0])
							{
								par = Lagrange(lp[i] + k - 2, 4, 8, tsyn[nxy][0]);
								for (int ki = 0; ki < 3; ki++)
								{
									rs[nxy][0][ki] = par[ki + 2];
									re[nxy][0][ki] = par[ki + 5];
								}
								free(par);
								break;
							}
						}
#endif // La_inter


						/*par = Lagrange(lp[i], nlp[i][0], 8, tsyn[nxy][0]);
						for (int k = 0; k < 3; k++)
						{
							rs[nxy][0][k] = par[k + 2];
							re[nxy][0][k] = par[k + 5];
						}
						free(par);*/


#ifdef Li_inter
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
#endif // Li_inter


#ifdef La_inter
						for (int k = 0; k < nlp[j][0] - 1; k++)
						{
							if (lp[j][k][0] <= tsyn[nxy][1] && lp[j][k + 1][0] >= tsyn[nxy][1])
							{
								par = Lagrange(lp[j] + k - 2, 4, 8, tsyn[nxy][1]);
								for (int ki = 0; ki < 3; ki++)
								{
									rs[nxy][1][ki] = par[ki + 2];
									re[nxy][1][ki] = par[ki + 5];
								}
								free(par);
								break;
							}
						}
#endif // La_inter


						/*par = Lagrange(lp[j], nlp[j][0], 8, tsyn[nxy][1]);
						for (int k = 0; k < 3; k++)
						{
							rs[nxy][1][k] = par[k + 2];
							re[nxy][1][k] = par[k + 5];
						}
						free(par);*/


#ifdef Li_inter
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
#endif // Li_inter


						nxy++;

					}
				}
			}
		}
	}
	printf("方程数量：%d\n", nxy);
	if (nxy < 5)
	{
		printf("未找到足够数量的方程。");
		exit(20);
	}


#ifdef print_EQ
	for (int i = 0; i < nxy; i++)
	{
		printf("%d %d %f %f %f\n", xy[i][0], xy[i][1], lag[i], tsyn[i][0], tsyn[i][1]);
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
#endif	//print_EQ


#ifdef eq_out
	FILE* eq;
	eq = fopen("E:\\eq.txt", "w");
	fprintf(eq, "%d\n", nxy);
	for (int i = 0; i < nxy; i++)
	{
		fprintf(eq, "%.6f ", tsyn[i][0]);
		for (int j = 0; j < 3; j++)
		{
			fprintf(eq, "%.8f ", rs[i][0][j]);
		}
		for (int j = 0; j < 3; j++)
		{
			fprintf(eq, "%.8f ", re[i][0][j]);
		}
		fprintf(eq, "%.6f ", tsyn[i][1]);
		for (int j = 0; j < 3; j++)
		{
			fprintf(eq, "%.8f ", rs[i][1][j]);
		}
		for (int j = 0; j < 3; j++)
		{
			fprintf(eq, "%.8f ", re[i][1][j]);
		}
		fprintf(eq, "\n");
	}
	fclose(eq);
#endif	//eqout


#ifdef dcf_test
	FILE* f1, * f2;
	for (int j = 0; j < nxy; j++)
	{
		char tname1[50] = "E:\\dcf_test\\";
		char tname2[50] = "E:\\dcf_test\\";
		char s[10];
		char s1[10] = "1.txt";
		char s2[10] = "2.txt";
		itoa(j, s, 10);
		strcat(tname1, s);
		strcat(tname1, s1);
		strcat(tname2, s);
		strcat(tname2, s2);

		f1 = fopen(tname1, "w");
		f2 = fopen(tname2, "w");

		for (int i = 0; i < nlp[xy[j][0]][0]; i++)
		{
			fprintf(f1, "%.4f %.4f\n", lp[xy[j][0]][i][0], lp[xy[j][0]][i][1]);
		}
		for (int i = 0; i < nlp[xy[j][1]][0]; i++)
		{
			fprintf(f2, "%.4f %.4f\n", lp[xy[j][1]][i][0] - lag[j], lp[xy[j][1]][i][1]);
		}
		fclose(f1);
		fclose(f2);
	}
#endif	//dcftest


#ifdef Newton_iter
	/*自转轴指向初值*/
	double chi2 = 0;
	double temp = 1e26;
	double Temp = 1e26;
	double Y[3] = { 0 };
	double A = 0;
	double dT = 0.5 * Tsyn * Tsyn / (jdmax - jdmin) * 0.8;

	int fsh = 0;

	double la, be;
	double dx[3] = { 0 };
	int in = 0;
	double psid = Tsyn;
	double X[3] = { 0 };
	for (double T = Tsyn; T < Tsyn + dT; T = T + dT)
	{
		Temp = 1e26;
		for (double a = -1.0; a < 2.0; a = a + 2.0)
		{
			//double** F, ** J, ** Jt, ** JtJ, ** JtF, ** JtJin, ** epsi;
			for (double lai = 0; lai < 2 * pi; lai += pi / 3)
			{
				for (double bei = -pi / 2; bei < pi / 2; bei += pi / 6)
				{
					la = lai;
					be = bei;
					psid = T;
					printf("%f %f\n", r2d(la), r2d(be));
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
						psid = psid - alpha * dx[0];
						la = la - alpha * dx[1];
						be = be - alpha * dx[2];
						chi2 = 0;
						for (int i = 0; i < nxy; i++)
						{
							F[i][0] = f(a, tsyn[i][0], tsyn[i][1], psid, rs[i][0], re[i][0], rs[i][1], re[i][1], la, be, Tsyn);
							chi2 = chi2 + pow(F[i][0], 2);
							J[i][0] = dfdPsid(tsyn[i][0], tsyn[i][1], psid);
							J[i][1] = dfdlp(a, rs[i][0], re[i][0], rs[i][1], re[i][1], la, be);
							J[i][2] = dfdbp(a, rs[i][0], re[i][0], rs[i][1], re[i][1], la, be);
						}
						Jt = TA(J, nxy, 3);
						JtJ = AB(Jt, 3, nxy, J, nxy, 3);
						JtJin = inv(JtJ, 3);
						JtF = AB(Jt, 3, nxy, F, nxy, 1);
						epsi = AB(JtJin, 3, 3, JtF, 3, 1);
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
						//printf("%15.10f %15.10f %15.10f\n", dx[0], dx[1], dx[2]);
						//printf("%f %f %f\n", psid * 24, r2d(la), r2d(be));
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
						if (in > MAX_ITER)
							break;
					} while (fsh);
					printf("x残差小：%d %f %f %f %f\n", (int)a, X[0] * 24, r2d(fmod(X[1], 2 * pi)), r2d(asin(sin(X[2]))), temp);
					printf("x迭代后：%d %f %f %f %f\n", (int)a, psid * 24, r2d(fmod(la, 2 * pi)), r2d(asin(sin(be))), chi2);
					printf("%d\n\n", in);
					if (Temp > temp)
					{
						Temp = temp;
						Y[0] = X[0];
						Y[1] = X[1];
						Y[2] = X[2];
						A = a;
					}
				}
			}
		}
		if (Y[0]>Tsyn)
		{
			Y[1] += 3 * pi;
			Y[2] = -Y[2];
		}
		printf("%f %f %f %f\n", Y[0] * 24, r2d(fmod(Y[1], 2 * pi)), r2d(asin(sin(Y[2]))), Temp);
	}
#endif //Newton_iter


#ifdef Tra
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
		for (double psid = T_asteroid / (T_asteroid / (Tsyn / 24) + atan(1 / 3) / (2 * pi) + 1); psid < T_asteroid / (T_asteroid / (Tsyn / 24) - atan(1 / 3) / (2 * pi) - 1); psid = psid + 1e-7)
		{
			for (double la = 0; la < 2 * pi; la = la + 0.017)
			{
				for (double be = -pi / 2; be < pi / 2; be = be + 0.017)
				{
					chi2 = 0;
					for (int i = 0; i < nxy; i++)
					{
						chi2 += pow(f(a, tsyn[i][0], tsyn[i][1], psid, rs[i][0], re[i][0], rs[i][1], re[i][1], la, be, Tsyn / 24), 2);
					}
					if (temp > chi2)
					{
						temp = chi2;
						X[0] = psid * 24;
						X[1] = r2d(la);
						X[2] = r2d(be);
					}
				}
			}
		}
		printf("%f %.6f %f %f\n", a, X[0], X[1], X[2]);
	}
#endif // Tra


#ifdef Tra_int
	double temp = 1e26;
	double X[3] = { 0 };
	double chi2 = 0;
	double* r1, * r2;
	double L1 = 0, L2 = 0;
	for (double psid = Tsyn - 0.001; psid < Tsyn + 0.001; psid = psid + tstep * 10)
	{
		for (double la = 0; la < 2 * pi; la = la + 0.05)
		{
			for (double be = -pi / 2; be < pi / 2; be = be + 0.05)
			{
				chi2 = 0;
				for (int i = 0; i < nxy; i++)
				{
					r1 = rc(re[i][0], rs[i][0], la, be);
					r2 = rc(re[i][1], rs[i][1], la, be);
					L1 = LL(r1);
					L2 = LL(r2);
					chi2 += pow((tsyn[i][0] - tsyn[i][1]) / psid - (L1 - L2) / (2 * pi) - (double)((int)((tsyn[i][0] - tsyn[i][1]) / psid - (L1 - L2) / (2 * pi) + 0.5)), 2);
					free(r1);
					free(r2);
				}
				if (temp > chi2)
				{
					temp = chi2;
					X[0] = psid * 24;
					X[1] = r2d(la);
					X[2] = r2d(be);
				}
			}
		}
	}
	printf("%f %f %f\n", X[0], X[1], X[2]);
#endif // Tra_int


#ifdef siga
	FILE* uix = fopen("E:\\uix.txt", "w");
	double T = 0;
	double temp = 1e26;
	double lambda = 0, beta = 0;
	for (double lambdap = 0; lambdap < 2*pi; lambdap = lambdap + 0.08)
	{
		for (double betap = -pi/2; betap <= pi/2; betap = betap + 0.08)
		{
			int n = 0;
			double s = 0;
			double* psid = (double*)malloc(nxy * sizeof(double));	//恒星周期
			double*** Rc = (double***)malloc(nxy * sizeof(double));	//赤道坐标
			double* mark = (double*)malloc(nxy * sizeof(double));	//去除标准差之外的恒星周期
			for (int i = 0; i < nxy; i++)
			{
				Rc[i] = (double**)malloc(2 * sizeof(double));
				Rc[i][0] = rc(re[i][0], rs[i][0], lambdap, betap);
				Rc[i][1] = rc(re[i][1], rs[i][1], lambdap, betap);
			}
			for (int i = 0; i < nxy; i++)
			{
				double ll1 = LL(Rc[i][0]);
				double ll2 = LL(Rc[i][1]);
				double dl = ll1 - ll2;	//赤经变化量
				if (ll1 - ll2 > pi)
				{
					dl = dl - 2 * pi;
				}
				if (ll1 - ll2 < -pi)
				{
					dl = dl + 2 * pi;
				}
				psid[i] = (tsyn[i][0] - tsyn[i][1]) / (dN(tsyn[i][0], tsyn[i][1], Tsyn) + (dl)/2/pi+(double)((int)((tsyn[i][0]-tsyn[i][1])/T_asteroid)));
			}
			double XXX = stdd(psid, nxy);	//标准差
			double me = mean(psid, nxy);	//平均值
			for (int i = 0; i < nxy; i++)
			{
				if (Abs(psid[i] - me) < XXX)	//标准差之内
				{
					mark[n] = psid[i];
					n++;
				}
			}
			XXX = stdd(mark, n);
			if (temp > XXX)
			{
				temp = XXX;
				T = mean(mark, n) * 24;
				lambda = r2d(lambdap);
				beta = r2d(betap);
			}
			fprintf(uix, "%f %f %f %.5e lp, bp, T, 方差\n", r2d(lambdap), r2d(betap), mean(mark, n), XXX*XXX);
			for (int i = 0; i < n; i++)
			{
				fprintf(uix, "%f\n", mark[i]);
			}
			for (int i = 0; i < nxy; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					free(Rc[i][j]);
				}
				free(Rc[i]);
			}
			free(psid);
			free(Rc);
			free(mark);
		}
	}
	printf("%f %f %f\n", T, lambda, beta);
	fprintf(uix, "%f %f %f\n", T, lambda, beta);
	fclose(uix);
#endif	//siga


	for (int i = 0; i < nlc; i++)
	{
		for (int j = 0; j < nlp[i][0]; j++)
		{
			free(lp[i][j]);
		}
		free(lp[i]);
		free(nlp[i]);
	}
	free(lp);
	free(nlp);
	free(qs);

	free(lmax);
	free(lmin);
	free(delt);
#endif	//spin_axis


	return 0;
}