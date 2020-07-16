#define pi 3.1415926535897932	//圆周率
double hhlx(double** arr, int na);	//行列式
double** inv(double** a, int n);	//逆矩阵
double** AB(double** a, int ma, int na, double** b, int mb, int nb);	//矩阵乘法
double** TA(double** a, int ma, int na);	//矩阵转置
double mod(double a, double b);	//取余
double mean(double* xulx, int na);	//数组平均值
double amax(double* xulx, int na);	//数组最大值
double amin(double* xulx, int na);	//数组最小值
double stdd(double* xulx, int na);	//数组标准差
double dotpro(double* xulx1, double* xulx2, int na);	//点乘
double* crosspro(double* xulx1, double* xulx2);	//叉乘
double* nor(double* xulx, int na);	//向量归一化
double r2d(double rad);