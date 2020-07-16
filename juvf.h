#define pi 3.1415926535897932	//圆周率
/*行列式*/
double hhlx(double** arr, int na);
/*矩阵求逆*/
double** inv(double** a, int n);
/*矩阵相乘*/
double** AB(double** a, int ma, int na, double** b, int mb, int nb);
/*矩阵转置*/
double** TA(double** a, int ma, int na);
/*取余*/
double mod(double a, double b);
/*数组平均值*/
double mean(double* xulx, int na);
/*数组最大值*/
double amax(double* xulx, int na);
/*数组最小值*/
double amin(double* xulx, int na);
/*数组标准差*/
double stdd(double* xulx, int na);
/*点乘*/
double dotpro(double* xulx1, double* xulx2, int na);
/*叉乘*/
double* crosspro(double* xulx1, double* xulx2);
/*向量归一化*/
double* nor(double* xulx, int na);
/*弧度转角度*/
double r2d(double rad);