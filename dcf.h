#include "juvf.h"

/*离散相关函数*/
double DCF(double** arr1, double** arr2, int na1, int na2, double dt, double dtau);

/*DCF求时间延迟*/
double timedelay_DCF(double** arr1, double** arr2, int n1, int n2, double tsyn);

/*离散相关函数_Jurkevich*/
double DCF_jur(double** arr1, double** arr2, int na1, int na2, double dt, double dtau);

/*结构函数*/
double strf(double** arr1, double** arr2, int na1, int na2, double tau, double dt);

/*结构函数求时间延迟*/
double timedelay_strf(double** arr1, double** arr2, int n1, int n2, double tsyn);