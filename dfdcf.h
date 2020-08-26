#include "juvf.h"
#include "dcf.h"

/*离散傅里叶变换*/
double** DFT(double** arr, int n);

/*z-DCF*/
double DFDCF(double** arr1, double** arr2, int na1, int na2, double dt, double tsyn);

/*DFDCF求时间延迟*/
double timedelay_DFDCF(double** arr1, double** arr2, int n1, int n2, double tsyn);