#include "myMath.h"

///return the maximum of x and y
double dMax(double x, double y)
{
	if (x >= y)
		return x;
	else
		return y;
}

///return the absolute value of x
double dAbs(double x)
{
	if (x >= 0)
		return x;
	else
		return -x;
}


///if y >=0, return dAbs(x), else return -dAbs(x)
double dSign(double x, double y)
{
	if (y >= 0)
		return dAbs(x);
	else
		return -dAbs(x);
}

///return x-y 
double dDiff(double x, double y)
{
	return x - y;
}
