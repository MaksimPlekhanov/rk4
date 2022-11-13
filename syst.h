#pragma once
#include <math.h>
#include <iostream>

struct Syst
{
	int param;
	double f(double x, double u)
	{
		if (param == 1)
			return 3 * u;
		if (param == 2)
			return (u * u) / (2 * x + x * x) + u - u * u * u * sin(10 * x);
	}
};

struct Syst2
{
	double F(double a, double ul)
	{
		return a * sqrt(ul * ul + 1.0);
	}
};