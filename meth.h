#pragma once
#include "syst.h"
#include <vector>

class Meth
{
protected:
	double x0, v0, x, v, h;
	int s;
public:
	Syst syst1;
	std::vector <std::vector<double>> solve;
	Meth(double _x0, double _v0, double _h, int _n);
	~Meth() {};
};

Meth::Meth(double _x0, double _v0, double _h, int _s)
{
	x0 = _x0;
	x = x0;
	v0 = _v0;
	v = v0;
	h = _h;
	s = _s;
}