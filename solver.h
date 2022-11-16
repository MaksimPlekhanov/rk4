#pragma once
#include "syst.h"
#include <vector>

class Solver
{
protected:
	double x0, v0, x, v, h;
	double v_v2, eps_e, eps_b, b;
	int n_max, n_counter, c1, c2;
	bool flag_cl;
	virtual void step(double &x, double &v, double &h) = 0;
	void clear();
	Syst syst1;
public:
	std::vector <std::vector<double>> solve;
	Solver(int p, double _x0, double _v0, double _h, int _n_max, double _eps_e, double _eps_b, double _b);
	~Solver() {};
};

Solver::Solver(int p, double _x0, double _v0, double _h, int _n_max, double _eps_e, double _eps_b, double _b)
{
	syst1.param = p;
	x0 = _x0;
	x = x0;
	v0 = _v0;
	v = v0;
	h = _h;
	v_v2 = 0;
	eps_e = _eps_e;
	eps_b = _eps_b;
	b = _b;
	n_max = _n_max;
	n_counter = 0;
	flag_cl = 1;
}

void Solver::clear()
{
	solve.clear();
	flag_cl = 1;
}