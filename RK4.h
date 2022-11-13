#pragma once
#include "meth.h"

class RK4 : public Meth
{
public:
	double k1, k2, k3, k4;
	RK4(double _x0, double _v0, double _h, int _n, int _p) : Meth(_x0, _v0, _h, _n) 
	{
		syst1.param = _p;
		std::vector<double> temp;
		for (size_t i = 0; i < s; i++)
		{
			x0 = x;
			v0 = v;
			x = x + h;
			k1 = syst1.f(x0, v0);
			k2 = syst1.f(x0 + h / 2, v0 + h / 2 * k1);
			k3 = syst1.f(x0 + h / 2, v0 + h / 2 * k2);
			k4 = syst1.f(x0 + h, v0 + h * k3);
			v = v0 + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
			temp.push_back(i);
			temp.push_back(x0);
			temp.push_back(v0);
			temp.push_back(h);
			solve.push_back(temp);
			temp.clear();
		}
	};
	friend std::ostream& operator<<(std::ostream& out, const RK4& sol)
	{
		out << "n " << "xi " << "vi " << "h\n";
		for (size_t i = 0; i < sol.s; i++)
		{
			out << sol.solve[i][0] << " " << sol.solve[i][1] << " " << sol.solve[i][2] << " " << sol.solve[i][3] << "\n";
		}
		return out;
	};
	~RK4() {};
};