#pragma once
#include "solver.h"
#include <iomanip>

class RK4 : public Solver
{
public:
	RK4(int p, double _x0, double _v0, double _h, int _n_max, double _eps_e, double _eps_b, double _b) : 
		Solver(p, _x0, _v0, _h, _n_max, _eps_e, _eps_b, _b) {};
	void main_solve();													// решение с половинным шагом
	friend std::ostream& operator<<(std::ostream& out, const RK4& sol);
	~RK4() {};
private:
	double k1, k2, k3, k4;												// коэффиценты метода
	virtual void step(double &x, double &v, double &h);	// шаг метода

	// служебные методы и переменные
	double vhalf, xhalf, vprev, x_c, v_c, h_c, v_accur, C;	
	int  c1_glob, c2_glob;
	double max_h, min_h, max_err;
	std::vector<double> temp;
	void double_step();
	void insert_step();
	void control_e();
};

void RK4::main_solve()
{
	C = v0 / pow(2.7218721872, 3 * x0);
	v_accur = v0;
	for (n_counter; n_counter < n_max; n_counter++)
	{
		vprev = v;
		if (n_counter == 0)
			v_accur = C * pow(2.7218721872, 3 * x);
	    insert_step();
		double_step();
		c1 = 0; c2 = 0;
		control_e();
		if (x + eps_b > b)
		{
			n_counter++;
			insert_step();
			for (size_t i = 0; i < n_counter - 1; i++)
			{
				solve[i][7] = solve[i + 1][7];
				solve[i][8] = solve[i + 1][8];
			}
			break;
		}
	}
	for (size_t i = 0; i < n_counter - 1; i++)
	{
		solve[i][7] = solve[i + 1][7];
		solve[i][8] = solve[i + 1][8];
	}
}

void RK4::step(double &x, double &v, double &h)
{
	k1 = syst1.f(x, v);
	k2 = syst1.f(x + h / 2, v + h / 2 * k1);
	k3 = syst1.f(x + h / 2, v + h / 2 * k2);
	k4 = syst1.f(x + h, v + h * k3);
	x = x + h;
	v = v + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
}

void RK4::double_step()
{
	h_c = h;
	h = h / 2;
	vhalf = v;
	v_c = v;
	xhalf = x;
	x_c = x;
	step(x_c, v_c, h_c);
	v_accur = C * pow(2.7218721872, 3 * x_c);
	step(xhalf, vhalf, h);
	step(xhalf, vhalf, h);
	h = h * 2;
	v = vhalf;
	x = xhalf;
	v_v2 = v_c - v;
}

void RK4::insert_step()
{
	temp.push_back(n_counter);
	temp.push_back(x);
	temp.push_back(v_c);
	temp.push_back(v);
	temp.push_back(v_v2);
	temp.push_back(abs(v_v2 / (pow(2, 4) - 1)));
	temp.push_back(h);
	temp.push_back(c1);
	temp.push_back(c2);
	temp.push_back(v_accur);
	temp.push_back(abs(v_accur - v));
	solve.push_back(temp);
	temp.clear();
}

void RK4::control_e()
{
	if (abs(v_v2 / (pow(2, 4) - 1)) > eps_e)
	{
		v = vprev;
		x = x - h;
		h = h / 2;
		double_step();
		c1++;
		c1_glob++;
	}
	if (abs(v_v2 / (pow(2, 4) - 1)) < eps_e / pow(2, 5))
	{
		h = h * 2;
		c2++;
		c2_glob++;
	}
}

std::ostream& operator<<(std::ostream& out, const RK4& sol)
{
	out << "i\t" << "xi\t"<< "vi\t"<< "v2i\t" << "vi-vi2\t" << "OLP\t" << "hi\t" << "c1\t" << "c2\n";
	for (size_t i = 0; i < sol.n_counter; i++)
	{
		out <<  sol.solve[i][0] << "\t";
		out << sol.solve[i][1] << "\t" << sol.solve[i][2]
			  << "\t" << sol.solve[i][3] << "\t" << sol.solve[i][4] << "\t" << sol.solve[i][5] << "\t" << sol.solve[i][6] 
			  << "\t" << sol.solve[i][7] << "\t" << sol.solve[i][8] << "\t" << sol.solve[i][9] << "\t" << sol.solve[i][10] << "\n";
	}
	return out;
};