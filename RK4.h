#pragma once
#include "solver.h"

class RK4 : public Solver
{
public:
	RK4(int p, double _x0, double _v0, double _h, int _n_max, double _eps_e, double _eps_b, double _b) : 
		Solver(p, _x0, _v0, _h, _n_max, _eps_e, _eps_b, _b) {};
	void first_solve();														// решение обычное
	void second_solve();													// решение с половинным шагом
	void third_solve();
	friend std::ostream& operator<<(std::ostream& out, const RK4& sol);
	~RK4() {};
private:
	double k1, k2, k3, k4;												// коэффиценты метода
	virtual void step(double &x, double &v, double &h);	// шаг метода

	// служебные методы и переменные
	double vhalf, xhalf, vprev, x_c, v_c, h_c;	
	std::vector<double> temp;
	void double_step();
	void insert_step();
	void control_e();
};

void RK4::first_solve()
{
	if (flag_cl == 1)
	{
		clear();
	}
	for (n_counter; n_counter < n_max; n_counter++)
	{
		if (v < 40 && v > -40)//v < Ylim
			vprev = v;
		else
			v = vprev;
		insert_step();
		step(x, v, h);
		flag_cl = 0;
		if (x + eps_b > b)
		{
			n_counter++;
			insert_step();
			break;
		}
	}
}

void RK4::second_solve()
{
	if (flag_cl == 1)
	{
		clear();
	}
	for (n_counter; n_counter < n_max; n_counter++)
	{
		if (v < 40 && v > -40)//v < Ylim
			vprev = v;
		else
			v = vprev;
		insert_step();
		double_step();
		control_e();
		flag_cl = 0;
		if (x + eps_b > b)
		{
			n_counter++;
			insert_step();
			break;
		}
	}
}

void RK4::third_solve()
{

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
	temp.push_back(v);
	temp.push_back(h);
	temp.push_back(v_v2);
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
	}
	if (abs(v_v2 / 15) < eps_e / pow(2, 5))
	{
		h = h * 2;
		c2++;
	}
}

std::ostream& operator<<(std::ostream& out, const RK4& sol)
{
	out << "n " << "\t" << "xi " << "\t" << "vi " << "\t" << "h\n";
	for (size_t i = 0; i < sol.n_counter; i++)
	{
		out << sol.solve[i][0] << "\t" << sol.solve[i][1] << "\t" << sol.solve[i][2] 
			  << "\t" << sol.solve[i][3] << "\t" << sol.solve[i][4] << "\n";
	}
	return out;
};