#include "RK4.h"

int main()
{
	RK4 table(1, 0.0, 1.0, 0.01, 10000, 0.001, 0.00001, 1.0);
	table.main_solve();
	std::cout << table;
	return 0;
}