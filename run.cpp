#include "RK4.h"

int main()
{
	RK4 table(0.0, 1.0, 0.01, 10, 1);
	std::cout << table;
	return 0;
}