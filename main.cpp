#include<iostream>
#include<fstream>

const double x_0 = 0.1;
const double y_0 = 0.1;
const long int max_time = 1000000;
const double dt = 0.01;

double x = x_0;
double y = y_0;
double y_prev = y_0;
long int t = 0;
long int check = 0;

double f1(const double x, const double y)
{
	return x - x * y - x * x;
}

double f2(const double x, const double y)
{
	return y * x - y * y;
}

void CountDxAndDy(const double x, const double y, const double dt, double &dx, double &dy)
{
	double k1_x, k2_x, k3_x, k4_x;
	double k1_y, k2_y, k3_y, k4_y;
	k1_x = f1(x, y);
	k1_y = f2(x, y);
	k2_x = f1(x + dt / 2 * k1_x, y + dt / 2 * k1_y);
	k2_y = f2(x + dt / 2 * k1_x, y + dt / 2 * k1_y);
	k3_x = f1(x + dt / 2 * k2_x, y + dt / 2 * k2_y);
	k3_y = f2(x + dt / 2 * k2_x, y + dt / 2 * k2_y);
	k4_x = f1(x + dt * k3_x, y + dt * k3_y);
	k4_y = f2(x + dt * k3_x, y + dt * k3_y);
	dx = dt / 6 * (k1_x + 2 * k2_x + 2 * k3_x + k4_x);
	dy = dt / 6 * (k1_y + 2 * k2_y + 2 * k3_y + k4_y);
}

int main()
{
	double dx, dy;
	std::ofstream f1out;
	std::ofstream f2out;
	f1out.open("C:\\Users\\stas2\\Desktop\\result\\result0.txt");//Enter your path
	f2out.open("C:\\Users\\stas2\\Desktop\\result\\result1.txt");//Enter your path
	while(t < max_time)
	{
		CountDxAndDy(x, y, dt, dx, dy);
		x += dx;
		y += dy;
		if (x - dx < x_0 && x >= x_0)
			check++;
		if (!check % 2)
		{
			f2out << y_prev << '\t' << y <<  '\t' << check/2 << std::endl;
			y_prev = y;
		}
		f1out << x << '\t' << y << '\t' << t << std::endl;
		t++;
	}
}