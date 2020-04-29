#include<iostream>
#include<fstream>
#include<algorithm>

const double x_0 = 0.1;
const double y_0 = 0.1;
const long int max_time = 10000;
const double dt = 0.01;
const double dh = 0; //Шаг векторного поля
const double E1 = 0.1;
const double E2 = 0.1;

double x = x_0;
double y = y_0;
//double y_prev = y_0;
double h = 0;
long int t = 0;
//long int check = 0;

double F1(const double x, const double y)
{
	return -2*x;
}

double DifF1x(const double x, const double y, const double eps = 0.0001)
{
	return (F1(x + eps, y) - F1(x - eps, y)) / (2 * eps);
}

double DifF1y(const double x, const double y, const double eps = 0.0001)
{
	return (F1(x, y + eps) - F1(x, y - eps)) / (2 * eps);
}

double F2(const double x, const double y)
{
	return -y;
}

double DifF2x(const double x, const double y, const double eps = 0.0001)
{
	return (F2(x + eps, y) - F2(x - eps, y)) / (2 * eps);
}

double DifF2y(const double x, const double y, const double eps = 0.0001)
{
	return (F2(x, y + eps) - F2(x, y - eps)) / (2 * eps);
}

bool IsHard(const double x, const double y)
{
	//Matrix A:
	//||a b||
	//||c d||
	double a = DifF1x(x, y);
	double b = DifF1y(x, y);
	double c = DifF2x(x, y);
	double d = DifF2y(x, y);
	//det(A-lE)=0:  //l-lambda
	//Discriminant:
	double D = (a + d) * (a + d) - 4 * (a * d - c * b);
	double l1, l2;
	if (D > 0)
		l1 = ((a + d) + sqrt(D))/2,
		l2 = ((a + d) - sqrt(D))/2;
	else
		l1 = (a + d) / 2,
		l2 = (a + d) / 2;
	if (l1 < 0 && l2 < 0 && std::min(l1, l2) / std::max(l1, l2) > 100)
		return true;
	return false;
}

void CountDxAndDy(const double x, const double y, const double dt, double &dx, double &dy, bool implicit = false)
{
	if (implicit)
	{
		dx = x / (1 + dt * 2) - x;
		dy = y / (1 + dt) - y;
		return;
	}

	double k1_x, k2_x, k3_x, k4_x;
	double k1_y, k2_y, k3_y, k4_y;
	k1_x = F1(x, y);
	k1_y = F2(x, y);
	k2_x = F1(x + dt / 2 * k1_x, y + dt / 2 * k1_y);
	k2_y = F2(x + dt / 2 * k1_x, y + dt / 2 * k1_y);
	k3_x = F1(x + dt / 2 * k2_x, y + dt / 2 * k2_y);
	k3_y = F2(x + dt / 2 * k2_x, y + dt / 2 * k2_y);
	k4_x = F1(x + dt * k3_x, y + dt * k3_y);
	k4_y = F2(x + dt * k3_x, y + dt * k3_y);
	dx = dt / 6 * (k1_x + 2 * k2_x + 2 * k3_x + k4_x);
	dy = dt / 6 * (k1_y + 2 * k2_y + 2 * k3_y + k4_y);
}

int main()
{
	double dx, dy;
	std::ofstream f1out;
	std::ofstream f2out;
	f1out.open("C:\\Users\\stas2\\Desktop\\result\\result0.txt");//Введите свой путь
	f2out.open("C:\\Users\\stas2\\Desktop\\result\\result1.txt");//Enter your path
	while(t < max_time)
	{
		if (t % 10 == 0)
		{
			f1out << x << '\t' << y << '\t' << t / 10 << std::endl;
		}
		if (!(t % 10000))
		{
			if (h >= 0.)
				x = x_0 + h,
				y = y_0 + h,
				h += dh;
		}
		CountDxAndDy(x, y, dt, dx, dy, IsHard(x,y));
		x += dx;
		y += dy;
		std::cout << "Time:" << t << ", (x;y)=(" << x  << ';' << y << "), Hard: " <<  IsHard(x,y) << std::endl;
		/*if (x - dx < x_0 && x >= x_0)
		{
			check++;
			if (!(check % 2))
			{
				f2out << y_prev << '\t' << y << '\t' << check / 2 << std::endl;
				y_prev = y;
			}
		}*/
		t++;
	}
}