#include<iostream>
#include<fstream>
#include<algorithm>
#include<vector>
#include<complex>
#include"matrix.h"
#include"fparser.hh"
#include"FindComplexRoots.h"

typedef std::vector<std::vector<long double>> matrix_double;
typedef std::complex<long double> complex;

const double x_0 = 0.1;
const double y_0 = 0.1;
const long int max_time = 0;
const double dt = 0.0001;
const double dh = 0; //Шаг векторного поля
const double E1 = 0.1;
const double E2 = 0.1;

std::vector<std::string> functions = { "-x1","-100*x2"};


double x = x_0;
double y = y_0;
//double y_prev = y_0;
double h = 0;
long int t = 0;
//long int check = 0;

std::string GetVariables(const std::vector<std::string>& functions)
{
	std::string variables = "";
	for (size_t i = 1; i <= functions.size(); i++)
		variables += 'x' + std::to_string(i) + ',';
	variables.pop_back();
	return variables;
}

matrix_double GetJacobianMatrix(const std::vector<long double>& coor, const std::vector<std::string>& functions, const long double eps = 0.0001)
{
	FunctionParser fp;
	matrix_double return_matrix{ functions.size(), std::vector<long double>(functions.size()) };
	std::vector<long double> copy_coor = coor;
	for (size_t i = 0; i < functions.size(); i++)
	{
		for (size_t j = 0; j < functions.size(); j++)
		{
			fp.Parse(functions[i], GetVariables(functions));
			copy_coor[j] += eps;
			return_matrix[i][j] = fp.Eval((double*)copy_coor.data());
			copy_coor[j] -= 2 * eps;
			return_matrix[i][j] = (return_matrix[i][j]-fp.Eval((double*)copy_coor.data()))/(2*eps);
			copy_coor[j] += eps;
		}
	}
	return return_matrix;
}

long double GetNorm(const matrix_double& curr_matrix)
{
	long double norm = 0;
	for (auto vec : curr_matrix)
		for (auto el : vec)
			norm += el * el;
	return norm;
}

complex Determinant(const std::vector<std::vector<complex>>& curr_matrix)
{
	if (curr_matrix.size() == 1)
		return curr_matrix[0][0];
	complex determinant = 0;
	int k = 1;
	for (size_t i = 0; i < curr_matrix.size(); i++)
	{
		if (std::abs(curr_matrix[0][i]) != 0)
		{
			std::vector<std::vector<complex>> new_matrix = curr_matrix;
			new_matrix.erase(new_matrix.begin());
			for (size_t j = 0; j < new_matrix.size(); j++)
			{
				new_matrix[j].erase(new_matrix[j].begin() + i);
			}
			determinant += complex(k,0) * curr_matrix[0][i] * Determinant(new_matrix);
		}
		k *= -1;
	}
	return determinant;
}

complex CharacteristicEquation(const complex& z, const matrix_double& curr_matrix)
{
	std::vector<std::vector<complex>> new_matrix{ curr_matrix.size(), std::vector<complex>(curr_matrix.size()) };
	for (size_t i = 0; i < new_matrix.size(); i++)
		new_matrix[i][i] = curr_matrix[i][i]-z;
	return Determinant(new_matrix);
}

bool IsHard(const std::vector<long double>& coor, const long double eps = 0.1)
{
	matrix_double jacobian_matrix = GetJacobianMatrix(coor, functions);
	std::vector<complex> eigenvalues;
	FindCR::FindRoots(sqrt(GetNorm(jacobian_matrix)), eps, 0, 0, eigenvalues, [jacobian_matrix](complex z) { return CharacteristicEquation(z, jacobian_matrix); });
	std::vector<long double> abs_eigenvalues;
	for (auto el : eigenvalues)
		if (std::abs(el) == 0)
			abs_eigenvalues.push_back(eps);
		else
			abs_eigenvalues.push_back(std::abs(el));
	if (*std::max_element(abs_eigenvalues.begin(), abs_eigenvalues.end()) / *std::min_element(abs_eigenvalues.begin(), abs_eigenvalues.end()) > 100)
		return true;
	return false;
}

std::vector<long double> CountNextCoor(const std::vector<long double>& coor, const std::vector<std::string>& functions, const double dt, bool implicit = false)
{
	if (implicit)
	{
		double *val = new double[coor.size()*coor.size()];
		matrix E;
		E.set(coor.size(), coor.size(), val);
		E.power(0);
		matrix_double jacob_mat = GetJacobianMatrix(coor, functions);
		matrix A;
		for (size_t i = 0; i < jacob_mat.size(); i++)
		{
			for (size_t j = 0; j < jacob_mat.size(); j++)
			{
				val[j + i * jacob_mat.size()] = jacob_mat[i][j];
			}
		}
		A.set(jacob_mat.size(), jacob_mat.size(), val);
		matrix B = E - (A * dt);
		matrix obB = B.power(-1);
		matrix X0;
		X0.set(coor.size(), 1, (double*)coor.data());
		matrix F0;
		FunctionParser fp;
		for (size_t i = 0; i < jacob_mat.size(); i++)
		{
			fp.Parse(functions[i], GetVariables(functions));
			val[i] = fp.Eval((double*)coor.data());
		}
		F0.set(coor.size(), 1, val);
		matrix X = X0 + obB * F0 * dt;
		std::vector<long double> result;
		for (size_t i = 0; i < coor.size(); i++)
			result.push_back(X.get_val(i + 1, 1));
		return result;
	}

	/*double k1_x, k2_x, k3_x, k4_x;
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
	dy = dt / 6 * (k1_y + 2 * k2_y + 2 * k3_y + k4_y);*/
}

int main()
{
	//double dx, dy;
	std::ofstream f1out;
	std::ofstream f2out;
	f1out.open("C:\\Users\\stas2\\Desktop\\result\\result0.txt");//Введите свой путь
	f2out.open("C:\\Users\\stas2\\Desktop\\result\\result1.txt");//Enter your path
	std::vector<long double> var(functions.size(), 0.1);
	matrix_double jacobian_matrix = GetJacobianMatrix(var, functions);
	std::vector<complex> complex_vec;
	FindCR::FindRoots(sqrt(GetNorm(jacobian_matrix)), 0.1, 0, 0, complex_vec, [jacobian_matrix](complex z) { return CharacteristicEquation(z, jacobian_matrix); });
	for (auto el : complex_vec)
		std::cout << el << '\t';
	std::cout << '\n';
	for (size_t i = 0; i < jacobian_matrix.size(); i++)
	{
		for (size_t j = 0; j < jacobian_matrix[i].size(); j++)
		{
			std::cout << jacobian_matrix[i][j] << '\t';
		}
		std::cout << '\n';
	}
	std::cout << "Hard: " << bool(IsHard(var)) << '\n';

	std::vector<long double> result = CountNextCoor(var, functions, dt, true);

	for (auto el : result)
		std::cout << el << '\t';

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
		//CountDxAndDy(x, y, dt, dx, dy, IsHard(x,y));
		//x += dx;
		//y += dy;
		//std::cout << "Time:" << t << ", (x;y)=(" << x  << ';' << y << "), Hard: " <<  IsHard(x,y) << std::endl;
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