#include<iostream>
#include<fstream>
#include<algorithm>
#include<vector>
#include<complex>
#include"fparser.hh"
#include"FindComplexRoots.h"

typedef std::vector<std::vector<long double>> matrix_double;
typedef std::complex<long double> complex;

const double x_0 = 0.1;
const double y_0 = 0.1;
const long int max_time = 100;
const double dt = 0.0001;
const double dh = 0; //Шаг векторного поля
const double E1 = 0.1;
const double E2 = 0.1;

std::vector<std::string> functions = { "-x1*x1","-150*x2*x2"};


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

template<typename T>
T Determinant(const std::vector<std::vector<T>>& curr_matrix)
{
	if (curr_matrix.size() == 1)
		return curr_matrix[0][0];
	T determinant = 0;
	int k = 1;
	for (size_t i = 0; i < curr_matrix.size(); i++)
	{
		if (curr_matrix[0][i] != 0)
		{
			std::vector<std::vector<T>> new_matrix = curr_matrix;
			new_matrix.erase(new_matrix.begin());
			for (size_t j = 0; j < new_matrix.size(); j++)
			{
				new_matrix[j].erase(new_matrix[j].begin() + i);
			}
			determinant += k * curr_matrix[0][i] * Determinant(new_matrix);
		}
		k *= -1;
	}
	return determinant;
}

matrix_double AdjugateT(const matrix_double& curr_matrix)
{
	matrix_double adjugateT{ curr_matrix.size() , std::vector<long double>(curr_matrix.size()) };
	int sign;
	for (size_t i = 0; i < curr_matrix.size(); i++)
	{
		for (size_t j = 0; j < curr_matrix.size(); j++)
		{
			matrix_double new_matrix = curr_matrix;
			new_matrix.erase(new_matrix.begin() + i);
			for (size_t k = 0; k < new_matrix.size(); k++)
			{
				new_matrix[k].erase(new_matrix[k].begin() + j);
			}
			sign = 1 - 2 * ((i + j) % 2);
			adjugateT[j][i] = sign * Determinant(new_matrix);
		}
	}
	return adjugateT;
}

matrix_double MatrixInverse(const matrix_double& curr_matrix)
{
	matrix_double matrix_inverse{ curr_matrix.size(),  std::vector<long double>(curr_matrix.size()) };
	long double det = Determinant(curr_matrix);
	int sign;
	for (size_t i = 0; i < curr_matrix.size(); i++)
	{
		for (size_t j = 0; j < curr_matrix.size(); j++)
		{
			matrix_double new_matrix = curr_matrix;
			new_matrix.erase(new_matrix.begin() + i);
			for (size_t k = 0; k < new_matrix.size(); k++)
			{
				new_matrix[k].erase(new_matrix[k].begin() + j);
			}
			sign = 1 - 2 * ((i + j) % 2);//  1 or -1
			matrix_inverse[j][i] = sign * Determinant(new_matrix) / det;
		}
	}
	return matrix_inverse;
}

complex CharacteristicEquation(const complex& z, const matrix_double& curr_matrix)
{
	std::vector<std::vector<complex>> new_matrix{ curr_matrix.size(), std::vector<complex>(curr_matrix.size()) };
	for (size_t i = 0; i < new_matrix.size(); i++)
		new_matrix[i][i] = curr_matrix[i][i]-z;
	return Determinant(new_matrix);
}

bool IsHard(const std::vector<long double>& coor, const long double eps = 0.0001)
{
	matrix_double jacobian_matrix = GetJacobianMatrix(coor, functions);
	std::vector<complex> eigenvalues;
	FindCR::FindRoots(sqrt(GetNorm(jacobian_matrix)), eps, 0, 0, eigenvalues, [jacobian_matrix](complex z) { return CharacteristicEquation(z, jacobian_matrix); });
	std::vector<long double> abs_eigenvalues;
	for (auto el : eigenvalues)
	{
		if (std::abs(el) < eps)
			return true;
		abs_eigenvalues.push_back(std::abs(el));
	}
	if (*std::max_element(abs_eigenvalues.begin(), abs_eigenvalues.end()) / *std::min_element(abs_eigenvalues.begin(), abs_eigenvalues.end()) > 100)
		return true;
	return false;
}

size_t KroneckerSymbol(size_t i, size_t j)
{
	return i == j ? 1 : 0;
}

template<typename T>
T Dot(const std::vector<T>& vec1, const std::vector<T>& vec2)
{
	T dot=0;
	if (vec1.size() != vec2.size())
		throw std::exception("Vectors aren't equal");
	for (size_t i = 0; i < vec1.size(); i++)
	{
		dot += vec1[i] * vec2[i];
	}
	return dot;
}

std::vector<long double> f(const std::vector<long double>& coor)
{
	std::vector<long double> result;
	FunctionParser fp;
	for (size_t i = 0; i < functions.size(); i++)
	{
		fp.Parse(functions[i], GetVariables(functions));
		result.push_back(fp.Eval((double*)coor.data()));
	}
	return result;
}

std::vector<long double> CountNextCoor(const std::vector<long double>& coor, const std::vector<std::string>& functions, const double dt)
{
	std::vector<long double> result;
	if (IsHard(coor))
	{
		const matrix_double jacobian_matrix = GetJacobianMatrix(coor, functions);
		matrix_double new_matrix{ jacobian_matrix.size(), std::vector<long double>(jacobian_matrix.size()) };
		for (size_t i = 0; i < new_matrix.size(); i++)
		{
			for (size_t j = 0; j < new_matrix.size(); j++)
			{
				new_matrix[i][j] = KroneckerSymbol(i, j) - jacobian_matrix[i][j] * dt;
			}
		}
		new_matrix = MatrixInverse(new_matrix);
		std::vector<long double> functions_coor_dt;
		FunctionParser fp;
		for (size_t i = 0; i < functions.size(); i++)
		{
			fp.Parse(functions[i], GetVariables(functions));
			functions_coor_dt.push_back(fp.Eval((double*)coor.data())*dt);
		}
		for (size_t i = 0; i < coor.size(); i++)
			result.push_back(coor[i]+Dot(new_matrix[i], functions_coor_dt));
		return result;
	}
	std::vector<long double> k1, k2, k3, k4, current_coor;
	current_coor = coor;
	k1 = f(current_coor);
	for (size_t i = 0; i < current_coor.size(); i++)
	{
		current_coor[i] = coor[i] + k1[i] * dt / 2;
	}
	k2 = f(current_coor);
	for (size_t i = 0; i < current_coor.size(); i++)
	{
		current_coor[i] = coor[i] + k2[i] * dt / 2;
	}
	k3 = f(current_coor);
	for (size_t i = 0; i < current_coor.size(); i++)
	{
		current_coor[i] = coor[i] + k3[i] * dt;
	}
	k4 = f(current_coor);
	for (size_t i = 0; i < coor.size(); i++)
	{
		result.push_back(coor[i] + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
	}
	return result;
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
	var[1] = 0.001;
	matrix_double jacobian_matrix = GetJacobianMatrix(var, functions);
	matrix_double inverse_matrix = MatrixInverse(jacobian_matrix);
	std::vector<complex> complex_vec;
	FindCR::FindRoots(sqrt(GetNorm(jacobian_matrix)), 0.0001, 0, 0, complex_vec, [jacobian_matrix](complex z) { return CharacteristicEquation(z, jacobian_matrix); });
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
	std::cout << "Inverse:\n";
	for (size_t i = 0; i < inverse_matrix.size(); i++)
	{
		for (size_t j = 0; j < inverse_matrix[i].size(); j++)
		{
			std::cout << inverse_matrix[i][j] << '\t';
		}
		std::cout << '\n';
	}
	std::cout << "Hard: " << bool(IsHard(var)) << '\n';
	std::cout << "Dot: " << Dot(jacobian_matrix[0], jacobian_matrix[1]) << '\n';

	std::vector<long double> result = CountNextCoor(var, functions, dt);

	for (auto el : result)
		std::cout << el << '\t';

	std::cout << "\nDeterminant: " << Determinant(jacobian_matrix) << '\n';

	while(t < max_time)
	{
		/*if (t % 10 == 0)
		{
			f1out << x << '\t' << y << '\t' << t / 10 << std::endl;
		}
		if (!(t % 10000))
		{
			if (h >= 0.)
				x = x_0 + h,
				y = y_0 + h,
				h += dh;
		}*/
		//CountDxAndDy(x, y, dt, dx, dy, IsHard(x,y));
		//x += dx;
		//y += dy;
		std::cout << "Time:" << t << ", (x;y)=(" << var[0]  << ';' << var[1] << "), Hard: " <<  IsHard(var) << std::endl;
		var = CountNextCoor(var, functions, dt);
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