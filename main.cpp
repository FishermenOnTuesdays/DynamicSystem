﻿#include<iostream>
#include<fstream>
#include<algorithm>
#include<vector>
#include<complex>
#include<cmath>
#include "fparser.hh"
#include "Eigen/Dense"

typedef std::vector<std::vector<long double>> matrix_double;
typedef std::complex<long double> complex;

const double x_0 = 0.1;
const double y_0 = 0.1;
long int max_time = 10;
double dt = 0.01;
const double dh = 0; //Шаг векторного поля
const double E1 = 0.1;
const double E2 = 0.1;
const std::vector<std::string> Lorenz_attractor = { "10*(x2-x1)", "x1*(28-x3)-x2", "x1*x2-(8/3)*x3" };
//const std::vector<std::string>  = {"10*(x2-x1)", "28*x1-x2-x1*x3-x1*x3"};
const std::vector<std::string> fixed_point = { "-x3", "-x1", "-x2" };

size_t N;
std::vector<std::string> functions;
std::vector<long double> lyapunov_exponents;
bool exact;
bool is_attractor = false;
Eigen::MatrixXd evolution_matrix;


double x = x_0;
double y = y_0;
double h = 0;
long int t = 0;

std::vector<long double> Sub(const std::vector<long double>& vec1, const std::vector<long double>& vec2)
{
	if (vec1.size() != vec2.size())
		throw std::exception("Vectors aren't equal");
	std::vector<long double> sub_vec(vec1.size());
	for (size_t i = 0; i < vec1.size(); i++)
		sub_vec[i] = vec1[i] - vec2[i];
	return sub_vec;
}

std::vector<long double> Add(const std::vector<long double>& vec1, const std::vector<long double>& vec2)
{
	if (vec1.size() != vec2.size())
		throw std::exception("Vectors aren't equal");
	std::vector<long double> add_vec(vec1.size());
	for (size_t i = 0; i < vec1.size(); i++)
		add_vec[i] = vec1[i] + vec2[i];
	return add_vec;
}

std::vector<long double> Mult(const std::vector<long double>& vec, long double scalar)
{
	std::vector<long double> mult_vec(vec.size());
	for (size_t i = 0; i < vec.size(); i++)
		mult_vec[i] = vec[i] * scalar;
	return mult_vec;
}

long double LengthVector(const std::vector<long double>& vec)
{
	long double length = 0;
	for (size_t i = 0; i < vec.size(); i++)
		length += vec[i] * vec[i];
	return sqrtl(length);
}

bool Equal(const std::vector<long double>& vec1, const std::vector<long double>& vec2, const long double eps = 0.000001)
{
	if (vec1.size() != vec2.size())
		return false;
	for (size_t i = 0; i < vec1.size(); i++)
		if (fabsl(vec1[i] - vec2[i]) > eps)
			return false;
	return true;
}

std::string GetVariables(const std::vector<std::string>& functions)
{
	std::string variables = "";
	for (size_t i = 1; i <= functions.size(); i++)
		variables += 'x' + std::to_string(i) + ',';
	variables.pop_back();
	return variables;
}

matrix_double GetJacobianMatrix(const std::vector<long double>& coor, const std::vector<std::string>& functions, const long double eps = 0.000001)
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

bool IsHard(const std::vector<long double>& coor, const long double eps = 0.000001)
{
	matrix_double jacobian_matrix = GetJacobianMatrix(coor, functions);
	Eigen::MatrixXd jacobian_matrix_new(jacobian_matrix.size(), jacobian_matrix.size());
	for (size_t i = 0; i < jacobian_matrix.size(); i++)
		for (size_t j = 0; j < jacobian_matrix.size(); j++)
			jacobian_matrix_new(i, j) = jacobian_matrix[i][j];
	//if(is_attractor)
		//evolution_matrix = jacobian_matrix_new * evolution_matrix;
	auto eigenvalues = jacobian_matrix_new.eigenvalues();
	//const double lyapunov_eps = 0.01;
	//if (!lyapunov_exponents.empty())
	//{
	//	exact = true;
	//	for (size_t i = 0; i < eigenvalues.size(); i++)
	//		if (fabs(eigenvalues[i].real() - lyapunov_exponents[i]) > lyapunov_eps)
	//		{
	//			exact = false;
	//			break;
	//		}
	//}
	//else
	//	exact = false;
	//lyapunov_exponents.clear();
	//for (size_t i = 0; i < eigenvalues.size(); i++)
	//{
	//	//std::cout << eigenvalues[i] << " ";
	//	lyapunov_exponents.push_back((evolution_matrix.transpose()*evolution_matrix).eigenvalues()[i].real());
	//}
	////std::cout << std::endl;
	std::vector<long double> abs_eigenvalues;
	for (size_t i = 0; i < eigenvalues.size(); i++)
	{
		if (std::abs(eigenvalues[i]) < eps)
			return true;
		abs_eigenvalues.push_back(std::abs(eigenvalues[i]));
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

void CountLogSum(std::vector<long double>& log_sum, const Eigen::MatrixXd& matrix_deviation)
{
	Eigen::HouseholderQR<Eigen::MatrixXd> qr(matrix_deviation);
	Eigen::MatrixXd orthonormal_matrix_deviation = qr.householderQ();
	//std::cout << "Deviation:\n" << matrix_deviation << std::endl;
	//std::cout << "Orthonormal:\n" << orthonormal_matrix_deviation << std::endl;
	for (size_t i = 0; i < N; i++)
	{
		Eigen::VectorXd current_row = matrix_deviation.row(i);
		//std::cout << "Row " << i << " : " << matrix_deviation.row(i) << std::endl;
		for (size_t j = 0; j < i; j++)
		{
			//std::cout << "Dot deviation_row " << i << " and orthonormal_row " << j << " : " << matrix_deviation.row(i).dot(orthonormal_matrix_deviation.row(j)) << std::endl;
			current_row -= matrix_deviation.row(i).dot(orthonormal_matrix_deviation.row(j)) * orthonormal_matrix_deviation.row(j);
			//std::cout << "Current row: " << current_row << std::endl;
		}
		//std::cout << "Current norm: " << current_row.norm() << std::endl;
		log_sum[i] += std::logl(current_row.norm());
	}
}

void CountNextCoorDeviation(matrix_double& var_deviation, const std::vector<std::string>& functions, const double dt)
{
	for (size_t i = 0; i < N; i++)
		var_deviation[i] = CountNextCoor(var_deviation[i], functions, dt);
	/*matrix_double var_deviation{ N, std::vector<long double>(N)};
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
			var_deviation[i][j] = var[j] + matrix_deviation(i,j);
		var_deviation[i] = CountNextCoor(var_deviation[i], functions, dt);
	}
	var = CountNextCoor(var, functions, dt);
	for (size_t i = 0; i < N; i++)
		for (size_t j = 0; j < N; j++)
			matrix_deviation(i, j) = var_deviation[i][j] - var[j];*/
}

void CountDeviation(Eigen::MatrixXd& deviation, const matrix_double& var_deviation, const std::vector<long double>& var)
{
	for (size_t i = 0; i < N; i++)
		for (size_t j = 0; j < N; j++)
			deviation(i, j) = var_deviation[i][j] - var[j];
}

int main()
{
	std::ofstream f1out, f2out;
	f1out.open("../wwwroot/output/result.csv");//Введите свой путь
	f2out.open("../wwwroot/output/laypunov.csv");//Введите свой путь
	std::vector<long double> var;
	#ifdef _DEBUG
		N = 3;
	#else
		std::cin >> N;
	#endif
	for (size_t i = 1; i <= N; i++)
		f1out << 'x' + std::to_string(i) + ',';
	for (size_t i = 1; i <= N; i++)
		f1out << 'l' + std::to_string(i) + ',';
	f1out << "t\n";
	functions.resize(N);
	#ifdef _DEBUG
		functions = fixed_point;
	#else
		for (size_t i = 0; i < N; i++)
			std::cin >> functions[i];
	#endif
	var.resize(N);
	#ifdef _DEBUG
		for (size_t i = 0; i < N; i++)
			var[i] = 0.1;
	#else
		for (size_t i = 0; i < N; i++)
			std::cin >> var[i];
	#endif
	std::vector<std::vector<long double>> dots;
	std::cin >> max_time;
	std::cin >> dt;
	while(t < max_time)
	{
		dots.push_back(var);
		for (size_t i = 0; i < var.size(); i++)
		{
			f1out << std::to_string(var[i])+',';
		}
		var = CountNextCoor(var, functions, dt);
		f1out << std::to_string(t * dt) + '\n';
		t++;
	}

	/*std::cout << "GOOD\n";
	for (size_t i = 0; i < N; i++)
		std::cout << var[i] << ", ";
	std::cout << "\nGOOD\n";*/

	long int M = 100000;
	long int T = 100;
	long double eps_deviation = 1;
	std::vector<long double> deviation_x(var.size(), 0);
	deviation_x[0] = eps_deviation;
	std::vector<long double> deviation_y(var.size(), 0);
	deviation_y[1] = eps_deviation;
	std::vector<long double> var_deviation_x = Add(var, deviation_x);
	std::vector<long double> var_deviation_y = Add(var, deviation_y);
	long double log_sum_x = 0, log_sum_y = 0;
	Eigen::MatrixXd deviation = Eigen::MatrixXd::Identity(N,N);
	matrix_double var_deviation{ N, std::vector<long double>(N) };
	for (size_t i = 0; i < N; i++)
		for (size_t j = 0; j < N; j++)
			var_deviation[i][j] = var[j] + deviation(i, j);
	std::vector<long double> log_sum(N, 0);
	
	for (long unsigned int i = 0; i < M; i++)
	{
		for (long unsigned int j = 0; j < T; j++)
		{
			CountNextCoorDeviation(var_deviation, functions, dt);
			var = CountNextCoor(var, functions, dt);
		}
		CountDeviation(deviation, var_deviation, var);
		CountLogSum(log_sum, deviation);
		Eigen::HouseholderQR<Eigen::MatrixXd> qr(deviation);
		deviation = qr.householderQ();
		for (size_t j = 0; j < N; j++)
			for (size_t k = 0; k < N; k++)
				var_deviation[j][k] = var[k] + deviation(j, k);
		for (size_t j = 0; j < N; j++)
			f2out << log_sum[j]/i/T << ",";
		f2out << std::endl;
	}
	f1out.close();
}