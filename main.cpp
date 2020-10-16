#include<iostream>
#include<fstream>
#include<algorithm>
#include<vector>
#include<complex>
#include<cmath>
#include "fparser.hh"
#include "Eigen/Dense"
#include "eigen-3.3.7/unsupported/Eigen/MatrixFunctions"
#include "DynamicSystem.h"
#include "symbolicc++/symbolicc++.h"

//typedef std::vector<std::vector<long double>> matrix_double;
//typedef std::complex<long double> complex;
//namespace Eigen
//{	
//	typedef Eigen::Matrix<long double, Dynamic, Dynamic> MatrixXld;
//	typedef Eigen::Matrix<long double, Dynamic, 1> VectorXld;
//}
//
//const double x_0 = 0.1;
//const double y_0 = 0.1;
//long int max_time = 10;
//long double dt = 0.01;
//const double dh = 0; //Шаг векторного поля
//const double E1 = 0.1;
//const double E2 = 0.1;
//const std::vector<std::string> Lorenz_attractor = { "10*(x2-x1)", "x1*(28-x3)-x2", "x1*x2-(8/3)*x3"};
//const std::vector<std::string> linear = {"-x1", "-2*x2", "-3*x3"};
//const std::vector<std::string> fixed_point = { "-x3", "-x1", "-x2" };
//const std::vector<std::string> variation = {"10*(x5-x4)", "28*x4-x5-x1*x6-x4*x3", "-8/3*x6+x1*x5+x4*x2"};
//
//std::vector<long double> var;
//size_t N;
//std::vector<long double> lyapunov_exponents;
//bool exact;
//bool is_attractor = false;
//Eigen::MatrixXld evolution_matrix;
//
//
//double x = x_0;
//double y = y_0;
//long int t = 0;
//
//std::string GetVariables(const std::vector<std::string>& functions)
//{
//	std::string variables = "";
//	for (size_t i = 1; i <= functions.size(); i++)
//		variables += 'x' + std::to_string(i) + ',';
//	variables.pop_back();
//	return variables;
//}
//
//std::vector<long double> f(const std::vector<long double>& coor, const std::vector<std::string>& functions, bool deviation = false)
//{
//	/*std::vector<long double> result{};
//	FunctionParser_ld fp;
//	for (size_t i = 0; i < functions.size(); i++)
//	{
//		fp.Parse(functions[i], GetVariables(functions));
//		result.push_back(fp.Eval((long double*)coor.data()));
//	}*/
//	std::vector<long double> result;
//	if (deviation)
//		result = { 10 * (coor[1] - coor[0]), 28*coor[0]  - coor[1] - var[0]*coor[2] - coor[0]*var[2], var[1]*coor[0] + var[0]*coor[1] - (8 / 3) * coor[2] };
//	else
//		result = { 10 * (coor[1] - coor[0]), coor[0] * (28 - coor[2]) - coor[1], coor[0] * coor[1] - (8 / 3) * coor[2] };
//	return result;
//}
//
//std::vector<long double> f(const Eigen::VectorXld& coor, const std::vector<std::string>& functions, bool deviation = false)
//{
//	/*std::vector<long double> result{};
//	FunctionParser_ld fp;
//	for (size_t i = 0; i < functions.size(); i++)
//	{
//		fp.Parse(functions[i], GetVariables(functions));
//		result.push_back(fp.Eval((long double*)coor.data()));
//	}*/
//	std::vector<long double> result;
//	if (deviation)
//		result = { 10 * (coor[1] - coor[0]), 28 * coor[0] - coor[1] - var[0] * coor[2] - coor[0] * var[2], var[1] * coor[0] + var[0] * coor[1] - (8 / 3) * coor[2] };
//	else
//		result = { 10 * (coor[1] - coor[0]), coor[0] * (28 - coor[2]) - coor[1], coor[0] * coor[1] - (8 / 3) * coor[2] };
//	return result;
//}
//
//std::vector<long double> Sub(const std::vector<long double>& vec1, const std::vector<long double>& vec2)
//{
//	if (vec1.size() != vec2.size())
//		throw std::exception("Vectors aren't equal");
//	std::vector<long double> sub_vec(vec1.size());
//	for (size_t i = 0; i < vec1.size(); i++)
//		sub_vec[i] = vec1[i] - vec2[i];
//	return sub_vec;
//}
//
//std::vector<long double> Add(const std::vector<long double>& vec1, const std::vector<long double>& vec2)
//{
//	if (vec1.size() != vec2.size())
//		throw std::exception("Vectors aren't equal");
//	std::vector<long double> add_vec(vec1.size());
//	for (size_t i = 0; i < vec1.size(); i++)
//		add_vec[i] = vec1[i] + vec2[i];
//	return add_vec;
//}
//
//std::vector<long double> Mult(const std::vector<long double>& vec, long double scalar)
//{
//	std::vector<long double> mult_vec(vec.size());
//	for (size_t i = 0; i < vec.size(); i++)
//		mult_vec[i] = vec[i] * scalar;
//	return mult_vec;
//}
//
//long double LengthVector(const std::vector<long double>& vec)
//{
//	long double length = 0;
//	for (size_t i = 0; i < vec.size(); i++)
//		length += vec[i] * vec[i];
//	return sqrtl(length);
//}
//
//bool Equal(const std::vector<long double>& vec1, const std::vector<long double>& vec2, const long double eps = 0.0000001)
//{
//	if (vec1.size() != vec2.size())
//		return false;
//	for (size_t i = 0; i < vec1.size(); i++)
//		if (fabsl(vec1[i] - vec2[i]) > eps)
//			return false;
//	return true;
//}
//
//matrix_double GetJacobianMatrix(const std::vector<long double>& coor, const std::vector<std::string>& functions, const long double eps = 0.0000001)
//{
//	matrix_double return_matrix{ functions.size(), std::vector<long double>(functions.size()) };
//	std::vector<long double> copy_coor_plus = coor;
//	std::vector<long double> copy_coor_minus = coor;
//	for (size_t i = 0; i < functions.size(); i++)
//	{
//		for (size_t j = 0; j < functions.size(); j++)
//		{
//			copy_coor_plus[j] += eps;
//			copy_coor_minus[j] -= eps;
//			return_matrix[i][j] = (f(copy_coor_plus, functions)[i] - f(copy_coor_minus, functions)[i]) / (2 * eps);
//			copy_coor_plus[j] = coor[j];
//			copy_coor_minus[j] = coor[j];
//		}
//	}
//	return return_matrix;
//}
//
//namespace Eigen
//{
//	Eigen::MatrixXld GetJacobianMatrix(const std::vector<long double>& coor, const std::vector<std::string>& functions, const long double eps = 0.0000001)
//	{
//		auto jacobian_matrix = ::GetJacobianMatrix(coor, functions);
//		Eigen::MatrixXld jacobian_matrix_eigen(N,N);
//		for (size_t i = 0; i < N; i++)
//			for (size_t j = 0; j < N; j++)
//				jacobian_matrix_eigen(i, j) = jacobian_matrix[i][j];
//		return jacobian_matrix_eigen;
//	}
//}
//
//long double GetNorm(const matrix_double& curr_matrix)
//{
//	long double norm = 0;
//	for (auto vec : curr_matrix)
//		for (auto el : vec)
//			norm += el * el;
//	return norm;
//}
//
//complex Determinant(const std::vector<std::vector<complex>>& curr_matrix)
//{
//	if (curr_matrix.size() == 1)
//		return curr_matrix[0][0];
//	complex determinant = 0;
//	int k = 1;
//	for (size_t i = 0; i < curr_matrix.size(); i++)
//	{
//		if (std::abs(curr_matrix[0][i]) != 0)
//		{
//			std::vector<std::vector<complex>> new_matrix = curr_matrix;
//			new_matrix.erase(new_matrix.begin());
//			for (size_t j = 0; j < new_matrix.size(); j++)
//			{
//				new_matrix[j].erase(new_matrix[j].begin() + i);
//			}
//			determinant += complex(k,0) * curr_matrix[0][i] * Determinant(new_matrix);
//		}
//		k *= -1;
//	}
//	return determinant;
//}
//
//template<typename T>
//T Determinant(const std::vector<std::vector<T>>& curr_matrix)
//{
//	if (curr_matrix.size() == 1)
//		return curr_matrix[0][0];
//	T determinant = 0;
//	int k = 1;
//	for (size_t i = 0; i < curr_matrix.size(); i++)
//	{
//		if (curr_matrix[0][i] != 0)
//		{
//			std::vector<std::vector<T>> new_matrix = curr_matrix;
//			new_matrix.erase(new_matrix.begin());
//			for (size_t j = 0; j < new_matrix.size(); j++)
//			{
//				new_matrix[j].erase(new_matrix[j].begin() + i);
//			}
//			determinant += k * curr_matrix[0][i] * Determinant(new_matrix);
//		}
//		k *= -1;
//	}
//	return determinant;
//}
//
//matrix_double AdjugateT(const matrix_double& curr_matrix)
//{
//	matrix_double adjugateT{ curr_matrix.size() , std::vector<long double>(curr_matrix.size()) };
//	int sign;
//	for (size_t i = 0; i < curr_matrix.size(); i++)
//	{
//		for (size_t j = 0; j < curr_matrix.size(); j++)
//		{
//			matrix_double new_matrix = curr_matrix;
//			new_matrix.erase(new_matrix.begin() + i);
//			for (size_t k = 0; k < new_matrix.size(); k++)
//			{
//				new_matrix[k].erase(new_matrix[k].begin() + j);
//			}
//			sign = 1 - 2 * ((i + j) % 2);
//			adjugateT[j][i] = sign * Determinant(new_matrix);
//		}
//	}
//	return adjugateT;
//}
//
//matrix_double MatrixInverse(const matrix_double& curr_matrix)
//{
//	matrix_double matrix_inverse{ curr_matrix.size(),  std::vector<long double>(curr_matrix.size()) };
//	long double det = Determinant(curr_matrix);
//	int sign;
//	for (size_t i = 0; i < curr_matrix.size(); i++)
//	{
//		for (size_t j = 0; j < curr_matrix.size(); j++)
//		{
//			matrix_double new_matrix = curr_matrix;
//			new_matrix.erase(new_matrix.begin() + i);
//			for (size_t k = 0; k < new_matrix.size(); k++)
//			{
//				new_matrix[k].erase(new_matrix[k].begin() + j);
//			}
//			sign = 1 - 2 * ((i + j) % 2);//  1 or -1
//			matrix_inverse[j][i] = sign * Determinant(new_matrix) / det;
//		}
//	}
//	return matrix_inverse;
//}
//
//complex CharacteristicEquation(const complex& z, const matrix_double& curr_matrix)
//{
//	std::vector<std::vector<complex>> new_matrix{ curr_matrix.size(), std::vector<complex>(curr_matrix.size()) };
//	for (size_t i = 0; i < new_matrix.size(); i++)
//		new_matrix[i][i] = curr_matrix[i][i]-z;
//	return Determinant(new_matrix);
//}
//
//bool IsHard(const std::vector<long double>& coor, const std::vector<std::string>& functions, const long double eps = 0.0000001)
//{
//	matrix_double jacobian_matrix = GetJacobianMatrix(coor, functions);
//	Eigen::MatrixXld jacobian_matrix_new(jacobian_matrix.size(), jacobian_matrix.size());
//	for (size_t i = 0; i < jacobian_matrix.size(); i++)
//		for (size_t j = 0; j < jacobian_matrix.size(); j++)
//			jacobian_matrix_new(i, j) = jacobian_matrix[i][j];
//	//if(is_attractor)
//		//evolution_matrix = jacobian_matrix_new * evolution_matrix;
//	auto eigenvalues = jacobian_matrix_new.eigenvalues();
//	//const double lyapunov_eps = 0.01;
//	//if (!lyapunov_exponents.empty())
//	//{
//	//	exact = true;
//	//	for (size_t i = 0; i < eigenvalues.size(); i++)
//	//		if (fabs(eigenvalues[i].real() - lyapunov_exponents[i]) > lyapunov_eps)
//	//		{
//	//			exact = false;
//	//			break;
//	//		}
//	//}
//	//else
//	//	exact = false;
//	//lyapunov_exponents.clear();
//	//for (size_t i = 0; i < eigenvalues.size(); i++)
//	//{
//	//	//std::cout << eigenvalues[i] << " ";
//	//	lyapunov_exponents.push_back((evolution_matrix.transpose()*evolution_matrix).eigenvalues()[i].real());
//	//}
//	////std::cout << std::endl;
//	std::vector<long double> abs_eigenvalues;
//	for (size_t i = 0; i < eigenvalues.size(); i++)
//	{
//		if (std::abs(eigenvalues[i]) < eps)
//			return true;
//		abs_eigenvalues.push_back(std::abs(eigenvalues[i]));
//	}
//	if (*std::max_element(abs_eigenvalues.begin(), abs_eigenvalues.end()) / *std::min_element(abs_eigenvalues.begin(), abs_eigenvalues.end()) > 100)
//		return true;
//	return false;
//}
//
//size_t KroneckerSymbol(size_t i, size_t j)
//{
//	return i == j ? 1 : 0;
//}
//
//template<typename T>
//T Dot(const std::vector<T>& vec1, const std::vector<T>& vec2)
//{
//	T dot=0;
//	if (vec1.size() != vec2.size())
//		throw std::exception("Vectors aren't equal");
//	for (size_t i = 0; i < vec1.size(); i++)
//	{
//		dot += vec1[i] * vec2[i];
//	}
//	return dot;
//}
//
//std::vector<long double> CountNextCoor(const std::vector<long double>& coor, const std::vector<std::string>& functions, const long double dt, bool deviation = false)
//{
//	std::vector<long double> result;
//	/*if (IsHard(coor, functions))
//	{
//		const matrix_double jacobian_matrix = GetJacobianMatrix(coor, functions);
//		matrix_double new_matrix{ jacobian_matrix.size(), std::vector<long double>(jacobian_matrix.size()) };
//		for (size_t i = 0; i < new_matrix.size(); i++)
//		{
//			for (size_t j = 0; j < new_matrix.size(); j++)
//			{
//				new_matrix[i][j] = KroneckerSymbol(i, j) - jacobian_matrix[i][j] * dt;
//			}
//		}
//		new_matrix = MatrixInverse(new_matrix);
//		std::vector<long double> functions_coor_dt = f(coor, functions);
//		for (size_t i = 0; i < coor.size(); i++)
//			result.push_back(coor[i]+Dot(new_matrix[i], functions_coor_dt));
//		return result;
//	}*/
//	std::vector<long double> k1, k2, k3, k4, current_coor;
//	current_coor = coor;
//	k1 = f(current_coor, functions, deviation);
//	for (size_t i = 0; i < current_coor.size(); i++)
//	{
//		current_coor[i] = coor[i] + k1[i] * dt / 2;
//	}
//	k2 = f(current_coor, functions, deviation);
//	for (size_t i = 0; i < current_coor.size(); i++)
//	{
//		current_coor[i] = coor[i] + k2[i] * dt / 2;
//	}
//	k3 = f(current_coor, functions, deviation);
//	for (size_t i = 0; i < current_coor.size(); i++)
//	{
//		current_coor[i] = coor[i] + k3[i] * dt;
//	}
//	k4 = f(current_coor, functions, deviation);
//	for (size_t i = 0; i < coor.size(); i++)
//	{
//		result.push_back(coor[i] + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
//	}
//	return result;
//	/*double k1_x, k2_x, k3_x, k4_x;
//	double k1_y, k2_y, k3_y, k4_y;
//	k1_x = F1(x, y);
//	k1_y = F2(x, y);
//	k2_x = F1(x + dt / 2 * k1_x, y + dt / 2 * k1_y);
//	k2_y = F2(x + dt / 2 * k1_x, y + dt / 2 * k1_y);
//	k3_x = F1(x + dt / 2 * k2_x, y + dt / 2 * k2_y);
//	k3_y = F2(x + dt / 2 * k2_x, y + dt / 2 * k2_y);
//	k4_x = F1(x + dt * k3_x, y + dt * k3_y);
//	k4_y = F2(x + dt * k3_x, y + dt * k3_y);
//	dx = dt / 6 * (k1_x + 2 * k2_x + 2 * k3_x + k4_x);
//	dy = dt / 6 * (k1_y + 2 * k2_y + 2 * k3_y + k4_y);*/
//}
//
//Eigen::VectorXld CountNextCoorTest(const Eigen::VectorXld& coor, const std::vector<std::string>& functions, const long double dt, bool deviation = false)
//{
//	//std::vector<long double> result;
//	/*if (IsHard(coor, functions))
//	{
//		const matrix_double jacobian_matrix = GetJacobianMatrix(coor, functions);
//		matrix_double new_matrix{ jacobian_matrix.size(), std::vector<long double>(jacobian_matrix.size()) };
//		for (size_t i = 0; i < new_matrix.size(); i++)
//		{
//			for (size_t j = 0; j < new_matrix.size(); j++)
//			{
//				new_matrix[i][j] = KroneckerSymbol(i, j) - jacobian_matrix[i][j] * dt;
//			}
//		}
//		new_matrix = MatrixInverse(new_matrix);
//		std::vector<long double> functions_coor_dt = f(coor, functions);
//		for (size_t i = 0; i < coor.size(); i++)
//			result.push_back(coor[i]+Dot(new_matrix[i], functions_coor_dt));
//		return result;
//	}*/
//	std::vector<long double> k1, k2, k3, k4;
//	Eigen::VectorXld current_coor, result(N);
//	current_coor = coor;
//	k1 = f(current_coor, functions, deviation);
//	for (size_t i = 0; i < current_coor.size(); i++)
//	{
//		current_coor[i] = coor[i] + k1[i] * dt / 2;
//	}
//	k2 = f(current_coor, functions, deviation);
//	for (size_t i = 0; i < current_coor.size(); i++)
//	{
//		current_coor[i] = coor[i] + k2[i] * dt / 2;
//	}
//	k3 = f(current_coor, functions, deviation);
//	for (size_t i = 0; i < current_coor.size(); i++)
//	{
//		current_coor[i] = coor[i] + k3[i] * dt;
//	}
//	k4 = f(current_coor, functions, deviation);
//	for (size_t i = 0; i < coor.size(); i++)
//	{
//		result[i] = coor[i] + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
//	}
//	return result;
//	/*double k1_x, k2_x, k3_x, k4_x;
//	double k1_y, k2_y, k3_y, k4_y;
//	k1_x = F1(x, y);
//	k1_y = F2(x, y);
//	k2_x = F1(x + dt / 2 * k1_x, y + dt / 2 * k1_y);
//	k2_y = F2(x + dt / 2 * k1_x, y + dt / 2 * k1_y);
//	k3_x = F1(x + dt / 2 * k2_x, y + dt / 2 * k2_y);
//	k3_y = F2(x + dt / 2 * k2_x, y + dt / 2 * k2_y);
//	k4_x = F1(x + dt * k3_x, y + dt * k3_y);
//	k4_y = F2(x + dt * k3_x, y + dt * k3_y);
//	dx = dt / 6 * (k1_x + 2 * k2_x + 2 * k3_x + k4_x);
//	dy = dt / 6 * (k1_y + 2 * k2_y + 2 * k3_y + k4_y);*/
//}
//
//void CountLogSum(std::vector<long double>& log_sum, const Eigen::MatrixXld& matrix_deviation, const long double eps_deviation)
//{
//	Eigen::HouseholderQR<Eigen::MatrixXld> qr(matrix_deviation);
//	Eigen::MatrixXld orthonormal_matrix_deviation = qr.householderQ();
//	//orthonormal_matrix_deviation *= eps_deviation;
//	for (size_t i = 0; i < N; i++)
//	{
//		Eigen::VectorXld current_row = matrix_deviation.row(i);
//		for (size_t j = 0; j < i; j++)
//		{
//			current_row -= matrix_deviation.row(i).dot(orthonormal_matrix_deviation.row(j)) * orthonormal_matrix_deviation.row(j);
//		}
//		log_sum[i] += std::logl(current_row.norm());
//	}
//}
//
//void CountNextCoorDeviation(matrix_double& var_deviation, const std::vector<std::string>& functions, const double dt)
//{
//	for (size_t i = 0; i < N; i++)
//		var_deviation[i] = CountNextCoor(var_deviation[i], functions, dt);
//	/*matrix_double var_deviation{ N, std::vector<long double>(N)};
//	for (size_t i = 0; i < N; i++)
//	{
//		for (size_t j = 0; j < N; j++)
//			var_deviation[i][j] = var[j] + matrix_deviation(i,j);
//		var_deviation[i] = CountNextCoor(var_deviation[i], functions, dt);
//	}
//	var = CountNextCoor(var, functions, dt);
//	for (size_t i = 0; i < N; i++)
//		for (size_t j = 0; j < N; j++)
//			matrix_deviation(i, j) = var_deviation[i][j] - var[j];*/
//}
//
//void CountDeviation(Eigen::MatrixXld& deviation, const matrix_double& var_deviation, const std::vector<long double>& var)
//{
//	for (size_t i = 0; i < N; i++)
//		for (size_t j = 0; j < N; j++)
//			deviation(i, j) = var_deviation[i][j] - var[j];
//}
//
//std::vector<std::string> GetLinearFunctions(Eigen::MatrixXld matrix)
//{
//	std::vector<std::string> functions{};
//	for (size_t i = 0; i < matrix.rows(); i++)
//	{
//		std::string current_function = "";
//		for (size_t j = 0; j < matrix.cols(); j++)
//		{
//			current_function += '+' + std::to_string(matrix(i, j)) + "*x" + std::to_string(j+1);
//		}
//		current_function.erase(current_function.begin());
//		functions.push_back(current_function);
//	}
//	return functions;
//}
//
//namespace Eigen
//{
//	MatrixXld Orthogonalize(MatrixXld matrix)
//	{
//		for (size_t i = 0; i < matrix.cols(); i++)
//		{
//			for (size_t j = 0; j < i; j++)
//			{
//				matrix.col(i) -= matrix.col(i).dot(matrix.col(j).normalized()) * matrix.col(j).normalized();
//			}
//		}
//		return matrix;
//	}
//
//	MatrixXld Normalize(MatrixXld matrix)
//	{
//		for (size_t i = 0; i < matrix.cols(); i++)
//		{
//			matrix.col(i).normalize();
//		}
//		return matrix;
//	}
//}
//
//int main()
//{
//	std::ofstream f1out, f2out;
//	f1out.open("../wwwroot/output/result.csv");//Введите свой путь
//	f2out.open("../wwwroot/output/laypunov.csv");//Введите свой путь
//	std::vector<std::string> functions;
//	#ifndef _DEBUG
//		N = 3;
//	#else
//		std::cin >> N;
//	#endif
//	for (size_t i = 1; i <= N; i++)
//		f1out << 'x' + std::to_string(i) + ',';
//	f1out << "t,";
//	for (size_t i = 1; i <= N-1; i++)
//		f1out << 'l' + std::to_string(i) + ',';
//	f1out << 'l' + std::to_string(N) + '\n';
//	functions.resize(N);
//	#ifndef _DEBUG
//		functions = Lorenz_attractor;
//	#else
//		for (size_t i = 0; i < N; i++)
//			std::cin >> functions[i];
//	#endif
//	var.resize(N);
//	#ifndef _DEBUG
//		for (size_t i = 0; i < N; i++)
//			var[i] = 0.1;
//	#else
//		for (size_t i = 0; i < N; i++)
//			std::cin >> var[i];
//	#endif
//	std::vector<std::vector<long double>> dots;
//	std::cin >> max_time;
//	std::cin >> dt;
//	while(t < max_time)
//	{
//		dots.push_back(var);
//		for (size_t i = 0; i < var.size(); i++)
//		{
//			f1out << std::to_string(var[i])+',';
//		}
//		var = CountNextCoor(var, functions, dt);
//		f1out << std::to_string(t * dt);
//		for (size_t i = 0; i < N; i++)
//			f1out << ',';
//		f1out << '\n';
//		t++;
//	}
//
//	std::cout << "End var\n";
//
//	/*std::cout << "GOOD\n";
//	for (size_t i = 0; i < N; i++)
//		std::cout << var[i] << ", ";
//	std::cout << "\nGOOD\n";*/
//	long int M = 10000;
//	long int T = 1;
//	long double eps_deviation = 1;
//	//std::vector<long double> deviation_x(var.size(), 0);
//	//deviation_x[0] = eps_deviation;
//	//std::vector<long double> deviation_y(var.size(), 0);
//	//deviation_y[1] = eps_deviation;
//	//std::vector<long double> var_deviation_x = Add(var, deviation_x);
//	//std::vector<long double> var_deviation_y = Add(var, deviation_y);
//	//long double log_sum_x = 0, log_sum_y = 0, log_sum_z = 0;
//	//Eigen::MatrixXld deviation = Eigen::MatrixXld::Identity(N,N);
//	//deviation *= eps_deviation;
//	//matrix_double var_deviation{ N, std::vector<long double>(N) };
//	/*for (size_t i = 0; i < N; i++)
//		for (size_t j = 0; j < N; j++)
//			var_deviation[i][j] = var[j] + deviation(i, j);*/
//	std::vector<long double> log_sum(N, 0);
//
//	//Eigen::MatrixXld Y = Eigen::MatrixXld::Identity(N,N);
//	//Eigen::VectorXld x_deviation(N);
//	/*for (size_t i = N / 2; i < N; i++)
//		var[i] = 1;*/
//	Eigen::MatrixXld deviation = Eigen::MatrixXld::Zero(N,N);
//	/*deviation(0, 0) = 10 * (var[1] - var[0]);
//	deviation(1, 0) = 28 * var[0] - var[1] - var[0] * var[2];
//	deviation(2, 0) = -8 / 3 * var[2] + var[0] * var[1];
//	deviation(1, 1) = 1;
//	deviation(2, 2) = 1;*/
//	//deviation = Eigen::Orthogonalize(deviation);
//	//deviation = Eigen::Normalize(deviation);
//	/*std::cout << Eigen::Orthogonalize(deviation).col(0).norm() << std::endl;*/
//	//for (size_t index_swaps = 0; index_swaps < N; index_swaps++)
//	//{
//		for (size_t i = 0; i < deviation.rows(); i++)
//			deviation(i, 0) = f(var, functions)[i];
//		for (size_t i = 1; i < N; i++)
//			deviation(i, i) = 1;
//
//		//deviation.col(0).swap(deviation.col(index_swaps));
//		//Eigen::HouseholderQR<Eigen::MatrixXld> qr(deviation);
//		/*Eigen::MatrixXld Q = qr.householderQ();
//		deviation = Q;*/
//		//std::cout << "Q: " << deviation << std::endl;
//		//deviation = Eigen::Orthogonalize(deviation);
//		//deviation = Eigen::Normalize(deviation);
//		deviation = deviation.householderQr().householderQ();
//		deviation *= eps_deviation;
//		for (long unsigned int i = 0; i < M; i++)
//		{
//			//Y += Eigen::GetJacobianMatrix(var, functions) * Y * dt;
//			//deviation += Eigen::GetJacobianMatrix(var, functions) * deviation * dt;
//			//std::cout << deviation.row(0).dot(deviation.row(1)) << deviation.row(1).dot(deviation.row(2)) << deviation.row(2).dot(deviation.row(0)) << std::endl;
//			//Eigen::MatrixXld Jacobian_Matrix = Eigen::GetJacobianMatrix(var,functions);
//			for (long unsigned int j = 0; j < T; j++)
//			{
//				//CountNextCoorDeviation(var_deviation, functions, dt);
//				/*Eigen::MatrixXld j_m = Eigen::GetJacobianMatrix(var, functions);*/
//				for (size_t v = 0; v < N; v++)
//				{
//					/*std::vector<long double> curr(N);*/
//					/*for (size_t k = 0; k < N; k++)
//					{
//						curr[k] = deviation(k, v);
//					}*/
//					deviation.col(v) = CountNextCoorTest(deviation.col(v), GetLinearFunctions(Eigen::GetJacobianMatrix(var, functions)), dt, true);
//					CountNextCoor(var, functions, dt);
//					/*for (auto el : GetLinearFunctions(Eigen::GetJacobianMatrix(var, functions)))
//					{
//						std::cout << el << std::endl;
//					}
//					std::cout << "----------------------\n";*/
//					/*var = CountNextCoor(var, functions, dt);*/
//					/*for (size_t k = 0; k < N; k++)
//					{
//						deviation(k, v) = curr[k];
//					}*/
//				}
//				//deviation += Eigen::GetJacobianMatrix(var, functions) * deviation * dt;
//			}
//
//			deviation = Eigen::Orthogonalize(deviation);
//			//std::cout << "R:\n" << (deviation.householderQr().householderQ().transpose()*deviation).diagonal()(0) << std::endl;
//			//deviation = deviation.householderQr().householderQ();
//			for (size_t j = 0; j < N; j++)
//				log_sum[j] += logl(fabsl(deviation.col(j).norm()));
//			/*for (size_t i = 0; i < deviation.rows(); i++)
//				deviation(i, 0) = f(var, functions)[i];
//			for (size_t i = 1; i < N; i++)
//				deviation(i, i) = 1;*/
//			//deviation = deviation.householderQr().householderQ();
//			//deviation *= eps_deviation;
//			//std::cout << log_sum[j] << std::endl;
//			//deviation = qr.householderQ();
//			deviation = Eigen::Normalize(deviation);
//			/*for (size_t i = 0; i < deviation.rows(); i++)
//				deviation(i, index_swaps) = f(var, functions)[i];
//			deviation = Eigen::Orthogonalize(deviation);*/
//			/*std::cout << deviation << std::endl;
//			Eigen::HouseholderQR<Eigen::MatrixXld> qr(deviation);
//			Eigen::MatrixXld Q = qr.householderQ();
//			std::cout << "Norms:\n";
//			std::cout << Q.col(0).norm() << '\n' << Q.col(1).norm() << '\n' << Q.col(2).norm() << '\n';*/
//
//			/*log_sum_x += logl(deviation.col(0).norm() / eps_deviation);
//			deviation.col(0) = deviation.col(0) * eps_deviation / deviation.col(0).norm();
//			deviation.col(1) = deviation.col(1) - deviation.col(1).dot(deviation.col(0)) * deviation.col(0);
//			log_sum_y += logl(deviation.col(1).norm() / eps_deviation);
//			deviation.col(1) = deviation.col(1) * eps_deviation / deviation.col(1).norm();
//			deviation.col(2) = deviation.col(2) - deviation.col(2).dot(deviation.col(0)) * deviation.col(0) - deviation.col(2).dot(deviation.col(1)) * deviation.col(1);
//			log_sum_z += logl(deviation.col(2).norm() / eps_deviation);
//			deviation.col(2) = deviation.col(2) * eps_deviation / deviation.col(2).norm();*/
//
//
//			//CountDeviation(deviation, var_deviation, var);
//			//CountLogSum(log_sum, deviation, eps_deviation);
//			//Eigen::HouseholderQR<Eigen::MatrixXld> qr(deviation);
//			//deviation = qr.householderQ();
//			//deviation *= eps_deviation;
//			/*for (size_t j = 0; j < N; j++)
//				for (size_t k = 0; k < N; k++)
//					var_deviation[j][k] = var[k] + deviation(j, k);*/
//					if (i != 0)
//					{
//						for (size_t j = 0; j < N; j++)
//						{
//							f1out << var[j] << ',';
//						}
//						f1out << i << ',';
//						for (size_t j = 0; j < N; j++)
//							if (j != N - 1)
//								f1out << log_sum[j] / i / T << ",";
//							else
//								f1out << log_sum[j] / i / T / dt;
//						f1out << std::endl;
//					}
//					//f2out << log_sum_x / i / T / dt << ',' << log_sum_y / i / T / dt << ',' << log_sum_z / i / T / dt << std::endl;
//					//std::cout << i << std::endl;
//		}
//	//}
//	/*Eigen::MatrixXld lymbda = (Y * Y.transpose()).log()/(2*M);
//	auto eigenvalues_lymbda = lymbda.eigenvalues();
//	for (size_t i = 0; i < N; i++)
//	{
//		std::cout << eigenvalues_lymbda(i) << " | ";
//	}*/
//	std::cout << "Lyapunov exponents: " << std::endl;
//	for (size_t i = 0; i < N; i++)
//		std::cout << log_sum[i] / M / T / dt << std::endl;
//	//std::cout << Eigen::GetJacobianMatrix({ 1,1,0 }, functions).row(0);
//	f1out.close();
//}

int main()
{
	std::string additional_variables = "a:= 12; b:=sin(2); c:=cos(3);";
	std::vector<std::string> linear = {"-x1","-2*x2","-3*x3"};
	std::vector<std::string> lorenz = {"10*(y-x)","28*x-y-x*z","-8/3*z+x*y"};
	Eigen::VectorXld initial_point(3);
	initial_point << 0.1, 0.1, 0.1;
	DynamicSystem dyns{ initial_point, lorenz, additional_variables};
	std::vector<Eigen::VectorXld> points = dyns.GetTrajectory(10);
	/*for (auto vector : points)
		std::cout << vector << std::endl << std::endl;*/
	auto spectrum = dyns.GetSpectrumLyapunov();
	for (auto exponent : spectrum)
		std::cout << exponent << " ";
	std::cout << std::endl;
}