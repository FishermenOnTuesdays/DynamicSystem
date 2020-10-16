#include "DynamicSystem.h"

//Other functions:

/*Other functions*/

//Public methods:

DynamicSystem::DynamicSystem(const Eigen::VectorXld& starting_point, const std::vector<std::string>& strings_functions, std::string additional_variables)
	: dimension(strings_functions.size())
{
	std::string variables;
	std::string sum_all_functions;
	FunctionParser_ld parser_for_variables;
	for (auto function : strings_functions)
		sum_all_functions += function + '+';
	sum_all_functions.pop_back();
	parser_for_variables.ParseAndDeduceVariables(additional_variables + sum_all_functions, variables);
	for (auto function : strings_functions)
	{
		FunctionParser_ld function_parser;
		function_parser.Parse(additional_variables+function, variables);
		this->functions.push_back(function_parser);
	}
	this->point_of_trajectory = starting_point;
	CalculateJacobianMatrix();
}

std::vector<Eigen::VectorXld> DynamicSystem::GetTrajectory(long double time)
{
	std::vector<Eigen::VectorXld> points_of_trajectory;
	for (size_t i = 0; i < static_cast<size_t>(time/dt); i++)
	{
		points_of_trajectory.push_back(this->point_of_trajectory);
		NextPointOfTrajectory();
	}
	return points_of_trajectory;
}

std::vector<long double> DynamicSystem::GetSpectrumLyapunov()
{
	size_t M = 10000;
	size_t T = 1;
	std::vector<long double> spectrum_of_lyapunov_exponents;
	std::vector<long double> sums_of_logarithms(this->dimension);
	Eigen::MatrixXld variation_matrix = Eigen::MatrixXld::Identity(this->dimension, this->dimension);
	for (size_t i = 0; i < M; i++)
	{
		for (size_t j = 0; j < T; j++)
		{
			Eigen::MatrixXld k1, k2, k3, k4, buffer_variation;
			buffer_variation = variation_matrix;
			k1 = this->jacobian_matrix*buffer_variation;
			buffer_variation = variation_matrix + k1 * this->dt / 2;
			k2 = this->jacobian_matrix * buffer_variation;
			buffer_variation = variation_matrix + k2 * this->dt / 2;
			k3 = this->jacobian_matrix * buffer_variation;
			buffer_variation = variation_matrix + k3 * this->dt;
			k4 = this->jacobian_matrix * buffer_variation;
			variation_matrix += this->dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
			NextPointOfTrajectory();
		}
		auto QR = variation_matrix.householderQr();
		Eigen::VectorXld diagonal = QR.matrixQR().diagonal();
		for (size_t j = 0; j < this->dimension; j++)
		{
			sums_of_logarithms[j] += logl(fabsl(diagonal(j)));
		}
		variation_matrix = QR.householderQ();
	}
	for (size_t i = 0; i < this->dimension; i++)
		spectrum_of_lyapunov_exponents.push_back(sums_of_logarithms[i] / M / T / this->dt);
	return spectrum_of_lyapunov_exponents;
}

//Private methods:

Eigen::VectorXld DynamicSystem::f(const Eigen::VectorXld& vector)
{
	Eigen::VectorXld result_vector(this->dimension);
	size_t i = 0;
	for (auto& function : this->functions)
		result_vector[i++] = function.Eval(vector.data());
	return result_vector;
}

void DynamicSystem::ExplicitRungeKuttaFourthOrder()
{
	Eigen::VectorXld k1, k2, k3, k4, buffer_point_of_trajectory;
	buffer_point_of_trajectory = this->point_of_trajectory;
	k1 = f(buffer_point_of_trajectory);
	for (size_t i = 0; i < this->dimension; i++)
	{
		buffer_point_of_trajectory[i] = this->point_of_trajectory[i] + k1[i] * this->dt / 2;
	}
	k2 = f(buffer_point_of_trajectory);
	for (size_t i = 0; i < this->dimension; i++)
	{
		buffer_point_of_trajectory[i] = this->point_of_trajectory[i] + k2[i] * this->dt / 2;
	}
	k3 = f(buffer_point_of_trajectory);
	for (size_t i = 0; i < this->dimension; i++)
	{
		buffer_point_of_trajectory[i] = this->point_of_trajectory[i] + k3[i] * this->dt;
	}
	k4 = f(buffer_point_of_trajectory);
	for (size_t i = 0; i < this->dimension; i++)
	{
		this->point_of_trajectory[i] += this->dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
	}
}

void DynamicSystem::NextPointOfTrajectory()
{
	if (false/*Dynamic system is hard?*/)
	{
		//Make implementation
		/*switch (this->implicit_method)
		{
		case ImplicitNumericalMethod::EulerImplicit:
			ImplicitEuler();
			break;
		}*/
	}
	else
	{
		switch (this->explicit_method)
		{
		//Make implementation
		/*case ExplicitNumericalMethod::EulerExplicit:
			ExplicitEuler();
			break;*/
		case ExplicitNumericalMethod::RungeKuttaFourthOrder:
			ExplicitRungeKuttaFourthOrder();
			break;
		}
	}
	CalculateJacobianMatrix();
}

void DynamicSystem::CalculateJacobianMatrix()
{
	this->jacobian_matrix = Eigen::MatrixXld::Zero(this->dimension, this->dimension);
	for (size_t i = 0; i < this->dimension; i++)
	{
		for (size_t j = 0; j < this->dimension; j++)
		{
			Eigen::VectorXld left_point = this->point_of_trajectory;
			left_point(j) -= this->epsilon;
			Eigen::VectorXld right_point = this->point_of_trajectory;
			right_point(j) += this->epsilon;
			this->jacobian_matrix(i, j) =
				(this->functions[i].Eval(right_point.data()) -
					this->functions[i].Eval(left_point.data())) / (2 * this->epsilon);
		}
	}
}
