#include "DynamicSystem.h"

namespace DynS
{

	//Other functions:

	std::vector<std::pair<std::pair<long double, long double>, long double>> GetMapLyapunovExponents(
		const Eigen::VectorXld& starting_point,
		const std::vector<std::string>& strings_functions,
		std::string variables, std::string additional_equations,
		const std::pair<std::string, std::string>& parameters,
		const std::pair<std::pair<long double, long double>, std::pair<long double, long double>>& ranges,
		const std::pair<long double, long double>& steps,
		long double time_access_to_attractor,
		long double time_calculation_lyapunov_spectrum,
		long double dt
	)
	{
		size_t number_of_dots = (ranges.first.second - ranges.first.first) * (ranges.second.second - ranges.second.first) / steps.first / steps.second;
		size_t number_of_y = (ranges.second.second - ranges.second.first) / steps.second;
		std::vector<std::pair<std::pair<long double, long double>, long double>> map_lyapunov_spectrum{ number_of_dots };
		#pragma omp parallel for
		for (int i = 0; i < number_of_dots; i++)
		{
			long double first_parameter = ranges.first.first + (i / number_of_y) * steps.first;
			long double second_parameter = ranges.second.first + (i % number_of_y) * steps.second;
			DynamicSystem dynamic_system{
					starting_point,
					strings_functions,
					variables,
					additional_equations +
					parameters.first + ":=" + std::to_string(first_parameter) + ";" +
					parameters.second + ":=" + std::to_string(second_parameter) + ";"
			};
			dynamic_system.SetDt(dt);
			dynamic_system.GetTrajectory(time_access_to_attractor);
			auto spectrum = dynamic_system.GetSpectrumLyapunov(time_calculation_lyapunov_spectrum);
			map_lyapunov_spectrum[i] = { {first_parameter, second_parameter}, *std::max_element(spectrum.begin(), spectrum.end()) };
		}
		/*for (long double first_parameter = ranges.first.first; first_parameter < ranges.first.second; first_parameter+=steps.first)
			for (long double second_parameter = ranges.second.first; second_parameter < ranges.second.second; second_parameter+=steps.second)
			{
				DynamicSystem dynamic_system{
					starting_point,
					strings_functions,
					variables,
					additional_equations +
					parameters.first  + ":=" + std::to_string(first_parameter)  + ";" +
					parameters.second + ":=" + std::to_string(second_parameter) + ";"
				};
				dynamic_system.SetDt(dt);
				dynamic_system.GetTrajectory(time_access_to_attractor);
				auto spectrum = dynamic_system.GetSpectrumLyapunov(time_calculation_lyapunov_spectrum);
				map_lyapunov_spectrum[i] = { {first_parameter, second_parameter}, *std::max_element(spectrum.begin(), spectrum.end()) };
				i++;
			}*/
		return map_lyapunov_spectrum;
	}

	//Public methods:

	DynamicSystem::DynamicSystem(const Eigen::VectorXld& starting_point, const std::vector<std::string>& strings_functions, std::string variables, std::string additional_variables)
		: dimension(strings_functions.size())
	{
		if (variables == "")//Defined variables
		{
			std::string variables;
			std::string sum_all_functions;
			FunctionParser_ld parser_for_variables;
			for (auto function : strings_functions)
				sum_all_functions += function + '+';
			sum_all_functions.pop_back();
			parser_for_variables.ParseAndDeduceVariables(additional_variables + sum_all_functions, variables);
		}
		for (auto function : strings_functions)
		{
			FunctionParser_ld function_parser;
			function_parser.Parse(additional_variables + function, variables);
			this->functions.push_back(function_parser);
		}
		this->point_of_trajectory = starting_point;
		CalculateJacobianMatrix();
	}

	std::vector<Eigen::VectorXld> DynamicSystem::GetTrajectory(long double time)
	{
		std::vector<Eigen::VectorXld> points_of_trajectory;
		for (size_t i = 0; i < static_cast<size_t>(time / dt); i++)
		{
			points_of_trajectory.push_back(this->point_of_trajectory);
			NextPointOfTrajectory();
		}
		return points_of_trajectory;
	}

	std::vector<long double> DynamicSystem::GetSpectrumLyapunov(long double time)
	{
		size_t M = time / this->dt;
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
				k1 = this->jacobian_matrix * buffer_variation;
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

	std::map<std::string, std::vector<long double>> DynamicSystem::GetTimeSeriesSpectrumLyapunov(long double time)
	{
		std::map<std::string, std::vector<long double>> series_spectrum_lyapunov;
		for (size_t i = 0; i < this->dimension; i++)
			series_spectrum_lyapunov.emplace("lambda" + std::to_string(i+1), std::vector<long double>{});
		series_spectrum_lyapunov.emplace("t", std::vector<long double>{});
		size_t M = time / this->dt;
		size_t T = 1;
		std::vector<long double> sums_of_logarithms(this->dimension);
		Eigen::MatrixXld variation_matrix = Eigen::MatrixXld::Identity(this->dimension, this->dimension);
		for (size_t i = 0; i < M; i++)
		{
			for (size_t j = 0; j < T; j++)
			{
				Eigen::MatrixXld k1, k2, k3, k4, buffer_variation;
				buffer_variation = variation_matrix;
				k1 = this->jacobian_matrix * buffer_variation;
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
				sums_of_logarithms[j] += logl(fabsl(diagonal(j)));
			variation_matrix = QR.householderQ();
			if (i != 0)
			{
				for (size_t j = 0; j < this->dimension; j++)
					series_spectrum_lyapunov["lambda" + std::to_string(j + 1)].push_back(sums_of_logarithms[j] / i / T / this->dt);
				series_spectrum_lyapunov["t"].push_back(i*this->dt);
			}
		}
		return series_spectrum_lyapunov;
	}

	void DynamicSystem::SetDt(long double dt)
	{
		this->dt = dt;
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
		buffer_point_of_trajectory = this->point_of_trajectory + k1 * this->dt / 2;
		k2 = f(buffer_point_of_trajectory);
		buffer_point_of_trajectory = this->point_of_trajectory + k2 * this->dt / 2;
		k3 = f(buffer_point_of_trajectory);
		buffer_point_of_trajectory = this->point_of_trajectory + k3 * this->dt;
		k4 = f(buffer_point_of_trajectory);
		this->point_of_trajectory += this->dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
	}

	void DynamicSystem::NextPointOfTrajectory()
	{
		if (this->point_of_trajectory.norm() > 1e30)
			return;
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
}