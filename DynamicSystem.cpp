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
			try
			{
				dynamic_system.GetTrajectory(time_access_to_attractor);
			}
			catch (std::exception& ex)
			{
				if (ex.what() == "Infinity trajectory")
				{
					dynamic_system.SetCurrentPointOfTrajectory(starting_point);
				}
			}
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

	// Poincare

	int Sign(long double x) {
		return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
	}

	int SideSign(PlaneEquation equation, Eigen::Vector3ld point) {
		return Sign(equation.A * point[0] + equation.B * point[1] + equation.C * point[2] + equation.D);
	}

	using Basis3ld = std::vector<Eigen::Vector3ld>;

	Eigen::Vector2ld applyBasis(Basis3ld basis, Eigen::Vector3ld point) {
		return Eigen::Vector2ld(basis[0].dot(point), basis[1].dot(point));
	}

	Basis3ld transformBasis(Eigen::Vector3ld vector) {
		Basis3ld basis = { Eigen::Vector3ld(1, 0, 0), Eigen::Vector3ld(0, 1, 0), Eigen::Vector3ld(0, 0, 1) };
		if (vector == basis[0] || vector == basis[1] || vector == basis[2])
			return basis;
		else {
			Eigen::Vector3ld z = vector.normalized();
			Eigen::Vector3ld x = (basis[1] - z * basis[1].dot(z)).normalized();
			Eigen::Vector3ld y = x.cross(z);
			return Basis3ld({ x, y, z });
		}
	}

	Eigen::Vector3ld intersectionCalc(PlaneEquation planeEquation, Eigen::Vector3ld pointA, Eigen::Vector3ld pointB) {
		long double t = (planeEquation.A * pointA[0] + planeEquation.B * pointA[1] + planeEquation.C * pointA[2] + planeEquation.D) / (planeEquation.A * (pointA[0] - pointB[0]) + planeEquation.B * (pointA[1] - pointB[1]) + planeEquation.C * (pointA[2] - pointB[2]));
		return Eigen::Vector3ld(pointA[0] + (pointB[0] - pointA[0]) * t, pointA[1] + (pointB[1] - pointA[1]) * t, pointA[2] + (pointB[2] - pointA[2]) * t);
	}

	PoincareMapData GetPoincareMap(PlaneEquation planeEquation, std::vector<Eigen::VectorXld> trajectory)
	{
		/*Use Eigen library for vector-matrix computation*/
		//Also you have this->trajectory for this method

		// assume 3d trajectory

		Basis3ld basis = transformBasis(Eigen::Vector3ld(planeEquation.A, planeEquation.B, planeEquation.C));

		std::vector<Eigen::Vector3ld> intersections3;
		std::vector<Eigen::Vector2ld> intersections2;

		Eigen::Vector3ld prevpoint;
		Eigen::Vector3ld point;
		Eigen::Vector3ld intersectionPoint;
		int prevsign;
		int sign;
		prevsign = SideSign(planeEquation, point);

		int N = trajectory.size();
		for (int i = 1; i < N; i++) {
			point = trajectory[i];
			sign = SideSign(planeEquation, point);
			if (sign == 0) {
				intersectionPoint = point;
				intersections3.push_back(intersectionPoint);
				intersections2.push_back(applyBasis(basis, intersectionPoint));
			}
			else if (sign != prevsign) {
				intersectionPoint = point;
				//intersections.push_back(intersectionPoint);
				///*
				intersectionPoint = intersectionCalc(planeEquation, prevpoint, point);
				//if (IsOnInterval(prevpoint, point, intersectionPoint))
				intersections3.push_back(intersectionPoint);
				intersections2.push_back(applyBasis(basis, intersectionPoint));
				//*/
			}
			prevpoint = point;
			prevsign = sign;
		}

		PoincareMapData result = PoincareMapData();
		result.intersections2D = intersections2;
		result.intersections3D = intersections3;
		return result;
	}

	// Bifurcation map
	std::vector<long double> GetBifurcationMap(std::vector<Eigen::VectorXld> trajectory)
	{
		/*Use Eigen library for vector-matrix computation*/
		//Also you have this->trajectory for this method

		// assume 3d trajectory
		PlaneEquation planeEquation;
		planeEquation.A = -1;
		planeEquation.B = -1;
		planeEquation.C = -1;
		planeEquation.D = 0;

		Basis3ld basis = transformBasis(Eigen::Vector3ld(planeEquation.A, planeEquation.B, planeEquation.C));

		//std::vector<Eigen::Vector3ld> intersections3;
		//std::vector<Eigen::Vector2ld> intersections2;
		std::vector<long double> intersections1;

		Eigen::Vector3ld prevpoint;
		Eigen::Vector3ld point;
		Eigen::Vector3ld intersectionPoint;
		Eigen::Vector2ld intersectionPoint2;
		//long double intersectionPoint1;
		int prevsign;
		int sign;
		prevsign = SideSign(planeEquation, point);

		int N = trajectory.size();
		for (int i = 1; i < N; i++) {
			point = trajectory[i];
			sign = SideSign(planeEquation, point);
			if (sign == 0) {
				intersectionPoint = point;
				//intersections3.push_back(intersectionPoint);
				intersectionPoint2 = applyBasis(basis, intersectionPoint);
				//intersections2.push_back(applyBasis(basis, intersectionPoint));
				intersections1.push_back(intersectionPoint2[0]);
			}
			else if (sign != prevsign) {
				intersectionPoint = point;
				//intersections.push_back(intersectionPoint);
				///*
				intersectionPoint = intersectionCalc(planeEquation, prevpoint, point);
				//if (IsOnInterval(prevpoint, point, intersectionPoint))
				//intersections3.push_back(intersectionPoint);
				intersectionPoint2 = applyBasis(basis, intersectionPoint);
				//intersections2.push_back(applyBasis(basis, intersectionPoint));
				intersections1.push_back(intersectionPoint2[1]);
				//*/
			}
			prevpoint = point;
			prevsign = sign;
		}

		//BifurcationMapData result = BifurcationMapData();
		//result.intersections2D = intersections2;
		//result.intersections3D = intersections3;
		//result.intersections1D = intersections1;
		return intersections1;
	}

	void TrajectoriesToObj(const std::vector<std::vector<Eigen::VectorXld>>& trajectories, Eigen::Vector3i axis_indexes, std::string name, std::string path)
	{
		Eigen::MatrixXd vertices;
		Eigen::MatrixXi faces;
		Eigen::VectorXi list_of_indexes;
		if (trajectories.size() == 0)
			throw std::exception("No trajectories");
		auto trajectory_with_min_size = *std::min_element(trajectories.begin(), trajectories.end(),
			[](const std::vector<Eigen::VectorXld>& first_trajectory, const std::vector<Eigen::VectorXld>& second_trajectory)
			{
				return first_trajectory.size() < second_trajectory.size();
			});
		size_t min_size = trajectory_with_min_size.size();
		if (min_size == 0)
			throw std::exception("There is an empty trajectory");
		if (trajectory_with_min_size[0].size() <= std::max({ axis_indexes[0], axis_indexes[1], axis_indexes[2] }))
			throw std::exception("There is an inappropriate index");
		vertices.resize(min_size * trajectories.size(), axis_indexes.size());
		for (size_t i = 0; i < trajectories.size(); i++)
			for (size_t j = 0; j < min_size; j++)
				vertices.row(i * min_size + j) = Eigen::Vector3d(trajectories[i][j][axis_indexes[0]], trajectories[i][j][axis_indexes[1]], trajectories[i][j][axis_indexes[2]]);
		faces.resize((min_size - 1) * (trajectories.size() - 1) * 2, axis_indexes.size());
		size_t triangle_number = 0;
		for (size_t i = 0; i < trajectories.size() - 1; i++)
		{
			for (size_t j = 0; j < min_size - 1; j++)
			{
				faces.row(triangle_number) = Eigen::Vector3i(j + min_size * i, j + min_size * (i + 1), j + 1 + min_size * i);
				faces.row(triangle_number + 1) = Eigen::Vector3i(j + 1 + min_size * i, j + min_size * (i + 1), j + 1 + min_size * (i + 1));
				triangle_number += 2;
			}
		}
		igl::decimate(vertices, faces, 10000, vertices, faces, list_of_indexes);
		igl::writeOBJ(path + name + ".obj", vertices, faces);
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
		variables += ",t";
		for (auto function : strings_functions)
		{
			FunctionParser_ld function_parser;
			function_parser.Parse(additional_variables + function, variables);
			this->functions.push_back(function_parser);
		}
		this->point_of_trajectory = starting_point;
		this->trajectory.push_back(starting_point);
		CalculateJacobianMatrix();
	}

	std::vector<Eigen::VectorXld> DynamicSystem::GetTrajectory(long double time)
	{
		bool is_infinity_trajectory = false;
		while (this->t <= time)
		{
			try
			{
				NextPointOfTrajectory();
				this->t += this->dt;
				this->timeSequence.push_back(this->t);
			}
			catch (InfinityTrajectoryException& ex)
			{
				this->comment.append(ex.what());
				break;
			}
		}
		return this->trajectory;
	}

	std::vector<long double> DynamicSystem::GetTimeSequence() {
		return this->timeSequence;
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
				try
				{
					NextPointOfTrajectory(true); // forced to be static dt
				}
				catch (InfinityTrajectoryException& ex)
				{
					for (size_t k = 0; k < this->dimension; k++)
						spectrum_of_lyapunov_exponents.push_back(sums_of_logarithms[k] / i / T / this->dt);
					return spectrum_of_lyapunov_exponents;
				}
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
			try {
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
					NextPointOfTrajectory(true); // forced to be static dt
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
					series_spectrum_lyapunov["t"].push_back(i * this->dt);
				}
			}
			catch (InfinityTrajectoryException& ex)
			{
				break;
			}
		}
		return series_spectrum_lyapunov;
	}

	/*For Rouol*/
	PoincareMapData DynamicSystem::GetPoincareMap(PlaneEquation planeEquation)
	{
		/*Use Eigen library for vector-matrix computation*/
		//Also you have this->trajectory for this method

		// assume 3d trajectory
		std::vector<Eigen::VectorXld> data = this->trajectory;

		Basis3ld basis = transformBasis(Eigen::Vector3ld(planeEquation.A, planeEquation.B, planeEquation.C));

		std::vector<Eigen::Vector3ld> intersections3;
		std::vector<Eigen::Vector2ld> intersections2;

		Eigen::Vector3ld prevpoint;
		Eigen::Vector3ld point;
		Eigen::Vector3ld intersectionPoint;
		int prevsign;
		int sign;
		prevsign = SideSign(planeEquation, point);

		int N = data.size();
		for (int i = 1; i < N; i++) {  
			point = trajectory[i];
			sign = SideSign(planeEquation, point);
			if (sign == 0) {
				intersectionPoint = point;
				intersections3.push_back(intersectionPoint);
				intersections2.push_back(applyBasis(basis, intersectionPoint));
			}
			else if (sign != prevsign) {
				intersectionPoint = point;
				//intersections.push_back(intersectionPoint);
				///*
				intersectionPoint = intersectionCalc(planeEquation, prevpoint, point);
				//if (IsOnInterval(prevpoint, point, intersectionPoint))
				intersections3.push_back(intersectionPoint);
				intersections2.push_back(applyBasis(basis, intersectionPoint));
				//*/
			}
			prevpoint = point;
			prevsign = sign;
		}
		
		PoincareMapData result = PoincareMapData();
		result.intersections2D = intersections2;
		result.intersections3D = intersections3;
		return result;
	}

	void DynamicSystem::SetDt(long double dt)
	{
		this->dt = dt;
	}

	void DynamicSystem::SetTime(long double time)
	{
		this->t = time;
	}

	void DynamicSystem::Reset(Eigen::VectorXld current_point)
	{
		if (current_point.size() != this->dimension)
			throw std::exception("current_point has invalid size");
		this->trajectory.clear();
		this->timeSequence.clear();
		SetTime(0);
		SetCurrentPointOfTrajectory(current_point);
	}

	void DynamicSystem::ResetWithTime(Eigen::VectorXld current_point_with_time)
	{
		if (current_point_with_time.size() != this->dimension + 1)
			throw std::exception("current_point_with_time has invalid size");
		this->trajectory.clear();
		this->timeSequence.clear();
		SetTime(current_point_with_time[0]);
		Eigen::VectorXld current_point(this->dimension);
		for (size_t i = 0; i < current_point.size(); i++)
			current_point[i] = current_point_with_time[i + 1];
		SetCurrentPointOfTrajectory(current_point);
	}

	void DynamicSystem::SetCurrentPointOfTrajectory(Eigen::VectorXld current_point)
	{
		this->point_of_trajectory = current_point;
		this->trajectory.clear();
		this->trajectory.push_back(this->point_of_trajectory);
	}

	std::string DynamicSystem::GetErrorComment()
	{
		return this->comment;
	}

	PartialDifferentialEquation::PartialDifferentialEquation(
		const std::vector<std::string>& strings_boundary_functions, 
		long double first_value_parameter, 
		long double second_value_parameter, 
		long double step_along_border, 
		const std::vector<std::string>& strings_functions_coefficients, 
		std::string variables, 
		std::string additional_variables, 
		long double dt
	)
		: first_value_parameter(first_value_parameter), second_value_parameter(second_value_parameter), step_along_border(step_along_border)
	{
		if (variables == "")//Defined variables
		{
			std::string sum_all_functions;
			FunctionParser_ld parser_for_variables;
			for (auto function : strings_functions_coefficients)
				sum_all_functions += function + '+';
			sum_all_functions.pop_back();
			parser_for_variables.ParseAndDeduceVariables(additional_variables + sum_all_functions, variables);
		}
		std::string parameter;
		std::string sum_all_functions;
		FunctionParser_ld parser_for_variables;
		for (auto function : strings_boundary_functions)
			sum_all_functions += function + '+';
		sum_all_functions.pop_back();
		parser_for_variables.ParseAndDeduceVariables(additional_variables + sum_all_functions, parameter);
		for (auto function : strings_boundary_functions)
		{
			FunctionParser_ld function_parser;
			function_parser.Parse(additional_variables + function, parameter);
			this->boundary_functions.push_back(function_parser);
		}
		Eigen::VectorXld starting_point;
		std::vector<std::string> strings_functions_dynamic_system;
		this->with_time = strings_functions_coefficients[0] == "0" ? false : true;
		if (!this->with_time)
		{
			starting_point = BoundaryFunctionWithTime(first_value_parameter);
			for (auto iterator_coefficient = strings_functions_coefficients.begin(); iterator_coefficient != strings_functions_coefficients.end(); iterator_coefficient++)
				strings_functions_dynamic_system.push_back(*iterator_coefficient);
		}
		else
		{
			starting_point = BoundaryFunctionWithoutTime(first_value_parameter);
			for (auto iterator_coefficient = strings_functions_coefficients.begin() + 1; iterator_coefficient != strings_functions_coefficients.end(); iterator_coefficient++)
				strings_functions_dynamic_system.push_back('(' + *iterator_coefficient + ")/(" + strings_functions_coefficients[0] + ')');
		}
		this->dynamic_system = new DynamicSystem(starting_point, strings_functions_dynamic_system, variables, additional_variables);
		this->dynamic_system->SetDt(dt);
	}

	std::vector<std::vector<Eigen::VectorXld>> PartialDifferentialEquation::GetSolution(long double time)
	{
		std::vector<std::vector<Eigen::VectorXld>> solution_surface;
		for (long double parameter = this->first_value_parameter; parameter < this->second_value_parameter; parameter += this->step_along_border)
		{
			if (this->with_time)
			{
				this->dynamic_system->ResetWithTime(BoundaryFunctionWithTime(parameter));
				solution_surface.push_back(this->dynamic_system->GetTrajectory(time));
			}
			else
			{
				this->dynamic_system->Reset(BoundaryFunctionWithTime(parameter));
				solution_surface.push_back(this->dynamic_system->GetTrajectory(time));
			}
		}
		return solution_surface;
	}

	std::vector<long double> PartialDifferentialEquation::GetTimeSequence()
	{
		return this->dynamic_system->GetTimeSequence();
	}

	//Private methods:

	Eigen::VectorXld DynamicSystem::f(const Eigen::VectorXld& vector)
	{
		Eigen::VectorXld vector_with_time = vector;
		vector_with_time.conservativeResize(vector.size() + 1);
		vector_with_time[vector.size()] = this->t;
		Eigen::VectorXld result_vector(this->dimension);
		size_t i = 0;
		for (auto& function : this->functions)
			result_vector[i++] = function.Eval(vector_with_time.data());
		return result_vector;
	}

	Eigen::VectorXld DynamicSystem::variableExplicitRungeKuttaFourthOrder(
		long double dt,
		Eigen::VectorXld point_of_trajectory
	)
	{
		Eigen::VectorXld k1, k2, k3, k4, buffer_point_of_trajectory;
		buffer_point_of_trajectory = point_of_trajectory;
		k1 = f(buffer_point_of_trajectory);
		buffer_point_of_trajectory = point_of_trajectory + k1 * dt / 2;
		k2 = f(buffer_point_of_trajectory);
		buffer_point_of_trajectory = point_of_trajectory + k2 * dt / 2;
		k3 = f(buffer_point_of_trajectory);
		buffer_point_of_trajectory = point_of_trajectory + k3 * dt;
		k4 = f(buffer_point_of_trajectory);
		return point_of_trajectory + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
	}

	Eigen::VectorXld DynamicSystem::IncrementExplicitRungeKuttaFourthOrder(
		long double dt,
		Eigen::VectorXld point_of_trajectory
	)
	{
		Eigen::VectorXld k1, k2, k3, k4, buffer_point_of_trajectory;
		buffer_point_of_trajectory = point_of_trajectory;
		k1 = f(buffer_point_of_trajectory);
		buffer_point_of_trajectory = point_of_trajectory + k1 * dt / 2;
		k2 = f(buffer_point_of_trajectory);
		buffer_point_of_trajectory = point_of_trajectory + k2 * dt / 2;
		k3 = f(buffer_point_of_trajectory);
		buffer_point_of_trajectory = point_of_trajectory + k3 * dt;
		k4 = f(buffer_point_of_trajectory);
		return dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
	}

	/*
	void DynamicSystem::AdaptiveExplicitRungeKuttaFourthOrder()
	{
		long double eps = 1e-3;
		long double adaptive_dt = this->dt;
		Eigen::VectorXld intactStep = this->variableExplicitRungeKuttaFourthOrder(adaptive_dt, this->point_of_trajectory);
		long double intactStepNorm = intactStep.norm();
		Eigen::VectorXld newStep = this->variableExplicitRungeKuttaFourthOrder(adaptive_dt * 2, this->point_of_trajectory);
		if (abs(newStep.norm() - intactStepNorm) / abs(intactStepNorm) < eps) {
			adaptive_dt *= 2;
			// try increase dt
			while (abs(newStep.norm() - intactStepNorm) / abs(intactStepNorm) < eps) {
				adaptive_dt *= 2; // Error is ok; increase step size.
				newStep = this->variableExplicitRungeKuttaFourthOrder(adaptive_dt, this->point_of_trajectory);
			}
		}
		else {
			newStep = this->variableExplicitRungeKuttaFourthOrder(adaptive_dt / 2, this->point_of_trajectory);
			while (abs(newStep.norm() - intactStepNorm) / abs(intactStepNorm) > eps) {
				adaptive_dt /= 2; // Error is too large; decrease step size.
				newStep = this->variableExplicitRungeKuttaFourthOrder(adaptive_dt, this->point_of_trajectory);
			}
		}
		this->point_of_trajectory = newStep;
		this->trajectory.push_back(this->point_of_trajectory);
	}
	*/
	/*
	void DynamicSystem::AdaptiveExplicitRungeKuttaFourthOrder()
	{
		//long double adaptive_dt = this->dt;
		long double intactStepNorm = 0;
		Eigen::VectorXld intactStep = this->variableExplicitRungeKuttaFourthOrder(this->dt, this->point_of_trajectory);
		while (abs(intactStep.norm() - intactStepNorm) > 1e-2) {
			this->dt /= 2; // Error is too large; decrease step size.
			intactStepNorm = intactStep.norm();
			intactStep = this->variableExplicitRungeKuttaFourthOrder(this->dt, this->point_of_trajectory);
		}
		//this->dt = adaptive_dt;
		this->point_of_trajectory = intactStep;
		this->trajectory.push_back(this->point_of_trajectory);
	}
	*/
	void DynamicSystem::FixedVExplicitRungeKuttaFourthOrder()
	{
		long double fixedStep = 1e-2;
		Eigen::VectorXld intactStep = this->IncrementExplicitRungeKuttaFourthOrder(this->dt, this->point_of_trajectory);
		this->dt *= fixedStep / intactStep.norm();
		this->ExplicitRungeKuttaFourthOrder();
	}

	const long double powl24 = powl(2, 4);

	long double RichardsonExtrapolation4Error(long double smallerStepNorm, long double largerStepNorm) {
		return abs(abs(powl24 * smallerStepNorm - largerStepNorm) / (powl24 - 1) - abs(smallerStepNorm));
	}
	
	void DynamicSystem::AdaptiveExplicitRungeKuttaFourthOrder()
	{
		Eigen::VectorXld halfStep, intactStep; // , doubleStep;
		long double halfStepNorm, intactStepNorm; //, doubleStepNorm;
		halfStep = this->variableExplicitRungeKuttaFourthOrder(this->dt / 2, this->point_of_trajectory);
		intactStep = this->variableExplicitRungeKuttaFourthOrder(this->dt, this->point_of_trajectory);
		//doubleStep = this->variableExplicitRungeKuttaFourthOrder(this->dt * 2, this->point_of_trajectory);

		halfStepNorm = halfStep.norm();
		intactStepNorm = intactStep.norm();

		long double richardsonExtrapolation4error = RichardsonExtrapolation4Error(halfStepNorm, intactStepNorm);
		this->dt = 0.9 * (this->dt / 2) * powl((epsilon * 10000) / richardsonExtrapolation4error, 1. / 4.);
		this->ExplicitRungeKuttaFourthOrder();
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
		this->trajectory.push_back(this->point_of_trajectory);
	}

	void DynamicSystem::ImplicitEuler()
	{
		this->point_of_trajectory += (Eigen::MatrixXld::Identity(this->dimension, this->dimension) - this->jacobian_matrix * this->dt).inverse() * this->dt * f(this->point_of_trajectory);
		this->trajectory.push_back(this->point_of_trajectory);
	}

	bool DynamicSystem::IsHard()
	{
		Eigen::VectorXcld eigenvalues = this->jacobian_matrix.eigenvalues();
		long double max_eigenvalue = fabsl(eigenvalues(0).real());
		//long double min_eigenvalue = fabsl(eigenvalues(0).real());
		for (size_t i = 0; i < eigenvalues.size(); i++)
		{
			max_eigenvalue = fabsl(eigenvalues(i).real()) > max_eigenvalue ? fabsl(eigenvalues(i).real()) : max_eigenvalue;
			//min_eigenvalue = fabsl(eigenvalues(i).real()) < min_eigenvalue ? fabsl(eigenvalues(i).real()) : min_eigenvalue;
		}
		return this->dt > 1. / max_eigenvalue;
		//return max_eigenvalue / min_eigenvalue > hard_number ? true : false;
	}

	void DynamicSystem::NextPointOfTrajectory(bool ForceStaticDt)
	{
		if (this->point_of_trajectory.norm() > 1e100)
			throw InfinityTrajectoryException("Infinity trajectory");
		if (this->IsHard()/*Dynamic system is hard?*/)
		{
			//Make implementation
			switch (this->implicit_method)
			{
			case ImplicitNumericalMethod::EulerImplicit:
				ImplicitEuler();
				break;
			}
		}
		else
		{
			if (ForceStaticDt) {
				ExplicitRungeKuttaFourthOrder();
			}
			else {
				switch (this->explicit_method)
				{
					//Make implementation
					/*case ExplicitNumericalMethod::EulerExplicit:
						ExplicitEuler();
						break;*/
				case ExplicitNumericalMethod::RungeKuttaFourthOrder:
					ExplicitRungeKuttaFourthOrder();
					break;
				case ExplicitNumericalMethod::AdaptiveRungeKuttaFourthOrder:
					AdaptiveExplicitRungeKuttaFourthOrder();
					break;
				case ExplicitNumericalMethod::FixedVRungeKuttaFourthOrder:
					FixedVExplicitRungeKuttaFourthOrder();
					break;
				}
			}
		}
		CalculateJacobianMatrix();
	}

	void DynamicSystem::CalculateJacobianMatrix()
	{
		this->jacobian_matrix = Eigen::MatrixXld::Zero(this->dimension, this->dimension);
		Eigen::VectorXld point_of_trajectory_with_time(this->dimension + 1);
		for (size_t i = 0; i < this->dimension; i++)
			point_of_trajectory_with_time[i] = this->point_of_trajectory[i];
		point_of_trajectory_with_time[this->dimension] = this->t;
		for (size_t i = 0; i < this->dimension; i++)
		{
			for (size_t j = 0; j < this->dimension; j++)
			{
				Eigen::VectorXld left_point = point_of_trajectory_with_time;
				left_point(j) -= this->epsilon;
				Eigen::VectorXld right_point = point_of_trajectory_with_time;
				right_point(j) += this->epsilon;
				this->jacobian_matrix(i, j) =
					(this->functions[i].Eval(right_point.data()) -
						this->functions[i].Eval(left_point.data())) / (2 * this->epsilon);
			}
		}
	}

	Eigen::VectorXld PartialDifferentialEquation::BoundaryFunctionWithTime(long double parameter)
	{
		Eigen::VectorXld result_vector(this->boundary_functions.size());
		size_t i = 0;
		for (auto& function : this->boundary_functions)
			result_vector[i++] = function.Eval(&parameter);
		return result_vector;
	}

	Eigen::VectorXld PartialDifferentialEquation::BoundaryFunctionWithoutTime(long double parameter)
	{
		Eigen::VectorXld result_vector(this->boundary_functions.size());
		size_t i = -1;
		for (auto& function : this->boundary_functions)
		{
			if (i == -1)
			{
				i++;
				continue;
			}
			result_vector[i++] = function.Eval(&parameter);
		}
		return result_vector;
	}

	//HyperbolicPartialDifferentialEquation

	HyperbolicPartialDifferentialEquation::HyperbolicPartialDifferentialEquation(
		std::string f, std::string g, std::string phi, std::string psi,
		std::tuple<long double, long double, long double> left_coefficients, 
		std::tuple<long double, long double, long double> right_coefficients, 
		std::pair<long double, long double> space_interval, 
		long double T, long double h, long double tau) :
		left_coefficients(left_coefficients), h(h), tau(tau),
		right_coefficients(right_coefficients), space_interval(space_interval), T(T)
	{
		// parse function strings
		this->f.Parse(f, "x");
		this->g.Parse(g, "x");
		this->phi.Parse(phi, "x");
		this->psi.Parse(psi, "x");


		if (2 * this->h * std::get<1>(this->left_coefficients) 
			- 
			3 * std::get<0>(this->left_coefficients) < DBL_EPSILON)
			this->h /= 2;
		if (2 * this->h * std::get<1>(this->right_coefficients) 
			+ 
			3 * std::get<0>(this->right_coefficients) < DBL_EPSILON)
			this->h /= 2;

		if (this->h < DBL_EPSILON)
			throw(std::exception("Very small step"));

		long double first_x = this->space_interval.first + this->h;
		long double next_first_x = this->space_interval.first + 2*this->h;
		long double previous_first_x = this->space_interval.first;
		long double minimum = this->f.Eval(&first_x) * this->g.Eval(&first_x);
		long double maximum = std::powl(this->f.Eval(&first_x)*
			(this->g.Eval(&next_first_x) - this->g.Eval(&previous_first_x) / 4), 2)
			+
			std::powl(this->f.Eval(&first_x) * this->g.Eval(&first_x), 2);
		for (long double x = this->space_interval.first + 2*this->h;
			x < this->space_interval.second;
			x += this->h)
		{
			long double next_x = x + this->h;
			long double previous_x = x - this->h;
			long double right_expression = this->f.Eval(&x) * this->g.Eval(&x);
			long double left_expression = std::powl(this->f.Eval(&x) *
				(this->g.Eval(&next_x) - this->g.Eval(&previous_x) / 4), 2)
				+
				std::powl(this->f.Eval(&x) * this->g.Eval(&x), 2);
			minimum = right_expression < minimum ? right_expression : minimum;
			maximum = left_expression > maximum ? left_expression : maximum;
		}

		if (maximum < DBL_EPSILON)
			this->tau = this->h;
		else
		{
			if (this->tau * this->tau >= minimum / maximum * this->h * this->h)
			{
				if (minimum / maximum > 0)
				{
					this->tau = sqrtl(minimum / maximum) * this->h / 2;
				}
				else
				{
					throw(std::exception("This problem is not solved by this scheme."));
				}
			}
		}

		if (this->tau < DBL_EPSILON)
			throw(std::exception("Very small step"));

		this->u = Eigen::MatrixXld::Zero(this->T / this->tau + 1,
			(this->space_interval.second - this->space_interval.first) / this->h + 1);

		this->is_solved = false;
	}

	Eigen::MatrixXld HyperbolicPartialDifferentialEquation::Solution()
	{
		if(!this->is_solved)
			Solve();
		return this->u;
	}

	void HyperbolicPartialDifferentialEquation::Solve()
	{
		long double x = this->space_interval.first;
		long double next_x = this->space_interval.first + this->h;
		long double previous_x = this->space_interval.first - this->h;
		for (size_t n = 0; n < this->u.cols(); n++)
		{
			this->u(0, n) = this->phi.Eval(&x);
			this->u(1, n) = this->phi.Eval(&x) + this->tau * this->psi.Eval(&x) +
				this->f.Eval(&x) *
				((this->g.Eval(&next_x) - this->g.Eval(&previous_x)) /
					(2 * this->h) *
					(this->phi.Eval(&next_x) - this->phi.Eval(&previous_x)) /
					(2 * this->h)
					+
					this->g.Eval(&x) * (this->phi.Eval(&next_x) -
						2 * this->phi.Eval(&x) +
						this->phi.Eval(&previous_x)) / (this->h * this->h));
			x += this->h;
			next_x += this->h;
			previous_x += this->h;
		}
		for (size_t m = 2; m < this->u.rows(); m++)
		{
			x = this->space_interval.first + this->h;
			next_x = this->space_interval.first + 2 * this->h;
			previous_x = this->space_interval.first;
			for (size_t n = 1; n < this->u.cols() - 1; n++)
			{
				this->u(m, n) = 2 * this->u(m - 1, n) - this->u(m - 2, n) +
					this->tau * this->tau * this->f.Eval(&x) * (
						(this->g.Eval(&next_x) - this->g.Eval(&previous_x)) *
						(this->u(m - 1, n + 1) - this->u(m - 1, n - 1)) /
						(4 * this->h * this->h)
						+
						this->g.Eval(&x) *
						(this->u(m - 1, n + 1) -
							2 * this->u(m - 1, n) +
							this->u(m - 1, n - 1)) /
						(this->h * this->h));
				x += this->h;
				next_x += this->h;
				previous_x += this->h;
			}
			this->u(m, 0) = (std::get<0>(this->left_coefficients) *
				(this->u(m, 2) - 4 * this->u(m, 1)) +
				2 * this->h * std::get<2>(this->left_coefficients)) /
				(2 * this->h * std::get<1>(this->left_coefficients) -
					3 * std::get<0>(this->left_coefficients));
			this->u(m, this->u.cols() - 1) = (std::get<0>(this->right_coefficients) *
				(4 * this->u(m, this->u.cols() - 2) - this->u(m, this->u.cols() - 3)) +
				2 * this->h * std::get<2>(this->right_coefficients)) /
				(2 * this->h * std::get<1>(this->right_coefficients) +
					3 * std::get<0>(this->right_coefficients));
		}
		this->is_solved = true;
	}

	//SecondOrderODESolver
	SecondOrderODESolver::SecondOrderODESolver(std::vector<std::pair<std::string, std::string>> functions_string_pairs, Eigen::Matrix<double, 2, 3> boundary_coefficients, std::pair<double, double> border, size_t N)
	{
		for (auto function_pair : functions_string_pairs) {
			FunctionParser parser;
			if (parser.Parse(function_pair.first, function_pair.second) >= 0)
				throw std::exception(parser.ErrorMsg());
			this->functions_coefficients.push_back(parser);
		}
		if (border.first > border.second)
			throw std::exception("Invalid \'border\' values");
		if (N < 3)
			throw std::exception("Invalid \'N\' value");
		if ((std::pow(boundary_coefficients(0, 0), 2) + std::pow(boundary_coefficients(0, 1), 2)) < 0.001
			&& (std::pow(boundary_coefficients(1, 0), 2) + std::pow(boundary_coefficients(1, 1), 2)) < 0.001)
			throw std::exception("Invalid \'boundary_coefficients\' values");
		this->A = Eigen::MatrixXd::Zero(N, N);
		this->Y = Eigen::VectorXd::Zero(N);
		this->PHI = Eigen::VectorXd::Zero(N);
		this->a = border.first;
		this->b = border.second;
		this->h = (this->b - this->a) / (N - 1);
		this->boundary_coefficients = std::move(boundary_coefficients);
		this->FillMatrixForSolving();
		this->FillVectorPHI();
		this->is_solved = false;
	}

	Eigen::VectorXd SecondOrderODESolver::GetSolution()
	{
		if (!this->is_solved)
			this->Solve();
		return this->Y;
	}

	double SecondOrderODESolver::GetConditionNumber()
	{
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(this->A);
		return svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
	}

	Eigen::MatrixXd SecondOrderODESolver::GetMatrixA()
	{
		return this->A;
	}

	void SecondOrderODESolver::FillMatrixForSolving()
	{
		//Filling 1st row of the matrix
		this->A(0, 0) = this->boundary_coefficients(0, 1) - 3 * this->boundary_coefficients(0, 0) / (2 * this->h);
		this->A(0, 1) = 2 * this->boundary_coefficients(0, 0) / this->h;
		this->A(0, 2) = -this->boundary_coefficients(0, 0) / (2 * this->h);

		//Filling the matrix from 1st to penultimate rows
		for (size_t i = 1; i < this->A.rows() - 1; i++)
		{
			double x = i * this->h + this->a;
			this->A(i, i - 1) = 1 / std::pow(this->h, 2) - this->functions_coefficients[0].Eval(&x) / (2 * this->h);
			if (this->functions_coefficients[0].EvalError())
				throw std::exception("Error \'Eval\' call");
			this->A(i, i) = this->functions_coefficients[1].Eval(&x) - 2 / std::pow(this->h, 2);
			if (this->functions_coefficients[1].EvalError())
				throw std::exception("Error \'Eval\' call");
			this->A(i, i + 1) = 1 / std::pow(this->h, 2) + this->functions_coefficients[0].Eval(&x) / (2 * this->h);
			if (this->functions_coefficients[0].EvalError())
				throw std::exception("Error \'Eval\' call");
		}

		//Filling last row of the matrix
		this->A(this->A.rows() - 1, this->A.rows() - 3) = this->boundary_coefficients(1, 0) / (2 * this->h);
		this->A(this->A.rows() - 1, this->A.rows() - 2) = -2 * this->boundary_coefficients(1, 0) / this->h;
		this->A(this->A.rows() - 1, this->A.rows() - 1) = this->boundary_coefficients(1, 1) + 3 * this->boundary_coefficients(1, 0) / (2 * this->h);
	}

	void SecondOrderODESolver::FillVectorPHI()
	{
		//Set 1st element of the vector
		this->PHI(0) = this->boundary_coefficients(0, 2);

		//Set the elements of the vector from 1st to penultimate
		for (size_t i = 1; i < this->PHI.size() - 1; i++)
		{
			double x = i * this->h + this->a;
			this->PHI(i) = this->functions_coefficients[2].Eval(&x);
			if (this->functions_coefficients[2].EvalError())
				throw std::exception("Error \'Eval\' call");
		}

		//Set last element of the vector
		this->PHI(this->PHI.size() - 1) = this->boundary_coefficients(1, 2);
	}

	void SecondOrderODESolver::Solve()
	{
		//Solve with colPivHouseholderQr decomposition
		this->Y = this->A.colPivHouseholderQr().solve(this->PHI);
		this->is_solved = true;
	}

}

/*
//#ifndef _DEBUG
// Needed for export to Python
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(dyns, module_handle) {
	module_handle.doc() = "Numerical Solvers library. Created in MEPhI in 2022";
	
	py::class_<DynS::HyperbolicPartialDifferentialEquation>(
		module_handle, "HyperbolicPartialDifferentialEquation"
		).def(py::init<std::string, std::string, std::string, std::string,
			std::tuple<long double, long double, long double>,
			std::tuple<long double, long double, long double>,
			std::pair<long double, long double>,
			long double, long double, long double>())
		.def("Solution", &DynS::HyperbolicPartialDifferentialEquation::Solution, "Returns solution by an explicit second-order method");
	py::class_<DynS::SecondOrderODESolver>(
		module_handle, "SecondOrderODESolver"
		).def(py::init<std::vector<std::pair<std::string, std::string>>, Eigen::Matrix<double, 2, 3>, std::pair<double, double>, size_t>())
		.def("GetConditionNumber", &DynS::SecondOrderODESolver::GetConditionNumber, "Returns condition number of matrix for solving ODE")
		.def("GetSolution", &DynS::SecondOrderODESolver::GetSolution, "Returns Y - approximation of the solution on a grid")
		.def("GetMatrixForSolution", &DynS::SecondOrderODESolver::GetMatrixA, "Returns matrix for solving ODE");

}
//#endif // !_DEBUG
*/