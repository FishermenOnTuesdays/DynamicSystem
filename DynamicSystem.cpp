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
		this->trajectory.clear();
		this->timeSequence.clear();
		SetTime(0);
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
		Eigen::VectorXld staring_point = BoundaryFunction(first_value_parameter);
		std::vector<std::string> strings_functions_dynamic_system;
		for (auto iterator_coefficient = strings_functions_coefficients.begin() + 1; iterator_coefficient != strings_functions_coefficients.end(); iterator_coefficient++)
			strings_functions_dynamic_system.push_back('(' + *iterator_coefficient + ")/(" + strings_functions_coefficients[0] + ')');
		this->dynamic_system = new DynamicSystem(staring_point, strings_functions_dynamic_system, variables, additional_variables);
		this->dynamic_system->SetDt(dt);
	}

	std::vector<std::vector<Eigen::VectorXld>> PartialDifferentialEquation::GetSolution(long double time)
	{
		std::vector<std::vector<Eigen::VectorXld>> solution_surface;
		for (long double parameter = this->first_value_parameter; parameter < this->second_value_parameter; parameter += this->step_along_border)
		{
			this->dynamic_system->Reset(BoundaryFunction(parameter));
			solution_surface.push_back(this->dynamic_system->GetTrajectory(time));
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
		Eigen::VectorXld point_of_trajectory_with_time = this->point_of_trajectory;
		point_of_trajectory_with_time.conservativeResize(this->point_of_trajectory.size() + 1);
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

	Eigen::VectorXld PartialDifferentialEquation::BoundaryFunction(long double parameter)
	{
		Eigen::VectorXld result_vector(this->boundary_functions.size());
		size_t i = 0;
		for (auto& function : this->boundary_functions)
			result_vector[i++] = function.Eval(&parameter);
		return result_vector;
	}
}