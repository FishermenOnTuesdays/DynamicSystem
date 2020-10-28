#pragma once
#include <iostream>
#include <functional>
#include <vector>
#include <map>
#include <omp.h>
#include <string>
#include "fparser.hh"
#include "Eigen/Dense"
#include "eigen-3.3.7/unsupported/Eigen/MatrixFunctions"

namespace Eigen
{
	using MatrixXld = Eigen::Matrix<long double, Dynamic, Dynamic>;
	using VectorXld = Eigen::Matrix<long double, Dynamic, 1>;
}

namespace DynS
{
	//Other functions

	//Returns a map of Lyapunov exponents this dynamic system
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
	);

	class DynamicSystem
	{
	public:
		//Enums:
		enum class ExplicitNumericalMethod
		{
			RungeKuttaFourthOrder,
			EulerExplicit
		};
		enum class ImplicitNumericalMethod
		{
			EulerImplicit
		};
		//Public methods:

		//Create a dynamic system that has specified initial values and defined by the specified functions
		DynamicSystem(const Eigen::VectorXld& starting_point, const std::vector<std::string>& strings_functions, std::string variables = "", std::string additional_variables = "");

		//Returns a sequence of trajectory's points at given time
		std::vector<Eigen::VectorXld> GetTrajectory(long double time);

		//Returns a spectrum of Lyapunov exponents this dynamic system
		std::vector<long double> GetSpectrumLyapunov(long double time);

		//Returns a series of Lypunov exponents spectrum at every step
		std::map<std::string, std::vector<long double>> GetTimeSeriesSpectrumLyapunov(long double time);

		/*For Rouol*/
		//Returns Poincare map
		std::vector<Eigen::VectorXld> GetPoincareMap(Eigen::VectorXld normal_vector, Eigen::VectorXld point_on_plane);

		//Set dt for this dynamic system
		void SetDt(long double dt);

		//Private methods:
	private:
		//Vector function defining a dynamic system
		Eigen::VectorXld f(const Eigen::VectorXld& vector);

		//Calculate next point of trajectory dynamic system by explicit Runge-Kutta fourth-order method
		void ExplicitRungeKuttaFourthOrder();

		//Calculate next point of dynamic system trajectory 
		void NextPointOfTrajectory();

		//Calculate Jacobian matrix in current point of dynamic system trajectory 
		void CalculateJacobianMatrix();

		//Private variables:
	private:
		//Dimension of space
		const size_t dimension;

		//Functions that define a dynamic system
		std::vector<FunctionParser_ld> functions;

		//Time integration step
		long double dt = 0.01;

		//Current point of dynamic system trajectory
		Eigen::VectorXld point_of_trajectory;

		//Trajectory of dynamic system
		std::vector<Eigen::VectorXld> trajectory;

		//Explicit method currently used
		ExplicitNumericalMethod explicit_method = ExplicitNumericalMethod::RungeKuttaFourthOrder;

		//Implicit method currently used
		ImplicitNumericalMethod implicit_method = ImplicitNumericalMethod::EulerImplicit;

		//Jacobian matrix in the current point of dynamic system trajectory
		Eigen::MatrixXld jacobian_matrix;

		//Accuracy
		const long double epsilon = 0.0000001;
	};

}
