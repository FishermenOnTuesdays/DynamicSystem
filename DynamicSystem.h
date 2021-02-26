#pragma once
#include <iostream>
#include <functional>
#include <vector>
#include <map>
#include <omp.h>
#include <string>
#include "fparser.hh"
#include "Eigen/Dense"
#include "unsupported/Eigen/MatrixFunctions"

namespace Eigen
{
	using MatrixXld = Eigen::Matrix<long double, Dynamic, Dynamic>;
	using VectorXld = Eigen::Matrix<long double, Dynamic, 1>;
	using VectorXcld = Eigen::Matrix<std::complex<long double>, Dynamic, 1>;
	using Vector2ld = Eigen::Matrix<long double, 2, 1>;
	using Vector3ld = Eigen::Matrix<long double, 3, 1>;
}

// needed for Poincare
typedef struct plane_equation {
	long double A, B, C, D;
} PlaneEquation;

typedef struct poincare_result {
	std::vector<Eigen::Vector3ld> intersections3D;
	std::vector<Eigen::Vector2ld> intersections2D;
} PoincareMapData;

/*
typedef struct bifurcation_result {
	std::vector<long double> intersections1D;
} BifurcationMapData;
*/

namespace DynS
{
	//Exceptions:
	class InfinityTrajectoryException : public std::exception
	{
	private:
		std::string m_error;
	public:
		InfinityTrajectoryException(std::string error) : m_error(error) {}
		const char* what() const noexcept { return m_error.c_str(); }
	};


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

	//Returns Poincare map from input trajectory
	PoincareMapData GetPoincareMap(PlaneEquation planeEquation, std::vector<Eigen::VectorXld> trajectory);

	//Returns Bifurcation map from input trajectory
	std::vector<long double> GetBifurcationMap(std::vector<Eigen::VectorXld> trajectory);

	class DynamicSystem
	{
	public:
		//Enums:
		enum class ExplicitNumericalMethod
		{
			RungeKuttaFourthOrder,
			AdaptiveRungeKuttaFourthOrder,
			FixedVRungeKuttaFourthOrder,
			EulerExplicit
		};
		enum class ImplicitNumericalMethod
		{
			EulerImplicit
		};
		//Public methods:

		//Create a dynamic system that has specified initial values and defined by the specified functions
		DynamicSystem(const Eigen::VectorXld& starting_point, const std::vector<std::string>& strings_functions, std::string variables = "", std::string additional_variables = "");

		/*
		//Create a dynamic system that has specified trajectory
		DynamicSystem(const std::vector<Eigen::VectorXld>& trajectory);
		*/

		//Returns a sequence of trajectory's points at given time
		std::vector<Eigen::VectorXld> GetTrajectory(long double time);

		//Returns a Time sequence of calculated trajectory
		std::vector<long double> GetTimeSequence();

		//Returns a spectrum of Lyapunov exponents this dynamic system
		std::vector<long double> GetSpectrumLyapunov(long double time);

		//Returns a series of Lyapunov exponents spectrum at every step
		std::map<std::string, std::vector<long double>> GetTimeSeriesSpectrumLyapunov(long double time);

		/*For Rouol*/
		//Returns Poincare map
		PoincareMapData GetPoincareMap(PlaneEquation planeEquation);

		//Set dt for this dynamic system
		void SetDt(long double dt);

		//Set current point of dynamic system trajectory
		void SetCurrentPointOfTrajectory(Eigen::VectorXld current_point);

		//Return error comment
		std::string GetErrorComment();

		//Explicit method currently used
		ExplicitNumericalMethod explicit_method = ExplicitNumericalMethod::RungeKuttaFourthOrder;

		//Implicit method currently used
		ImplicitNumericalMethod implicit_method = ImplicitNumericalMethod::EulerImplicit;

		//Private methods:
	private:
		//Vector function defining a dynamic system
		Eigen::VectorXld f(const Eigen::VectorXld& vector);

		//Calculate and return next point of trajectory dynamic system by explicit Runge-Kutta fourth-order method
		Eigen::VectorXld variableExplicitRungeKuttaFourthOrder(long double dt,
			Eigen::VectorXld point_of_trajectory
		);

		//Calculate and return increment of trajectory dynamic system by explicit Runge-Kutta fourth-order method
		Eigen::VectorXld IncrementExplicitRungeKuttaFourthOrder(
			long double dt,
			Eigen::VectorXld point_of_trajectory
		);

		//Calculate next point of trajectory dynamic system by explicit Runge-Kutta fourth-order method with Fixed Velocity step
		void FixedVExplicitRungeKuttaFourthOrder();

		//Calculate next point of trajectory dynamic system by explicit Runge-Kutta fourth-order method with Adaptive step
		void AdaptiveExplicitRungeKuttaFourthOrder();

		//Calculate next point of trajectory dynamic system by explicit Runge-Kutta fourth-order method
		void ExplicitRungeKuttaFourthOrder();

		//Calculate next point of trajectory dynamic system by explicit Runge-Kutta fourth-order method
		void ImplicitEuler();

		//Determines whether the system is hard
		bool IsHard(long double hard_number);

		//Calculate next point of dynamic system trajectory 
		void NextPointOfTrajectory(bool FORCESTATICDT = false);

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

		//Time
		long double t = 0;
		
		//Current point of dynamic system trajectory
		Eigen::VectorXld point_of_trajectory;

		//Trajectory of dynamic system
		std::vector<Eigen::VectorXld> trajectory;

		//Time Sequence
		std::vector<long double> timeSequence;

		//Jacobian matrix in the current point of dynamic system trajectory
		Eigen::MatrixXld jacobian_matrix;

		//Accuracy
		const long double epsilon = 0.0000001;

		//Error comment:
		std::string comment = "";
	};

}
