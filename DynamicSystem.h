//                                  
//                                   
//                                   ██████╗ ██╗   ██╗███╗   ██╗███████╗
//                                   ██╔══██╗╚██╗ ██╔╝████╗  ██║██╔════╝
//                                   ██║  ██║ ╚████╔╝ ██╔██╗ ██║███████╗
//                                   ██║  ██║  ╚██╔╝  ██║╚██╗██║╚════██║
//                                   ██████╔╝   ██║   ██║ ╚████║███████║
//                                   ╚═════╝    ╚═╝   ╚═╝  ╚═══╝╚══════╝
//                                                     

#pragma once

#include <iostream>
#include <functional>
#include <algorithm>
#include <vector>
#include <utility>
#include <map>
#include <omp.h>
#include <string>
#include <cmath>
#include "fparser.hh"
#include "Eigen/Dense"
#include "unsupported/Eigen/MatrixFunctions"
#include <igl/readOFF.h>
#include <igl/decimate.h>
#include <igl/writeOBJ.h>


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

	//Convert surface from trajectories to obj
	//trajectories - array of trajectories with n-dimensional points
	//axis_indexes - indexes of the x, y, z axes at the n-dimensional point of the trajectory
	//name - the name of the .obj file
	//path - the path to save the .obj file
	void TrajectoriesToObj(const std::vector<std::vector<Eigen::VectorXld>>& trajectories, Eigen::Vector3i axis_indexes, std::string name = "surface", std::string path = "");

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

		//Set time for this dynamic system
		void SetTime(long double time);

		//Reset dynamic system (clear trajectory and time sequence, set time to zero and set current point of trajectory)
		void Reset(Eigen::VectorXld current_point);

		//Reset dynamic system with time point on first position (clear trajectory and time sequence, set time to current time point and set current point of trajectory)
		void ResetWithTime(Eigen::VectorXld current_point_with_time);

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
		Eigen::VectorXld variableExplicitRungeKuttaFourthOrder(
			long double dt,
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
		bool IsHard();

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
	
	class PartialDifferentialEquation
	{
		//Public methods
	public:
		//Create a partial differential equation with cooficient functions and boundary equations
		PartialDifferentialEquation(
			const std::vector<std::string>& strings_boundary_functions, 
			long double first_value_parameter, 
			long double second_value_parameter, 
			long double step_along_border, 
			const std::vector<std::string>& strings_functions_coefficients, 
			std::string variables = "", 
			std::string additional_variables = "", 
			long double dt = 0.01
		);

		//Calculate a solution of current partial differential equation through desired time
		std::vector<std::vector<Eigen::VectorXld>> GetSolution(long double time);

		//Returns a Time sequence of solution
		std::vector<long double> GetTimeSequence();

		//Private methods
	private:
		//Return a value of boundary function in this coordinates with time on first position
		Eigen::VectorXld BoundaryFunctionWithTime(long double coordinates);

		//Return a value of boundary function in this coordinates without time
		Eigen::VectorXld BoundaryFunctionWithoutTime(long double coordinates);

		//Variables
	private:
		//Equavelent characteristic ODE system
		DynamicSystem* dynamic_system;
		
		//Left side of boundary curve
		long double first_value_parameter;

		//Right side of a boundary curve
		long double second_value_parameter;

		//Step along a boundary curve
		long double step_along_border;

		//System of boundary functions in parametric form
		std::vector<FunctionParser_ld> boundary_functions;

		//If the partial differential equation has time derivative it's true else it's false
		bool with_time;
	};

	class HyperbolicPartialDifferentialEquation
	{
		//Public methods
	public:
		//Create a hyperbolic partial differential equation
		//f - Function in front of compound derivative
		//g - Function in front of derivative u
		//q - Heterogeneity function
		//phi - Initial offset
		//psi - Initial velocity
		//left_coefficients - Coefficients at the left end
		//right_coefficients - Coefficients at the right end
		//space_interval - Borders of space
		//T - Simulation time
		//h - Step in space
		//tau - Step in time
		HyperbolicPartialDifferentialEquation(
			std::string f,
			std::string g,
			std::string q,
			std::string phi,	
			std::string psi,
			std::tuple<long double, long double, long double> left_coefficients,
			std::tuple<long double, long double, long double> right_coefficients,
			std::pair<long double, long double> space_interval,
			long double T,
			long double h,
			long double tau
		);

		//Get solution of this hyperbolic partial differential equation by an explicit second-order method
		Eigen::MatrixXld Solution();

		//Get x coordinates of matrix
		std::vector<long double> GetXs();

		//Get t coordinates of matrix
		std::vector<long double> GetTs();

		//Private methods
	private:
		//Solve this hyperbolic partial differential equation by an explicit second-order method
		void Solve();

		//Variables
	private:
		//Step in space
		long double h;

		//Step in time
		long double tau;

		//Borders of space
		std::pair<long double, long double> space_interval;

		//Initial offset
		FunctionParser_ld phi;

		//Initial velocity
		FunctionParser_ld psi;

		//Coefficients at the left end
		std::tuple<long double, long double, long double> left_coefficients;

		//Coefficients at the right end
		std::tuple<long double, long double, long double> right_coefficients;

		//Simulation time
		long double T;

		//Function in front of compound derivative
		FunctionParser_ld f;

		//Function in front of derivative u
		FunctionParser_ld g;

		//Heterogeneity function
		FunctionParser_ld q;

		//Matrix of solution
		Eigen::MatrixXld u;

		//Indicates the presence of a solution
		bool is_solved;

		//Offset along the matrix columns to store the value
		size_t offset_h;

		//Offset along the matrix rows to store the value
		size_t offset_tau;

		//x value vector
		std::vector<long double> xs;

		//t value vector
		std::vector<long double> ts;
	};

	class SecondOrderODESolver
	{
	public:
		//Constructor
		//Example:
		//ODE is y'' + p(x)·y' + q(x)·y = φ(x)
		//The boundary conditions are:
		//α1·y'(a) + β1·y(a) = γ1
		//α2·y'(b) + β2·y(b) = γ2 ,
		//where x∈[a,b]
		//functions_string_pairs is {{"p(x)", "x"}, {"q(x)", "x"}, {"φ(x)", "x"}}
		//boundary_coefficients is:
		//α1 β1 γ1
		//α2 β2 γ2
		//border is {a,b}
		//N is number of points of the partition for approximation
		SecondOrderODESolver(std::vector<std::pair<std::string, std::string>> functions_string_pairs, Eigen::Matrix<double, 2, 3> boundary_coefficients, std::pair<double, double> border, size_t N);

		//Return Y - approximation of the solution on a grid
		Eigen::VectorXd GetSolution();

		//Return condition number of matrix A (the description is in the comments to the private variable A)
		double GetConditionNumber();

		//Return matrix for solving ODE - A (the description is in the comments to the private variable A)
		Eigen::MatrixXd GetMatrixA();

	private:
		//Filling the matrix A (the description is in the comments to the private variable A)
		void FillMatrixForSolving();

		//Filling the vector PHI (the description is in the comments to the private variable PHI)
		void FillVectorPHI();

		//Find Y in equation A·Y = Φ
		void Solve();

	private:
		//Functions that stand as coefficients in ODE
		//Example:
		//ODE is y'' + p(x)·y' + q(x)·y = φ(x)
		//functions_coefficients[0] is p(x)
		//functions_coefficients[1] is q(x)
		//functions_coefficients[2] is φ(x)
		std::vector<FunctionParser> functions_coefficients;

		//Matrix that stores the coefficients of the boundary conditions
		//Example:
		//The boundary conditions are:
		//α1·y'(a) + β1·y(a) = γ1
		//α2·y'(b) + β2·y(b) = γ2 ,
		//where x∈[a,b]
		//boundary_coefficients is:
		//α1 β1 γ1
		//α2 β2 γ2
		Eigen::Matrix<double, 2, 3> boundary_coefficients;

		//Left border by x
		//x∈[a,b]
		double a;

		//Right border by x
		//x∈[a,b]
		double b;

		//Grid step
		double h;

		//Matrix for solving ODE
		//A·Y = Φ ,
		//where Y is solution approximation vector, Φ is approximation vector of φ(x)
		Eigen::MatrixXd A;

		//Solution approximation vector
		//A·Y = Φ ,
		//where A is matrix for solving ODE, Φ is approximation vector of φ(x)
		Eigen::VectorXd Y;

		//Approximation vector of φ(x)
		//A·Y = Φ ,
		//where A is matrix for solving ODE, Y is solution approximation vector
		Eigen::VectorXd PHI;

		//A value indicating whether the ODE has been resolved
		bool is_solved;
	};

	class ParabolicPartialDifferentialEquation
	{
		//Public methods
	public:
		//Create a parabolic partial differential equation
		//q - Function in front of compound derivative (function of the variables u, x, t)
		//k - Function in front of derivative u (function of the variables u, x, t)
		//f - Heterogeneity function (function of the variables u, x, t)
		//phi - Initial offset (function of the variable x)
		//left_coefficients - Coefficients at the left end (three functions of the variable t)
		//right_coefficients - Coefficients at the right end (three functions of the variable t)
		//space_interval - Borders of space (two real numbers)
		//T - Simulation time (positive real number)
		//h - Step in space (positive real number less then 1)
		//tau - Step in time (positive real number less then 1)
		//rarefaction_ratio_x - How many times to increase the step along the x-axis to save to the u-matrix
		//rarefaction_ratio_t - How many times to increase the step along the t-axis to save to the u-matrix
		ParabolicPartialDifferentialEquation(
			const std::string& q, 
			const std::string& k, 
			const std::string& f, 
			const std::string& phi,
			const std::vector<std::string>& left_coefficients,
			const std::vector<std::string>& right_coefficients,
			std::pair<long double, long double> space_interval,
			long double T,
			long double h,
			long double tau,
			size_t rarefaction_ratio_x,
			size_t rarefaction_ratio_t
		);

		//Get solution of this hyperbolic partial differential equation by an explicit second-order method
		const Eigen::MatrixXld& Solution();

		//Get x coordinates of matrix
		const std::vector<long double>& GetXs();

		//Get t coordinates of matrix
		const std::vector<long double>& GetTs();

		//Private methods
	private:
		//Solve this hyperbolic partial differential equation by an explicit second-order method
		void Solve();

		//Variables
	private:
		//Step in space
		long double h;

		//Step in time
		long double tau;

		//Borders of space
		std::pair<long double, long double> space_interval;

		//Initial offset
		FunctionParser_ld phi;

		//Coefficients at the left end
		std::vector<FunctionParser_ld> left_coefficients;

		//Coefficients at the right end
		std::vector<FunctionParser_ld> right_coefficients;

		//Simulation time
		long double T;

		//Function in front of compound derivative
		FunctionParser_ld q;

		//Function in front of derivative u
		FunctionParser_ld k;

		//Heterogeneity function
		FunctionParser_ld f;

		//Matrix of solution
		Eigen::MatrixXld u;

		//Indicates the presence of a solution
		bool is_solved;

		//Offset along the matrix columns to store the value
		size_t offset_h;

		//Offset along the matrix rows to store the value
		size_t offset_tau;

		//x value vector
		std::vector<long double> xs;

		//t value vector
		std::vector<long double> ts;
	};

}
