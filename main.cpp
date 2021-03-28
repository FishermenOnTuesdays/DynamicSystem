#include<iostream>
#include<fstream>
#include<algorithm>
#include<vector>
#include<map>
#include<complex>
#include<cmath>
#include "fparser.hh"
#include "Eigen/Dense"
#include "unsupported/Eigen/MatrixFunctions"
#include "DynamicSystem.h"
#include <nlohmann/json.hpp>

enum class InputData
{
	Main,
	LyapunovMap,
	Bifurcation,
	PoincareMap
};

struct InputDataMain
{
	Eigen::VectorXld starting_values;
	std::vector<std::string> functions;
	std::string variables;
	std::string additional_equations;
	std::pair<std::string, std::string> parameters;
	std::pair<std::pair<long double, long double>, std::pair<long double, long double>> ranges;
	std::pair<long double, long double> steps;
	std::string parameter;
	std::pair<long double, long double> range;
	long double step;
	long double time;
	long double dt;
	std::vector<Eigen::VectorXld> trajectory;
	std::vector<long double> planeEquation;
	int ExplicitNumericalMethodCode;
};

struct OutputDataMain
{
	std::vector<Eigen::VectorXld> trajectory;
	std::map<std::string, std::vector<long double>> series_of_spectrum_lyapunov_exponents;
	std::string variables;
	std::vector<std::pair<std::pair<long double, long double>, long double>> map_lyapunov_exponents;
	std::vector<Eigen::Vector3ld> intersections3D;
	std::vector<Eigen::Vector2ld> intersections2D;
	long double dt;
	std::string comment;
	std::vector<long double> timeSequence;
};

void from_json(const nlohmann::json& json, InputDataMain& input_data)
{
	try {
		std::vector<long double> starting_values = json.at("start values[]").get<std::vector<long double>>();
		input_data.starting_values.resize(starting_values.size());
		for (size_t i = 0; i < starting_values.size(); i++)
			input_data.starting_values(i) = starting_values[i];
		json.at("functions[]").get_to(input_data.functions);
		json.at("variables").get_to(input_data.variables);
		json.at("additional equations").get_to(input_data.additional_equations);
		json.at("time").get_to(input_data.time);
		json.at("dt").get_to(input_data.dt);
		json.at("ExplicitNumericalMethodCode").get_to(input_data.ExplicitNumericalMethodCode);
	}
	catch (nlohmann::json::out_of_range & ex) {}
	try {
		// trajectory
		std::vector<std::vector<long double>> trajectory = json.at("trajectory[]").get<std::vector<std::vector<long double>>>();
		input_data.trajectory.resize(trajectory.size());
		for (size_t i = 0; i < trajectory.size(); i++)
		{
			input_data.trajectory[i].resize(trajectory[i].size());
			for (size_t j = 0; j < trajectory[i].size(); j++)
				input_data.trajectory[i](j) = trajectory[i][j];
		}
	}
	catch (nlohmann::json::out_of_range & ex) {}
	try { json.at("plane equation[]").get_to(input_data.planeEquation); }
	catch (nlohmann::json::out_of_range& ex) {}
	try { json.at("parameters[]").get_to(input_data.parameters); } 
	catch (nlohmann::json::out_of_range& ex) {}
	try { json.at("ranges[]").get_to(input_data.ranges); }
	catch (nlohmann::json::out_of_range& ex) {}
	try { json.at("steps[]").get_to(input_data.steps); }
	catch (nlohmann::json::out_of_range& ex) {}
	try { json.at("parameter").get_to(input_data.parameter); }
	catch (nlohmann::json::out_of_range& ex) {}
	try { json.at("range[]").get_to(input_data.range); }
	catch (nlohmann::json::out_of_range& ex) {}
	try { json.at("step").get_to(input_data.step); }
	catch (nlohmann::json::out_of_range& ex) {}
}

void to_json(nlohmann::json& json, const OutputDataMain& output_data)
{
	std::map<std::string, std::vector<long double>> trajectory;
	std::string temp_variables = output_data.variables;
	std::replace(temp_variables.begin(), temp_variables.end(), ',', ' ');
	std::istringstream iss(temp_variables);
	std::vector<std::string> variables(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
	for (auto variable : variables)
		trajectory.emplace(variable, std::vector<long double>{});
	trajectory.emplace("t", std::vector<long double>{});
	//long double time = 0;
	int trajlen = output_data.trajectory.size();
	if (trajlen > 0) {
		int counter = 0;
		const int maxtrajlen = 100000;
		int step = floor(trajlen / maxtrajlen);
		if (step < 1) step = 1;
		for (int counter = 0; counter < output_data.trajectory.size() - 1; counter += step) {
			Eigen::VectorXld point = output_data.trajectory[counter];
			for (size_t i = 0; i < variables.size(); i++)
				trajectory[variables[i]].push_back(point(i));
			trajectory["t"].push_back(output_data.timeSequence[counter]);
		}
	}
	// intersections3D
	std::vector<std::vector<long double>> intersections3D;
	for (size_t i = 0; i < output_data.intersections3D.size(); i++)
	{
		intersections3D.push_back({});
		for (size_t j = 0; j < output_data.intersections3D[i].size(); j++)
			intersections3D[i].push_back(output_data.intersections3D[i][j]);
	}
	// intersections2D
	std::vector<std::vector<long double>> intersections2D;
	for (size_t i = 0; i < output_data.intersections2D.size(); i++)
	{
		intersections2D.push_back({});
		for (size_t j = 0; j < output_data.intersections2D[i].size(); j++)
			intersections2D[i].push_back(output_data.intersections2D[i][j]);
	}
	json = nlohmann::json
	{ 
		{"trajectory", trajectory},
		{"intersections3D", intersections3D},
		{"intersections2D", intersections2D},
		{"series of spectrum lyapunov exponents", output_data.series_of_spectrum_lyapunov_exponents}, 
		{"map_lyapunov_exponents", output_data.map_lyapunov_exponents},
		{"comment", output_data.comment}
	};
}

nlohmann::json Main(nlohmann::json& input_json)
{
	InputDataMain input_data = input_json;
	OutputDataMain output_data{};
	DynS::DynamicSystem dynamic_system{ input_data.starting_values, input_data.functions, input_data.variables, input_data.additional_equations };
	dynamic_system.SetDt(input_data.dt);
	switch (input_data.ExplicitNumericalMethodCode)
	{
		case 0:
			dynamic_system.explicit_method = DynS::DynamicSystem::ExplicitNumericalMethod::RungeKuttaFourthOrder;
			break;
		case 1:
			dynamic_system.explicit_method = DynS::DynamicSystem::ExplicitNumericalMethod::AdaptiveRungeKuttaFourthOrder;
			break;
		case 2:
			dynamic_system.explicit_method = DynS::DynamicSystem::ExplicitNumericalMethod::EulerExplicit;
			break;
		default:
			break;
	}
	output_data.trajectory = dynamic_system.GetTrajectory(input_data.time);
	output_data.timeSequence = dynamic_system.GetTimeSequence();
	output_data.comment = dynamic_system.GetErrorComment();
	if(output_data.comment == "Infinity trajectory")
		dynamic_system.SetCurrentPointOfTrajectory(input_data.starting_values);
	output_data.series_of_spectrum_lyapunov_exponents = dynamic_system.GetTimeSeriesSpectrumLyapunov(input_data.time);
	output_data.variables = input_data.variables;
	output_data.dt = input_data.dt;
	return nlohmann::json{ output_data };
}

nlohmann::json PoincareMap(nlohmann::json& input_json)
{
	InputDataMain input_data = input_json;
	PlaneEquation planeEquation;
	planeEquation.A = input_data.planeEquation[0];
	planeEquation.B = input_data.planeEquation[1];
	planeEquation.C = input_data.planeEquation[2];
	planeEquation.D = input_data.planeEquation[3];
	std::vector<Eigen::VectorXld> trajectory = input_data.trajectory;
	PoincareMapData result = DynS::GetPoincareMap(planeEquation, trajectory);
	OutputDataMain output_data{};
	output_data.intersections2D = result.intersections2D;
	output_data.intersections3D = result.intersections3D;
	return nlohmann::json{ output_data };
}

nlohmann::json LyapunovMap(nlohmann::json& input_json)
{
	InputDataMain input_data = input_json;
	return nlohmann::json{ DynS::GetMapLyapunovExponents(input_data.starting_values, input_data.functions, input_data.variables, input_data.additional_equations, input_data.parameters, input_data.ranges, input_data.steps, input_data.time, input_data.time, input_data.dt) };
}

nlohmann::json Bifurcation(nlohmann::json& input_json)
{
	InputDataMain input_data = input_json;
	size_t number_of_trajectories = (input_data.range.second - input_data.range.first) / input_data.step;
	//std::vector<std::map<std::string, std::vector<long double>>> trajectories(number_of_trajectories);
	std::vector<std::vector<long double>> BifurcationMap(number_of_trajectories);
	//std::vector<long double> parameter_values(number_of_trajectories);
	#pragma omp parallel for
	for (int i = 0; i < number_of_trajectories; i++)
	{
		long double parameter = input_data.range.first + i * input_data.step;
		DynS::DynamicSystem dynamic_system{
				input_data.starting_values,
				input_data.functions,
				input_data.variables,
				input_data.additional_equations +
				input_data.parameter + ":=" + std::to_string(parameter) + ";"
		};
		dynamic_system.SetDt(input_data.dt);
		auto trajectory = dynamic_system.GetTrajectory(input_data.time);
		/*
		std::string temp_variables = input_data.variables;
		std::replace(temp_variables.begin(), temp_variables.end(), ',', ' ');
		std::istringstream iss(temp_variables);
		std::vector<std::string> variables(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
		for (auto variable : variables)
			trajectories[i].emplace(variable, std::vector<long double>{});
		trajectories[i].emplace("t", std::vector<long double>{});
		*/

		// make Bifurcation Map
		BifurcationMap[i] = DynS::GetBifurcationMap(trajectory);

		/*
		long double time = 0;
		for (const auto& point : trajectory)
		{
			for (size_t j = 0; j < variables.size(); j++)
				trajectories[i][variables[j]].push_back(point(j));
			trajectories[i]["t"].push_back(time);
			time += input_data.dt;
		}
		parameter_values[i] = parameter;
		*/
	}
	//OutputDataMain output_data{};
	/*output_data.map_lyapunov_exponents =*/ //return nlohmann::json{ {"parameter", input_data.parameter}, {"parameter_values", parameter_values}, {"trajectories", trajectories} };
	//return nlohmann::json{ output_data };
	return nlohmann::json{ {"BifurcationMap", BifurcationMap} };
}

int main()
{
	try
	{
		nlohmann::json input_json{};
		std::cin >> input_json;
		nlohmann::json output_json{};
		switch (input_json.at("request type").get<InputData>())
		{
		case InputData::Main:
			output_json = Main(input_json);
			break;
		case InputData::LyapunovMap:
			output_json = LyapunovMap(input_json);
			break;
		case InputData::Bifurcation:
			output_json = Bifurcation(input_json);
			break;
		case InputData::PoincareMap:
			output_json = PoincareMap(input_json);
			break;
		}
		std::cout << output_json;
	}
	catch(std::exception& ex)
	{
		std::cout << "Error:" << ex.what();
	}
}