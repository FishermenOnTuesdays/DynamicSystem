#include<iostream>
#include<fstream>
#include<algorithm>
#include<vector>
#include<map>
#include<complex>
#include<cmath>
#include "fparser.hh"
#include "Eigen/Dense"
#include "eigen-3.3.7/unsupported/Eigen/MatrixFunctions"
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
};

struct OutputDataMain
{
	std::vector<Eigen::VectorXld> trajectory;
	std::map<std::string, std::vector<long double>> series_of_spectrum_lyapunov_exponents;
	std::string variables;
	std::vector<std::pair<std::pair<long double, long double>, long double>> map_lyapunov_exponents;
	long double dt;
};

void from_json(const nlohmann::json& json, InputDataMain& input_data)
{
	std::vector<long double> starting_values = json.at("start values[]").get<std::vector<long double>>();
	input_data.starting_values.resize(starting_values.size());
	for (size_t i = 0; i < starting_values.size(); i++)
		input_data.starting_values(i) = starting_values[i];
	json.at("functions[]").get_to(input_data.functions);
	json.at("variables").get_to(input_data.variables);
	json.at("additional equations").get_to(input_data.additional_equations);
	json.at("time").get_to(input_data.time);
	json.at("dt").get_to(input_data.dt);
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
	long double time = 0;
	for (const auto& point : output_data.trajectory)
	{
		for (size_t i = 0; i < variables.size(); i++)
			trajectory[variables[i]].push_back(point(i));
		trajectory["t"].push_back(time);
		time += output_data.dt;
	}
	json = nlohmann::json
	{ 
		{"trajectory", trajectory}, 
		{"series of spectrum lyapunov exponents", output_data.series_of_spectrum_lyapunov_exponents}, 
		{"map_lyapunov_exponents", output_data.map_lyapunov_exponents}
	};
}

nlohmann::json Main(nlohmann::json& input_json)
{
	InputDataMain input_data = input_json;
	OutputDataMain output_data{};
	DynS::DynamicSystem dynamic_system{ input_data.starting_values, input_data.functions, input_data.variables, input_data.additional_equations };
	dynamic_system.SetDt(input_data.dt);
	output_data.trajectory = dynamic_system.GetTrajectory(input_data.time);
	output_data.series_of_spectrum_lyapunov_exponents = dynamic_system.GetTimeSeriesSpectrumLyapunov(input_data.time);
	output_data.variables = input_data.variables;
	output_data.dt = input_data.dt;
	return nlohmann::json{ output_data };
}

nlohmann::json LyapunovMap(nlohmann::json& input_json)
{
	InputDataMain input_data = input_json;
	//OutputDataMain output_data{};
	/*output_data.map_lyapunov_exponents =*/ return nlohmann::json{ DynS::GetMapLyapunovExponents(input_data.starting_values, input_data.functions, input_data.variables, input_data.additional_equations, input_data.parameters, input_data.ranges, input_data.steps, input_data.time, input_data.time, input_data.dt) };
	//return nlohmann::json{ output_data };
}

nlohmann::json PoincareMap(nlohmann::json& input_json)
{
	InputDataMain input_data = input_json;
	return nlohmann::json{ DynS::GetMapLyapunovExponents(input_data.starting_values, input_data.functions, input_data.variables, input_data.additional_equations, input_data.parameters, input_data.ranges, input_data.steps, input_data.time, input_data.time, input_data.dt) };
}

nlohmann::json Bifurcation(nlohmann::json& input_json)
{
	InputDataMain input_data = input_json;
	size_t number_of_trajectories = (input_data.range.second - input_data.range.first) / input_data.step;
	std::vector<std::map<std::string, std::vector<long double>>> trajectories(number_of_trajectories);
	std::vector<long double> parameter_values(number_of_trajectories);
	//#pragma omp parallel for
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
		std::string temp_variables = input_data.variables;
		std::replace(temp_variables.begin(), temp_variables.end(), ',', ' ');
		std::istringstream iss(temp_variables);
		std::vector<std::string> variables(std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>());
		for (auto variable : variables)
			trajectories[i].emplace(variable, std::vector<long double>{});
		trajectories[i].emplace("t", std::vector<long double>{});
		long double time = 0;
		for (const auto& point : trajectory)
		{
			for (size_t j = 0; j < variables.size(); j++)
				trajectories[i][variables[j]].push_back(point(j));
			trajectories[i]["t"].push_back(time);
			time += input_data.dt;
		}
		parameter_values[i] = parameter;
	}
	//OutputDataMain output_data{};
	/*output_data.map_lyapunov_exponents =*/ return nlohmann::json{ {"parameter", input_data.parameter}, {"parameter_values", parameter_values}, {"trajectories", trajectories} };
	//return nlohmann::json{ output_data };
}

int main()
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
	}
	std::cout << output_json;
}