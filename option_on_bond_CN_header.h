#pragma once
#ifndef Header_h
#define Header_h

#include <iostream>
#include <iomanip> 
#include <cmath>
#include <vector>

void SOR(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d,
	std::vector<double>& B_new, int iter_max, int jmax, double tol, double omega);
	// SOR algorithm that solves matrix equations of the form Ax = b, where A is a known matrix, b is a known vector
	// and x is an unknown vector. omega is a parameter that can be varied to increase the speed of the SOR algorithm.

void american_SOR(const std::vector<double>& aB, const std::vector<double>& bB, const std::vector<double>& cB, const std::vector<double>& dB,
	const std::vector<double>& aV, const std::vector<double>& bV, const std::vector<double>& cV, const std::vector<double>& dV,
	std::vector<double>& B_new, std::vector<double>& V_new, int iter_max, int jmax, double tol, double omega, std::vector<double> params);
	// SOR algorithm that solves matrix equations of the form Ax = b, where A is a known matrix, b is a known vector
	// and x is an unknown vector. omega is a parameter that can be varied to increase the speed of the SOR algorithm.
	// Incorporates an extra condition into the algorithm due to the ability to exercise an American option whenever.

std::vector<double> thomas_solver(const std::vector<double>& a, const std::vector<double>& b_,
	const std::vector<double>& c, std::vector<double>& d);
	// Thomas algorithm that solves matrix equations of the form Ax = b, where A is a known matrix, b is a known vector and x is an unknown vector

void penalty_method(const std::vector<double>& aB, const std::vector<double>& bB, const std::vector<double>& cB, std::vector<double>& dB,
	const std::vector<double>& aV, const std::vector<double>& bV, const std::vector<double>& cV, const std::vector<double>& dV,
	std::vector<double>& B_new, std::vector<double>& V_new, double rho, int jmax, int iter_max, double tol, std::vector<double> params);
	// Penalty Method solves the matrix equation using a Thomas solver while also simulating the ability to exercise an option at any time

std::vector<std::vector<double>> Crank_Nicolson_bond_setup(std::vector<double> params, int jmax, bool deriv_BC,
	std::vector<double> B_old, double dt, double dr, double i);
	// Sets up matrix equation values for the specific PDEs of the bond

std::vector<std::vector<double>> Crank_Nicolson_european_option_setup(std::vector<double> params, int jmax, bool deriv_BC,
	std::vector<double> V_old, double dt, double dr, double i);
	// Sets up matrix equation values for the specific PDEs of the European option

std::vector<std::vector<double>> Crank_Nicolson_american_option_setup(std::vector<double> params, int jmax, bool deriv_BC,
	std::vector<double> B_old, std::vector<double> V_old, double dt, double dr, double i);
	// Sets up matrix equation values for the specific PDEs of the American option

double interpolator(std::vector<double> input, std::vector<double> r, std::vector<double> params, double dr);
	// linearly interpolates a vector to give an approximation for a value that could be between data points

#endif