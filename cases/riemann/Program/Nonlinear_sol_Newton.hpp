#ifndef NONLINEAR_SHOCK_HPP
#define NONLINEAR_SHOCK_HPP
#pragma once

#include <iostream>
#include <cmath>
#include <iomanip>
#include <tuple>

// Function declarations
double solve_shock_pressure(); // solves for the pressure ratio across the shock
double shock_velocity(double p2); // returns the shock velocity, using p2
double rarefaction_velocity(double p2); // velocity across the rarefaction, given p2

// Function to solve for u and p using Newton's method for sods
std::tuple<double, double> solve_nonlinear_shock();

#endif // NONLINEAR_SHOCK_HPP