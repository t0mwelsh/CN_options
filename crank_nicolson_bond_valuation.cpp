/* Program uses Crank-Nicolson methods to find the value of a coupon bond. Matrix equations created by the
Crank-Nicolson method are solved by both SOR and Thomas algorithms In this problem,
the bond paid out a continuous coupon at the rate of C^{e−αt} for constants C = 2.91 and α = 0.01.
The risk-neutral process followed by a stochastic interest rate r is given by: dr = κ(θeµt − r)dt + σrβ dW,
where W(t) denotes the Wiener process and the following constants are set at κ = 0.04227, θ = 0.0319,
µ = 0.0036, σ = 0.241 and β = 0.527. It then follows that the market value of the coupon bond B(r, t; T) 
satisfies the following PDE: ∂B/∂t + 1/2σ^2r^{2β} ∂^2B/∂r^2 + κ(θe^{µt} − r)∂B/∂r − rB + Ce{−αt} = 0;
the domain of this problem is r ∈ [0, ∞) and 0 ≤ t < T. Furthermore, B(r, t; T) then satisfies boundary
conditions B(r, t = T; T) = F; ∂B/∂t + κθ^e{µt}∂B/∂r + Ce{−αt} = 0 at r = 0; B(r, t; T) → 0 as r → ∞,
with the latter known as a Dirichlet condition. The following boundary condition, known as a Neumann condition, 
is also appropriate: ∂B/∂r → as r → ∞. */

#include <iostream>
#include <iomanip> 
#include <cmath>
#include <vector>

void SOR(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d,
	std::vector<double>& B_new, int iter_max, int jmax, double tol, double omega)
// SOR algorithm that solves matrix equations of the form Ax = b, where A is a known matrix, b is a known vector
// and x is an unknown vector. omega is a parameter that can be varied to increase the speed of the SOR algorithm.
{
	int sor;
	for (sor = 0; sor < iter_max; sor++)
	{
		double error = 0.; // calculate residual
		// SOR equations in here
		{
			double y = (d[0] - c[0] * B_new[1]) / b[0];
			y = B_new[0] + omega * (y - B_new[0]);
			error += fabs(B_new[0] - y);
			B_new[0] = y;
		}
		for (int j = 1; j < jmax; j++)
		{
			double y = (d[j] - a[j] * B_new[j - 1] - c[j] * B_new[j + 1]) / b[j];
			y = B_new[j] + omega * (y - B_new[j]);
			error += fabs(B_new[j] - y);
			B_new[j] = y;
		}
		{
			double y = (d[jmax] - a[jmax] * B_new[jmax - 1]) / b[jmax];
			y = B_new[jmax] + omega * (y - B_new[jmax]);
			error += fabs(B_new[jmax] - y);
			B_new[jmax] = y;
		}

		// check for convergence and exit loop if converged
		if (error < tol)
			break;
	}
	if (sor == iter_max)
		std::cout << "\n NOT CONVERGED \n";
}

std::vector<double> thomas_solver(const std::vector<double>& a, const std::vector<double>& b_,
	const std::vector<double>& c, std::vector<double>& d)
// Thomas algorithm that solves matrix equations of the form Ax = b, where A is a known matrix, b is a known vector and x is an unknown vector
{
	int n = a.size();
	std::vector<double> b(n), temp(n);
	// initial first value of b
	b[0] = b_[0];
	for (int j = 1; j < n; j++)
	{
		b[j] = b_[j] - c[j - 1] * a[j] / b[j - 1];
		d[j] = d[j] - d[j - 1] * a[j] / b[j - 1];
	}
	// calculate solution
	temp[n - 1] = d[n - 1] / b[n - 1];
	for (int j = n - 2; j >= 0; j--)
		temp[j] = (d[j] - c[j] * temp[j + 1]) / b[j];
	return temp;
}

std::vector<std::vector<double>> Crank_Nicolson_setup(std::vector<double> params, int jmax, bool deriv_BC, 
	std::vector<double> B_old, double dt, double dr, int i)
// Sets up matrix equation values for the specific PDEs
{
	// declare vectors for matrix equations
	std::vector<double> a(jmax + 1), b(jmax + 1), c(jmax + 1), d(jmax + 1);
	// set up matrix equations
	a[0] = 0.;
	b[0] = 1. / dt + 1. / dr * params[4] * params[2] * exp(params[5] * (i + 0.5) * dt);
	c[0] = -1. / dr * params[4] * params[2] * exp(params[5] * (i + 0.5) * dt);
	d[0] = 1. / dt * B_old[0] + params[6] * exp(-params[7] * (i + 0.5) * dt);
	for (int j = 1; j <= jmax - 1; j++)
	{
		a[j] = 0.25 * (params[9] * params[9] * j * j * pow(j * dr, 2 * params[8] - 2) -
			1. / dr * params[4] * params[2] * exp(params[5] * (i + 0.5) * dt) + params[4] * j);
		b[j] = -1. / dt - 0.5 * params[9] * params[9] * j * j * pow(j * dr, 2 * params[8] - 2) - j * dr / 2.;
		c[j] = 0.25 * (params[9] * params[9] * j * j * pow(j * dr, 2 * params[8] - 2) +
			1. / dr * params[4] * params[2] * exp(params[5] * (i + 0.5) * dt) - params[4] * j);
		d[j] = -a[j] * B_old[j - 1] - (b[j] + 2. / dt) * B_old[j] - c[j] * B_old[j + 1] -
			params[6] * exp(-params[7] * (i + 0.5) * dt);
	}
	if (deriv_BC) {
		a[jmax] = -1. / dr;
		b[jmax] = 1. / dr;
		d[jmax] = 0.;
	}
	else {
		a[jmax] = 0.;
		b[jmax] = 1.;
		d[jmax] = 0.;
	}
	
	std::vector<std::vector<double>> result{ a, b, c, d };
	return result;
}

double interpolator(std::vector<double> B_final, std::vector<double> r, std::vector<double> params, double dr)
// linearly interpolates a vector to give an approximation for a value that could be between data points
{
	double bond_value;
	{
		int jStar = params[3] / dr;
		double sum = 0.;
		sum += (params[3] - r[jStar]) / dr * B_final[jStar + 1];
		sum += (r[jStar + 1] - params[3]) / dr * B_final[jStar];
		bond_value = sum;
	}
	return bond_value;
}

double Crank_Nicolson_method(std::vector<double> params, int iMax, int jmax,
	double rmax, double omega, int iter_max, double tol, bool deriv_BC, bool use_SOR)
// Crank Nicolson scheme that utilises SOR to solve the matrix equation and gives a value of the bond for a chosen r value
{
	// declare and initialise local variables (ds,dt) 
	double dr = rmax / jmax;
	double dt = params[0] / iMax;
	// create storage for the stock price and bond price (old and new)
	std::vector<double> r(jmax + 1), B_old(jmax + 1), B_new(jmax + 1);
	// setup and initialise the stock price 
	for (int j = 0; j <= jmax; j++)
	{
		r[j] = j * dr;
	}
	// setup and initialise the final conditions on the bond price 
	for (int j = 0; j <= jmax; j++)
	{
		B_old[j] = params[1];
		B_new[j] = params[1];
	}
	// start looping through time levels
	for (int i = iMax - 1; i >= 0; i--)
	{
		// declare vectors for matrix equations
		std::vector<std::vector<double>> matrix_vectors = Crank_Nicolson_setup(params, jmax, deriv_BC, B_old, dt, dr, i);

		// Solve matrix equations with SOR or thomas_solver using the individual vectors
		if (use_SOR) {
			SOR(matrix_vectors[0], matrix_vectors[1], matrix_vectors[2], matrix_vectors[3], B_new, iter_max, jmax, tol, omega);
		}
		else {
			B_new = thomas_solver(matrix_vectors[0], matrix_vectors[1], matrix_vectors[2], matrix_vectors[3]);
		}
		
		// set old=new 
		B_old = B_new;
	}
	// finish looping through time levels

	// output the estimated bond price
	return interpolator(B_new, r, params, dr);
}


int main()
{
	// declare and initialise Black Scholes parameters
	double T = 5, F = 86, theta = 0.0319, r0 = 0.0412, kappa = 0.04227,
		mu = 0.0036, C = 2.91, alpha = 0.01, beta = 0.527, sigma = 0.241;
	std::vector<double> params = { T, F, theta, r0, kappa, mu, C, alpha, beta, sigma };
	std::cout << std::setprecision(10);

	std::cout << Crank_Nicolson_method(params, 100, 100, 2, 1.43, 10000, 1.e-8, true, true) << std::endl;
	std::cout << Crank_Nicolson_method(params, 100, 100, 2, 1.43, 10000, 1.e-8, true, false) << std::endl;

	return 0;
}