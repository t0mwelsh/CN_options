/* Program uses Crank-Nicolson methods to find the value of options on a coupon bond. The option is either
a European or American call option and is bought d at time T1 ∈ [0, T] with strike price X. Matrix equations 
created by the Crank-Nicolson method are solved by both SOR and Thomas algorithms. In this problem,
the bond paid out a continuous coupon at the rate of C^{e−αt} for constants C = 2.91 and α = 0.01.
The risk-neutral process followed by a stochastic interest rate r is given by: dr = κ(θeµt − r)dt + σrβ dW,
where W(t) denotes the Wiener process and the following constants are set at κ = 0.04227, θ = 0.0319,
µ = 0.0036, σ = 0.241 and β = 0.527. It then follows that the market value of the coupon bond B(r, t; T)
satisfies the following PDE: ∂B/∂t + 1/2σ^2r^{2β} ∂^2B/∂r^2 + κ(θe^{µt} − r)∂B/∂r − rB + Ce{−αt} = 0;
the domain of this problem is r ∈ [0, ∞) and 0 ≤ t < T. Furthermore, B(r, t; T) then satisfies boundary
conditions B(r, t = T; T) = F; ∂B/∂t + κθ^e{µt}∂B/∂r + Ce{−αt} = 0 at r = 0; B(r, t; T) → 0 as r → ∞,
with the latter known as a Dirichlet condition. The following boundary condition, known as a Neumann condition,
is also appropriate: ∂B/∂r → as r → ∞. 

Moreover,  it can be shown that V(r, t; T1, T) satisfies the following PDE:
∂V/∂t + 1/2σ^2r{2β} ∂^2V/∂r^2 + κ(θe^{µt} − r)∂V/∂r − rV = 0. If this is a European call option then V(r, t; T1, T) 
satisfies the boundary conditions V(r, t = T1; T1, T) = max{B(r, T1; T) − X, 0}; ∂V/∂t + κθe^{µt} ∂V/ ∂r = 0 at r = 0; 
V(r, t; T1, T) → 0 as r → ∞. While if this is an American call option the boundary conditions are given by
V(r, t; T1, T) ≥ max{B(r, t; T) − X, 0}, t ≤ T1; V (r, t; T1, T) = B(r, t; T) − X at r = 0;
V(r, t; T1, T) → 0 as r → ∞. */

#include "option_on_bond_CN_header.h"

std::vector<double> bond_Crank_Nicolson(std::vector<double> params, int imax, int jmax,
	double rmax, bool deriv_BC, bool fix_grid)
	// function calculates value of a bond up to T1 using a Crank-Nicolson scheme
{
	// declare and initialise local variables (dr,dt) 
	double dr = rmax / jmax;
	double dt = params[0] / imax;

	// create storage for the stock price and bond price (old and new)
	std::vector<double> r(jmax + 1), B_old(jmax + 1), B_new(jmax + 1);
	// setup and initialise the stock price 
	for (int j = 0; j <= jmax; j++)
	{
		r[j] = j * dr;
	}
	// setup and initialise the final conditions on the option price 
	for (int j = 0; j <= jmax; j++)
	{
		B_old[j] = params[1];
		B_new[j] = params[1];
	}
	// start looping through time levels
	double di = params[10] / dt - floor(params[10] / dt);
	//remainder corresponds to where T1 is between grid points
	int iT1; di < 0.5 ? iT1 = floor(params[10] / dt) : iT1 = ceil(params[10] / dt);
	//iT1 will be the closest grid point to T1
	for (int itemp = imax - 1; itemp >= iT1; itemp--)
	{
		double i;
		if (fix_grid) { di < 0.5 ? i = itemp + di : i = itemp - (1 - di); }
		// shifts all inner grid points so that T1 is on grid point. Careful to make sure arguments of following functions then take i as a double
		else { i = itemp; }

		// declare vectors for matrix equations
		std::vector<std::vector<double>> bond_matrix_vectors = Crank_Nicolson_bond_setup(params, jmax, deriv_BC, B_old, dt, dr, i);

		// solve matrix equations with thomas_solverr
		B_new = thomas_solver(bond_matrix_vectors[0], bond_matrix_vectors[1], bond_matrix_vectors[2], bond_matrix_vectors[3]);

		// set old=new 
		B_old = B_new;
	}
	// finish looping through time levels
	
	return B_new;
}

double European_Call_Crank_Nicolson(std::vector<double> params, int imax, int jmax,
	double rmax, bool deriv_BC, bool fix_grid, bool early_return)
	// Function calculates the value of a European call option on a specific coupon bond using a Crank-Nicolson scheme
{
	// declare and initialise local variables (dr,dt) 
	double dr = rmax / jmax;
	double dt = params[0] / imax;

	// create storage for the stock price and option price (old and new)
	std::vector<double> r(jmax + 1), V_old(jmax + 1), V_new(jmax + 1);
	// setup and initialise the stock price 
	for (int j = 0; j <= jmax; j++)
	{
		r[j] = j * dr;
	}
	std::vector<double> B_new = bond_Crank_Nicolson(params, imax, jmax, rmax, deriv_BC, fix_grid);

	if (early_return) {
		double bond_value = interpolator(B_new, r, params, dr);
		return bond_value - params[11];
	}

	else {
		//now use bond_value for BC of option
		double di = params[10] / dt - floor(params[10] / dt);
		//remainder corresponds to where T1 is between grid points
		int iT1; di < 0.5 ? iT1 = floor(params[10] / dt) : iT1 = ceil(params[10] / dt);

		for (int j = 0; j <= jmax; j++)
		{
			V_old[j] = std::max(B_new[j] - params[11], 0.);
			V_new[j] = std::max(B_new[j] - params[11], 0.);
		}

		// start looping through time levels
		for (int itemp = iT1 - 1; itemp >= 0; itemp--)
		{
			double i;
			if (fix_grid && itemp != 0) { di < 0.5 ? i = itemp + di : i = itemp - (1 - di); }
			// shifts all inner grid points so that T1 is on grid point except t=0
			else { i = itemp; }
			
			// declare vectors for matrix equations
			std::vector<std::vector<double>> option_matrix_vectors = Crank_Nicolson_european_option_setup(params, jmax, deriv_BC, V_old, dt, dr, i);
			
			// solve matrix equations with thomas_solverr
			V_new = thomas_solver(option_matrix_vectors[0], option_matrix_vectors[1], option_matrix_vectors[2], option_matrix_vectors[3]);
			// set old=new 
			V_old = V_new;
		}
		// finish looping through time levels

		// output the estimated option price
		return interpolator(V_new, r, params, dr);;
	}
}

double American_Call_Crank_Nicolson(std::vector<double> params, int imax, int jmax, double rmax, 
	double omega, int iter_max, double tol, bool deriv_BC, bool fix_grid, bool use_SOR)
	// Function calculates the value of a European call option on a specific coupon bond using a Crank-Nicolson scheme.
	// It gives the option to solve the ensuing matrix equations with SOR or a penalty method.
{
	// declare and initialise local variables (dr,dt) 
	double dr = rmax / jmax;
	double dt = params[0] / imax;

	// create storage for the stock price and option price (old and new)
	std::vector<double> r(jmax + 1), V_old(jmax + 1), V_new(jmax + 1), B_old(jmax + 1), B_new(jmax + 1);

	// setup and initialise the stock price 
	for (int j = 0; j <= jmax; j++)
	{
		r[j] = j * dr;
	}

	B_new = bond_Crank_Nicolson(params, imax, jmax, rmax, deriv_BC, fix_grid);
	B_old = B_new;
	
	//now use bond_value for BC of option
	for (int j = 0; j <= jmax; j++)
	{
		V_old[j] = std::max(B_old[j] - params[11], 0.);
		V_new[j] = std::max(B_new[j] - params[11], 0.);
	}

	// start looping through time levels
	double di = params[10] / dt - floor(params[10] / dt);
	int iT1; di < 0.5 ? iT1 = floor(params[10] / dt) : iT1 = ceil(params[10] / dt);

	for (int itemp = iT1 - 1; itemp >= 0; itemp--)
	{
		double i;
		if (fix_grid && itemp != 0) { di < 0.5 ? i = itemp + di : i = itemp - (1 - di); }
		// shifts all inner grid points so that T1 is on grid point
		else { i = itemp; }
		
		// declare vectors for matrix equations
		std::vector<std::vector<double>> bond_matrix_vectors = Crank_Nicolson_bond_setup(params, jmax, deriv_BC, B_old, dt, dr, i);
		std::vector<std::vector<double>> option_matrix_vectors = Crank_Nicolson_american_option_setup(params, jmax, deriv_BC, B_old, V_old, dt, dr, i);
		
		// Solve matrix equations with SOR or Penalty Method using the individual vectors
		if (use_SOR) {
			american_SOR(bond_matrix_vectors[0], bond_matrix_vectors[1], bond_matrix_vectors[2], bond_matrix_vectors[3],
				option_matrix_vectors[0], option_matrix_vectors[1], option_matrix_vectors[2], option_matrix_vectors[3], B_new, V_new,
				iter_max, jmax, tol, omega, params);
		}
		else {
			penalty_method(bond_matrix_vectors[0], bond_matrix_vectors[1], bond_matrix_vectors[2], bond_matrix_vectors[3],
				option_matrix_vectors[0], option_matrix_vectors[1], option_matrix_vectors[2], option_matrix_vectors[3], B_new,
				V_new, 1/tol, jmax, iter_max, tol, params);
		}

		// set old=new 
		V_old = V_new;
		B_old = B_new;
	}
	// finish looping through time levels

	// output the estimated option price
	return interpolator(V_new, r, params, dr);
}

int main()
{
	// declare and initialise Black Scholes parameters
	double T = 5, F = 86, theta = 0.0319, r0 = 0.0412, kappa = 0.04227, mu = 0.0036, C = 2.91,
		alpha = 0.01, beta = 0.527, sigma = 0.241, T1 = 1.5111, X = 86.4;

	std::vector<double> params = { T, F, theta, r0, kappa, mu, C, alpha, beta, sigma, T1, X };
	
	std::cout << std::setprecision(10) << std::endl;

	std::cout << European_Call_Crank_Nicolson(params, 100, 100, 1, true, true, false) << std::endl;

	std::cout << American_Call_Crank_Nicolson(params, 100, 100, 1, 1, 10000, 1.e-8, true,
	  true, true) << std::endl;
	std::cout << American_Call_Crank_Nicolson(params, 100, 100, 1, 1, 10000, 1.e-8, true,
		true, false) << std::endl;
	
	return 0;
}
