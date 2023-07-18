#include "option_on_bond_CN_header.h"

void SOR(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d,
	std::vector<double>& B_new, int iter_max, int jmax, double tol, double omega)
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

void american_SOR(const std::vector<double>& aB, const std::vector<double>& bB, const std::vector<double>& cB, const std::vector<double>& dB,
	const std::vector<double>& aV, const std::vector<double>& bV, const std::vector<double>& cV, const std::vector<double>& dV,
	std::vector<double>& B_new, std::vector<double>& V_new, int iter_max, int jmax, double tol, double omega, std::vector<double> params)
{
	int sor;
	for (sor = 0; sor < iter_max; sor++)
	{
		double error = 0.; // calculate residual
		// SOR equations in here
		{
			double yB = (dB[0] - cB[0] * B_new[1]) / bB[0];
			B_new[0] = B_new[0] + omega * (yB - B_new[0]);
		}
		for (int j = 1; j < jmax; j++)
		{
			double yB = (dB[j] - aB[j] * B_new[j - 1] - cB[j] * B_new[j + 1]) / bB[j];
			B_new[j] = B_new[j] + omega * (yB - B_new[j]);
		}
		{
			double yB = (dB[jmax] - aB[jmax] * B_new[jmax - 1]) / bB[jmax];
			B_new[jmax] = B_new[jmax] + omega * (yB - B_new[jmax]);
		}

		{
			double yV = (dV[0] - cV[0] * V_new[1]) / bV[0];
			yV = std::max(B_new[0] - params[11], V_new[0] + omega * (yV - V_new[0]));
			//not B_old as trying to update V_new
			error += fabs(V_new[0] - yV);
			V_new[0] = yV;
		}
		for (int j = 1; j < jmax; j++)
		{
			double yV = (dV[j] - aV[j] * V_new[j - 1] - cV[j] * V_new[j + 1]) / bV[j];
			yV = std::max(B_new[j] - params[11], V_new[j] + omega * (yV - V_new[j]));
			error += fabs(V_new[j] - yV);
			V_new[j] = yV;
		}
		{
			double yV = (dV[jmax] - aV[jmax] * V_new[jmax - 1]) / bV[jmax];
			yV = std::max(B_new[jmax] - params[11], V_new[jmax] + omega * (yV - V_new[jmax]));
			error += fabs(V_new[jmax] - yV);
			V_new[jmax] = yV;
		}
		// check for convergence and exit loop if converged
		if (error < tol)
			break;
		if (sor == iter_max - 1)
			std::cout << error << std::endl;
	}
	if (sor == iter_max)
		std::cout << "\n NOT CONVERGED \n";
}

std::vector<double> thomas_solver(const std::vector<double>& a, const std::vector<double>& b_,
	const std::vector<double>& c, std::vector<double>& d)
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

void penalty_method(const std::vector<double>& aB, const std::vector<double>& bB, const std::vector<double>& cB, std::vector<double>& dB,
	const std::vector<double>& aV, const std::vector<double>& bV, const std::vector<double>& cV, const std::vector<double>& dV,
	std::vector<double>& B_new, std::vector<double>& V_new, double rho, int jmax, int iter_max, double tol, std::vector<double> params)
{
	//calculate new B with thomas_solver
	B_new = thomas_solver(aB, bB, cB, dB);

	int penaltyIt; //penalty scheme to calculate V
	for (penaltyIt = 0; penaltyIt < iter_max; penaltyIt++)
	{
		// create new vectors containing a copy of the FD approximations
		std::vector<double> aHat(aV), bHat(bV), cHat(cV), dHat(dV);

		// apply penalty here to finite difference scheme
		for (int j = 1; j < jmax; j++)
		{
			// if current value suggests apply penalty, adjust matrix equations
			if (V_new[j] < B_new[j] - params[11])
			{
				bHat[j] = bV[j] - rho; dHat[j] = dV[j] - rho * (B_new[j] - params[11]);
			}
		}

		// solve with thomas method
		std::vector<double> y = thomas_solver(aHat, bHat, cHat, dHat);
		// y now contains next guess at solution
		// check for differences between V_new and y
		double error = 0.;
		for (int j = 0; j <= jmax; j++)
			error += (V_new[j] - y[j]) * (V_new[j] - y[j]);

		// update value of V_new
		V_new = y;
		// make an exit condition when solution is converged
		if (error < tol * tol)
		{
			break;
		}
	}
	if (penaltyIt >= iter_max)
	{
		std::cout << " Error NOT converging within required iterations" << std::endl;
		throw;
	}
}

std::vector<std::vector<double>> Crank_Nicolson_bond_setup(std::vector<double> params, int jmax, bool deriv_BC,
	std::vector<double> B_old, double dt, double dr, double i)
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
		b[j] = -1. / dt - 0.5 * params[9] * params[9] * j * j * pow(j * dr, 2 * params[8] - 2)
			- j * dr / 2.;
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

std::vector<std::vector<double>> Crank_Nicolson_european_option_setup(std::vector<double> params, int jmax, bool deriv_BC,
	std::vector<double> V_old, double dt, double dr, double i)
{
	// declare vectors for matrix equations
	std::vector<double> aV(jmax + 1), bV(jmax + 1), cV(jmax + 1), dV(jmax + 1);
	// set up matrix equations
	aV[0] = 0.;
	bV[0] = 1. / dt + 1. / dr * params[4] * params[2] * exp(params[5] * (i + 0.5) * dt);
	cV[0] = -1. / dr * params[4] * params[2] * exp(params[5] * (i + 0.5) * dt);
	dV[0] = 1. / dt * V_old[0];
	for (int j = 1; j <= jmax - 1; j++)
	{
		aV[j] = 0.25 * (params[9] * params[9] * j * j * pow(j * dr, 2 * params[8] - 2) -
			1. / dr * params[4] * params[2] * exp(params[5] * (i + 0.5) * dt) + params[4] * j);
		bV[j] = -1. / dt - 0.5 * params[9] * params[9] * j * j * pow(j * dr, 2 * params[8] - 2)
			- j * dr / 2.;
		cV[j] = 0.25 * (params[9] * params[9] * j * j * pow(j * dr, 2 * params[8] - 2) +
			1. / dr * params[4] * params[2] * exp(params[5] * (i + 0.5) * dt) - params[4] * j);
		dV[j] = -aV[j] * V_old[j - 1] - (bV[j] + 2. / dt) * V_old[j] - cV[j] * V_old[j + 1];
	}
	if (deriv_BC) {
		aV[jmax] = -1. / dr;
		bV[jmax] = 1. / dr;
		dV[jmax] = 0.;
	}
	else {
		aV[jmax] = 0.;
		bV[jmax] = 1.;
		dV[jmax] = 0.;
	}
	std::vector<std::vector<double>> result{ aV, bV, cV, dV };
	return result;
}

std::vector<std::vector<double>> Crank_Nicolson_american_option_setup(std::vector<double> params, int jmax, bool deriv_BC,
	std::vector<double> B_old, std::vector<double> V_old, double dt, double dr, double i)
{
	// declare vectors for matrix equations
	std::vector<double> aV(jmax + 1), bV(jmax + 1), cV(jmax + 1), dV(jmax + 1);
	// set up matrix equations
	aV[0] = 0.;
	bV[0] = 1.;
	cV[0] = 0;
	dV[0] = B_old[0] - params[11]; //doesn't matter if B_new as both equal at this point
	for (int j = 1; j <= jmax - 1; j++)
	{
		aV[j] = 0.25 * (params[9] * params[9] * j * j * pow(j * dr, 2 * params[8] - 2) -
			1. / dr * params[4] * params[2] * exp(params[5] * (i + 0.5) * dt) + params[4] * j);
		bV[j] = -1. / dt - 0.5 * params[9] * params[9] * j * j * pow(j * dr, 2 * params[8] - 2)
			- j * dr / 2.;
		cV[j] = 0.25 * (params[9] * params[9] * j * j * pow(j * dr, 2 * params[8] - 2) +
			1. / dr * params[4] * params[2] * exp(params[5] * (i + 0.5) * dt) - params[4] * j);
		dV[j] = -aV[j] * V_old[j - 1] - (bV[j] + 2. / dt) * V_old[j] - cV[j] * V_old[j + 1];
	}
	if (deriv_BC) {
		aV[jmax] = -1. / dr;
		bV[jmax] = 1. / dr;
		dV[jmax] = 0.;
	}
	else {
		aV[jmax] = 0.;
		bV[jmax] = 1.;
		dV[jmax] = 0.;
	}

	std::vector<std::vector<double>> result{ aV, bV, cV, dV };
	return result;
}

double interpolator(std::vector<double> input, std::vector<double> r, std::vector<double> params, double dr)
{
	double bond_value;
	{
		int jStar = params[3] / dr;
		double sum = 0.;
		sum += (params[3] - r[jStar]) / dr * input[jStar + 1];
		sum += (r[jStar + 1] - params[3]) / dr * input[jStar];
		bond_value = sum;
	}
	return bond_value;
}
