#include <fstream>
#include <functional>
#include <cmath>

const double GAMMA = 0.1;
const double BETA = 0.001;
const double DELTA = 0.1;
const double TMAX = 100.0;
const double TOL = 1e-6;
const double MI0 = 1.0;
const int MIMAX = 20;
const int N = 500;

const double ALPHA = BETA*N - GAMMA;

inline double f(double x)
{
	return ALPHA*x - BETA*x*x; 
}

double iterPicard(double res)
{
	const double oldRes = res;
	for(int i = 0; i < MIMAX; ++i)
	{
		double newRes = oldRes + (DELTA/2.0)*(f(oldRes) + f(res));
		if(fabs(newRes - res) < TOL)
			i = MIMAX;
		res = newRes;
	}
	return res;
}

double iterNewton(double res)
{
	const double oldRes = res;
	for(int i = 0; i < MIMAX; ++i)
	{
		double newRes = res-(res-oldRes-(DELTA/2.0)*(f(oldRes)+f(res)))/(1.0-(DELTA/2.0)*(ALPHA-2*BETA*res));
		if(fabs(newRes - res) < TOL)
			i = MIMAX;
		res = newRes;
	}
	return res;
}

void solve(std::function<double(double)> method, std::string fileName)
{
	std::ofstream fout;
	fout.open(fileName);
	double res = MI0;
	fout << 0.0 << " " << res << " " << N - res << "\n";
	for(double t = DELTA; t <= TMAX; t += DELTA)
	{
		res = method(res);
		fout << t << " " << res << " " << N - res << "\n";
	}
	fout << "\n\n";
	fout.close();
}

//////////////////////////////////////

const double a11 = 0.25;
const double a12 = 0.25 - sqrt(3.)/6.0;
const double a21 = 0.25 + sqrt(3.)/6.0;
const double a22 = 0.25;

const double b1 = 0.5;
const double b2 = 0.5;

const double c1 = 0.5 - sqrt(3.)/6.0;
const double c2 = 0.5 + sqrt(3.)/6.0;

double f1(double U1, double U2, double res)
{
	return U1 - res - DELTA*(a11*f(U1) + a12*f(U2));
}

double f2(double U1, double U2, double res)
{
	return U2 - res - DELTA*(a21*f(U1) + a22*f(U2));
}

void findU(double& U1, double& U2)
{
	double u1 = U1, u2 = U2, res = U1;
	for(int i = 0; i < MIMAX; ++i)
	{
		double F1 = f1(u1, u2, res);
		double F2 = f2(u1, u2, res);
		
		double m11 = 1.0 - DELTA*a11*(ALPHA - 2*BETA*u1);
		double m12 = -DELTA*a12*(ALPHA - 2*BETA*u2);
		double m21 = -DELTA*a21*(ALPHA - 2*BETA*u1);
		double m22 = 1.0 - DELTA*a22*(ALPHA - 2*BETA*u2);
		double denominator = m11*m22 - m12*m21;

		double deltaU1 = (F2*m12 - F1*m22)/denominator;
		double deltaU2 = (F1*m21 - F2*m11)/denominator;
		
		double newU1 = u1 + deltaU1;
		double newU2 = u2 + deltaU2;
		if(fabs(newU1 - u1) < TOL && fabs(newU2 - u2) < TOL)
			i = MIMAX;
			
		u1 = newU1;
		u2 = newU2;
	}
	U1 = u1;
	U2 = u2;
}

void rk2(std::string fileName)
{
	std::ofstream fout;
	fout.open(fileName);
	double res = MI0;
	double U1 = res, U2 = res;
	fout << 0.0 << " " << res << " " << N - res << "\n";
	for(double t = DELTA; t <= TMAX; t += DELTA)
	{
		findU(U1, U2);
		res = res + DELTA*(b1*f(U1) + b2*f(U2));
		U1 = U2 = res;
		fout << t << " " << res << " " << N - res << "\n";
	}
	fout << "\n\n";
	fout.close();
}

int main()
{
	solve(iterPicard, "out1.dat");
	solve(iterNewton, "out2.dat");
	rk2("out3.dat");
	return 0;
}
