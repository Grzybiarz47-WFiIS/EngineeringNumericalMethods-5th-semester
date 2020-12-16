#include <fstream>
#include <cmath>

const double LAMBDA = -1.;

inline double y(double t)
{
	return exp(LAMBDA*t);
}
void euler(double tMin, double tMax, double deltaT, std::ofstream& fout)
{
	double res = y(tMin);
	fout << tMin << " " << res << " " << 0.0 << std::endl;
	for(double i = tMin + deltaT; i <= tMax; i += deltaT)
	{
		res = res*(LAMBDA*deltaT + 1.0);
		fout << i << " " << res << " " << fabs(res - y(i)) << std::endl;
	}
	fout << std::endl << std::endl;
}
void rk2(double tMin, double tMax, double deltaT, std::ofstream& fout)
{
	double res = y(tMin);
	fout << tMin << " " << res << " " << 0.0 << std::endl;
	for(double i = tMin + deltaT; i <= tMax; i += deltaT)
	{
		double k1 = LAMBDA*res;
		double k2 = LAMBDA*(res + deltaT*k1);
		res = res + (deltaT/2.0)*(k1 + k2);
		fout << i << " " << res << " " << fabs(res - y(i)) << std::endl;
	}
	fout << std::endl << std::endl;
}
void rk4(double tMin, double tMax, double deltaT, std::ofstream& fout)
{
	double res = y(tMin);
	fout << tMin << " " << res << " " << 0.0 << std::endl;
	for(double i = tMin + deltaT; i <= tMax; i += deltaT)
	{
		double k1 = LAMBDA*res;
		double k2 = LAMBDA*(res + (deltaT/2.0)*k1);
		double k3 = LAMBDA*(res + (deltaT/2.0)*k2);
		double k4 = LAMBDA*(res + deltaT*k3);
		res = res + (deltaT/6.0)*(k1 + 2*k2 + 2*k3 + k4);
		fout << i << " " << res << " " << fabs(res - y(i)) << std::endl;
	}
	fout << std::endl << std::endl;
}

////RLC
const double R = 100, L = 0.1, C = 0.001;
const double w0 = 1.0/(sqrt(L*C));
const double T0 = 2*M_PI/w0;

inline double V(double x)
{
	return 10*sin(x);
}
void rrz2(double factor, double deltaT, std::ofstream& fout)
{
	double wv = factor*w0;	
	double Qres = 0.0, Ires = 0.0;
	double tMax = 4*T0;
	fout << 0.0 << " " << Qres << " " << Ires << std::endl;
	for(double t = deltaT; t <= tMax; t += deltaT)
	{
		double Qk1 = Ires;
		double Ik1 = V(wv*(t-deltaT)) - (1.0/(L*C))*Qres - (R/L)*Ires;
		double Qk2 = Ires + (deltaT/2.0)*Qk1;
		double Ik2 = V(wv*(t-0.5*deltaT)) - (1.0/(L*C))*(Qres+(deltaT/2.0)*Qk1) - (R/L)*(Ires+(deltaT/2.0)*Ik1);
		double Qk3 = Ires + (deltaT/2.0)*Qk2;
		double Ik3 = V(wv*(t-0.5*deltaT)) - (1.0/(L*C))*(Qres+(deltaT/2.0)*Qk2) - (R/L)*(Ires+(deltaT/2.0)*Ik2);
		double Qk4 = Ires + deltaT*Qk3;
		double Ik4 = V(wv*t) - (1.0/(L*C))*(Qres+deltaT*Qk3) - (R/L)*(Ires+deltaT*Ik3);
		Qres = Qres + (deltaT/6.0)*(Qk1+2*Qk2+2*Qk3+Qk4);
		Ires = Ires + (deltaT/6.0)*(Ik1+2*Ik2+2*Ik3+Ik4);
		fout << t << " " << Qres << " " << Ires << std::endl;
	}
	fout << std::endl << std::endl;
}

int main()
{
	//1 - Metoda jawna Eulera
	std::ofstream fout;
	fout.open("out1.dat");
	euler(0, 5, 0.01, fout);
	euler(0, 5, 0.1, fout);
	euler(0, 5, 1., fout);
	fout.close();
	//2 - Metoda RK2
	fout.open("out2.dat");
	rk2(0, 5, 0.01, fout);
	rk2(0, 5, 0.1, fout);
	rk2(0, 5, 1., fout);
	fout.close();
	//3 - Metoda RK4
	fout.open("out3.dat");
	rk4(0, 5, 0.01, fout);
	rk4(0, 5, 0.1, fout);
	rk4(0, 5, 1., fout);
	fout.close();
	//4 - RRZ rzedu II
	fout.open("out4.dat");
	rrz2(0.5, 1e-4, fout);
	rrz2(0.8, 1e-4, fout);
	rrz2(1.0, 1e-4, fout);
	rrz2(1.2, 1e-4, fout);
	fout.close();
	return 0;
}
