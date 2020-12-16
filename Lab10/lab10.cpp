#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <cmath>

typedef std::function<double(double)> func;
typedef std::vector <double> vecD;

const int nx = 150, nt = 1000;
const double delta = 0.1, dt = 0.05;
const double sigma = 0.5, xA = 7.5, xF = 2.5;
std::ofstream Eout, fout;

inline double gaussian(double x)
{
    return exp(-pow(x-xA, 2)/(2*pow(sigma, 2)));
}
inline double aF(double x, double t)
{
    return cos((50*t)/(nt*dt))*(x == xF ? 1 : 0);
}
inline void wp(vecD& u, func f)
{
    u.at(0) = u.at(nx) = 0.0;
    for(int i = 1; i < nx; ++i)
        u.at(i) = f(i*delta);
}
void findA(vecD& a, vecD& u, vecD& u0, double alpha, double beta, double t)
{
    for(int i = 1; i < nx; ++i)
    {
        double expr1 = (u.at(i+1)-2*u.at(i)+u.at(i-1))/pow(delta, 2);
        double expr2 = beta*(u.at(i)-u0.at(i))/dt;
        double expr3 = alpha*aF(i*delta, t);
        a.at(i) = expr1-expr2+expr3;
    }
}
void findVorVp(vecD& vp, vecD& v, vecD& a)
{
    for(int i = 0; i <= nx; ++i)
        vp.at(i) = v.at(i)+(dt/2)*a.at(i);
}
void findU(vecD& u, vecD& vp)
{
    for(int i = 0; i <= nx; ++i)
        u.at(i) += dt*vp.at(i);
}
double E(vecD& u, vecD& v)
{
    double expr = (delta/4)*(pow((u.at(1)-u.at(0))/delta, 2)+pow((u.at(nx)-u.at(nx-1))/delta, 2));
    double sum = 0;
    for(int i = 1; i < nx; ++i)
        sum += pow(v.at(i), 2)+pow((u.at(i+1)-u.at(i-1))/(2*delta), 2);
    return (delta/2)*sum + expr;
}
void solve(double alpha, double beta, func f)
{
    vecD v(nx+1, 0.0), vp(nx+1, 0.0);
    vecD u(nx+1), u0(nx+1), a(nx+1);

    wp(u, f);
    u0 = u;
    findA(a, u, u0, alpha, beta, 0.0);

    for(int n = 1; n <= nt; ++n)
    {
        findVorVp(vp, v, a);
        u0 = u;
        findU(u, vp);
        findA(a, u, u0, alpha, beta, dt*n);
        findVorVp(v, vp, a);

        Eout << n*dt << " " << E(u, v) << std::endl;
        for(int i = 0; i <= nx; ++i)
            fout << n*dt << " " << i*delta << " " << u.at(i) << std::endl;
        fout << std::endl;
    }
    fout << "\n\n";
    Eout << "\n\n";
}
int main()
{
    fout.open("out.dat");
    Eout.open("Eout.dat");
    solve(0.0, 0.0, gaussian);
    solve(0.0, 0.1, gaussian);
    solve(0.0, 1.0, gaussian);
    solve(1.0, 1.0, [](double x){return 0.0;});
    Eout.close();
    fout.close();
    return 0;
}

