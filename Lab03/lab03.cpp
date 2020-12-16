#include <fstream>
#include <utility>
#include <functional>
#include <cmath>

using pair = std::pair<double,double>; 

const double x0 = 0.01;
const double v0 = 0.0;
const double t0 = 0.0, tMax = 40.0;
const double deltaT0 = 1.0;

const double S = 0.75;
const double alpha = 5.0;
const double p = 2.0;


inline void print(double t, double deltaT, double x, double v, std::ofstream& fout)
{
    fout << t << " " << deltaT << " " << x << " " << v << std::endl;
}

inline double f(double x, double v)
{
    return v;
}

inline double g(double x, double v)
{
    return alpha*(1-x*x)*v-x;
}

///////////////

pair trapeziumM(double x, double v, double deltaT, const double alpha)
{
    double sigma = 1e-10;
    double newX = x, newV = v;
    double deltaX, deltaV;
    do{
        double F = newX-x-(deltaT/2.0)*(f(x,v)+f(newX, newV));
        double G = newV-v-(deltaT/2.0)*(g(x,v)+g(newX, newV));

        double a11 = 1.0;
        double a12 = -deltaT/2.0;
        double a21 = -(deltaT/2.0)*(-2*alpha*newX*newV-1);
        double a22 = 1-(deltaT/2.0)*alpha*(1-newX*newX);

        deltaX = ((-F)*a22-(-G)*a12)/(a11*a22-a12*a21);
        deltaV = ((-G)*a11-(-F)*a21)/(a11*a22-a12*a21);

        newX += deltaX;
        newV += deltaV;

        deltaX = fabs(deltaX);
        deltaV = fabs(deltaV);
    }while(deltaX >= sigma || deltaV >= sigma);
    return pair(newX, newV);
}

pair rk2M(double x, double v, double deltaT, const double alpha)
{
    double k1x = f(x, v);
    double k1v = g(x, v);

    double k2x = v+deltaT*k1v;
    double k2v = alpha*(1-(x+deltaT*k1x)*(x+deltaT*k1x))*k2x-(x+deltaT*k1x); 

    return pair(x+(deltaT/2)*(k1x+k2x), v+(deltaT/2)*(k1v+k2v));
}

void timeStepControl(const double TOL, std::function<pair(double,double,double,const double)> method, std::ofstream& fout)
{
    double x = x0;
    double v = v0;
    double deltaT = deltaT0;
    for(double t = t0; t < tMax;)
    {
        pair p1 = method(x, v, deltaT, alpha);
        p1 = method(p1.first, p1.second, deltaT, alpha);

        pair p2 = method(x, v, 2.0*deltaT, alpha);

        double Ex = fabs((p1.first-p2.first)/(pow(2.0, p)-1));
        double Ev = fabs((p1.second-p2.second)/(pow(2.0, p)-1));
        if(fmax(Ex, Ev) < TOL)
        {
            t += 2*deltaT;
            x = p1.first;
            v = p1.second;
            print(t, deltaT, x, v, fout);
        }
        deltaT = pow((S*TOL)/fmax(Ex, Ev), 1.0/(p+1))*deltaT;
    }
    fout << "\n\n";
}

int main()
{
    std::ofstream fout; 
    fout.open("out1.dat");
    timeStepControl(1e-2, trapeziumM, fout);
    timeStepControl(1e-5, trapeziumM, fout);
    fout.close(); 
 
    fout.open("out2.dat");
    timeStepControl(1e-2, rk2M, fout);
    timeStepControl(1e-5, rk2M, fout);
    fout.close(); 

    return 0;
}
