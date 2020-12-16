#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <ctime>

typedef std::vector<double> vec;

const double eps = 1, delta = 0.1;
const double V1 = 10, V2 = 0;
const double TOL = 1e-8;
const int nx = 150, ny = 100;

const double xmax = delta*nx, ymax = delta*ny;
const double sigmax = xmax*0.1, sigmay = ymax*0.1;

std::vector<vec> density(nx+1, vec(ny+1));

double inline p1(double x, double y)
{
    double a = (x-0.35*xmax)/(sigmax);
    double b = (y-0.5*ymax)/(sigmay);
    return exp(-a*a-b*b);
}
double inline p2(double x, double y)
{
    double a = (x-0.65*xmax)/(sigmax);
    double b = (y-0.5*ymax)/(sigmay);
    return -exp(-a*a-b*b);
}
double inline p(int i, int j)
{
    return p1(delta*i, delta*j) + p2(delta*i, delta*j);
}
double inline discreetResult(std::vector<vec>& v, int i, int j)
{
    double sum = v[i+1][j]+v[i-1][j]+v[i][j+1]+v[i][j-1];
    return 0.25*(sum+(pow(delta,2)*density[i][j])/eps);
}
double discreetIntegral(std::vector<vec>& v)
{
    double res = 0.0;
    for(int i = 0; i < nx; ++i)
        for(int j = 0; j < ny; ++j)
        {
            double par1 = 0.5*pow((v[i+1][j]-v[i][j])/delta,2);
            double par2 = 0.5*pow((v[i][j+1]-v[i][j])/delta,2);
            res += pow(delta,2)*(par1+par2-density[i][j]*v[i][j]);
        }
    return res;
}
void printVec(std::vector<vec>& v, std::ofstream& fout)
{
    for(int i = 0; i <= nx; ++i)
    {
        for(int j = 0; j <= ny; ++j)
            fout << i*delta << " " << j*delta << " " << v[i][j] << std::endl;
        fout << std::endl;
    }
}
void printErr(std::vector<vec>& v, std::ofstream& Errout)
{
    for(int i = 1; i < nx; ++i)
    {
        for(int j = 1; j < ny; ++j)
        {
            double err = (v[i+1][j] - 2*v[i][j] + v[i-1][j])/pow(delta,2) +
            (v[i][j+1] - 2*v[i][j] + v[i][j-1])/pow(delta,2) + density[i][j]/eps;
            Errout << i*delta << " " << j*delta << " " << err << std::endl;
        }
        Errout << std::endl;
    }
}
void global(double wG, std::ofstream& fout, std::ofstream& Sout, std::ofstream& Errout)
{
    clock_t t = clock();
    std::vector<vec> oldVec(nx+1,(vec(ny+1, 0.0))), newVec(nx+1,(vec(ny+1, 0.0)));
    for(int i = 0; i <= nx; ++i)
    {
        oldVec[i][0] = V1;
        oldVec[i][ny] = V2;
        newVec[i][0] = V1;
        newVec[i][ny] = V2;
    }
    double Sprev, Snext = 1;
    int count = 0;

    do{
        for(int i = 1; i < nx; ++i)
            for(int j = 1; j < ny; ++j)
                newVec[i][j] = discreetResult(oldVec, i, j);
        for(int j = 1; j < ny; ++j)
        {
            newVec[0][j] = newVec[1][j];
            newVec[nx][j] = newVec[nx-1][j];
        }
        for(int i = 0; i <= nx; ++i)
            for(int j = 0; j <= ny; ++j)
                oldVec[i][j] = (1.0-wG)*oldVec[i][j]+wG*newVec[i][j];

        Sprev = Snext;
        Snext = discreetIntegral(oldVec);
        Sout << count << " " << Snext << std::endl;
        count++;

    }while(fabs((Snext-Sprev)/(Sprev)) > TOL);

    printVec(oldVec, fout);
    printErr(oldVec, Errout);

    fout << std::endl << std::endl;
    Sout << std::endl << std::endl;
    Errout << std::endl << std::endl;
    std::cout << "iter: " << count << " czas: " << (clock()-(double)t)/CLOCKS_PER_SEC << std::endl;
}
void local(double wL, std::ofstream& Sout)
{
    clock_t t = clock();
    std::vector<vec> oldVec(nx+1,(vec(ny+1, 0.0)));
    for(int i = 0; i <= nx; ++i)
    {
        oldVec[i][0] = V1;
        oldVec[i][ny] = V2;
    }
    double Sprev, Snext = 1;
    int count = 0;
    do{
        for(int i = 1; i < nx; ++i)
            for(int j = 1; j < ny; ++j)
                oldVec[i][j] = (1-wL)*oldVec[i][j]+wL*discreetResult(oldVec, i, j);
        for(int j = 1; j < ny; ++j)
        {
            oldVec[0][j] = oldVec[1][j];
            oldVec[nx][j] = oldVec[nx-1][j];
        }

        Sprev = Snext;
        Snext = discreetIntegral(oldVec);
        Sout << count << " " << Snext << std::endl;
        count++;

    }while(fabs((Snext-Sprev)/(Sprev)) > TOL);

    Sout << std::endl << std::endl;
    std::cout << "iter: " << count << " czas: " << (clock()-(double)t)/CLOCKS_PER_SEC << std::endl;
}
int main()
{
    //prepare density function
    for(int i = 0; i <= nx; ++i)
        for(int j = 0; j <= ny; ++j)
            density[i][j] = p(i,j);
    //preparing output
    std::ofstream fout, Sout, Errout;
    fout.open("out.dat");
    Sout.open("S.dat");
    Errout.open("err.dat");
    //global 
    global(0.6, fout, Sout, Errout);
    global(1.0, fout, Sout, Errout);
    Errout.close();
    fout.close();
    //local
    local(1.0, Sout);
    local(1.4, Sout);
    local(1.8, Sout);
    local(1.9, Sout);
    Sout.close();
    return 0;
}

