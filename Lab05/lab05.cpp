#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

typedef std::vector<double> vec;

const double delta = 0.2; 
const double TOL = 1e-8;
const int nx = 128, ny = 128, kmax = 16;
const double xmax = delta*nx, ymax = delta*ny;

int iterCounter = 0;

double inline wb(int i, int j, bool mode)
{
    double x = delta*i;
    double y = delta*j;
    if(mode)
        return sin(2*M_PI*(x/xmax));
    else
        return sin(M_PI*(y/ymax));
}
double stopCondition(std::vector<vec>& V, const int k)
{
    double res = 0.0;
    for(int i = 0; i < nx; i += k)
        for(int j = 0; j < ny; j += k)
        {
            double den = 2*k*delta;
            double par1 = (V[i+k][j]-V[i][j])/den;
            double par2 = (V[i+k][j+k]-V[i][j+k])/den;
            double par3 = (V[i][j+k]-V[i][j])/den;
            double par4 = (V[i+k][j+k]-V[i+k][j])/den;
            res += (pow(delta*k, 2)/2)*(pow(par1+par2, 2)+pow(par3+par4, 2));
        }
    return res;
}
void print(std::vector<vec>& V, const int k, std::ofstream& fout)
{
    for(int i = 0; i <= nx; i += k)
    {
        for(int j = 0; j <= ny; j += k)
            fout << delta*i << " " << delta*j << " " << V[i][j] << "\n";
        fout << "\n";
    }
    fout << "\n\n";
}
void compaction(std::vector<vec>& V, const int k)
{
    for(int i = 0; i < nx; i += k)
        for(int j = 0; j < ny; j += k)
        {
            V[i+k/2][j+k/2] = 0.25*(V[i][j]+V[i+k][j]+V[i][j+k]+V[i+k][j+k]);
            if(i < nx-k)
                V[i+k][j+k/2] = 0.5*(V[i+k][j]+V[i+k][j+k]);
            if(j < ny-k)
                V[i+k/2][j+k] = 0.5*(V[i][j+k]+V[i+k][j+k]);
            if(j > 0)
                V[i+k/2][j] = 0.5*(V[i][j]+V[i+k][j]);
            if(i > 0)
                V[i][j+k/2] = 0.5*(V[i][j]+V[i][j+k]);
        }
}
void relax(std::vector<vec>& V, const int k, std::ofstream& Sout)
{
    double Sprev, Snext = 1;
    do
    {
        for(int i = k; i < nx; i += k)
            for(int j = k; j < ny; j += k)
                V[i][j] = 0.25*(V[i+k][j]+V[i-k][j]+V[i][j+k]+V[i][j-k]);
        Sprev = Snext;
        Snext = stopCondition(V, k);
        Sout << iterCounter << " " << Snext << "\n"; 
        ++iterCounter;
    } while(fabs((Snext-Sprev)/Sprev) > TOL);
    Sout << "\n\n";
    std::cout << "k = " << k << " iter: " << iterCounter << std::endl;
}
int main()
{
    //preparing grid
    std::vector<vec> V(nx+1, vec(ny+1, 0.0));
    for(int i = 0; i <= nx; ++i)
    {
        V[i][0] = wb(i, 0, true);
        V[i][ny] = -V[i][0];
    }
    for(int j = 0; j <= ny; ++j)
    {
        V[0][j] = wb(0, j, false);
        V[nx][j] = V[0][j];
    }
    //multigrid relaxation
    std::ofstream fout, Sout;
    Sout.open("S.dat");
    fout.open("out.dat");
    for(int i = kmax; i; i /= 2)
    {
        relax(V, i, Sout);
        print(V, i, fout);
        if(i > 1)
            compaction(V, i);
    }
    fout.close();
    Sout.close();
    return 0;
}