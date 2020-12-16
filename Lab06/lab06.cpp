#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <cmath>
#include "mgmres.h"

typedef std::function<double(double, double, const int, const int)> f;
typedef std::ofstream stream;
typedef std::vector <double> vecD;
typedef std::vector <int> vecI;
const double delta = 0.1;

//inline
inline double defaultp(double x, double y, const int nx, const int ny)
{
    return 0;
}
inline double p1(double x, double y, const int nx, const int ny)
{
    double xmax = delta*nx, ymax = delta*ny;
    double sigma = xmax/10;
    return exp(-pow((x-0.25*xmax)/(sigma*sigma), 2)-pow((y-0.5*ymax)/(sigma*sigma), 2));
}
inline double p2(double x, double y, const int nx, const int ny)
{
    double xmax = delta*nx, ymax = delta*ny;
    double sigma = xmax/10;
    return -exp(-pow((x-0.75*xmax)/(sigma*sigma), 2)-pow((y-0.5*ymax)/(sigma*sigma), 2));
}
inline int findL(int i, int j, const int nx)
{
    return i + j*(nx+1);
}
inline int findJ(int l, const int nx)
{
    return l/(nx+1);
}
inline int findI(int l, int j, const int nx)
{
    return l-j*(nx+1);
}
inline const double eps(int l, const int nx, const double eps1, const double eps2)
{
    int j = findJ(l, nx);
    int i = findI(l, j, nx);
    return i > nx/2 ? eps2 : eps1;
}
//printing
void printTestA(vecD& v, const int N, stream& fout)
{
    for(int l = 0; l < N; ++l)
        fout << l << " " << v[l] << std::endl;
    fout << std::endl << std::endl;
}
void printTestB(vecD& v, const int N, const int nx, stream& fout)
{
    for(int l = 0; l < N; ++l)
    {
        int j = findJ(l, nx);
        int i = findI(l, j, nx);
        fout << l << " " << i << " " << j << " " << v[l] << std::endl;
    }
    fout << std::endl << std::endl;
}
void print(vecD& v, const int N, const int nx, stream& fout)
{
    for(int l = 0; l < N; ++l)
    {
        int j = findJ(l, nx);
        int i = findI(l, j, nx);
        fout << i*delta << " " << j*delta << " " << v[l] << std::endl;
        if(i % nx == 0 && i > 0)
            fout << std::endl;
    }
    fout << std::endl << std::endl;
}
//solving
void solve(vecD params, stream& fout, bool test = false, f p1 = defaultp, f p2 = defaultp)
{
    //const
    const double eps1 = params[0], eps2 = params[1];
    const double V1 = params[2], V2 = params[3], V3 = params[4], V4 = params[5];
    const int nx = params[6], ny = params[7];
    const int N = (nx+1)*(ny+1);
    
    //vectors
    vecD a(5*N), b(N), V(N);
    vecI ja(5*N), ia(N+1, -1);

    //algorithm
    int k = -1;
    for(int l = 0; l < N; ++l)
    {
        int j = findJ(l, nx);
        int i = findI(l, j, nx);
        int edge = 0;
        double Vb = 0;
        if(i == 0)
        {
            edge = 1;
            Vb = V1;
        }
        if(j == ny)
        {
            edge = 1;
            Vb = V2;
        }
        if(i == nx)
        {
            edge = 1;
            Vb = V3;
        }
        if(j == 0)
        {
            edge = 1;
            Vb = V4;
        }
        b[l] = -(p1(i*delta, j*delta, nx, ny) + p2(i*delta, j*delta, nx, ny));
        if(edge == 1)
            b[l] = Vb;
        ia[l] = -1;
        //1
        if(l-nx-1 >= 0 && edge == 0)
        {
            ++k;
            if(ia[l] < 0)
                ia[l] = k;
            a[k] = eps(l, nx, eps1, eps2)/(delta*delta);
            ja[k] = l-nx-1;
        }
        //2
        if(l-1 >= 0 && edge == 0)
        {
            ++k;
            if(ia[l] < 0)
                ia[l] = k;
            a[k] = eps(l, nx, eps1, eps2)/(delta*delta);
            ja[k] = l-1;
        }
        //3
        ++k;
        if(ia[l] < 0)
            ia[l] = k;
        if(edge == 0)
            a[k] = -(2*eps(l, nx, eps1, eps2)+eps(l+1, nx, eps1, eps2)+eps(l+nx+1, nx, eps1, eps2))/(delta*delta); 
        else 
            a[k] = 1;
        ja[k] = l;
        //4
        if(l < N && edge == 0)
        {
            ++k;
            a[k] = eps(l+1, nx, eps1, eps2)/(delta*delta);
            ja[k] = l+1;
        }
        //5
        if(l < N-nx-1 && edge == 0)
        {
            ++k;
            a[k] = eps(l+nx+1, nx, eps1, eps2)/(delta*delta);
            ja[k] = l+nx+1;
        }
    }
    int nz_num = k+1;
    ia[N] = nz_num;

    //test
    if(test == true)
    {
        stream testout;
        testout.open("test.dat");
        printTestA(a, 5*N, testout);
        printTestB(b, N, nx, testout);
        testout.close();
    }

    //function pmgmres_ilu_cr
    pmgmres_ilu_cr(N, nz_num, ia.data(), ja.data(), a.data(), V.data(), b.data(), 500, 500, 1e-8, 1e-8);

    //printResult
    print(V, N, nx, fout);
}
int main()
{
    stream fout; 
    fout.open("out.dat");   
    solve({1,1,10,-10,10,-10,4,4}, fout, true);
    solve({1,1,10,-10,10,-10,50,50}, fout);
    solve({1,1,10,-10,10,-10,100,100}, fout);
    solve({1,1,10,-10,10,-10,200,200}, fout);
    solve({1,1,0,0,0,0,100,100}, fout, false, p1, p2);
    solve({1,2,0,0,0,0,100,100}, fout, false, p1, p2);
    solve({1,10,0,0,0,0,100,100}, fout, false, p1, p2);
    fout.close();
    return 0;
}