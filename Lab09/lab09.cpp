#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>

//typedef
typedef std::ofstream stream;
typedef gsl_matrix matrix;
typedef gsl_vector vec;
auto& mset = gsl_matrix_set;
auto& vset = gsl_vector_set;
auto& mget = gsl_matrix_get;
auto& vget = gsl_vector_get;

//const
const int nx = 40, ny = 40;
const int N = (nx+1)*(ny+1);
const int IT_MAX = 2000;
const double delta = 1, deltaT = 1;
const double Ta = 40, Tb = 0, Tc = 30, Td = 0;
const double kb = 0.1, kd = 0.6;

//inline
inline int findL(int i, int j)
{
    return i + j*(nx+1);
}
inline int findJ(int l)
{
    return l/(nx+1);
}
inline int findI(int l, int j)
{
    return l-j*(nx+1);
}

//printing
void print(vec* T, stream& fout)
{
    for(int l = 0; l < N; ++l)
    {
        int j = findJ(l);
        int i = findI(l, j);
        fout << i*delta << " " << j*delta << " " << vget(T, l) << std::endl;
        if(i % nx == 0 && i > 0)
            fout << std::endl;
    }
    fout << std::endl << std::endl;
}

//preparations
void prepare(matrix* A, matrix* B, vec* c, vec* T)
{
    for(int l = 0; l < N; ++l)
    {
        int j = findJ(l);
        int i = findI(l, j);
        if(i == 0 || i == nx) //WB Dirichlet
        {
            mset(A, l, l, 1);
            mset(B, l, l, 1);
            vset(c, l, 0);
            vset(T, l, i == 0 ? Ta : Tc);
        }
        else if(j == 0 || j == ny) //WB von Neumann
        {
            if(j == 0)
            {
                double temp = 1 + 1/(kd*delta);
                mset(A, l, l, temp);
                temp = -1/(kd*delta);
                mset(A, l, l+nx+1, temp);
                vset(c, l, Td);
            }
            else
            {
                double temp = -1/(kb*delta);
                mset(A, l, l-nx-1, temp);
                temp = 1 + 1/(kb*delta);
                mset(A, l, l, temp);
                vset(c, l, Tb);
            }
            for(int k = 0; k < N; ++k)
                mset(B, l, k, 0); 
        }
        else
        {
            double temp = deltaT/(2*delta*delta);
            mset(A, l, l-nx-1, temp);
            mset(A, l, l-1, temp);
            mset(A, l, l+1, temp);
            mset(A, l, l+nx+1, temp);
            temp = -(2*deltaT)/(delta*delta)-1;
            mset(A, l, l, temp);
            temp = -deltaT/(2*delta*delta);
            mset(B, l, l-nx-1, temp);
            mset(B, l, l-1, temp);
            mset(B, l, l+1, temp);
            mset(B, l, l+nx+1, temp);
            temp = (2*deltaT)/(delta*delta)-1;
            mset(B, l, l, temp);
        }
    }
}

//solving
void solve(std::vector<int> iter, stream& fout)
{
    //allocation
    matrix* A = gsl_matrix_calloc(N, N);
    matrix* B = gsl_matrix_calloc(N, N);
    vec* c = gsl_vector_calloc(N);
    vec* d = gsl_vector_calloc(N);
    vec* T = gsl_vector_calloc(N);
    vec* Tcopy = gsl_vector_calloc(N);
    gsl_permutation* perm = gsl_permutation_alloc(N);
    prepare(A, B, c, T);
    int signum = 0;

    //first LU
    gsl_linalg_LU_decomp(A, perm, &signum);

    //CN
    int index = 0;
    for(int it = 0; it < IT_MAX; ++it)
    {
        for(int l = 0; l < N; ++l)
        {
            vset(d, l, vget(c, l));
            if(it+1 == iter.at(index))
                vset(Tcopy, l, vget(T, l));
        }
        gsl_blas_dgemv(CblasNoTrans, 1, B, T, 1, d);
        gsl_linalg_LU_solve(A, perm, d, T);
        if(it + 1 == iter.at(index))
        {
            print(T, fout);
            for(int l = 0; l < N; ++l)
                vset(Tcopy, l, vget(T, l)-vget(Tcopy, l));
            print(Tcopy, fout);
            index++;
        }
    }

    //deallocation
    gsl_permutation_free(perm);
    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_vector_free(c);
    gsl_vector_free(d);
    gsl_vector_free(T);
    gsl_vector_free(Tcopy);
}

int main()
{
    //START
    clock_t t = clock();
    std::cout << "START" << std::endl;

    stream fout;
    fout.open("out.dat");
    solve({100,200,500,1000,IT_MAX}, fout);
    fout.close();

    //STOP
    std::cout << "STOP" << std::endl;
    std::cout << static_cast<double>(clock()-t)/CLOCKS_PER_SEC << std::endl;
    return 0;
}