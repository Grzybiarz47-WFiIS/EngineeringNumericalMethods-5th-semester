#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

typedef std::vector<double> vecD;
typedef std::vector<vecD> grid;

const int nx = 200, ny = 90, I1 = 50, J1 = 55;
const int IT_MAX = 20000;

const double delta = 0.01, p = 1.0, mi = 1.0;

bool checkEdge(int i, int j)
{
    if(i == 0 || i == nx || j == 0 || j == ny)
        return false;
    if(i <= I1 && j <= J1)
        return false;
    return true;
}
double findQout(double Qin)
{
    double yny = delta*ny;
    double yJ1 = delta*J1;
    return Qin*(pow(yny, 3)-pow(yJ1, 3)-3*yny*yny*yJ1+3*yny*yJ1*yJ1)/pow(yny, 3);
}
double findErr(grid& psi, grid& dzeta)
{
    double res = 0.0;
    int j2 = J1+2;
    for(int i = 1; i < nx; ++i)
        res += psi[i+1][j2]+psi[i-1][j2]+psi[i][j2+1]+psi[i][j2-1]-4*psi[i][j2]-pow(delta, 2)*dzeta[i][j2];
    return res;
}
void print(grid& v, std::ostream& fout)
{
    for(int i = 0; i <= nx; ++i)
    {
        for(int j = 0; j <= ny; ++j)
            fout << i*delta << " " << j*delta << " " << v[i][j] << std::endl;
        fout << std::endl;
    }
    fout << "\n\n"; 
}
void printU(grid& v, std::ostream& fout)
{
    for(int i = 1; i < nx; ++i)
    {
        for(int j = 1; j < ny; ++j)
        {
            if(i > I1 || j > J1)
                fout << i*delta << " " << j*delta << " " << (v[i][j-1]-v[i][j+1])/(-2*delta) << std::endl;
            else
                fout << i*delta << " " << j*delta << " " << 0.0 << std::endl; 
        }
        fout << std::endl;
    }
    fout << "\n\n"; 
}
void printV(grid& v, std::ostream& fout)
{
    for(int i = 1; i < nx; ++i)
    {
        for(int j = 1; j < ny; ++j)
        {
            if(i > I1 || j > J1)
                fout << i*delta << " " << j*delta << " " << (v[i+1][j]-v[i-1][j])/(-2*delta) << std::endl;
            else
                fout << i*delta << " " << j*delta << " " << 0.0 << std::endl; 
        }
        fout << std::endl;
    }
    fout << "\n\n"; 
}
void wbPsi(grid& v, double Qout, double Qin)
{
    //A - in
    for(int j = J1; j <= ny; ++j)
    {
        double y = j*delta;
        double yJ1 = J1*delta;
        double yny = ny*delta;
        v[0][j] = (Qin/(2*mi))*((pow(y, 3)/3)-(pow(y, 2)/2)*(yJ1+yny)+y*yJ1*yny);
    }
    //C - out
    for(int j = 0; j <= ny; ++j)
    {
        double y = j*delta;
        double yJ1 = J1*delta;
        double yny = ny*delta;
        double temp = (Qin*pow(yJ1, 2)*(3*yny-yJ1))/(12*mi);
        v[nx][j] = (Qout/(2*mi))*((pow(y, 3)/3)-(pow(y, 2)/2)*yny)+temp;
    }
    //B 
    for(int i = 1; i < nx; ++i)
        v[i][ny] = v[0][ny];
    //D
    for(int i = I1; i < nx; ++i)
        v[i][0] = v[0][J1];
    //E
    for(int j = 1; j <= J1; ++j)
        v[I1][j] = v[0][J1];
    //F
    for(int i = 1; i <= I1; ++i)
        v[i][J1] = v[0][J1];
}
void wbDzeta(grid& v, grid& psi, double Qout, double Qin)
{
    //A - in
    for(int j = J1; j <= ny; ++j)
    {
        double y = j*delta;
        double yJ1 = J1*delta;
        double yny = ny*delta;
        v[0][j] = (Qin/(2*mi))*(2*y-yJ1-yny);
    }
    //C - out
    for(int j = 0; j <= ny; ++j)
    {
        double y = j*delta;
        double yny = ny*delta;
        v[nx][j] = (Qout/(2*mi))*(2*y-yny);
    }
    //B 
    for(int i = 1; i < nx; ++i)
        v[i][ny] = (2/pow(delta, 2))*(psi[i][ny-1]-psi[i][ny]);
    //D
    for(int i = I1+1; i < nx; ++i)
        v[i][0] = (2/pow(delta, 2))*(psi[i][1]-psi[i][0]);
    //E
    for(int j = 1; j < J1; ++j)
        v[I1][j] = (2/pow(delta, 2))*(psi[I1+1][j]-psi[I1][j]);
    //F
    for(int i = 1; i <= I1; ++i)
        v[i][J1] = (2/pow(delta, 2))*(psi[i][J1+1]-psi[i][J1]);
    //E/F
    v[I1][J1] = 0.5*(v[I1-1][J1]+v[I1][J1-1]);
}
void solve(double Q)
{
    //preparing
    std::ofstream fout(std::string("out") + std::to_string(int(Q)) + std::string(".dat"));
    const double Qin = Q;
    const double Qout = findQout(Qin);
    grid psi(nx+1, vecD(ny+1, 0.0)), dzeta(nx+1, vecD(ny+1, 0.0));
    grid psi_copy;
    wbPsi(psi, Qout, Qin);

    //algorithm
    for(int it = 1; it <= IT_MAX; ++it)
    {
        double omega;
        if(it < 2000)
            omega = 0.0;
        else
            omega = 1.0;
        for(int i = 1; i < nx; ++i)
            for(int j = 1; j < ny; ++j)
                if(checkEdge(i, j))
                {
                    psi[i][j] = 0.25*(psi[i+1][j]+psi[i-1][j]+psi[i][j+1]+psi[i][j-1]-pow(delta, 2)*dzeta[i][j]);
                    double temp = (omega*p)/(16*mi);
                    double temp1 = (psi[i][j+1]-psi[i][j-1])*(dzeta[i+1][j]-dzeta[i-1][j]);
                    double temp2 = (psi[i+1][j]-psi[i-1][j])*(dzeta[i][j+1]-dzeta[i][j-1]);
                    dzeta[i][j] = 0.25*(dzeta[i+1][j]+dzeta[i-1][j]+dzeta[i][j+1]+dzeta[i][j-1]);
                    dzeta[i][j] -= temp*(temp1-temp2);
                }
        wbDzeta(dzeta, psi, Qout, Qin);

        if(it % (IT_MAX/400) == 0)
            std::cout << "Iter: " << it << " Error: " << findErr(psi, dzeta) << std::endl;
    }

    //printing
    print(psi, fout);
    print(dzeta, fout);
    printU(psi, fout);
    printV(psi, fout);

    fout.close();
}
int main()
{
    solve(-1000);
    solve(-4000);
    solve(4000);
    return 0;
}