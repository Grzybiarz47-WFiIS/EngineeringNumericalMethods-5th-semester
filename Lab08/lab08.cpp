#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>

typedef std::vector <double> vecD;
typedef std::vector<vecD> grid;

const int nx = 400, ny = 90, I1 = 200, I2 = 210, J1 = 50, IT_PICARD = 20, IT_MAX = 10000;
const double delta = 0.01, xA = 0.45, yA = 0.45;
const double sigma = 10*delta;
std::string inputFile = "psi.dat";
std::string outputVFile = "out1.dat";
std::string outputUFile = "out2.dat";
std::string outputCXFile = "out3.dat";

inline double u0(int i, int j)
{
    double x = i*delta, y = j*delta;
    double temp1 = 1/(2*M_PI*pow(sigma, 2));
    double temp2 = exp(-(pow(x-xA, 2)+pow(y-yA, 2))/(2*pow(sigma, 2)));
    return temp1*temp2;
}
inline double expr(double a, double b)
{
    return (a-b)/(2*delta);
}
bool checkEdge(int i, int j)
{
    if((i >= I1 && i <= I2) && (j <= J1))
        return false;
    return true;
}
void readFile(grid& psi)
{
    std::ifstream fin(inputFile);
    int i, j;
    double x;
    while(fin >> i >> j >> x)
        psi[i][j] = x;
    fin.close();
}
void print(grid& psi, std::ostream& fout)
{
    for(int i = 0; i <= nx; ++i)
    {
        for(int j = 0; j <= ny; ++j)
            fout << i*delta << " " << j*delta << " " << psi[i][j] << std::endl;
        fout << std::endl;
    }
    fout << "\n\n";
}
void prepareU(grid& u)
{
    for(int i = 0; i <= nx; ++i)
        for(int j = 0; j <= ny; ++j)
            u[i][j] = u0(i, j);
}
double findV(grid& vx, grid& vy, grid& psi)
{
    for(int i = 1; i < nx; ++i)
        for(int j = 1; j < ny; ++j)
        {
            vx[i][j] = (psi[i][j+1]-psi[i][j-1])/(2*delta);
            vy[i][j] = -(psi[i+1][j]-psi[i-1][j])/(2*delta);
        }
    for(int i = I1; i <= I2; ++i)
        for(int j = 0; j <= J1; ++j)
            vx[i][j] = vy[i][j] = 0.0;
    for(int i = 1; i < nx; ++i)
        vx[i][0] = vy[i][ny] = 0.0;
    for(int j = 0; j <= ny; ++j)
    {
        vx[0][j] = vx[1][j];
        vx[nx][j] = vx[nx-1][j];
    }
    double vmax = sqrt(pow(vx[0][0], 2)+pow(vy[0][0], 2));    
    for(int i = 0; i <= nx; ++i)
        for(int j = 0; j <= ny; ++j)
        {
            double temp = sqrt(pow(vx[i][j], 2)+pow(vy[i][j], 2));
            vmax = temp > vmax ? temp : vmax;
        }
    return vmax;
}
std::pair <double, double> findCX(grid& u)
{
    std::pair <double, double> p(0, 0);
    for(int i = 0; i <= nx; ++i)
        for(int j = 0; j <= ny; ++j)
        {
            p.first += u[i][j];
            p.second += u[i][j]*i*delta;
        }
    return p;
}
void solve(grid& u, grid& vx, grid& vy, const double D, const double dt, std::ostream& fout2, std::ostream& fout3)
{
    grid newU(u);
    for(int it = 1; it <= IT_MAX; ++it)
    {
        for(int k = 1; k <= IT_PICARD; ++k)
            for(int i = 0; i <= nx; ++i)
                for(int j = 1; j < ny; ++j)
                    if(checkEdge(i, j))
                    {
                        double factor = 1/(1+(2*D*dt/pow(delta, 2)));
                        double first, second, third, fourth, fifth;
                        first = u[i][j];
                        third = 0.5*dt*vy[i][j]*(expr(u[i][j+1], u[i][j-1])+expr(newU[i][j+1], newU[i][j-1]));
                        if(i == 0)
                        {
                            second = 0.5*dt*vx[i][j]*(expr(u[i+1][j], u[nx][j])+expr(newU[i+1][j], newU[nx][j]));
                            fourth = 0.5*dt*D*(u[i+1][j]+u[nx][j]+u[i][j+1]+u[i][j-1]-4*u[i][j])/pow(delta, 2);
                            fifth = 0.5*dt*D*(newU[i+1][j]+newU[nx][j]+newU[i][j+1]+newU[i][j-1])/pow(delta, 2);
                        }
                        else if(i == nx)
                        {
                            second = 0.5*dt*vx[i][j]*(expr(u[0][j], u[i-1][j])+expr(newU[0][j], newU[i-1][j]));
                            fourth = 0.5*dt*D*(u[0][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-4*u[i][j])/pow(delta, 2);
                            fifth = 0.5*dt*D*(newU[0][j]+newU[i-1][j]+newU[i][j+1]+newU[i][j-1])/pow(delta, 2);
                        }
                        else
                        {
                            second = 0.5*dt*vx[i][j]*(expr(u[i+1][j], u[i-1][j])+expr(newU[i+1][j], newU[i-1][j]));
                            fourth = 0.5*dt*D*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-4*u[i][j])/pow(delta, 2);
                            fifth = 0.5*dt*D*(newU[i+1][j]+newU[i-1][j]+newU[i][j+1]+newU[i][j-1])/pow(delta, 2);
                        }
                        newU[i][j] = factor*(first-second-third+fourth+fifth);
                    }
        u = newU;
        std::pair<double, double> cxPair = findCX(u);
        fout3 << it*dt << " " << pow(delta, 2)*cxPair.first << " " << pow(delta, 2)*cxPair.second << std::endl;
        if(it % (IT_MAX/20) == 0)
        {
            std::cout << "Iteracja: " << it << std::endl;
            print(u, fout2);
        }
    }
    fout3 << "\n\n";
}
int main()
{
    clock_t t = clock();

    //v
    std::ofstream fout1(outputVFile);
    
    grid psi(nx+1, vecD(ny+1));
    grid vx(nx+1, vecD(ny+1)), vy(nx+1, vecD(ny+1));
    double vmax, deltaT;

    readFile(psi);
    vmax = findV(vx, vy, psi);
    deltaT = delta/(4*vmax);

    print(vx, fout1);
    print(vy, fout1);

    fout1.close();

    //u
    std::ofstream fout2(outputUFile);
    std::ofstream fout3(outputCXFile);
    grid u(nx+1, vecD(ny+1));

    prepareU(u);
    solve(u, vx, vy, 0.0, deltaT, fout2, fout3);

    prepareU(u);
    solve(u, vx, vy, 0.1, deltaT, fout2, fout3);

    fout2.close();
    fout3.close();

    std::cout << "Czas: " << double(clock()-t)/CLOCKS_PER_SEC << std::endl;
    return 0;
}