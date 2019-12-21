#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <conio.h>
using namespace std;


int X_coef(int i)
{
    return (i % 2);
}

int Y_coef(int i)
{
    return ((int)(i / 2)) % 2 ;
}

int Z_coef(int i)
{
    return (int)(i / 4);
}

std::vector<std::vector<double>> formxyz(const std::vector<double >& x, const std::vector<double >& y, const std::vector<double >& z) {
    std::ofstream fout("xyz.txt");
    std::vector<std::vector<double>> tmp;
    int count = 0;
    for (double iz : z) {
        for (double iy : y) {
            for (double ix : x) {
                std::vector<double> tmpxyz;
                tmpxyz.emplace_back(ix);
                tmpxyz.emplace_back(iy);
                tmpxyz.emplace_back(iz);
                tmp.emplace_back(tmpxyz);
                fout << count++ << ' ' << ix << ' ' << iy << ' ' << iz << std::endl;
            }
        }
    }
    fout.close();
    return tmp;
}

std::vector<std::vector<double>> formnvtr(const std::vector<double >& x, const std::vector<double >& y, const std::vector<double >& z) {
    std::vector<std::vector<double>> tmp;
    std::vector<double> a;
    std::ofstream fout("nvtr.txt");
    int xsize = x.size();
    int ysize = y.size();
    int zsize = z.size();
    int count = 0;
    int b = 0;
    for (size_t iz = 0; iz < zsize - 1; ++iz) {
        count = (int)iz * xsize * ysize;
        for (size_t iy = 0; iy < ysize - 1; ++iy) {
            count += (int)iy * xsize;
            for (size_t ix = 0; ix < xsize - 1; ++ix) {
                a.emplace_back(ix + count);
                a.emplace_back(ix + count + 1);
                a.emplace_back(ix + xsize + count);
                a.emplace_back(ix + xsize + count + 1);
                a.emplace_back(ix + count + xsize * ysize);
                a.emplace_back(ix + count + xsize * ysize + 1);
                a.emplace_back(ix + xsize + xsize * ysize + count);
                a.emplace_back(ix + xsize + xsize * ysize + count + 1);
                fout << ix + count << ' '
                     << ix + count + 1 << ' '
                     << ix + xsize + count << ' '
                     << ix + xsize + count + 1 << ' '
                     <<ix + count + xsize * ysize << ' '
                     << ix + count + xsize * ysize + 1 << ' '
                     <<ix + xsize + xsize * ysize + count << ' '
                     << ix + xsize + xsize * ysize + count + 1 << std::endl;
                tmp.push_back(a);
                a.clear();
            }
            count = (int)iz * xsize * ysize;
        }
    }
    fout.close();
    return tmp;
}



vector<vector<double>> locG(double koef, double x0, double x1, double y0, double y1, double z0, double z1) //локальная матрица жёсткости для конечного элемента
{
    const int g[2][2] = {
            {1, -1},
            {-1, 1}
    };
    const double m[2][2] = { // коэффициент 1/6 перенесён в место пользования
            {2, 1 },
            {1, 2}
    };
    double hx = x1 - x0;
    double hy = y1 - y0;
    double hz = z1 - z0;
    vector<vector<double>> G;
    G.resize(8);
    for (int i = 0, size = G.size(); i < size; ++i)
        G[i].resize(8);

    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            G[i][j] = (koef/36)*((hy*hz/hx) * g[X_coef(i)][X_coef(j)] * m[Y_coef(i)][Y_coef(j)] * m[Z_coef(i)][Z_coef(j)] +
                                 (hx*hz/hy) * m[X_coef(i)][X_coef(j)] * g[Y_coef(i)][Y_coef(j)] * m[Z_coef(i)][Z_coef(j)]+
                                 (hx*hy/hz) * m[X_coef(i)][X_coef(j)] * m[Y_coef(i)][Y_coef(j)] * g[Z_coef(i)][Z_coef(j)]);
        }
    }
    return G;
}

vector<vector<double>> locM(double koef, double x0, double x1, double y0, double y1, double z0, double z1) //локальная матрица массы для конечного элемента
{
    const double m[2][2] = { // коэффициент 1/6 перенесён в место пользования
            {2, 1 },
            {1, 2}
    };
    double hx = x1 - x0;
    double hy = y1 - y0;
    double hz = z1 - z0;
    vector<vector<double>> M;
    M.resize(8);
    for (int i = 0, size = M.size(); i < size; ++i)
        M[i].resize(8);

    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            M[i][j] = (koef/216) * ((hy*hz*hx) * m[X_coef(i)][X_coef(j)] * m[Y_coef(i)][Y_coef(j)] * m[Z_coef(i)][Z_coef(j)]);
        }
    }
    return M;
}

vector<vector<vector<double>>> GlobG3(double koef, const std::vector<double >& x, const std::vector<double >& y, const std::vector<double >& z) //глобальная матрица жёсткости для конечного элемента
{
    vector<vector<vector<double>>> G;
    int number = 0;
    double x0, x1, y0, y1, z0, z1;// 0 - предыдущий, 1 - текущий
    x0 = x[0];
    y0 = y[0];
    z0 = z[0];
    for (int i = 1; i < x.size(); i++)
    {
        for (int j = 1; j < y.size(); j++)
        {
            for (int k = 1; k < z.size(); k++)
            {
                G.emplace_back(locG(koef, x0, x[i], y0, y[j], z0, z[k]));
                z0 = z[k];
            }
            z0 = z[0];
            y0 = y[j];
        }
        y0 = y[0];
        x0 = x[i];
    }
    return G;
}

vector<vector<double>> GlobalG(double koef, const std::vector<double >& x, const std::vector<double >& y, const std::vector<double >& z, const std::vector<std::vector<double>>& nvtr)
{
    vector<vector<vector<double>>> G3 = GlobG3(koef, x, y, z);
    int _size = x.size()*y.size()*z.size();
    vector<vector<double>> G;
    for (int i = 0; i < _size; ++i) {
        vector<double> tmp;
        G.push_back(tmp);
    }
    for (auto &vect : G) {
        for (int i = 0; i < _size; ++i) {
            vect.push_back(0);
        }
    }
    for (int numElem = 0; numElem < G3.size(); numElem++) {
        for (int i = 0; i < G3[numElem].size(); i++) {
            for (int j = 0; j < G3[numElem].size(); j++) {
                if (G3[numElem][i][j] != 0) {
                    G[nvtr[numElem][i]][nvtr[numElem][j]] += G3[numElem][i][j];
                }
            }
        }
    }
    return G;
}

vector<vector<vector<double>>> GlobM3(double koef, const std::vector<double >& x, const std::vector<double >& y, const std::vector<double >& z) //глобальная матрица жёсткости для конечного элемента
{
    vector<vector<vector<double>>> M;
    double x0, x1, y0, y1, z0, z1;// 0 - предыдущий, 1 - текущий
    x0 = x[0];
    y0 = y[0];
    z0 = z[0];
    for (int i = 1; i < x.size(); i++)
    {
        for (int j = 1; j < y.size(); j++)
        {
            for (int k = 1; k < z.size(); k++)
            {
                M.emplace_back(locM(koef, x0, x[i], y0, y[j], z0, z[k]));
                z0 = z[k];
            }
            z0 = z[0];
            y0 = y[j];
        }
        y0 = y[0];
        x0 = x[i];
    }
    return M;
}
vector<vector<double>> GlobalM(double koef, const std::vector<double >& x, const std::vector<double >& y, const std::vector<double >& z, const std::vector<std::vector<double>>& nvtr)
{
    vector<vector<vector<double>>> M3 = GlobM3(koef, x, y, z);
    int _size = x.size()*y.size()*z.size();
    vector<vector<double>> M;
    for (int i = 0; i < _size; ++i) {
        vector<double> tmp;
        M.push_back(tmp);
    }
    for (auto &vect : M) {
        for (int i = 0; i < _size; ++i) {
            vect.push_back(0);
        }
    }
    for (int numElem = 0; numElem < M3.size(); numElem++) {
        for (int i = 0; i < M3[numElem].size(); i++) {
            for (int j = 0; j < M3[numElem].size(); j++) {
                if (M3[numElem][i][j] != 0) {
                    M[nvtr[numElem][i]][nvtr[numElem][j]] += M3[numElem][i][j];
                }
            }
        }
    }
    return M;
}

double f(int x, int y, int z)
{
    return x*x + y*y + z*z;
}

vector<double> MultMatrix(vector<vector<double>> matrix, vector<double> vector1)
{
    vector<double> out;
    out.resize(8);
    for (int ix = 0; ix < 8; ix++)
    {
        out[ix] = 0;
        for (int jx = 0; jx < 8; jx++)
            out[ix] += matrix[ix][jx] * vector1[jx];
    }
    return out;
}

vector<double> GlobF(double koef, vector<double> x, vector<double> y, vector<double> z, const vector<vector<double>>& nvtr, vector<vector<double>> xyz)
{
    vector<double> F;
    vector<double> F1;
    int size = x.size() * y.size() * z.size();
    F.resize(size);
    F1.resize(8);
    for (int i = 0; i < size; i++)
    {
        F[i] = 0;
    }
    vector<vector<double>> lM;
    lM.resize(8);
    for (int i = 0, size = lM.size(); i < size; ++i)
        lM[i].resize(8);

    for (const auto& element : nvtr)
    {
        for (size_t i = 0; i < 8; ++i)
        {
            F1[i] = f(xyz[element[i]][0], xyz[element[i]][1], xyz[element[i]][2]);
        }
        lM = locM(koef, xyz[element[0]][0], xyz[element[1]][0], xyz[element[0]][1], xyz[element[2]][1], xyz[element[0]][2], xyz[element[4]][2]);
        F1 = MultMatrix(lM, F1);
        for (size_t i = 0; i < element.size(); ++i)
        {
            F[element[i]] += F1[i];
        }
    }
    return F;
}


int main() {
    std::ifstream fin("x.txt");
    std::vector<double> x, y, z;
    while (!fin.eof()) {
        double tmp;
        fin >> tmp;
        x.emplace_back(tmp);
    }
    fin.close();
    fin.open("y.txt");
    while (!fin.eof()) {
        double tmp;
        fin >> tmp;
        y.emplace_back(tmp);
    }
    fin.close();
    fin.open("z.txt");
    while (!fin.eof()) {
        double tmp;
        fin >> tmp;
        z.emplace_back(tmp);
    }
    fin.close();
    int range = x.size()*y.size()*z.size();

    std::vector<std::vector<double>> xyz = formxyz(x, y, z);

    std::vector<std::vector<double>> nvtr = formnvtr(x, y, z);

    vector<vector<double>> G = GlobalG(216, x, y, z, nvtr);
    setlocale(LC_ALL, "Russian");
    printf_s("Матрица жёсткости:\n");
    for (int i = 0; i < range; i++)
    {
        for (int j = 0; j < range; j++)
        {
            printf_s("%3.0lf ", G[i][j]);
        }
        printf_s("\n");
    }
    printf_s("\n");

    vector<vector<double>> M = GlobalM(216, x, y, z, nvtr);
    printf_s("Матрица массы:\n");
    for (int i = 0; i < range; i++)
    {
        for (int j = 0; j < range; j++)
        {
            printf_s("%3.0lf ", M[i][j]);
        }
        printf_s("\n");
    }
    printf_s("\n");

    vector<double> F = GlobF(1, x, y, z, nvtr, xyz);
    for (int j = 0; j < range; j++)
    {
        printf_s("%3.7lf ", F[j]);
    }
    printf_s("\n");

    _getch();
    return 0;
}