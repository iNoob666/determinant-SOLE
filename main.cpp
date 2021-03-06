#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>

double fabs(double &d);

using namespace std;


double f(int x, int y, int z, int t)
{
    return 0;
}

double u(int x, int y, int z, double t){
    return x;
}

int X_coef(int i){
    return (i % 2);
}

int Y_coef(int i){
    return ((int)(i / 2)) % 2 ;
}

int Z_coef(int i){
    return (int)(i / 4);
}

std::vector<std::vector<double>> formxyz(const vector<double >& x, const vector<double >& y, const vector<double >& z) {
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

std::vector<std::vector<double>> formnvtr(const vector<double >& x, const vector<double >& y, const vector<double >& z) {
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
                     << ix + count + xsize * ysize << ' '
                     << ix + count + xsize * ysize + 1 << ' '
                     << ix + xsize + xsize * ysize + count << ' '
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


vector<vector<double>> locG(double koef, double x0, double x1, double y0, double y1, double z0, double z1){ //локальная матрица жёсткости для конечного элемента
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
    int size = G.size();
    for (int i = 0; i < size; ++i)
        G[i].resize(8);

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
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

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            M[i][j] = (koef/216) * ((hy*hz*hx) * m[X_coef(i)][X_coef(j)] * m[Y_coef(i)][Y_coef(j)] * m[Z_coef(i)][Z_coef(j)]);
        }
    }
    return M;
}


vector<vector<double>> GlobalG(vector<double> koef, const vector<double >& x, const vector<double >& y, const vector<double >& z, const vector<vector<double>>& nvtr, vector<vector<double>> xyz)
{
    int _size = x.size()*y.size()*z.size();
    int countOfElements = (x.size() - 1) * (y.size() - 1) * (z.size() - 1);
    vector<vector<double>> G;
    G.resize(_size);
    for (int i = 0; i < _size; ++i)
        G[i].resize(_size);

    for(int k = 0; k < countOfElements; ++k) {
        vector<vector<double>> tmp = locG(koef[k], xyz[nvtr[k][0]][0], xyz[nvtr[k][1]][0], xyz[nvtr[k][0]][1], xyz[nvtr[k][2]][1], xyz[nvtr[k][0]][2], xyz[nvtr[k][4]][2]);
        for(int i = 0; i < 8; ++i) {
            for(int j = 0; j < 8; ++j) {
                G[nvtr[k][i]][nvtr[k][j]] += tmp[i][j];
            }
        }
    }
    return G;
}


vector<double> MultMatrix(vector<vector<double>> matrix, vector<double> vect)
{
    if(matrix.size() == vect.size()) {
        vector<double> out;
        out.resize(matrix.size());
        for (int ix = 0; ix < matrix.size(); ix++) {
            out[ix] = 0;
            for (int jx = 0; jx < matrix.size(); jx++)
                out[ix] += matrix[ix][jx] * vect[jx];
        }
        return out;
    }
    else
    {
        cout << "Error in multiplying matrix by vector" << endl;
        exit(1);
    }
}


vector<vector<double>> SummMatrix(vector<vector<double>> B, const vector<vector<double>> C){
    vector<vector<double>> A;
    int n = B.size();
    A.resize(n);
    for(auto & i : A){
        i.resize(n);
    }

    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            A[i][j] = B[i][j] + C[i][j];
        }
    }
    return A;
}


vector<vector<double>> GlobalM(vector<double> koef, const vector<double >& x, const vector<double >& y, const vector<double >& z, const vector<vector<double>>& nvtr, vector<vector<double>> xyz)
{
    int _size = x.size()*y.size()*z.size();
    int countOfElements = (x.size() - 1) * (y.size() - 1) * (z.size() - 1);

    vector<vector<double>> M;
    M.resize(_size);
    for (int i = 0; i < _size; ++i)
        M[i].resize(_size);

    for(int k = 0; k < countOfElements; ++k) {
        vector<vector<double>> tmp = locM(koef[k], xyz[nvtr[k][0]][0], xyz[nvtr[k][1]][0], xyz[nvtr[k][0]][1], xyz[nvtr[k][2]][1], xyz[nvtr[k][0]][2], xyz[nvtr[k][4]][2]);
        for(int i = 0; i < 8; ++i) {
            for(int j = 0; j < 8; ++j) {
                M[nvtr[k][i]][nvtr[k][j]] += tmp[i][j];
            }
        }
    }
    return M;
}


vector<double> GlobalB(const vector<double >& x, const vector<double >& y, const vector<double >& z, const vector<vector<double>>& nvtr, vector<vector<double>> xyz, int t){
    int countOfElements = (x.size() - 1) * (y.size() - 1) * (z.size() - 1);
    vector<double> F;
    F.resize(x.size()*y.size()*z.size());

    for(int k = 0; k < countOfElements; ++k) {
        vector<vector<double>> tmp = locM(1, xyz[nvtr[k][0]][0], xyz[nvtr[k][1]][0], xyz[nvtr[k][0]][1], xyz[nvtr[k][2]][1], xyz[nvtr[k][0]][2], xyz[nvtr[k][4]][2]);
        vector<double> F1;
        F1.resize(8);
        for (int i = 0; i < 8; ++i) {
            F1[i] = f(xyz[nvtr[k][i]][0], xyz[nvtr[k][i]][1], xyz[nvtr[k][i]][2], t);
        }
        F1 = MultMatrix(tmp, F1);
        for (int i = 0; i < 8; ++i) {
            F[nvtr[k][i]] += F1[i];
        }
    }
    return F;
}


void FirstBoundaryConditions(vector<vector<double>> &A, vector<double> &F, const vector<vector<double>> &xyz, double t)
{
    int size = A.size();
    for (int i = 0; i < size; ++i) {
        if (xyz[i][0] == xyz[0][0] || xyz[i][1] == xyz[0][1] || xyz[i][2] == xyz[0][2] || xyz[i][0] == xyz[size - 1][0] || xyz[i][1] == xyz[size - 1][1] || xyz[i][2] == xyz[size - 1][2]) {
            F[i] = u(xyz[i][0], xyz[i][1], xyz[i][2], t);
            for (int j = 0; j < size; ++j) {
                if (i == j)
                    A[i][j] = 1;
                else
                    A[i][j] = 0;
            }
        }
    }
}


void luDecompose(vector<vector<double>> &A, vector<double> &F){
    int n = A.size();

    for(int i = 1; i < n; ++i){
        A[i][0] /= A[0][0];
    }

    for(int i = 1; i < n; ++i){
        for(int j = 1; j < n; ++j){
            if(i > j){
                double sum = 0;
                for(int k = 0; k < j; k++){
                    sum += A[i][k] * A[k][j];
                }
                A[i][j] = (A[i][j] - sum) / A[j][j];
            } else{
                double sum = 0;
                for(int k = 0; k < i; k++){
                    sum += A[i][k] * A[k][j];
                }
                A[i][j] -= sum;
            }
        }
    }

    for(int i = 1; i < n; ++i){
        double sum = 0;
        for(int k = 0; k < i; ++k){
            sum += F[k] * A[i][k];
        }
        F[i] -= sum;
    }

    F[n - 1] /= A[n - 1][n - 1];
    for(int i = n - 2; i >= 0; --i){
        double sum = 0;
        for(int k = i + 1; k < n; ++k){
            sum += A[i][k] * F[k];
        }
        F[i] = (F[i] - sum) / A[i][i];
    }
}

vector<vector<double>> multiplyMatrixByConst(vector<vector<double>> &A, double c){
    vector<vector<double>> tmp = A;
    for(auto & i : tmp){
        for(auto & j : i){
            j *= c;
        }
    }
    return tmp;
}

vector<double> summVectors(vector<double> &a, vector<double> &b){
    if(a.size() == b.size()) {
        vector<double> tmp;
        tmp.resize(a.size());
        for (int i = 0; i < a.size(); ++i){
            tmp[i] = a[i] + b[i];
        }
        return tmp;
    }
    else{
        cout << "Wrong size in summVectors" << endl;
        exit(1);
    }
}

vector<double> operator+(vector<double> a, vector<double> b){
    if(a.size() == b.size()) {
        for (int i = 0; i < a.size(); ++i) {
            a[i] += b[i];
        }
        return a;
    }
    else
    {
        cout << "Wrong size in operator+" << a.size() << ' ' << b.size() << endl;
        exit(1);
    }
}

vector<double> operator-(vector<double> a, vector<double> b){
    if(a.size() == b.size()) {
        for (int i = 0; i < a.size(); ++i) {
            a[i] -= b[i];
        }
        return a;
    }
    else
    {
        cout << "Wrong size in operator-" << endl;
        exit(1);
    }
}


int main() {
    vector<double> x, y, z, lyambda, sigma, xhi, t;
    vector<vector<double>> q;

    ifstream fin("x.txt");
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

    fin.open("lyambda.txt");
    while (!fin.eof()) {
        double tmp;
        fin >> tmp;
        lyambda.push_back(tmp);
    }
    fin.close();

    fin.open("sigma.txt");
    while (!fin.eof()) {
        double tmp;
        fin >> tmp;
        sigma.push_back(tmp);
    }
    fin.close();

    fin.open("xhi.txt");
    while (!fin.eof()) {
        double tmp;
        fin >> tmp;
        xhi.push_back(tmp);
    }
    fin.close();

    fin.open("t.txt");
    while (!fin.eof()) {
        double tmp;
        fin >> tmp;
        t.emplace_back(tmp);
    }
    fin.close();

    int n = (x.size() - 1) * (y.size() - 1) * (z.size() - 1);

    if(lyambda.size() == 1){
        for(int i = 1; i < n; ++i){
            lyambda.push_back(lyambda[0]);
        }
    }
    if(sigma.size() == 1){
        for(int i = 1; i < n; ++i){
            sigma.push_back(sigma[0]);
        }
    }
    if(xhi.size() == 1){
        for(int i = 1; i < n; ++i){
            xhi.push_back(xhi[0]);
        }
    }

    std::vector<std::vector<double>> xyz = formxyz(x, y, z);
    std::vector<std::vector<double>> nvtr = formnvtr(x, y, z);

    q.resize(t.size());
    int count = x.size() * y.size() * z.size();
    for(auto & vect : q){
        vect.resize(count);
    }

    for(int k = 0; k < n; ++k) {
        for (int i = 0; i < 8; ++i) {
            q[0][nvtr[k][i]] = u(xyz[nvtr[k][i]][0], xyz[nvtr[k][i]][1], xyz[nvtr[k][i]][2], t[0]);
            q[1][nvtr[k][i]] = u(xyz[nvtr[k][i]][0], xyz[nvtr[k][i]][1], xyz[nvtr[k][i]][2], t[1]);
        }
    }

    for(int i = 2; i < t.size(); ++i) {
        double dt = 1. / (t[i] - t[i - 1]);
        // Левая часть
        vector<vector<double>> Mxhi = GlobalM(xhi, x, y, z, nvtr, xyz);
        vector<vector<double>> Msigma = GlobalM(sigma, x, y, z, nvtr, xyz);
        vector<vector<double>> A = SummMatrix(multiplyMatrixByConst(Mxhi, dt * dt), multiplyMatrixByConst(Msigma, dt / 2));
        //-----------------------------------------------------
        // Правая часть
        vector<vector<double>> G = GlobalG(lyambda, x, y, z, nvtr, xyz);
        vector<double> b = GlobalB(x, y, z, nvtr, xyz, t[i - 1]);

        //FirstBoundaryConditions(A, b, xyz, t[i]);// учет первых краевых условий для формирования b
            // Формирование полной правой части d
            vector<double> d = b + MultMatrix(multiplyMatrixByConst(Mxhi, dt * dt * 2), q[i - 1]) -
                                    MultMatrix(multiplyMatrixByConst(Mxhi, dt * dt), q[i - 2]) +
                                    MultMatrix(multiplyMatrixByConst(Msigma, dt / 2), q[i - 2]) -
                                    MultMatrix(G, q[i - 1]);
            //-------------------------------------------------
        //-----------------------------------------------------
        // Решение
        luDecompose(A, d);
        //-----------------------------------------------------
        // Вывод ответа на временном отрезке t[i]
        cout << t[i] << ':' << endl;
        for (auto &elem : d) {
            cout << fixed << elem << endl;
        }
        cout << endl << endl;

        q[i] = d;// Сохранение данных о решении на этом временном отрезке
        //-----------------------------------------------------
    }
    return 0;
}