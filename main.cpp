#include <iostream>
#include <fstream>
#include <vector>


int main() {
    std::ifstream fin("x.txt");
    std::vector<double> x, y, z;
    while(!fin.eof()){
        double tmp;
        fin >> tmp;
        x.emplace_back(tmp);
    }
    fin.close();
    fin.open("y.txt");
    while(!fin.eof()){
        double tmp;
        fin >> tmp;
        y.emplace_back(tmp);
    }
    fin.close();
    fin.open("z.txt");
    while(!fin.eof()){
        double tmp;
        fin >> tmp;
        z.emplace_back(tmp);
    }

    std::ofstream fout("xyz.txt");

    int count = 0;
    for(double iz : z){
        for(double iy : y){
            for (double ix : x) {
                fout << count++ << ' ' << ix << ' ' << iy << ' ' << iz << std::endl;
            }
        }
    }
    fout.close();

    fout.open("nvtr.txt");

    for(size_t iz = 0; iz < z.size() - 1; ++iz) {
        for (size_t iy = 0; iy < y.size() - 1; ++iy) {
            for (size_t ix = 0; ix < x.size() - 1; ++ix) {
                fout << ix + iy * x.size() + iz * x.size() * y.size() << ' '
                << ix + iy * x.size() + iz * x.size() * y.size() + 1 << ' '
                << ix + (iy + 1) * x.size() + iz * x.size() * y.size() << ' '
                << ix + (iy + 1) * x.size() + iz * x.size() * y.size() + 1 << ' '
                << ix + iy * x.size() + (iz + 1) * x.size() * y.size() << ' '
                << ix + iy * x.size() + (iz + 1) * x.size() * y.size() + 1 << ' '
                << ix + (iy + 1) * x.size() + (iz + 1) * x.size() * y.size() << ' '
                << ix + (iy + 1) * x.size() + (iz + 1) * x.size() * y.size() + 1 << std::endl;
            }
        }
    }

    return 0;
}