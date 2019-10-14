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

    int xsize = x.size();
    int ysize = y.size();
    int zsize = z.size();
    int tmp = 0;
    for(size_t iz = 0; iz < zsize - 1; ++iz) {
        tmp = iz * xsize * ysize;
        for (size_t iy = 0; iy < ysize - 1; ++iy) {
            tmp += iy * xsize;
            for (size_t ix = 0; ix < xsize - 1; ++ix) {
                fout << ix + tmp << ' '
                << ix + tmp + 1 << ' '
                << ix + xsize + tmp << ' '
                << ix + xsize + tmp + 1 << ' '
                << ix + tmp + xsize * ysize << ' '
                << ix + tmp + xsize * ysize + 1 << ' '
                << ix + xsize + xsize * ysize + tmp << ' '
                << ix + xsize + xsize * ysize + tmp + 1 << std::endl;
            }
        }
    }

    return 0;
}