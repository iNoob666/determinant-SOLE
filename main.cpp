#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>


std::vector<std::vector<double>> formxyz(const std::vector<double >& x, const std::vector<double >& y, const std::vector<double >& z){
    std::ofstream fout("xyz.txt");
    std::vector<std::vector<double>> tmp;
    int count = 0;
    for(double iz : z){
        for(double iy : y){
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

std::vector<std::vector<double>> formnvtr(const std::vector<double >& x, const std::vector<double >& y, const std::vector<double >& z){
    std::vector<std::vector<double>> tmp;
    std::vector<double> a;
    std::ofstream fout("nvtr.txt");
    int xsize = x.size();
    int ysize = y.size();
    int zsize = z.size();
    int count = 0;
    for(size_t iz = 0; iz < zsize - 1; ++iz) {
        count = (int)iz * xsize * ysize;
        for (size_t iy = 0; iy < ysize - 1; ++iy) {
            count += (int)iy * xsize;
            for (size_t ix = 0; ix < xsize - 1; ++ix) {
                fout << a.emplace_back(ix + count) << ' '
                     << a.emplace_back(ix + count + 1) << ' '
                     << a.emplace_back(ix + xsize + count) << ' '
                     << a.emplace_back(ix + xsize + count + 1) << ' '
                     << a.emplace_back(ix + count + xsize * ysize) << ' '
                     << a.emplace_back(ix + count + xsize * ysize + 1) << ' '
                     << a.emplace_back(ix + xsize + xsize * ysize + count) << ' '
                     << a.emplace_back(ix + xsize + xsize * ysize + count + 1) << std::endl;
                tmp.push_back(a);
                a.clear();
            }
        }
    }
    fout.close();
    return tmp;
}


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
    fin.close();

    std::vector<std::vector<double>> xyz = formxyz(x, y, z);

    std::vector<std::vector<double>> nvtr = formnvtr(x, y, z);

    std::vector<std::vector<double>> listOfConnections;

    int count = x.size() * y.size() * z.size();
    for(size_t i = 0; i < count; ++i){
        std::vector<double> a;
        listOfConnections.push_back(a);
    }

    for(const auto& element : nvtr){
        for(size_t i = 0; i < element.size(); ++i){
            for(size_t j = 0; j < i; ++j){
                listOfConnections[element[i]].push_back(element[j]);
            }
        }
    }
    for (auto & item : listOfConnections) {
        std::sort(item.begin(), item.end());
    }
    for(auto & item : listOfConnections){
        for(auto i = ++item.begin(); i < item.end(); ++i){
            if(*i == *(i - 1)){
                item.erase(i);
            }
        }
    }


    return 0;
}