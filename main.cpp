#include <iostream>
#include <fstream>
#include <vector>


struct point{
    void addPoint(double x, double y, double z, int count){
        this->x = x;
        this->y = y;
        this->z = z;
        globalPos = count;
    }

    double getX(){
        return x;
    }

    double getY(){
        return y;
    }

    double getZ(){
        return z;
    }

    double getGlobalPos(){
        return globalPos;
    }

private:
    double x;
    double y;
    double z;
    int globalPos;
};

struct cube{
    cube(){
        for(size_t i = 0; i < 8; ++i){
            localPos[i] = 0;
            globalPos[i] = 0;
        }
        points = nullptr;
        size = 0;
    }

    ~cube(){
        for(size_t i = 0; i < 8; ++i){
            localPos[i] = 0;
            globalPos[i] = 0;
        }
        delete [] points;
        size = 0;
    }

    void addPoint(point a){
        if(size != 0){
            point *tmp = points;
            delete [] points;
            points = new point [size + 1];
            for(size_t i = 0; i < size; ++i){
                points[i] = tmp[i];
            }
            points[size] = a;
            size++;
        }
        else{

        }
    }
private:
    int size;
    int localPos [8];
    int globalPos [8];
    point *points;
};

struct web{
private:
    cube* elements;
};

int main() {
    std::ifstream fin("x.txt");
    std::vector<double> x, y, z;
    std::vector<point> a;
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
//    for(int i = 0; i < x.size(); ++i)
//        std::cout << x.at(i) << ' ';
//    std::cout << std::endl;
//    for(int i = 0; i < y.size(); ++i)
//        std::cout << y.at(i) << ' ';
//    std::cout << std::endl;
//    for(int i = 0; i < z.size(); ++i)
//        std::cout << z.at(i) << ' ';

    int count = 1;
    for(int iz = 0; iz < z.size(); ++iz){
        for(int iy = 0; iy < y.size(); ++iy){
            for (int ix = 0; ix < x.size(); ++ix) {
                point tmp;
                tmp.addPoint(x.at(ix), y.at(iy), z.at(iz), count++);
                a.emplace_back(tmp);
                std::cout << a.at(count - 2).getX() << ' ' << a.at(count - 2).getY() << ' ' << a.at(count - 2).getZ() << ' ' << a.at(count - 2).getGlobalPos();
                std::cout << std::endl;
            }
        }
    }
    std::cout << "___" << std::endl;
    unsigned int sum = (x.size() - 1) * (y.size() - 1) * (z.size() - 1);

    return 0;
}