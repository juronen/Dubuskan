#include "dubuskan.hpp"
#include <time.h>
#include <ctime>
#include <iostream>
#include <vector>
#include <deque>
#include <unordered_map>
#include <string>
using std::vector;

DataPoint<2, long> mkpt(long d, double x, double y)
{
    DataPoint<2, long> p(d, x, y);
    return p;
}

int main() 
{
    auto t = time(NULL);
    srand(0);
    std::cout << t << std::endl;
    Space<2, long> space(10, 100, 100, 1000);
    for (long l = 0; l < 30000; l++)
    {
        long x = rand() % 1000;
        long y = rand() % 1000;
        auto p = mkpt(l, x, y);
        space.add(p);
    } 
    int num = 0;
    for (auto &c : space.clusters)
    {
        std::cout << num++ << ": " << c.second->points.size() << std::endl;
    }

}
