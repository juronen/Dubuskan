#include "dubuskan.hpp"

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
    Space<2, long> space(10, 100, 5);
    for (long l = 0; l < 1000000; l++)
    {
        long x = rand() % 1000;
        long y = rand() % 1000;
        auto p = mkpt(l, x, y);
        space.add(p);
    } 
}
