#include "dubuskan.hpp"

#include <iostream>
#include <vector>
using std::vector;

int main() 
{
    Space<2, long> space(1, 100000, 100);
    DataPoint<2, long> p1(1, 1, 1);
    DataPoint<2, long> p2(2, 2, 2);
    DataPoint<2, long> p3(3, 3, 3);
    DataPoint<2, long> p4(4, 30, 30);
    auto v = {p1, p2, p3, p4};
    for (auto p : v)
    {
        space.add(p);
        std::cout << space.clusters.size() << std::endl;
    } 
}
