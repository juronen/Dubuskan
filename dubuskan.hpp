#ifndef DUBUSKAN
#define DUBUSKAN

#define SpaceType template <size_t dims, typename T>

#include <string>
#include <set>
#include <map>
#include <deque>
#include <algorithm>
#include <array>
#include <iostream>
#include <vector>
#include <stack>
#include <cmath>
#include <unordered_map>
#include <sstream>

using std::size_t;

SpaceType struct Cluster;

SpaceType struct DataPoint
{
    Cluster<dims, T>* cluster;
    T data;
    std::array<double, dims> vec;

    DataPoint(T d) : data(d) {}
       
    template<class... Vec>
    DataPoint(T d, Vec... v) : vec{v...}, data(d) {}

    double fast_distance(DataPoint* d)
    {
        double sum = 0;
        for (size_t i = 0; i < dims; i++)
        {
            double diff = d->vec[i] - vec[i];
            sum += diff * diff;
        }
        return sum;
    }

    bool can_reach(DataPoint* d, double& dist) 
    {
        return fast_distance(d) <= dist;     
    }    
};

SpaceType struct Cluster
{
    std::deque<DataPoint<dims, T>> points;

    Cluster() {}

    void add(DataPoint<dims, T> p)
    {
        points.push_back(p);
    }    
};

SpaceType class Cell
{
private:
    std::deque<DataPoint<dims, T>> data;

    bool can_host(DataPoint<dims, T> p, double dist, Cluster<dims, T>** c)
    {
        return std::any_of(data.begin(), data.end(),
            [&c, &p, &dist, this] (DataPoint<dims, T> cell_point)
            {
                if (p.can_reach(&cell_point, dist))
                {
                    *c = cell_point.cluster;
                    return true;                
                }
                return false;
            }
        );
    }

public:
    Cell() {}

    void add(DataPoint<dims, T> p)
    {
        data.push_back(p);
    }

    Cluster<dims, T>* find_cluster(DataPoint<dims, T> p, double dist)
    {
        Cluster<dims, T>* cluster;
        if (can_host(p, dist, &cluster))
        {
            return cluster;
        }
        return nullptr;
    }
};

typedef std::vector<std::vector<size_t>> size_2D;

SpaceType class Space
{

private:
    double max_dist;
    double cell_side_length;
    double average_pts;
    long total_pts;
    long cell_count_target;
    std::unordered_map<std::string, Cell<dims, T>> cells;
    
    void recompute_average()
    {
        average_pts = (double)total_pts / (double)cells.size();
    }

    std::array<size_t, dims> to_lattice(std::array<double, dims> &arr)
    {
        std::array<size_t, dims> floored;
        std::transform(arr.begin(), arr.end(), floored.begin(),
                [this] (size_t coord)
                {
                    return floor(coord / this->cell_side_length);
                }
        );
        return floored; 
    }

    std::string make_key(std::array<double, dims> &arr)
    {
        auto floored = to_lattice(arr);
        std::ostringstream oss;
        for (auto &coord : floored)
        {
            oss << " " << coord;
        }
        return oss.str();
    }

    void add_point(DataPoint<dims, T> &p, std::string key)
    {
        cells[key].add(p);
    }


public:
    std::vector<Cluster<dims, T>*> clusters;

    Space(long cell_count_target, double cell_side_length, double max_dist) 
    {
        this->max_dist = max_dist;
        this->cell_count_target = cell_count_target;
        this->cell_side_length = cell_side_length;
        std::array<double, dims> zero{};
        Cell<dims, T> cell;
        cells[make_key(zero)] = cell;
    }

    std::vector<std::array<double, dims>> cartesian_product(size_2D opts, long n)
    {
        long num_neighbors = std::pow(n, opts.size());
        // Avoid a resize by computing the needed number of elements
        std::vector<std::array<double, dims>> neighbors(num_neighbors);
        std::vector<std::vector<size_t>::iterator> iterators;
        for (auto &a : opts) 
        {
            iterators.push_back(a.begin());
        }
        bool flag = false;
        size_t idx = 0;
        while (!flag)
        {
            std::array<double, dims> neighbor;
            std::transform(iterators.begin(), iterators.end(), neighbor.begin(),
                    [] (std::vector<size_t>::iterator it)
                    {
                        return *it;
                    }
            );
            neighbors[idx++] = neighbor;
            for (size_t l = 0; l < iterators.size(); l++)
            {
                if (++iterators[l] != opts[l].end())
                    break;
                else 
                {
                    if (l < iterators.size() - 1)
                        iterators[l] = opts[l].begin();
                    else
                        flag = true;
                }
            }
        }
        return neighbors;
    }

    std::set<Cluster<dims, T>*> find_clusters(DataPoint<dims, T> &p)
    {
        long span = max_dist/cell_side_length;
        // Round the coordinates to get a cell location in terms of
        // the cell side lengths
        std::array<size_t, dims> cell_loc = to_lattice(p.vec);
        // Get all the neighboring cells. This means for each dimension,
        // we get the previous and next coordinate for a total of 3,
        // and all the possible choices form a square in 2 dimensions,
        // a cube in 3 dimensions etc... 
        std::vector<std::vector<size_t>> opts(dims);
        for (size_t dim = 0; dim < dims; dim++)
        {
            for (long offset = -span; offset <= span; offset++)
            {
                opts[dim].push_back(cell_loc[dim] + offset);
            }
        }
        auto neighbors = cartesian_product(opts, 2 * span + 1);
        std::set<Cluster<dims, T>*> result;
        for (auto &neighbor : neighbors)
        {
            std::string key = make_key(neighbor);
            if (cells.count(key))
            {
                auto cell = cells[key];
                Cluster<dims, T>* cluster;
                if ((cluster = cell.find_cluster(p, max_dist)) != nullptr)
                {
                    result.insert(cluster);
                }
            }     
        }
        return result;   
    }

    void add(DataPoint<dims, T> &p)
    {
        std::string key = make_key(p.vec);
        std::set<Cluster<dims, T>*> found = find_clusters(p);
        Cluster<dims, T>* cluster = new Cluster<dims, T>();
        switch (found.size())
        {
            case 0:
                cluster->add(p);
                p.cluster = cluster;
                clusters.push_back(cluster);
                break;
            case 1:
                delete cluster;
                cluster = *clusters.begin();
                cluster->add(p);
                p.cluster = cluster;
                break;
            default:
                for (auto it : found)
                {
                    for (auto &pt : it->points)
                    {
                        cluster->add(pt);
                        pt.cluster = cluster;
                    }
                    delete it;
                }
                clusters.push_back(cluster); 
        }
        cells[key].add(p);
    }

};

#endif
