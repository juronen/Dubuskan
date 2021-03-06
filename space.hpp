#ifndef DUBUSKAN_SPACE 
#define DUBUSKAN_SPACE

#define SpaceType template <size_t dims, typename T>

#include <unordered_set>
#include <list>
#include <memory>
#include <set>
#include <deque>
#include <algorithm>
#include <array>
#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <sstream>
#include <string.h>

using std::size_t;
using std::shared_ptr;
using std::make_shared;

SpaceType struct Cluster;
SpaceType struct DataPoint;

template<size_t dims, typename T> 
using SharedData = shared_ptr<DataPoint<dims, T>>;

SpaceType struct DataPoint
{
    Cluster<dims, T>* cluster;
    bool flag;
    T data;
    std::array<double, dims> vec;

    DataPoint(T d) : data(d) {}

    DataPoint(const DataPoint& p)
    {
        data = p.data;
        vec = p.vec;
    }
       
    template<class... Vec>
    DataPoint(T d, Vec... v) : vec{v...}, flag(false), data(d) {}

    double fast_distance(DataPoint<dims, T> *d)
    {
        double sum = 0;
        for (size_t i = 0; i < dims; i++)
        {
            double diff = d->vec[i] - vec[i];
            sum += diff * diff;
        }
        return sum;
    }

    bool can_reach(DataPoint<dims, T> *d, double& dist) 
    {
        return fast_distance(d) <= dist;     
    }    
};

SpaceType struct Cluster
{

    long id;

    std::list<SharedData<dims, T>> points;

    Cluster(long uid) : id(uid) {}

    void add(SharedData<dims, T> p)
    {
        points.push_front(p);
    }    
};

SpaceType struct Cell
{
    std::list<SharedData<dims, T>> data;

    bool can_host(DataPoint<dims, T>* p, 
                  double dist,
                  Cluster<dims, T>* &c,
                  long &count)
    {
        for (auto &cell_point : data)
        {
            count++;
            if (cell_point->can_reach(p, dist))
            {
                c = cell_point->cluster;
                return true;
            }
        }
        return false;
    }

    Cell() {}

    void add(SharedData<dims, T> p)
    {
        data.push_front(p);
    }

    Cluster<dims, T>* find_cluster(DataPoint<dims, T> *p, double dist, long &count)
    {
        Cluster<dims, T>* cluster;
        if (can_host(p, dist, cluster, count))
        {
            return cluster;
        }
        return nullptr;
    }
};

typedef std::vector<std::vector<size_t>> size_2D;

template<size_t dims> struct Key
{
    size_t operator()(const std::array<size_t, dims> &arr) const noexcept
    {
        size_t result = arr.size();
        for (size_t i = 0; i < arr.size(); i++)
        {
            result *= result * 63 + arr[i];
        }
        return result;
    }
};

SpaceType class Space
{

private:
    double max_dist;
    double cell_side_length;
    double actual_cell_side;
    long last_search_count;
    long search_count_target;
    long cluster_num;
    std::unordered_map<std::array<size_t, dims>, Cell<dims, T>, Key<dims>> cells;

    std::array<size_t, dims> to_lattice(std::array<double, dims> &arr)
    {
        std::array<size_t, dims> floored;
        std::transform(arr.begin(), arr.end(), floored.begin(),
                [this] (size_t coord)
                {
                    return floor(coord / this->actual_cell_side);
                }
        );
        return floored; 
    }

    void divide_grid()
    {
        double excess_ratio =
            (double)last_search_count / (double)search_count_target;
        while (excess_ratio > 1.0 && actual_cell_side >= 2*cell_side_length)
        {
            actual_cell_side /= 2.0;
            excess_ratio /= pow(2.0, dims);
        }
        std::cout << "New cell side length: " << actual_cell_side;
        std::cout << " - reassigning points..." << std::endl;
        std::unordered_map<std::array<size_t, dims>, Cell<dims, T>, Key<dims>> tmp;
        for (auto &cell : cells) 
        {
            for (auto &point : cell.second.data)
            {
                auto key = to_lattice(point->vec);
                if (!tmp.count(key))
                {
                    Cell<dims, T> c;
                    tmp[key] = c;
                }
                tmp[key].add(point);
            }   
        } 
        cells.clear();
        cells = tmp;
        std::cout << "Cell reassignment complete" << std::endl;
    }


public:
    std::unordered_map<long, Cluster<dims, T>*> clusters;

    Space(long search_count_target, 
          double cell_side_length,
          double max_dist,
          double upper_bound) 
    {
        actual_cell_side = upper_bound;
        cluster_num = 0;
        last_search_count = 0;
        this->max_dist = max_dist;
        this->search_count_target = search_count_target;
        this->cell_side_length = cell_side_length;
        std::array<size_t, dims> zero{};
        Cell<dims, T> cell;
        cells[zero] = cell;
    }

    ~Space()
    {
        for (auto it = clusters.begin(); it != clusters.end(); ++it)
        {
            delete it->second;
        }
    }

    std::vector<std::array<size_t, dims>> outer_product(size_2D opts, long n)
    {
        
        long num_neighbors = std::pow(n, opts.size());
        // Avoid a resize by computing the needed number of elements
        std::vector<std::array<size_t, dims>> neighbors(num_neighbors);
        std::vector<std::vector<size_t>::iterator> iterators;
        for (auto &a : opts) 
        {
            iterators.push_back(a.begin());
        }
        bool flag = false;
        size_t idx = 0;
        while (!flag)
        {
            std::array<size_t, dims> neighbor;
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

    std::unordered_set<Cluster<dims, T>*> find_clusters(DataPoint<dims, T> *p)
    {
        long span = max_dist/cell_side_length;
        // Round the coordinates to get a cell location in terms of
        // the cell side lengths
        std::array<size_t, dims> cell_loc = to_lattice(p->vec);
        // Get all the neighboring cells. This means for each dimension,
        // we get the previous and next coordinate for a total of 3,
        // and all the possible choices form a square in 2 dimensions,
        // a cube in 3 dimensions etc... 
        std::vector<std::vector<size_t>> opts(dims);
        for (size_t dim = 0; dim < dims; dim++)
        {
            for (long offset = -span; offset <= span; offset++)
            {
                long coord = cell_loc[dim] + offset;
                if (coord < 0)
                    coord = 0;
                opts[dim].push_back(coord);
            }
        }
        auto neighbors = outer_product(opts, 2 * span + 1);
        std::unordered_set<Cluster<dims, T>*> result;
        std::unordered_set<std::array<size_t, dims>, Key<dims>> searched;
        last_search_count = 0;
        for (auto &neighbor : neighbors)
        {
            if (cells.count(neighbor) && !searched.count(neighbor))
            {
                searched.insert(neighbor);
                Cluster<dims, T>* cluster;
                long count = 0;
                if ((cluster = cells[neighbor].find_cluster(p, max_dist, count)) 
                        != nullptr)
                {
                    result.insert(cluster);
                }
                last_search_count += count;
            }     
        }
        return result;   
    }

    /* Adds the data point to the space. Attempts to add the point to a
     * pre-existing cluster, and if it is nearby to multiple clusters, it
     * will merge those clusters. If it took more than the user defined
     * number of points to be searched, the cells will be divided into
     * smaller ones.  */
    long add(DataPoint<dims, T> &point)
    {
        auto p = make_shared<DataPoint<dims, T>>(point);
        auto key = to_lattice(point.vec);
        std::unordered_set<Cluster<dims, T>*> found = find_clusters(&point);
        Cluster<dims, T>* cluster = new Cluster<dims, T>(cluster_num++);
        switch (found.size())
        {
            case 0:
                cluster->add(p);
                p->cluster = cluster;
                clusters[cluster->id] = cluster;
                break;
            case 1:
                delete cluster;
                cluster_num--;
                (*found.begin())->add(p);
                p->cluster = (*found.begin());
                break;
            default:
                for (auto it = found.begin(); it != found.end(); it++)
                {
                    for (auto pt : (*it)->points)
                    {
                        cluster->add(pt);
                        pt->cluster = cluster;
                    }

                    clusters.erase(clusters.find((*it)->id));
                    delete *it;
                }
                p->cluster = cluster;
                cluster->add(p);
                clusters[cluster->id] = cluster; 
        }
        cells[key].add(p);
        if (last_search_count > search_count_target && actual_cell_side >= 2*cell_side_length)
        {
            std::cout << "Search count: " << last_search_count;
            std::cout << " - Target: " << search_count_target;
            std::cout << " - Side: " << actual_cell_side;
            std::cout << " - Dividing..." << std::endl;
            divide_grid();
        }
        return last_search_count;
    }

};

#endif
