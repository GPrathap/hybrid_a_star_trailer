#ifndef GRID_SEARCH_A_STAR
#define GRID_SEARCH_A_STAR

#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <initializer_list>
#include <Eigen/Dense>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <thread>
#include <iostream>
#include <cstdlib> 
#include <execution>
#include "math_utility.hpp"
#include <nanoflann.hpp>
#include <queue>

namespace grid_search
{
    using namespace math_utility;

    typedef struct Node
    {
        Eigen::Vector2i pose;
        double cost;
        int pind;
        
        Node(const Eigen::Vector2i position, double cost_, int id){
            pose = position;
            cost = cost_;
            pind = id;
        }

        Node(){
        }

    };

    typedef struct CostNode
    {
        double cost;
        int id;
        CostNode(double cost_, int index){
            cost = cost_;
            id =  index;
        }
    };

    struct CompareNode {
        bool operator()(const CostNode& a, const CostNode& b) {
            return a.cost >= b.cost;
        }
    };
    

    struct PointCloud {
        Eigen::MatrixXd points;

        inline size_t kdtree_get_point_count() const { return points.rows(); }

        inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
            return points(idx, dim); // accessing the coordinate
        }

        template <class BBOX> bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
    };

    using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double
                                                                            , PointCloud>, PointCloud, 2>;


    class GridAStar{
        public:
            GridAStar();

            void calc_obstacle_map(Eigen::MatrixXd& obses, double& reso, double& vr);
            void calc_dist_policy(Eigen::Vector2d g
                            , Eigen::MatrixXd obses, double reso, double vr
                            , Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& pmap);
            int calc_index(const Node& node) const;
            double h(int x, int y) const;
            double calc_cost(const Node& n, const Node& ngoal) const;
            Eigen::MatrixXd get_motion_model() const;
            bool verify_node(Node node) const;
            void calc_policy_map(std::unordered_map<int, Node>& closed
                    ,  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& pmap) const;
            void get_final_path(std::unordered_map<int, Node>& closed, const Node ngoal, const Node nstart
                            , const double reso);
            void calc_astar_path(Eigen::Vector2d s, Eigen::Vector2d g
                        , Eigen::MatrixXd obses, double reso, double vr);
            Node search_min_cost_node(const std::unordered_map<int, Node*> open, const Node ngoal);
            inline std::vector<Eigen::Vector2d> getPath(){
                return final_path_;
            }

        private:
            Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> obs_map;
            int minx_, miny_, maxx_, maxy_;
            int xwidth_, ywidth_;
            std::vector<Eigen::Vector2d> final_path_;
    };   
}


#endif 