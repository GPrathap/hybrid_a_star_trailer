#ifndef HYBRID_A_STAR_TRAILER_LIB
#define HYBRID_A_STAR_TRAILER_LIB

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
#include "lib_cpp_hybrid_a_star/grid_a_star.hpp"
#include "lib_cpp_hybrid_a_star/rs_paths.hpp"
#include "lib_cpp_hybrid_a_star/trailerlib.hpp"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

namespace planning
{
    using namespace math_utility; 

    typedef struct PlannerParams
    {
        static constexpr double XY_GRID_RESOLUTION = 2.0; //[m]
        static constexpr double YAW_GRID_RESOLUTION = math_utility::deg2rad(15.0); //[rad]
        static constexpr int N_PATHS_NEEDED = 1;
        static constexpr double GOAL_TYAW_TH = math_utility::deg2rad(5.0); //[rad]
        static constexpr double MOTION_RESOLUTION = 0.1; //[m] path interporate resolution
        static constexpr double N_STEER = 20.0; // number of steer command
        static constexpr double EXTEND_AREA= 5.0; //[m] map extend length
        static constexpr double SKIP_COLLISION_CHECK= 4; // skip number for collision check
        
        static constexpr double SB_COST = 100.0; // switch back penalty cost
        static constexpr double BACK_COST = 5.0; // backward penalty cost
        static constexpr double STEER_CHANGE_COST = 5.0; // steer angle change penalty cost
        static constexpr double STEER_COST = 1.0; // steer angle change penalty cost
        static constexpr double JACKKNIF_COST= 200.0; // Jackknif cost
        static constexpr double H_COST = 5.0; // Heuristic cost
        
        static constexpr double WB = VehicleParams::WB; //[m] Wheel base
        static constexpr double LT = VehicleParams::LT; //[m] length of trailer
        static constexpr double MAX_STEER = VehicleParams::MAX_STEER; //[rad] maximum steering angle
    };

    typedef struct HybridNode{
        int xind; //x index
        int yind; //y index
        int yawind; //yaw index
        bool direction;; // moving direction forword:true, backword:false
        Eigen::MatrixXd poses;
        // x::Array{Float64}; // x position [m]
        // y::Array{Float64}; // y position [m]
        // yaw::Array{Float64}; // yaw angle [rad]
        // yaw1::Array{Float64}; // trailer yaw angle [rad]
        // directions::Array{Bool}; // directions of each points forward: true, backward:false
        double steer; // steer input
        double cost; // cost
        int pind; // parent index
        HybridNode(){};
        HybridNode(int xind_, int yind_, int yawind_, bool direction_, 
                const Eigen::MatrixXd& poses_, double steer_, double cost_, int pind_)
                    : xind(xind_), yind(yind_), yawind(yawind_), 
                    direction(direction_), poses(poses_), steer(steer_), cost(cost_), pind(pind_) 
                {};

        friend std::ostream& operator<<(std::ostream& os, const HybridNode& node) {
            os << "HybridNode("
            << "xind: " << node.xind << ", "
            << "yind: " << node.yind << ", "
            << "yawind: " << node.yawind << ", "
            << "direction: " << (node.direction ? "forward" : "backward")
            << "steer: " << node.steer << ", "
            << "cost: " << node.cost << ", "
            << "pind: " << node.pind
            << " poses size : " << node.poses.rows()
            << ")";
            return os;
        }
    };

    template <typename T>
    void print_vec(const std::vector<T>& vec) {
        std::cout << "[ ";
        for (const auto& elem : vec) {
            std::cout << elem << " ";
        }
        std::cout << "]" << std::endl;
    }


    typedef struct Config{
        int minx;
        int miny;
        int minyaw;
        int minyawt;
        int maxx;
        int maxy;
        int maxyaw;
        int maxyawt;
        int xw;
        int yw;
        int yaww;
        int yawtw;
        double xyreso;
        double yawreso;
        Config(){};
        Config(int minx_, int miny_, int minyaw_, int minyawt_, int maxx_, int maxy_, int maxyaw_, int maxyawt_
                        , int xw_, int yw_, int yaww_, int yawtw_, double xyreso_, double yawreso_):
                        minx(minx_), miny(miny_), minyaw(minyaw_), minyawt(minyawt_), maxx(maxx_), maxy(maxy_), maxyaw(maxyaw_), maxyawt(maxyawt_), xw(xw_), yw(yw_), yaww(yaww_), yawtw(yawtw_),
                        xyreso(xyreso_), yawreso(yawreso_){};
    };

    typedef struct HybridPath{
        Eigen::MatrixXd poses;
        // x::Array{Float64} # x position [m]
        // y::Array{Float64} # y position [m]
        // yaw::Array{Float64} # yaw angle [rad]
        // yaw1::Array{Float64} # trailer angle [rad]
        // direction::Array{Bool} # direction forward: true, back false
        double cost; // cost
        HybridPath(const Eigen::MatrixXd& poses_init, double cost_init)
            : poses(poses_init), cost(cost_init) {}
        HybridPath() : poses(Eigen::MatrixXd()), cost(0.0) {}
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

    struct CompareCostNode {
        bool operator()(const CostNode& a, const CostNode& b) {
            return a.cost > b.cost;
        }
    };

    struct CompareNode {
        bool operator()(const HybridPath& a, const HybridPath& b) {
            return a.cost > b.cost;
        }
    };

    class TrailerHybridAStar{
        public:
            TrailerHybridAStar();
            void calc_config(Eigen::MatrixXd& obses, double xyreso, double yawreso);
            void calc_holonomic_with_obstacle_heuristic(HybridNode& gnode
                                    , Eigen::MatrixXd& obses, double xyreso);
            int calc_index(HybridNode& node);
            bool is_same_grid(HybridNode& node1, HybridNode& node2);
            HybridNode calc_next_node(HybridNode& current, int c_id, double u, double d);
            bool verify_index(HybridNode& node, Eigen::MatrixXd& obses
                            , double inityaw1, grid_search::KDTree& kdtree);
            void calc_motion_inputs(std::vector<double>& u, std::vector<double>& d);
            double calc_rs_path_cost(rs_paths::Path& rspath, Eigen::VectorXd& yaw1);
            bool analystic_expantion(HybridNode& node, HybridNode& ngoal
                        ,  Eigen::MatrixXd& obses, grid_search::KDTree& kdtree
                        , rs_paths::Path& selected_path);
            bool update_node_with_analystic_expantion(HybridNode& current
                , HybridNode& ngoal,  Eigen::MatrixXd& obses, grid_search::KDTree& kdtree
                , double gyaw1, HybridNode& updated_node);
            double calc_cost(HybridNode& n, HybridNode& ngoal);
            void get_final_path(std::unordered_map<int, HybridNode>& closed
                                                    , HybridNode& ngoal
                                                    , HybridNode& nstart, HybridPath& path);
            bool calc_hybrid_astar_path(Eigen::Vector4d s
                        , Eigen::Vector4d g, Eigen::MatrixXd& obses, HybridPath& path);
        private:
            Config config_;
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> h_dp;
            grid_search::GridAStar grid_a_star_;
            TrailerLib trailerlib_;
            rs_paths::RSPaths rs_path_;
    }; 
}


#endif 