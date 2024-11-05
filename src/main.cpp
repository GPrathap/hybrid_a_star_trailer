#include <iostream>
#include <Eigen/Dense>

#include "lib_cpp_hybrid_a_star/rs_paths.hpp"
#include "lib_cpp_hybrid_a_star/grid_a_star.hpp"
#include "lib_cpp_hybrid_a_star/trailerlib.hpp"
#include "lib_cpp_hybrid_a_star/trailer_hybrid_a_star.hpp"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main(int arc, char *argv[]){

    // Eigen::Vector4d s(-10.0, 6.0, math_utility::deg2rad(0.0), math_utility::deg2rad(0.0));
    Eigen::Vector4d s(20.0, 6.0, math_utility::deg2rad(0.0), math_utility::deg2rad(0.0));
    // Eigen::Vector4d s(14.0, 10.0, math_utility::deg2rad(0.0), math_utility::deg2rad(0.0));
    Eigen::Vector4d g(0.0, 0.0, math_utility::deg2rad(90.0), math_utility::deg2rad(90.0));
    
    Eigen::MatrixXd obss(144, 2);
    int obs_index = 0;
    for(int i=-25; i<= 25; i++){
        obss.row(obs_index) << i*1.0, 15.0;
        obs_index++;
    }
    for(int i=-25; i<= -4; i++){
        obss.row(obs_index) << i*1.0, 4.0;
        obs_index++;
    }
    for(int i=-15; i<= 4; i++){
        obss.row(obs_index) << -4.0, i*1.0;
        obs_index++;
    }
    for(int i=-15; i<= 4; i++){
        obss.row(obs_index) << 4.0, i*1.0;
        obs_index++;
    }
    for(int i=4; i<= 25; i++){
        obss.row(obs_index) << i*1.0, 4.0;
        obs_index++;
    }
    for(int i=-4; i<= 4; i++){
        obss.row(obs_index) << i*1.0, -15.0;
        obs_index++;
    }

    planning::HybridPath path;
    planning::TrailerHybridAStar trailer_hybrid_astar;
    planning::TrailerLib trailer_lib;
    planning::PlannerParams planner_params;
    bool find_path = trailer_hybrid_astar.calc_hybrid_astar_path(s, g, obss, path);

    plt::figure();
    std::vector<double> ox, oy;
    for(int i=0; i< obss.rows(); i++){
        ox.push_back(obss.row(i)[0]);
        oy.push_back(obss.row(i)[1]);
    }

    std::vector<double> path_x, path_y;
    for(int i=0; i< path.poses.rows(); i++){
        path_x.push_back(path.poses.row(i)[0]);
        path_y.push_back(path.poses.row(i)[1]);
    }
   
    Eigen::VectorXd x = path.poses.col(0);
    Eigen::VectorXd y = path.poses.col(1);
    Eigen::VectorXd yaw = path.poses.col(2);
    Eigen::VectorXd yaw1 = path.poses.col(3);
    Eigen::VectorXd direction = path.poses.col(4);
    double steer = 0.0;
    for (size_t ii = 0; ii < x.size()-1; ++ii) {
        plt::clf();
        plt::plot(ox, oy, ".r");
        plt::plot(path_x, path_y, ".y");
        if (ii < x.size() - 1) {
            double k = (yaw[ii + 1] - yaw[ii]) / planner_params.MOTION_RESOLUTION;
            if (direction[ii]< 0.0) {
                k *= -1.0;
            }
            steer = std::atan2(planner_params.WB * k, 1.0);  // Equivalent of Julia's `Base.atan(WB*k, 1.0)`
        } else {
            steer = 0.0;
        }
        trailer_lib.plot_trailer(x[ii], y[ii], yaw[ii], yaw1[ii], steer);
        plt::grid(true);
        plt::axis("equal");
        plt::pause(0.0001);  // Small pause for animation
    }
    return 0;
}
