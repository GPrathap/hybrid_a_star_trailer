#ifndef TRAILER_LIB
#define TRAILER_LIB

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
#include "matplotlibcpp.h"
#include <cstdlib>


namespace plt = matplotlibcpp;

namespace planning
{
    using namespace math_utility; 

    typedef struct VehicleParams
    {
        static constexpr double WB = 3.7;  //[m] wheel base: rear to front steer
        static constexpr double LT = 8.0; //[m] rear to trailer wheel
        static constexpr double W = 2.6; //[m] width of vehicle
        static constexpr double LF = 4.5; //[m] distance from rear to vehicle front end of vehicle
        static constexpr double LB = 1.0; //[m] distance from rear to vehicle back end of vehicle
        static constexpr double LTF = 1.0; //[m] distance from rear to vehicle front end of trailer
        static constexpr double LTB = 9.0; //[m] distance from rear to vehicle back end of trailer
        static constexpr double MAX_STEER = 0.6; //[rad] maximum steering angle 
        static constexpr double TR = 0.5; // Tyre radius [m] for plot
        static constexpr double TW = 1.0; // Tyre width [m] for plot
        static constexpr double WBUBBLE_DIST = 3.5; //distance from rear and the center of whole bubble
        static constexpr double WBUBBLE_R = 10.0; // whole bubble radius
        static constexpr double B = 4.45; // distance from rear to vehicle back end
        static constexpr double C = 11.54; // distance from rear to vehicle front end
        static constexpr double I = 8.55;// width of vehicle
        static const std::vector<double> VRX;
        static const std::vector<double> VRY;
    };

    class TrailerLib{
        public:
            TrailerLib();
            bool rect_check(const Eigen::Vector3d& current_pose, const Eigen::MatrixXd& obses
                                , const std::vector<double> vrx, const std::vector<double> vry);

            bool check_collision(Eigen::MatrixXd& poses, grid_search::KDTree& kd_tree
                        , const Eigen::MatrixXd& obses, double wbd, double wbr
                        , const std::vector<double> vrx, const std::vector<double> vry);

            bool calc_trailer_yaw_from_xyyaw(Eigen::MatrixXd& poses, const double init_tyaw
                                                        , Eigen::VectorXd& steps, Eigen::VectorXd& yaws);
            bool check_trailer_collision(const Eigen::MatrixXd& obses
                                    , Eigen::MatrixXd& poses, grid_search::KDTree& kd_tree);
            void plot_trailer(double x, double y, double yaw, double yaw1, double steer);
        private:
            
    }; 
}


#endif 