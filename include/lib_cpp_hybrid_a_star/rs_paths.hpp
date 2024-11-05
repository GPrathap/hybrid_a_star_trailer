#ifndef RS_PATHS
#define RS_PATHS

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

namespace rs_paths
{
    using namespace math_utility;
    constexpr double STEP_SIZE = 0.1;

    typedef struct 
    {
        std::vector<double> lengths;
        std::vector<std::string> ctypes;
        double L{0.0};
        Eigen::MatrixXd poses;
        double cost;
        // Path() : lengths{} {}
    } Path;

    struct CompareNode {
        bool operator()(const Path& a, const Path& b) {
            return a.cost > b.cost;
        }
    };

    template <typename T>
    void printVector(const std::vector<T>& vec) {
        std::cout << "[ ";
        for (const auto& elem : vec) {
            std::cout << elem << " ";
        }
        std::cout << "]" << std::endl;
    }


    class RSPaths{
        public:
            RSPaths();
            void interpolate(int ind, double l, std::string m, double maxc, double ox, double oy, double oyaw, Eigen::MatrixXd& poses);
            void set_path(std::vector<Path>& paths, std::vector<double> lengths, std::vector<std::string> ctypes);
            bool SLS(Eigen::Vector3d& p, Eigen::Vector3d& p_up);
            void SCS(Eigen::Vector3d& p, std::vector<Path>& paths);
            bool LSL(const Eigen::Vector3d& p, Eigen::Vector3d& p_up);
            bool LSR(const Eigen::Vector3d& p, Eigen::Vector3d& p_up);
            bool LRL(const Eigen::Vector3d& p, Eigen::Vector3d& p_up);
            void get_label(Path& path, std::string& label);
            void calc_tauOmega(double u, double v, double xi, double eta, double phi, double& tau, double& omega);
            bool LRSR(const Eigen::Vector3d& p, Eigen::Vector3d& p_up);
            bool LRSL(const Eigen::Vector3d& p, Eigen::Vector3d& p_up);
            bool LRLRn(const Eigen::Vector3d& p, Eigen::Vector3d& p_up);
            bool LRLRp(const Eigen::Vector3d& p, Eigen::Vector3d& p_up);
            bool LRSLR(const Eigen::Vector3d& p, Eigen::Vector3d& p_up);
            void CCCC(const Eigen::Vector3d& p,  std::vector<Path>& paths);
            void CCSC(const Eigen::Vector3d& p, std::vector<Path>& paths);
            void CCSCC(const Eigen::Vector3d& p, std::vector<Path>& paths);
            void CCC(const Eigen::Vector3d& p, std::vector<Path>& paths);
            void CSC(const Eigen::Vector3d& p, std::vector<Path>& paths);
            void generate_path(Eigen::Vector3d q0, Eigen::Vector3d q1, double maxc, std::vector<Path>& paths);
            void generate_local_course(double L, std::vector<double> lengths
                            , std::vector<std::string> mode, double maxc, double step_size, Eigen::MatrixXd& poses);
            void calc_paths(Eigen::Vector3d s, Eigen::Vector3d g, double maxc,  std::vector<Path>& paths, double step_size=STEP_SIZE);
            Path calc_shortest_path(Eigen::Vector3d s, Eigen::Vector3d g, double maxc, double step_size=STEP_SIZE);
            double calc_shortest_path_length(Eigen::Vector3d s, Eigen::Vector3d g, double maxc, double step_size=STEP_SIZE);
            void calc_curvature(Eigen::MatrixXd& poses, std::vector<double>& c, std::vector<double>& ds);
            void check_path(Eigen::Vector3d s, Eigen::Vector3d g, double max_curvature);
            void test();
    };   
}


#endif 