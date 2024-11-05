#include "lib_cpp_hybrid_a_star/trailerlib.hpp"

namespace planning
{
    const std::vector<double> VehicleParams::VRX = {VehicleParams::C, VehicleParams::C
                                    , -VehicleParams::B, -VehicleParams::B, VehicleParams::C };
    const std::vector<double> VehicleParams::VRY = {-VehicleParams::I/2.0, VehicleParams::I/2.0
                        , VehicleParams::I/2.0, -VehicleParams::I/2.0, -VehicleParams::I/2.0};

    TrailerLib::TrailerLib(){

    }

    bool TrailerLib::rect_check(const Eigen::Vector3d& current_pose, const Eigen::MatrixXd& obses
                                , const std::vector<double> vrx, const std::vector<double> vry){

        double c = cos(-current_pose.z());
        double s = sin(-current_pose.z());
        for(int j=0; j< obses.rows(); j++){
            Eigen::VectorXd obs = obses.row(j);
            double tx = obs.x() - current_pose.x();
            double ty = obs.y() - current_pose.y();
            double lx = (c*tx - s*ty);
            double ly = (s*tx + c*ty);
            double sumangle = 0.0;
            
            for(int i=0; i< vrx.size()-1; i++){
                double x1 = vrx[i] - lx;
                double y1 = vry[i] - ly;
                double x2 = vrx[i+1] - lx;
                double y2 = vry[i+1] - ly;
                double d1 = hypot(x1,y1);
                double d2 = hypot(x2,y2);
                double theta1 = atan2(y1,x1);
                double tty = (-sin(theta1)*x2 + cos(theta1)*y2);
                double tmp = (x1*x2+y1*y2)/(d1*d2);
                tmp = std::clamp(tmp, 0.0, 1.0);
                sumangle += (tty >= 0.0) ? acos(tmp) :  -acos(tmp);
            }
            if(sumangle >= M_PI){
                return false;
            }
        }
        return true;
    }

    bool TrailerLib::check_collision(Eigen::MatrixXd& poses, grid_search::KDTree& kd_tree
                        , const Eigen::MatrixXd& obses, double wbd, double wbr
                        , const std::vector<double> vrx, const std::vector<double> vry){
        
        double radius = wbr;
        for(int i=0; i<poses.rows(); i++){
            Eigen::VectorXd next_pose = poses.row(i);
            double cx = next_pose[0] + wbd*cos(next_pose[2]);
            double cy = next_pose[1] + wbd*sin(next_pose[2]);
            Eigen::Vector3d next_estimated_pose = next_pose.head(3);
            std::vector<nanoflann::ResultItem<size_t, double>> indices_dists;
            nanoflann::RadiusResultSet<double, size_t> result_set(radius, indices_dists);
            Eigen::Vector2d query_point(cx, cy);
            kd_tree.findNeighbors(result_set, query_point.data(), nanoflann::SearchParameters());
            if(indices_dists.size() < 1){
                continue;
            }
            Eigen::MatrixXd obs_subset(indices_dists.size(), 2);
            for (size_t i = 0; i < indices_dists.size(); ++i) {
                obs_subset.row(i) << obses(indices_dists[i].first, 0), obses(indices_dists[i].first, 1);
            }
            if (!rect_check(next_estimated_pose, obs_subset, vrx, vry)){
                return false;
            }
        }
        return true;
    }

    bool TrailerLib::calc_trailer_yaw_from_xyyaw(Eigen::MatrixXd& poses, const double init_tyaw
                                                        , Eigen::VectorXd& steps, Eigen::VectorXd& yaws){
        Eigen::VectorXd yaw = poses.col(2);
        yaws.resize(yaw.size());
        yaws.setZero();
        yaws[0] = init_tyaw;
        for(int i=1; i<yaws.size(); i++){
            double steps_d =  steps[i-1];
            double delta_theta = (steps_d/VehicleParams::LT)*sin(yaw[i-1] - yaws[i-1]);
            yaws[i] += yaws[i-1] + delta_theta;
        }
        return true;
    } 

    bool TrailerLib::check_trailer_collision(const Eigen::MatrixXd& obses
                                    , Eigen::MatrixXd& poses, grid_search::KDTree& kd_tree){
        std::vector<double> vrxt = {VehicleParams::LTF, VehicleParams::LTF, -VehicleParams::LTB
                                    , -VehicleParams::LTB, VehicleParams::LTF};
        std::vector<double> vryt = {-VehicleParams::W/2.0, VehicleParams::W/2.0
                                    , VehicleParams::W/2.0, -VehicleParams::W/2.0, -VehicleParams::W/2.0};

        double DT = (VehicleParams::LTF + VehicleParams::LTB)/2.0 - VehicleParams::LTB;
        double DTR = (VehicleParams::LTF + VehicleParams::LTB)/2.0 + 0.3;


        Eigen::MatrixXd back_path = poses;
        Eigen::MatrixXd front_path = poses;
        back_path.col(2) = poses.col(3);

        if (!check_collision(back_path, kd_tree, obses, DT, DTR, vrxt, vryt)){
            return false;
        }

        std::vector<double> vrxf = {VehicleParams::LF, VehicleParams::LF, -VehicleParams::LB, -VehicleParams::LB, VehicleParams::LF};
        std::vector<double> vryf = {-VehicleParams::W/2.0, VehicleParams::W/2.0
                                    , VehicleParams::W/2.0, -VehicleParams::W/2.0, -VehicleParams::W/2.0};
        double DF = (VehicleParams::LF + VehicleParams::LB)/2.0 - VehicleParams::LB;
        double DFR = (VehicleParams::LF + VehicleParams::LB)/2.0 + 0.3;

        if (!check_collision(front_path, kd_tree, obses, DF, DFR, vrxf, vryf)){
            return false;
        }
        return true;
    }

    void TrailerLib::plot_trailer(double x, double y, double yaw, double yaw1, double steer){

        std::string truckcolor = "-k";
        typedef VehicleParams VP;
        double LENGTH = VP::LB+VP::LF;
        double LENGTHt = VP::LTB + VP::LTF;
        Eigen::MatrixXd truckOutLine(2, 5);
        truckOutLine << -VP::LB, (LENGTH - VP::LB), (LENGTH - VP::LB), (-VP::LB), -VP::LB
                        , VP::W/2, VP::W/2, -VP::W/2, -VP::W/2, VP::W/2;
        Eigen::MatrixXd trailerOutLine(2, 5);
        trailerOutLine << -VP::LTB, (LENGTHt - VP::LTB), (LENGTHt - VP::LTB), (-VP::LTB), (-VP::LTB)
                            ,VP::W/2, VP::W/2, -VP::W/2, -VP::W/2, VP::W/2;

        Eigen::MatrixXd rr_wheel(2, 5);
        rr_wheel << VP::TR, -VP::TR, -VP::TR, VP::TR, VP::TR, -VP::W/12.0+VP::TW
                    ,  -VP::W/12.0+VP::TW, VP::W/12.0+VP::TW, VP::W/12.0+VP::TW, -VP::W/12.0+VP::TW;

        Eigen::MatrixXd rl_wheel(2, 5);
        rl_wheel << VP::TR, -VP::TR, -VP::TR, VP::TR, VP::TR, -VP::W/12.0-VP::TW
                    ,  -VP::W/12.0-VP::TW, VP::W/12.0-VP::TW, VP::W/12.0-VP::TW, -VP::W/12.0-VP::TW;
                    
        Eigen::MatrixXd fr_wheel(2, 5);
        fr_wheel << VP::TR, -VP::TR, -VP::TR, VP::TR, VP::TR,
                    -VP::W/12.0+VP::TW,  -VP::W/12.0+VP::TW, VP::W/12.0+VP::TW, VP::W/12.0+VP::TW, -VP::W/12.0+VP::TW;

        Eigen::MatrixXd fl_wheel(2, 5);
        fl_wheel <<  VP::TR, -VP::TR, -VP::TR, VP::TR, VP::TR,
                    -VP::W/12.0-VP::TW,  -VP::W/12.0-VP::TW, VP::W/12.0-VP::TW, VP::W/12.0-VP::TW, -VP::W/12.0-VP::TW; 

        Eigen::MatrixXd tr_wheel(2, 5);
        tr_wheel << VP::TR, -VP::TR, -VP::TR, VP::TR, VP::TR,
                    -VP::W/12.0+VP::TW,  -VP::W/12.0+VP::TW, VP::W/12.0+VP::TW, VP::W/12.0+VP::TW, -VP::W/12.0+VP::TW;
        
        Eigen::MatrixXd tl_wheel(2, 5);
        tl_wheel << VP::TR, -VP::TR, -VP::TR, VP::TR, VP::TR,
                    -VP::W/12.0-VP::TW,  -VP::W/12.0-VP::TW, VP::W/12.0-VP::TW, VP::W/12.0-VP::TW, -VP::W/12.0-VP::TW;
        
        Eigen::MatrixXd Rot1(2, 2);
        Rot1 << cos(yaw), sin(yaw), -sin(yaw), cos(yaw);

        Eigen::MatrixXd Rot2(2, 2);
        Rot2 << cos(steer), sin(steer), -sin(steer), cos(steer);

        Eigen::MatrixXd Rot3(2, 2);
        Rot3 << cos(yaw1), sin(yaw1), -sin(yaw1), cos(yaw1);

        fr_wheel = Rot2*fr_wheel;
        fl_wheel = Rot2*fl_wheel;
        fr_wheel.row(0) = fr_wheel.row(0).array() + VP::WB;
        fl_wheel.row(0) = fl_wheel.row(0).array() + VP::WB;
        fr_wheel = Rot1.transpose()*fr_wheel;
        fl_wheel = Rot1.transpose()*fl_wheel;
        tr_wheel.row(0) = tr_wheel.row(0).array() - VP::LT;
        tl_wheel.row(0) = tl_wheel.row(0).array() - VP::LT;
        tr_wheel = Rot3.transpose()*tr_wheel;
        tl_wheel = Rot3.transpose()*tl_wheel;

        truckOutLine = Rot1.transpose()*truckOutLine;
        trailerOutLine = Rot3.transpose()*trailerOutLine;

        rr_wheel = Rot1.transpose()*rr_wheel;
        rl_wheel = Rot1.transpose()*rl_wheel;

        truckOutLine.row(0) = truckOutLine.row(0).array() + x;
        truckOutLine.row(1) = truckOutLine.row(1).array() + y;
        trailerOutLine.row(0) = trailerOutLine.row(0).array() + x;
        trailerOutLine.row(1) = trailerOutLine.row(1).array() + y;

        fr_wheel.row(0) = fr_wheel.row(0).array() + x;
        fr_wheel.row(1) = fr_wheel.row(1).array() + y;
        rr_wheel.row(0) = rr_wheel.row(0).array() + x;
        rr_wheel.row(1) = rr_wheel.row(1).array() + y;
        fl_wheel.row(0) = fl_wheel.row(0).array() + x;
        fl_wheel.row(1) = fl_wheel.row(1).array() + y;
        rl_wheel.row(0) = rl_wheel.row(0).array() + x;
        rl_wheel.row(1) = rl_wheel.row(1).array() + y;

        tr_wheel.row(0) = tr_wheel.row(0).array() + x;
        tr_wheel.row(1) = tr_wheel.row(1).array() + y;
        tl_wheel.row(0) = tl_wheel.row(0).array() + x;
        tl_wheel.row(1) = tl_wheel.row(1).array() + y;


        auto plotOutline = [&](const Eigen::MatrixXd& outline) {
            std::vector<double> vec_x(outline.cols()), vec_y(outline.cols());
            for (int i = 0; i < outline.cols(); ++i) {
                vec_x[i] = outline(0, i);
                vec_y[i] = outline(1, i);
            }
            plt::plot(vec_x, vec_y, truckcolor);
        };

        // plt::figure();
        plotOutline(truckOutLine);
        plotOutline(trailerOutLine);
        plotOutline(fr_wheel);
        plotOutline(rr_wheel);
        plotOutline(rl_wheel);
        plotOutline(fl_wheel);
        plotOutline(tr_wheel);
        plotOutline(tl_wheel);
        plt::plot({x}, {y}, "*");
        // plt::axis("equal");
        // plt::show();

    }
}

// int main(int argc, char *argv[]) {

//     double x = 0.0;
//     double y = 0.0;
//     double yaw0 = math_utility::deg2rad(10.0);
//     double yaw1 = math_utility::deg2rad(-10.0);
//     planning::TrailerLib trailer_lib;
//     trailer_lib.plot_trailer(x, y, yaw0, yaw1, 0.0);
// }

