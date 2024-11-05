#include "lib_cpp_hybrid_a_star/rs_paths.hpp"


namespace rs_paths
{

    RSPaths::RSPaths(){

    }

    void RSPaths::interpolate(int ind, double l, std::string m, double maxc, double ox, double oy, double oyaw, Eigen::MatrixXd& poses){

            if(m == "S"){
                poses.row(ind).head(3) << ox + l / maxc * cos(oyaw), oy + l / maxc * sin(oyaw), oyaw;
            }else{
                double ldx = sin(l) / maxc;
                double ldy;
                if( m == "L" ){  // left turn
                    ldy = (1.0 - cos(l)) / maxc;
                } else if ( m == "R" ) {  // right turn
                    ldy = (1.0 - cos(l)) / -maxc ;
                }
                double gdx = cos(-oyaw) * ldx + sin(-oyaw) * ldy;
                double gdy = -sin(-oyaw) * ldx + cos(-oyaw) * ldy;
                poses.row(ind).head(2) << ox + gdx, oy + gdy;
            }
            if (m == "L"){  // left turn
                poses.row(ind)[2] = oyaw + l;
            } else if(m == "R"){  // right turn
                poses.row(ind)[2] = oyaw - l;
            }  
            if(l > 0.0){
                poses.row(ind)[3] = 1.0;
            }else{
                poses.row(ind)[3] = -1.0;
            }
    }

    void RSPaths::set_path(std::vector<Path>& paths, std::vector<double> lengths, std::vector<std::string> ctypes){
        Path path;
        path.ctypes = ctypes;
        path.lengths = lengths;
        // check same path exist
        std::for_each(paths.begin(), paths.end(), [&](Path &tpath) {
            bool typeissame = (tpath.ctypes == path.ctypes);
            if(typeissame){
                if (tpath.lengths.size() != path.lengths.size()) {
                    std::cerr << "Vectors must have the same size!" << std::endl;
                    return;
                }else{
                    std::vector<double> result_diff_path;
                    std::transform(tpath.lengths.begin(), tpath.lengths.end()
                            , path.lengths.begin(), std::back_inserter(result_diff_path), std::minus<>());
                    double sum = std::accumulate(result_diff_path.begin(), result_diff_path.end(), 0.0);
                    if(sum <= 0.01){
                        return;
                    }
                }
            }
        });
        path.L = std::transform_reduce(path.lengths.begin(), path.lengths.end(), 0.0, std::plus<double>(), 
                    [](double x) { return std::abs(x); });

        paths.push_back(path);
    }
    
    bool RSPaths::SLS(Eigen::Vector3d& p, Eigen::Vector3d& p_up){
        double x = p.x();
        double y = p.y();
        double phi = p.z();
        phi = mod2pi(phi);
        if (y > 0.0 && phi > 0.0 && phi < M_PI*0.99){
            double xd = - y/tan(phi) + x;
            double t =  xd - tan(phi/2.0);
            double u = phi;
            double v = sqrt(pow((x-xd),2)+pow(y, 2))-tan(phi/2.0);
            p_up << t, u, v; 
            return true;
        }else if(y < 0.0 && phi > 0.0 && phi < M_PI*0.99){
            double xd = - y/tan(phi) + x;
            double t =  xd - tan(phi/2.0);
            double u = phi;
            double v = -sqrt(pow((x-xd),2)+pow(y, 2))-tan(phi/2.0);
            p_up << t, u, v; 
            return true;
        }
        p_up << 0.0, 0.0, 0.0; 
        return false;
    }

    void RSPaths::SCS(Eigen::Vector3d& p, std::vector<Path>& paths){
        Eigen::Vector3d p_up;
        bool flag = SLS(p, p_up);
        double t = p_up[0];
        double u = p_up[1];
        double v = p_up[2];
        if (flag){
            std::vector<double> lengths = {t, u, v};
            std::vector<std::string> ctypes = {"S","L","S"};
            set_path(paths, lengths, ctypes);
        }
        double x = p.x();
        double y = p.y();
        double phi = p.z();
        Eigen::Vector3d p_est( x, -y, -phi);
        flag = SLS(p_est, p_up);
        if (flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {t, u, v};
            std::vector<std::string> ctypes = {"S","R","S"};
            set_path(paths, lengths, ctypes);
        }
    }

    bool RSPaths::LSL(const Eigen::Vector3d& p, Eigen::Vector3d& p_up){
        double x = p.x();
        double y = p.y();
        double phi = p.z();
        double u;
        double t;
        double v;
        polar(x - sin(phi), y - 1.0 + cos(phi), u, t);
        if(t >= 0.0){
            v = mod2pi(phi - t);
            if (v >= 0.0){
                p_up << t, u, v;
                return true;
            }
        }
        p_up << 0.0, 0.0, 0.0; 
        return false;
    }

    bool RSPaths::LSR(const Eigen::Vector3d& p, Eigen::Vector3d& p_up){
        double x = p.x();
        double y = p.y();
        double phi = p.z();
        double u1;
        double t1;
        double v;
        polar(x + sin(phi), y - 1.0 - cos(phi), u1, t1);
        u1 = std::pow(u1, 2);
        if(u1 >= 4.0){
            double u = sqrt(u1 - 4.0);
            double theta = atan2(2.0, u);
            double t = mod2pi(t1 + theta);
            double v = mod2pi(t - phi);
            if (t >= 0.0 && v >= 0.0){
                p_up << t, u, v;
                return true;
            }
        }
        p_up << 0.0, 0.0, 0.0; 
        return false;
    }

    bool RSPaths::LRL(const Eigen::Vector3d& p, Eigen::Vector3d& p_up){
        double x = p.x();
        double y = p.y();
        double phi = p.z();
        double u1;
        double t1;
        double v;
        polar(x - sin(phi), y - 1.0 + cos(phi), u1, t1);
        if(u1 <= 4.0){
            double u = -2.0*asin(0.25 * u1);
            double t = mod2pi(t1 + 0.5 * u + M_PI);
            double v = mod2pi(phi - t + u);
            if (t >= 0.0 && u <= 0.0){
                p_up << t, u, v;
                return true;
            }
        }
        p_up << 0.0, 0.0, 0.0; 
        return false;
    }

    void RSPaths::get_label(Path& path, std::string& label){
        label = "";
        for(int i=0; i< path.ctypes.size(); i++){
            label += path.ctypes[i];
            if(path.lengths[i] > 0.0){
                label += "+";
            }else{
                label += "-";
            }
        }
    }

    void RSPaths::calc_tauOmega(double u, double v, double xi, double eta, double phi, double& tau, double& omega){
        double delta = mod2pi(u-v);
        double A = sin(u) - sin(delta);
        double B = cos(u) - cos(delta) - 1.0;
        double t1 = atan2(eta*A - xi*B, xi*A + eta*B);
        double t2 = 2.0 * (cos(delta) - cos(v) - cos(u)) + 3.0;
        tau = (t2 < 0) ? mod2pi(t1 + M_PI): mod2pi(t1);
        omega = mod2pi(tau - u + v - phi);
        return;
    }

    bool RSPaths::LRSR(const Eigen::Vector3d& p, Eigen::Vector3d& p_up){
        double x = p.x();
        double y = p.y();
        double phi = p.z();
        double xi = x + sin(phi);
        double eta = y - 1.0 - cos(phi);
        double rho, theta;
        polar(-eta, xi, rho, theta);
        if (rho >= 2.0){
            double t = theta;
            double u = 2.0 - rho;
            double v = mod2pi(t + 0.5*M_PI - phi);
            if (t >= 0.0 && u <= 0.0 && v <=0.0){
                p_up << t, u, v;
                return true;
            }        
        }
        p_up << 0.0, 0.0, 0.0; 
        return false;
    }

    bool RSPaths::LRSL(const Eigen::Vector3d& p, Eigen::Vector3d& p_up){
        double x = p.x();
        double y = p.y();
        double phi = p.z();
        double xi = x - sin(phi);
        double eta = y - 1.0 + cos(phi);
        double rho, theta;
        polar(xi, eta, rho, theta);
        if (rho >= 2.0){
            double r = sqrt(rho*rho - 4.0);
            double u = 2.0 - r;
            double t = mod2pi(theta + atan2(r, -2.0));
            double v = mod2pi(phi - 0.5*M_PI - t);
            if (t >= 0.0 && u<=0.0 && v<=0.0){
                p_up << t, u, v;
                return true;
            }        
        }
        p_up << 0.0, 0.0, 0.0; 
        return false;
    }

    bool RSPaths::LRLRn(const Eigen::Vector3d& p, Eigen::Vector3d& p_up){
        double x = p.x();
        double y = p.y();
        double phi = p.z();
        double xi = x + sin(phi);
        double eta = y - 1.0 - cos(phi);
        double rho = 0.25 * (2.0 + sqrt(xi*xi + eta*eta));
        if(rho <= 1.0){
            double u = acos(rho);
            double t, v;
            calc_tauOmega(u, -u, xi, eta, phi, t, v);
            if (t >= 0.0 && v <= 0.0){
                p_up << t, u, v;
                return true;
            }
        }
        p_up << 0.0, 0.0, 0.0; 
        return false;
    }


    bool RSPaths::LRLRp(const Eigen::Vector3d& p, Eigen::Vector3d& p_up){
        double x = p.x();
        double y = p.y();
        double phi = p.z();
        double xi = x + sin(phi);
        double eta = y - 1.0 - cos(phi);
        double rho = (20.0 - xi*xi - eta*eta) / 16.0;
        if(rho>=0.0 && rho<=1.0){
            double u = -acos(rho);
            if (u >= -0.5 * M_PI){
                double t, v;
                calc_tauOmega(u, u, xi, eta, phi, t, v);
                if (t >= 0.0 && v >= 0.0){
                    p_up << t, u, v;
                    return true;
                }
            }
            
        }
        p_up << 0.0, 0.0, 0.0; 
        return false;
    }

    bool RSPaths::LRSLR(const Eigen::Vector3d& p, Eigen::Vector3d& p_up){
        double x = p.x();
        double y = p.y();
        double phi = p.z();
        double xi = x + sin(phi);
        double eta = y - 1.0 - cos(phi);
        double rho, theta;
        polar(xi, eta, rho, theta);
        if(rho >= 2.0){
            double u = 4.0 - sqrt(rho*rho - 4.0);
            if (u <= 0.0){
                double t, v;
                t = mod2pi(atan2((4.0-u)*xi -2.0*eta, -2.0*xi + (u-4.0)*eta));
                v = mod2pi(t - phi);
                if ( t >= 0.0 && v >=0.0){
                    p_up << t, u, v;
                    return true;
                }
            }
            
        }
        p_up << 0.0, 0.0, 0.0; 
        return false;
    }


    void RSPaths::CCCC(const Eigen::Vector3d& p,  std::vector<Path>& paths){

        double x = p.x();
        double y = p.y();
        double phi = p.z();
        Eigen::Vector3d p_up;
        bool flag = LRLRn(p, p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {t, u, -u, v};
            std::vector<std::string> ctypes = {"L","R","L","R"};
            set_path(paths, lengths, ctypes); 
        }

        
        flag = LRLRn(Eigen::Vector3d(-x, y, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-t, -u, u, -v};
            std::vector<std::string> ctypes = {"L","R","L","R"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRLRn(Eigen::Vector3d(x, -y, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {t, u, -u, v};
            std::vector<std::string> ctypes = {"R","L","R","L"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRLRn(Eigen::Vector3d(-x, -y, phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-t, -u, u, -v};
            std::vector<std::string> ctypes = {"R","L","R","L"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRLRp(Eigen::Vector3d(x, y, phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {t, u, u, v};
            std::vector<std::string> ctypes = {"L","R","L","R"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRLRp(Eigen::Vector3d(-x, y, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-t, -u, -u, -v};
            std::vector<std::string> ctypes = {"L","R","L","R"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRLRp(Eigen::Vector3d(x, -y, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {t, u, u, v};
            std::vector<std::string> ctypes = {"R","L","R","L"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRLRp(Eigen::Vector3d(-x, -y, phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-t, -u, -u, -v};
            std::vector<std::string> ctypes = {"R","L","R","L"};
            set_path(paths, lengths, ctypes);
        }
    }

    void RSPaths::CCSC(const Eigen::Vector3d& p, std::vector<Path>& paths){

        double x = p.x();
        double y = p.y();
        double phi = p.z();
        Eigen::Vector3d p_up;
        bool flag = LRSL(p, p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {t, -0.5*M_PI, u, v};
            std::vector<std::string> ctypes = {"L","R","S","L"};
            set_path(paths, lengths, ctypes); 
        }

        
        flag = LRSL(Eigen::Vector3d(-x, y, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-t, 0.5*M_PI, -u, -v};
            std::vector<std::string> ctypes = {"L","R","S","L"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRSL(Eigen::Vector3d(x, -y, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {t, -0.5*M_PI, u, v};
            std::vector<std::string> ctypes = {"R","L","S","R"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRSL(Eigen::Vector3d(-x, -y, phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-t, 0.5*M_PI, -u, -v};
            std::vector<std::string> ctypes = {"R","L","S","R"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRSR(Eigen::Vector3d(x, y, phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {t, -0.5*M_PI, u, v};
            std::vector<std::string> ctypes = {"L","R","S","R"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRSR(Eigen::Vector3d(-x, y, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-t, 0.5*M_PI, -u, -v};
            std::vector<std::string> ctypes = {"L","R","S","R"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRSR(Eigen::Vector3d(x, -y, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {t, -0.5*M_PI, u, v};
            std::vector<std::string> ctypes = {"R","L","S","L"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRSR(Eigen::Vector3d(-x, -y, phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-t, 0.5*M_PI, -u, -v};
            std::vector<std::string> ctypes = {"R","L","S","L"};
            set_path(paths, lengths, ctypes);
        }

        // backwards
        double xb = x*cos(phi) + y*sin(phi);
        double yb = x*sin(phi) - y*cos(phi);
        flag = LRSL(Eigen::Vector3d(xb, yb, phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {v, u, -0.5*M_PI, t};
            std::vector<std::string> ctypes = {"L","S","R","L"};
            set_path(paths, lengths, ctypes);
        }


        flag = LRSL(Eigen::Vector3d(-xb, yb, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-v, -u, 0.5*M_PI, -t};
            std::vector<std::string> ctypes = {"L","S","R","L"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRSL(Eigen::Vector3d(xb, -yb, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {v, u, -0.5*M_PI, t};
            std::vector<std::string> ctypes = {"R","S","L","R"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRSL(Eigen::Vector3d(-xb, -yb, phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-v, -u, 0.5*M_PI, -t};
            std::vector<std::string> ctypes = {"R","S","L","R"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRSR(Eigen::Vector3d(xb, yb, phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {v, u, -0.5*M_PI, t};
            std::vector<std::string> ctypes = {"R","S","R","L"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRSR(Eigen::Vector3d(-xb, yb, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-v, -u, 0.5*M_PI, -t};
            std::vector<std::string> ctypes = {"R","S","R","L"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRSR(Eigen::Vector3d(xb, -yb, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {v, u, -0.5*M_PI, t};
            std::vector<std::string> ctypes = {"L","S","L","R"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRSR(Eigen::Vector3d(-xb, -yb, phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-v, -u, 0.5*M_PI, -t};
            std::vector<std::string> ctypes = {"L","S","L","R"};
            set_path(paths, lengths, ctypes);
        }
    }

    void RSPaths::CCSCC(const Eigen::Vector3d& p, std::vector<Path>& paths){

        double x = p.x();
        double y = p.y();
        double phi = p.z();
        Eigen::Vector3d p_up;
        bool flag = LRSLR(p, p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {t, -0.5*M_PI, u, -0.5*M_PI, v};
            std::vector<std::string> ctypes = {"L","R","S","L","R"};
            set_path(paths, lengths, ctypes); 
        }

        flag = LRSLR(Eigen::Vector3d(-x, y, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-t, 0.5*M_PI, -u, 0.5*M_PI, -v};
            std::vector<std::string> ctypes = {"L","R","S","L","R"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRSLR(Eigen::Vector3d(x, -y, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {t, -0.5*M_PI, u, -0.5*M_PI, v};
            std::vector<std::string> ctypes = {"R","L","S","R","L"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRSLR(Eigen::Vector3d(-x, -y, phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-t, 0.5*M_PI, -u, 0.5*M_PI, -v};
            std::vector<std::string> ctypes = {"R","L","S","R","L"};
            set_path(paths, lengths, ctypes);
        }
    }

    void RSPaths::CCC(const Eigen::Vector3d& p, std::vector<Path>& paths){

        double x = p.x();
        double y = p.y();
        double phi = p.z();
        Eigen::Vector3d p_up;
        bool flag = LRL(p, p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {t, u, v};
            std::vector<std::string> ctypes = {"L","R","L"};
            set_path(paths, lengths, ctypes); 
        }

        
        flag = LRL(Eigen::Vector3d(-x, y, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-t, -u, -v};
            std::vector<std::string> ctypes = {"L","R","L"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRL(Eigen::Vector3d(x, -y, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {t, u, v};
            std::vector<std::string> ctypes = {"R","L","R"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRL(Eigen::Vector3d(-x, -y, phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-t, -u, -v};
            std::vector<std::string> ctypes = {"R","L","R"};
            set_path(paths, lengths, ctypes);
        }

        //  backwards
        double xb = x*cos(phi) + y*sin(phi);
        double yb = x*sin(phi) - y*cos(phi);

        flag = LRL(Eigen::Vector3d(xb, yb, phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {v, u, t};
            std::vector<std::string> ctypes = {"L","R","L"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRL(Eigen::Vector3d(-xb, yb, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-v, -u, -t};
            std::vector<std::string> ctypes = {"L","R","L"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRL(Eigen::Vector3d(xb, -yb, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {v, u, t};
            std::vector<std::string> ctypes = {"R","L","R"};
            set_path(paths, lengths, ctypes);
        }

        flag = LRL(Eigen::Vector3d(-xb, -yb, phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-v, -u, -t};
            std::vector<std::string> ctypes = {"R","L","R"};
            set_path(paths, lengths, ctypes);
        }
    }

    void RSPaths::CSC(const Eigen::Vector3d& p, std::vector<Path>& paths){
        double x = p.x();
        double y = p.y();
        double phi = p.z();
        Eigen::Vector3d p_up;
        bool flag = LSL(p, p_up);
        
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {t, u, v};
            std::vector<std::string> ctypes = {"L","S","L"};
            set_path(paths, lengths, ctypes); 
        }

        
        flag = LSL(Eigen::Vector3d(-x, y, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-t, -u, -v};
            std::vector<std::string> ctypes = {"L","S","L"};
            set_path(paths, lengths, ctypes);
        }

        flag = LSL(Eigen::Vector3d(x, -y, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {t, u, v};
            std::vector<std::string> ctypes = {"R","S","R"};
            set_path(paths, lengths, ctypes);
        }

        flag = LSL(Eigen::Vector3d(-x, -y, phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-t, -u, -v};
            std::vector<std::string> ctypes = {"R","S","R"};
            set_path(paths, lengths, ctypes);
        }
        flag = LSR(Eigen::Vector3d(x, y, phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {t, u, v};
            std::vector<std::string> ctypes = {"L","S","R"};
            set_path(paths, lengths, ctypes);
        }

        flag = LSR(Eigen::Vector3d(-x, y, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-t, -u, -v};
            std::vector<std::string> ctypes = {"L","S","R"};
            set_path(paths, lengths, ctypes);
        }

        flag = LSR(Eigen::Vector3d(x, -y, -phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {t, u, v};
            std::vector<std::string> ctypes = {"R","S","L"};
            set_path(paths, lengths, ctypes);
        }

        flag = LSR(Eigen::Vector3d(-x, -y, phi), p_up);
        if(flag){
            double t = p_up[0];
            double u = p_up[1];
            double v = p_up[2];
            std::vector<double> lengths = {-t, -u, -v};
            std::vector<std::string> ctypes = {"R","S","L"};
            set_path(paths, lengths, ctypes);
        }
    }

    void RSPaths::generate_path(Eigen::Vector3d q0, Eigen::Vector3d q1, double maxc, std::vector<Path>& paths){
        double dx = q1[0] - q0[0];
        double dy = q1[1] - q0[1];
        double dth = q1[2] - q0[2];
        double c = cos(q0[2]);
        double s = sin(q0[2]);
        double x = (c*dx + s*dy)*maxc;
        double y = (-s*dx + c*dy)*maxc;
        Eigen::Vector3d p_est(x, y, dth);
        SCS(p_est, paths);
        CSC(p_est, paths);
        CCC(p_est, paths);
        CCCC(p_est, paths);
        CCSC(p_est, paths);
        CCSCC(p_est, paths);
        
    }

    void RSPaths::generate_local_course(double L, std::vector<double> lengths
                    , std::vector<std::string> mode, double maxc, double step_size, 
                    Eigen::MatrixXd& poses){
                
                int npoint = std::abs(trunc(L/step_size)) + static_cast<int>(lengths.size()) + 3;
                poses = Eigen::MatrixXd::Zero(npoint, 4);
                int ind = 1;
                if (lengths[0] > 0.0){
                    poses.row(0)[3] = 1.0;
                } else{
                    poses.row(0)[3] = -1.0;
                } 
                double d = (lengths[0] > 0.0) ? step_size : -step_size;
                double pd = d;
                double ll = 0.0;
                for (size_t i = 0; i < mode.size(); i++)
                {
                    double l = lengths[i];
                    std::string m = mode[i];
                    d = (l > 0.0) ? step_size : -step_size;
                    // set prigin state
                    double ox = poses.row(ind)[0];
                    double oy = poses.row(ind)[1];
                    double oyaw = poses.row(ind)[2];
                    ind -= 1;
                    if ( i >= 1 && (lengths[i-1]*lengths[i])>0.0 ){
                        pd = -d - ll;
                    } else{
                        pd = d - ll;
                    }
                    while (std::abs(pd) <= std::abs(l) && ind < poses.rows()-1){
                        ind += 1;
                        interpolate(ind, pd, m, maxc, ox, oy, oyaw, poses);
                        pd += d;
                    }
                    ll = l - pd - d; // calc remain length
                    ind += 1;
                    ind = (ind >= poses.rows()) ? ind-1 : ind;
                    interpolate(ind, l, m, maxc, ox, oy, oyaw, poses);
                }
                // remove unused data
                int last_index = poses.rows()-1;
                while ( poses.row(last_index)[0] == 0.0 && last_index > 0){
                    last_index--;
                }
                poses.conservativeResize(last_index+1, poses.cols());
    }

    void RSPaths::calc_paths(Eigen::Vector3d s, Eigen::Vector3d g, double maxc,  std::vector<Path>& paths, double step_size){

        Eigen::Vector3d q0 = s;
        Eigen::Vector3d q1 = g;
        generate_path(q0, q1, maxc, paths);
        std::for_each(paths.begin(), paths.end(), [&](Path &path) {
            Eigen::MatrixXd poses;
            // printVector(path.lengths);
            generate_local_course(path.L, path.lengths, path.ctypes, maxc, step_size*maxc, poses);
             // convert global coordinate
            for(int i=0; i<poses.rows(); i++){
                double ix = poses.row(i)[0];
                double iy = poses.row(i)[1];
                double iyaw = poses.row(i)[2];
                double direction = poses.row(i)[3];
                double x = cos(-q0[2]) * ix + sin(-q0[2]) * iy + q0[0];
                double y = -sin(-q0[2]) * ix + cos(-q0[2]) * iy + q0[1];
                double yaw = pi_to_pi(iyaw + q0[2]);
                poses.row(i) << x, y, yaw, direction;
            }
            for (int i=0; i< path.lengths.size(); i++){
                path.lengths[i] = path.lengths[i]/maxc;
            }
            path.L = path.L/maxc;
            path.poses = poses;
        });
    }

    Path RSPaths::calc_shortest_path(Eigen::Vector3d s, Eigen::Vector3d g, double maxc, double step_size){
        std::vector<Path> paths;
        calc_paths(s, g, maxc, paths, step_size);
        double minL = std::numeric_limits<double>::max();
        int best_path_index = 0;
        int index = 0;
        std::for_each(paths.begin(), paths.end(), [&](Path &path) {
            if (path.L <= minL){
                minL = path.L;
                best_path_index = index;
            }
            ++index;
        });
        Path tmp;
        return (paths.size() > best_path_index ) ? paths[best_path_index] : tmp;
    }

    double RSPaths::calc_shortest_path_length(Eigen::Vector3d s, Eigen::Vector3d g, double maxc, double step_size){

        Eigen::Vector3d q0 = s;
        Eigen::Vector3d q1 = g;
        std::vector<Path> paths;
        generate_path(q0, q1, maxc, paths);
        double minL = std::numeric_limits<double>::max();
        std::for_each(paths.begin(), paths.end(), [&](Path &path) {
            double L = path.L/maxc;
            if (L <= minL){
                minL = path.L;
            }
        });
        return minL;
    }

    void RSPaths::calc_curvature(Eigen::MatrixXd& poses, std::vector<double>& c, std::vector<double>& ds){
        
        for(int i=1; i< poses.rows()-1; i++){
            double dxn = poses.row(i)[0]-poses.row(i-1)[0];
            double dxp = poses.row(i+1)[0]-poses.row(i)[0];
            double dyn = poses.row(i)[1]-poses.row(i-1)[1];
            double dyp = poses.row(i+1)[1]-poses.row(i)[1];
            double dn = sqrt(pow(dxn,2.0)+pow(dyn,2.0));
            double dp = sqrt(pow(dxp,2.0)+pow(dyp,2.0));
            double dx = 1.0/(dn+dp)*(dp/dn*dxn+dn/dp*dxp);
            double ddx = 2.0/(dn+dp)*(dxp/dp-dxn/dn);
            double dy = 1.0/(dn+dp)*(dp/dn*dyn+dn/dp*dyp);
            double ddy = 2.0/(dn+dp)*(dyp/dp-dyn/dn);
            double curvature = (ddy*dx-ddx*dy)/(pow(dx,2)+pow(dy,2.0));
            double d = (dn+dp)/2.0;
            if(std::isnan(curvature)){
                curvature = 0.0;
            }
            int direction = poses.row(i)[3];
            if(direction <= 0.0){
                curvature = -curvature;
            }
            if(c.size() == 0){
                ds.push_back(d);
                c.push_back(curvature);
            }
        }
        if(ds.size()>0){
            ds.push_back(ds.back());
            c.push_back(c.back());  
        }
    }

    void RSPaths::check_path(Eigen::Vector3d s, Eigen::Vector3d g, double max_curvature){
        std::vector<Path> paths;
        calc_paths(s, g, max_curvature, paths);
        std::for_each(paths.begin(), paths.end(), [&](Path &path) {
            // Matrix to store the differences (it will have one less row)
            Eigen::MatrixXd diffs(path.poses.rows() - 1, path.poses.cols());
            // Calculate column-wise differences
            for (int col = 0; col < path.poses.cols(); ++col) {
                diffs.col(col) = path.poses.col(col).tail(path.poses.rows() - 1) - path.poses.col(col).head(path.poses.rows() - 1);
            }
            // create a vector to store the result (size 20)
            Eigen::VectorXd result(diffs.rows());
            result = (diffs.col(0).array().square() + diffs.col(1).array().square()).sqrt(); 
            for(int i=0; i<result.size(); i++){
                if(std::abs(result[i] - STEP_SIZE) <= 0.001){
                    std::cerr<< "trajectory is not feasible " << std::endl;
                    return;
                }
            }
        });
    }


    void RSPaths::test(){
        Eigen::Vector3d s(0.0, 0.0, deg2rad(10.0));
        Eigen::Vector3d g(7.0, -8.0, deg2rad(50.0));
        double max_curvature = 2.0;
        check_path(s, g, max_curvature);

        s<< (0.0, 10.0, deg2rad(-10.0));
        g<< (-7.0, -8.0, deg2rad(-50.0));
        max_curvature = 2.0;
        check_path(s, g, max_curvature);

        s<< (0.0, 10.0, deg2rad(-10.0));
        g<< (-7.0, -8.0, deg2rad(150.0));
        max_curvature = 1.0;
        check_path(s, g, max_curvature);

        s<< (0.0, 10.0, deg2rad(-10.0));
        g<< (7.0, 8.0, deg2rad(150.0));
        max_curvature = 2.0;
        check_path(s, g, max_curvature);

        s<< (-40.0, 549.0, 2.443);
        g<< (36.0, 446.0, 0.698);
        max_curvature = 0.0589;
        check_path(s, g, max_curvature);

        for(int i=0; i< 100; i++){
            s<< (rand()*100.0 - 50.0, rand()*100.0 - 50.0, deg2rad(rand()*360.0 - 180.0));
            g<< (rand()*100.0 - 50.0, rand()*100.0 - 50.0, deg2rad(rand()*360.0 - 180.0));
            max_curvature = rand()/10.0;
            std::cout<< i << "," << s.transpose() << ", " << g.transpose() << std::endl;
            check_path(s, g, max_curvature);
        }
    }

}


// int main(int argc, char *argv[]) {

//     rs_paths::RSPaths rs_path;
//     Eigen::Vector3d s(3.0, 10.0, math_utility::deg2rad(40.0));
//     Eigen::Vector3d g(0.0, 1.0, math_utility::deg2rad(0.0));
//     double max_curvature = 0.1;
    
    
    
//     rs_paths::Path path = rs_path.calc_shortest_path(s, g, max_curvature);
//     std::vector<double> rc;
//     std::vector<double> rds;
//     std::cout<< "========1main 0" << path.poses.rows() << std::endl;
//     rs_path.calc_curvature(path.poses, rc, rds);
//     std::vector<double> path_short_x, path_short_y;
//     for(int i=0; i<path.poses.rows(); i++){
//         path_short_x.push_back(path.poses.row(i)[0]);
//         path_short_y.push_back(path.poses.row(i)[1]);
//     }
  

//     std::vector<std::vector<double>> bpath_x, bpath_y;
//     std::vector<rs_paths::Path> paths;
//     rs_path.calc_paths(s, g, max_curvature, paths);
//     for(auto path : paths){
//         std::vector<double> path_info_x, path_info_y;
//         for(int i=0; i<path.poses.rows(); i++){
//             path_info_x.push_back(path.poses.row(i)[0]);
//             path_info_y.push_back(path.poses.row(i)[1]);
//         }
//         bpath_x.push_back(path_info_x);
//         bpath_y.push_back(path_info_y);
//     }
//     // std::cout<< "========1main 1" << path.poses.rows() << std::endl;

//     // // First subplot
//     plt::figure();
//     for(int i=0; i< bpath_x.size(); i++){
//         plt::plot(bpath_x[i], bpath_y[i]);
//     }

//     plt::plot({s.x()}, {s.y()}, "bo"); // Start point
//     plt::plot({g.x()}, {g.y()}, "go");    // End point

//     plt::legend();
//     plt::grid(true);
//     plt::axis("equal");

//     // // // Second subplot for curvature
//     plt::figure();
//     plt::plot(path_short_x, path_short_y, "-r");
//     // plt::plot(rc, ".r");
//     plt::grid(true);
//     // plt::title("Curvature");

//     // // // Show all plots
//     plt::show();




//     return 0;
// }