#include "lib_cpp_hybrid_a_star/trailer_hybrid_a_star.hpp"

namespace planning
{
    TrailerHybridAStar::TrailerHybridAStar(){
    }

    void TrailerHybridAStar::calc_config(Eigen::MatrixXd& obses, double xyreso, double yawreso){
        
        Eigen::VectorXd col_min = obses.colwise().minCoeff(); 
        Eigen::VectorXd col_max = obses.colwise().maxCoeff();
        col_min = col_min.array() - PlannerParams::EXTEND_AREA;
        col_max = col_max.array() + PlannerParams::EXTEND_AREA;
        obses.row(obses.rows()-2) << col_min[0], col_min[1];
        obses.row(obses.rows()-1) << col_max[0], col_max[1];
        int minx = std::round(col_min[0]/xyreso);
        int miny = std::round(col_min[1]/xyreso);
        int maxx = std::round(col_max[0]/xyreso);
        int maxy = std::round(col_max[1]/xyreso);
        
        int xw = maxx - minx;
        int yw = maxy - miny;

        int minyaw = std::round(-M_PI/yawreso) - 1;
        int maxyaw = std::round(M_PI/yawreso);
        int yaww = std::round((maxyaw - minyaw));

        int minyawt = minyaw;
        int maxyawt = maxyaw;
        int yawtw = yaww;

        config_ = Config(minx, miny, minyaw, minyawt, maxx, maxy, maxyaw, maxyawt, xw, yw, yaww, yawtw, xyreso, yawreso);

    }

    void TrailerHybridAStar::calc_holonomic_with_obstacle_heuristic(HybridNode& gnode
                                    , Eigen::MatrixXd& obses, double xyreso){
        Eigen::Vector2d goal_pose = gnode.poses.row(gnode.poses.rows()-1).head(2);
        
        grid_a_star_.calc_dist_policy(goal_pose, obses, xyreso, 1.0, h_dp);
    }

    int TrailerHybridAStar::calc_index(HybridNode& node){
        
        int ind = (node.yawind - config_.minyaw)*config_.xw*config_.yw
                    +(node.yind - config_.miny)*config_.xw + (node.xind - config_.minx);
        // 4D grid
        int rows = node.poses.rows();
        if( rows>0 ){
            double yaw1 = node.poses.row(node.poses.rows()-1)[3];
            int yaw1ind = std::round(yaw1/config_.yawreso);
            ind += (yaw1ind - config_.minyawt) *config_.xw*config_.yw*config_.yaww;
        } else{
            std::cout<< " Error(calc_index): " << std::endl;
        }
        if(ind <= 0){
            std::cout<< " Error(calc_index): " << std::endl;
        }
        return ind;
    }

    void TrailerHybridAStar::calc_motion_inputs(std::vector<double>& u, std::vector<double>& d){
        std::vector<double> up;
        double increment = PlannerParams::MAX_STEER/PlannerParams::N_STEER;
        for (double i=increment; i <= PlannerParams::MAX_STEER+0.002; i+=increment) {
            up.push_back(i);
        }
        u.clear();
        d.clear();

        u.push_back(0.0);
        u.insert(u.end(), up.begin(), up.end());
        std::transform(up.begin(), up.end(), std::back_inserter(u), [](double i) { return -i; });

        int len_u = u.size();
        d.resize(2*len_u);
        double* d_ptr = d.data();
        double* mid_ptr = d_ptr + len_u;
        for (int i = 0; i < len_u; ++i) {
            *(d_ptr++) = 1.0;
            *(mid_ptr++) = -1.0;
        }
        u.reserve(2 * u.size());
        u.insert(u.end(), u.begin(), u.end());
    }

    double TrailerHybridAStar::calc_rs_path_cost(rs_paths::Path& rspath, Eigen::VectorXd& yaw1){
        double cost = 0.0;
        for(int i=0; i<rspath.lengths.size(); i++){
            double l = rspath.lengths[i];
            cost += (rspath.lengths[i] >= 0) ? l : std::abs(l) * PlannerParams::BACK_COST;
        }
        // std::cout<< " 1cost " << cost << std::endl;
        // swich back penalty
        for(int i=0; i<rspath.lengths.size()-1; i++){
            if(rspath.lengths[i] * rspath.lengths[i+1] < 0.0){
                cost += PlannerParams::SB_COST;
            }
        }
        // steer penalyty
        for(int i=0; i<rspath.ctypes.size(); i++){
            if(rspath.ctypes[i] != "S"){
                cost += PlannerParams::STEER_COST*std::abs(PlannerParams::MAX_STEER);
            }
        }
        // steer change penalty
        // calc steer profile
        int nctypes = rspath.ctypes.size();
        std::vector<double> ulist(nctypes, 0.0);
        for(int i=0; i<nctypes; i++){
            if(rspath.ctypes[i] == "R"){
                ulist[i] = - PlannerParams::MAX_STEER;
            }else if(rspath.ctypes[i] == "L"){
                ulist[i] = PlannerParams::MAX_STEER;
            }
        }
        for(int i=0; i<rspath.ctypes.size()-1; i++){
            cost += PlannerParams::STEER_CHANGE_COST*std::abs(ulist[i+1] - ulist[i]);
        }
        Eigen::VectorXd yaw_diff = rspath.poses.col(2) - yaw1;
        Eigen::VectorXd wrapped_yaw_diff = yaw_diff.unaryExpr([](double x) {
            return pi_to_pi(x);
        });
        // calculate the absolute values
        Eigen::VectorXd abs_yaw_diff = wrapped_yaw_diff.array().abs();
        // sum the absolute values
        double sum_abs_yaw_diff = abs_yaw_diff.sum();
        cost += PlannerParams::JACKKNIF_COST * sum_abs_yaw_diff;
        return cost;
    }

    bool TrailerHybridAStar::update_node_with_analystic_expantion(HybridNode& current
                , HybridNode& ngoal,  Eigen::MatrixXd& obses, grid_search::KDTree& kdtree
                , double gyaw1, HybridNode& updated_node){
            
            rs_paths::Path apath;
            
            bool find_path = analystic_expantion(current, ngoal, obses, kdtree, apath);
            if(find_path){
                Eigen::VectorXd steps = apath.poses.col(3).array()*PlannerParams::MOTION_RESOLUTION;
                Eigen::VectorXd current_pose = current.poses.row(current.poses.rows()-1);
                Eigen::VectorXd yaw1;
                bool can_estimate = trailerlib_.calc_trailer_yaw_from_xyyaw(apath.poses, current_pose[3], steps, yaw1);
                if (std::abs(math_utility::pi_to_pi(yaw1[yaw1.size()-1] - gyaw1)) >= PlannerParams::GOAL_TYAW_TH){
                    return false;
                }
                double fcost = current.cost + calc_rs_path_cost(apath, yaw1);
                int fpind = calc_index(current);
                // to get rows from index 1 to the end (2nd to last rows)
                Eigen::MatrixXd updated_poses(apath.poses.rows() - 1, apath.poses.cols()+1);
                updated_poses.block(0, 0, apath.poses.rows() - 1, apath.poses.cols()) = apath.poses.block(1, 0, apath.poses.rows() - 1, apath.poses.cols());
                updated_poses.col(4) = updated_poses.col(3);
                updated_poses.col(3) = yaw1.head(yaw1.size()-1);
                double  fsteer = 0.0;
                HybridNode est_node(current.xind, current.yind, current.yawind
                            , current.direction, updated_poses, fsteer, fcost, fpind);
                updated_node = est_node;
                return true;
            }
        return false;
    }


    void TrailerHybridAStar::get_final_path(std::unordered_map<int, HybridNode>& closed
                                                    , HybridNode& ngoal, HybridNode& nstart, HybridPath& path){
        
        Eigen::MatrixXd g_poses = ngoal.poses.colwise().reverse();
        int nid = ngoal.pind;
        int finalcost = ngoal.cost;
        std::vector<Eigen::MatrixXd> ref_path;
        ref_path.push_back(g_poses);
        int total_rows = g_poses.rows();
        int total_cols = 5;
        while(true){
            if (closed.find(nid) == closed.end()) {
                std::cout<< nid << " cant find the requested node" << std::endl;
                break;
            }
            HybridNode n = closed[nid];
            Eigen::MatrixXd n_poses = n.poses.colwise().reverse();
            ref_path.push_back(n_poses);
            total_rows += n_poses.rows();
            nid = n.pind;
            if (is_same_grid(n, nstart)){
                break;
            }
        }
        Eigen::MatrixXd final_path(total_rows, total_cols);
        int current_row = 0;
        for (const auto& mat : ref_path) {
            final_path.block(current_row, 0, mat.rows(), mat.cols()) = mat;
            current_row += mat.rows();
        }
        Eigen::MatrixXd final_path_up  = final_path.colwise().reverse();
        if(final_path_up.rows() > 2 && final_path_up.cols() > 3 ){
            final_path_up.row(0)[3] = final_path_up.row(1)[3];
            HybridPath path_final(final_path_up, finalcost);
            path = path_final;
        }else{
            std::cout << "Path can not be found..." << std::endl;
        }
    }

    double TrailerHybridAStar::calc_cost(HybridNode& n, HybridNode& ngoal){
        int index_x = (n.xind - config_.minx) - 1;
        int index_y = (n.yind - config_.miny) - 1;
        double total_cost = n.cost + PlannerParams::H_COST*h_dp(index_x, index_y);
        return total_cost;
    }


    bool TrailerHybridAStar::calc_hybrid_astar_path(Eigen::Vector4d s, Eigen::Vector4d g
                        , Eigen::MatrixXd& obses, HybridPath& path){
            double syaw = math_utility::pi_to_pi(s[2]);
            double gyaw = math_utility::pi_to_pi(g[2]);

            double xyreso = PlannerParams::XY_GRID_RESOLUTION;
            double yawreso = PlannerParams::YAW_GRID_RESOLUTION;
            grid_search::PointCloud cloud;
            cloud.points = obses;
            grid_search::KDTree kdtree(2 , cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
            kdtree.buildIndex();
            calc_config(obses, xyreso, yawreso);
            int xind_s = std::round(s[0]/xyreso);
            int yind_s = std::round(s[1]/xyreso);
            int yawind_s = std::round(s[2]/yawreso);

            int xind_g = std::round(g[0]/xyreso);
            int yind_g = std::round(g[1]/xyreso);
            int yawind_g = std::round(g[2]/yawreso);
            Eigen::MatrixXd poses_s(1,5);
            poses_s << s[0], s[1], s[2], s[3], 1.0;

            Eigen::MatrixXd poses_g(1,5);
            poses_g << g[0], g[1], g[2], g[3], 1.0;

            HybridNode nstart(xind_s, yind_s, yawind_s, true, poses_s, 0.0, 0.0, -1);
            HybridNode ngoal(xind_g, yind_g, yawind_g, true, poses_g, 0.0, 0.0, -1);
            calc_holonomic_with_obstacle_heuristic(ngoal, obses, xyreso);

            std::unordered_map<int, HybridNode> open;
            std::unordered_map<int, HybridNode> closed;

            open[calc_index(nstart)] = nstart;

            std::priority_queue<CostNode, std::vector<CostNode>, CompareCostNode> pq;
            CostNode cost_start(calc_cost(nstart, ngoal), calc_index(nstart));
            pq.push(cost_start);

            std::vector<double> u;
            std::vector<double> d;
            calc_motion_inputs(u, d);

            int nmotion = u.size();
            HybridNode fnode;
            int counter = 0;

            while (true){
                if(open.empty()){
                    std::cout<< "Error: Cannot find path, No open set" << std::endl;
                    break;
                }
                if(pq.empty()){
                    std::cerr << " Can not find a path " << std::endl;
                    break;
                }
                CostNode next_node = pq.top();
                pq.pop();
                int c_id = next_node.id;
                HybridNode current = open[c_id];
                open.erase(c_id);
                closed[c_id] = current;
                HybridNode fpath;
                bool isupdated = update_node_with_analystic_expantion(current, ngoal, obses, kdtree, g[3], fpath);
                if (isupdated){
                    fnode = fpath;
                    break;
                }
                double inityaw1 = current.poses.row(0)[3];
                for(int i=0; i<nmotion; i++){
                    HybridNode node = calc_next_node(current, c_id, u[i], d[i]);
                    if (!verify_index(node, obses, inityaw1, kdtree)){
                        continue;
                    }
                    int node_ind = calc_index(node);
                    //  If it is already in the closed set, skip it
                    if (closed.find(node_ind) != closed.end()) {
                        continue;
                    }
                    if (open.find(node_ind) == open.end()) {
                        open[node_ind] = node;
                        double cost = calc_cost(node, ngoal);
                        CostNode cost_node(cost, node_ind);
                        pq.push(cost_node);
                    }else{
                        if(open[node_ind].cost > node.cost){
                            open[node_ind] = node;
                        }
                    }
                }
            }
            // std::cout<< "final expand node:" << open.size() + closed.size() << std::endl;
            get_final_path(closed, fnode, nstart, path);
            return true;
    }


    bool TrailerHybridAStar::analystic_expantion(HybridNode& node, HybridNode& ngoal
                        ,  Eigen::MatrixXd& obses, grid_search::KDTree& kdtree
                        , rs_paths::Path& selected_path){

        Eigen::VectorXd current_pose = node.poses.row(node.poses.rows()-1);
        Eigen::VectorXd target_pose = ngoal.poses.row(ngoal.poses.rows()-1);
        double sx = current_pose[0];
        double sy = current_pose[1];
        double syaw = current_pose[2];
        double max_curvature = tan(PlannerParams::MAX_STEER)/PlannerParams::WB;

        std::vector<rs_paths::Path> paths;
        rs_path_.calc_paths(current_pose.head(3), target_pose.head(3)
                            , max_curvature, paths, PlannerParams::MOTION_RESOLUTION);

        if(paths.size() == 0){
            return false;
        }
        std::priority_queue<rs_paths::Path, std::vector<rs_paths::Path>, rs_paths::CompareNode> pathqueue;
        

        for(auto path : paths){
            Eigen::VectorXd steps = path.poses.col(3).array()*PlannerParams::MOTION_RESOLUTION;
            Eigen::VectorXd yaw1;
            bool can_estimate = trailerlib_.calc_trailer_yaw_from_xyyaw(path.poses, current_pose[3], steps, yaw1);
            path.cost = calc_rs_path_cost(path, yaw1);
            pathqueue.push(path);
        }

        while (!pathqueue.empty()) {
            rs_paths::Path path = pathqueue.top();
            pathqueue.pop();

            Eigen::VectorXd steps = path.poses.col(3).array()*PlannerParams::MOTION_RESOLUTION;
            Eigen::VectorXd yaws1;
            bool can_estimate = trailerlib_.calc_trailer_yaw_from_xyyaw(node.poses, current_pose[3], steps, yaws1);
            std::vector<int> indices;
            for(int i=0; i<path.poses.rows(); i+=PlannerParams::SKIP_COLLISION_CHECK){
                indices.push_back(i);
            }

            Eigen::MatrixXd selected_poses(indices.size(), path.poses.cols());
            for (size_t i = 0; i < indices.size(); ++i) {
                selected_poses.row(i) = path.poses.row(indices[i]);
            }

            if (trailerlib_.check_trailer_collision(obses, selected_poses, kdtree)){
                selected_path = path;
                return true;
            }
        } 
        return false;       
    }

    bool TrailerHybridAStar::is_same_grid(HybridNode& node1, HybridNode& node2){
        if(node1.xind != node2.xind){
            return false;
        }
        if(node1.yind != node2.yind){
            return false;
        }
        if(node1.yawind != node2.yawind){
            return false;
        }
        return true;
    }

    bool TrailerHybridAStar::verify_index(HybridNode& node, Eigen::MatrixXd& obses
                            , double inityaw1, grid_search::KDTree& kdtree){
        
        if ((node.xind - config_.minx) >= config_.xw){
            return false;
        }else if((node.xind - config_.minx) <= 0){
            return false;
        }
        if ((node.yind - config_.miny) >= config_.yw){
            return false;
        }else if((node.yind - config_.miny) <= 0){
            return false;
        }
         // check collisiton
        Eigen::VectorXd steps = node.poses.col(4).array()*PlannerParams::MOTION_RESOLUTION;
        Eigen::VectorXd yaws1;
        bool can_estimate = trailerlib_.calc_trailer_yaw_from_xyyaw(node.poses, inityaw1, steps, yaws1);
        std::vector<int> indices;
        for(int i=0; i<node.poses.rows(); i+=PlannerParams::SKIP_COLLISION_CHECK){
            indices.push_back(i);
        }
        Eigen::MatrixXd selected_poses(indices.size(), node.poses.cols());
        for (size_t i = 0; i < indices.size(); ++i) {
            selected_poses.row(i) = node.poses.row(indices[i]);
        }
        if (!trailerlib_.check_trailer_collision(obses,  selected_poses, kdtree)){
            return false;
        }
        return true;
    }


    HybridNode TrailerHybridAStar::calc_next_node(HybridNode& current, int c_id, double u, double d){

        double arc_l = PlannerParams::XY_GRID_RESOLUTION*1.5;
        int nlist = std::floor(arc_l/PlannerParams::MOTION_RESOLUTION) + 1;
        Eigen::MatrixXd poses(nlist, 5);
        Eigen::VectorXd current_pose = current.poses.row(current.poses.rows()-1);
        if( nlist >0 ){
            poses.row(0)[0] = current_pose.x() + d * PlannerParams::MOTION_RESOLUTION*cos(current_pose[2]);
            poses.row(0)[1] = current_pose.y() + d * PlannerParams::MOTION_RESOLUTION*sin(current_pose[2]);
            poses.row(0)[2] = math_utility::pi_to_pi(current_pose[2] 
                                + d*PlannerParams::MOTION_RESOLUTION/PlannerParams::WB * tan(u));
            poses.row(0)[3] = math_utility::pi_to_pi(current_pose[3] 
                        + d*PlannerParams::MOTION_RESOLUTION/PlannerParams::LT*sin(current_pose[2]-current_pose[3]));
        }

        for(int i=0; i<(nlist-1); i++){
            poses.row(i+1)[0] = poses.row(i)[0] 
                                + d * PlannerParams::MOTION_RESOLUTION*cos(poses.row(i)[2]);
            poses.row(i+1)[1] = poses.row(i)[1] 
                                + d * PlannerParams::MOTION_RESOLUTION*sin(poses.row(i)[2]);
            poses.row(i+1)[2] = math_utility::pi_to_pi(poses.row(i)[2] 
                                + d*PlannerParams::MOTION_RESOLUTION/PlannerParams::WB * tan(u));  
            poses.row(i+1)[3] = math_utility::pi_to_pi(poses.row(i)[3] 
                                + d*PlannerParams::MOTION_RESOLUTION/PlannerParams::LT*sin(poses.row(i)[2]-poses.row(i)[3]));
        }

        Eigen::VectorXd updated_pose = poses.row(poses.rows()-1);
        int xind = std::round(updated_pose[0]/config_.xyreso);
        int yind = std::round(updated_pose[1]/config_.xyreso);
        int yawind = std::round(updated_pose[2]/config_.yawreso);
        double addedcost = 0.0;
        double direction = 1.0;
        if (d > 0){
            direction = 1.0;
            addedcost += std::abs(arc_l);
        }else{
            direction = 0.0;
            addedcost += std::abs(arc_l) * PlannerParams::BACK_COST;
        }
        // swich back penalty
        if (direction != current.direction){ // switch back penalty
            addedcost += PlannerParams::SB_COST;
        }
        // steer penalyty
        addedcost += PlannerParams::STEER_COST*std::abs(u);
        // steer change penalty
        addedcost += PlannerParams::STEER_CHANGE_COST*std::abs(current.steer - u);
        double total_angle_diff = (poses.col(2) - poses.col(3))
                    .unaryExpr([](double angle) { return std::abs(math_utility::pi_to_pi(angle)); })
                    .sum();
        // jacknif cost
        addedcost += PlannerParams::JACKKNIF_COST*total_angle_diff;
        double cost = current.cost + addedcost;
        poses.col(4).setConstant(direction);
        HybridNode next_node(xind, yind, yawind, direction, poses, u, cost, c_id);
        return next_node;
    }
    
}


