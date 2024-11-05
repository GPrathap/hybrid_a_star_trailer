#include "lib_cpp_hybrid_a_star/grid_a_star.hpp"


namespace grid_search
{

    GridAStar::GridAStar(){

    }


    void GridAStar::calc_obstacle_map(Eigen::MatrixXd& obses, double& reso, double& vr){

        Eigen::VectorXd col_min = obses.colwise().minCoeff(); 
        Eigen::VectorXd col_max = obses.colwise().maxCoeff();

        minx_ = std::floor(col_min[0]);
        miny_ = std::floor(col_min[1]);
        maxx_ = std::ceil(col_max[0]);
        maxy_ = std::ceil(col_max[1]);
        
        xwidth_ = maxx_ - minx_;
        ywidth_ = maxy_ - miny_;

        std::cout<< "xwidth_: " << xwidth_ << " ywidth_: " << ywidth_ << std::endl;
        obs_map.resize(xwidth_, ywidth_);
        obs_map.setZero();

        PointCloud cloud;
        cloud.points = obses;

        KDTree kd_tree(2 , cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
        kd_tree.buildIndex();
        const size_t num_results = 1;
        std::vector<size_t> ret_indexes(num_results); 
        std::vector<double> out_dist_sqr(num_results);
        nanoflann::KNNResultSet<double> resultSet(num_results);
        for(int ix=0; ix<xwidth_; ix++ ){
            int x = ix + minx_;
            for(int iy=0; iy<ywidth_; iy++ ){
                int y = iy + miny_;
                const double query_pt[2] = {x*1.0, y*1.0};
                resultSet.init(ret_indexes.data(), out_dist_sqr.data());
                kd_tree.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParameters(10));
                double dis = std::sqrt(out_dist_sqr[0]);
                if (dis <= vr/reso){
                    obs_map(ix,iy) = 1;
                } 
            }
        }
    }

    int GridAStar::calc_index(const Node& node) const {
        return static_cast<int>((node.pose.y() - miny_)*xwidth_ + (node.pose.x() - minx_));
    }

    double GridAStar::calc_cost(const Node& n, const Node& ngoal) const {
        return (n.cost + h(n.pose.x() - ngoal.pose.x(), n.pose.y() - ngoal.pose.y()));
    }

    double GridAStar::h(int x, int y) const {
        return sqrt( pow(x,2) + pow(y,2));
    }

    Eigen::MatrixXd GridAStar::get_motion_model() const {
        Eigen::MatrixXd motion(8, 3);
        motion<< 1, 0, 1, 0, 1, 1, -1, 0, 1, 0, -1, 1, -1, -1, sqrt(2), -1, 1, sqrt(2)
        , 1, -1, sqrt(2), 1, 1, sqrt(2);
        return motion;
    }

    void GridAStar::calc_dist_policy(Eigen::Vector2d g
                        , Eigen::MatrixXd obses, double reso, double vr
                        , Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& pmap){
        
        Eigen::Vector2i goal_pose (g.x()/reso, g.y()/reso);
        Node ngoal(goal_pose, 0.0, -1);  
        for(int i=0; i<obses.rows(); i++){
            obses.row(i) = obses.row(i)/reso;
        }   
        calc_obstacle_map(obses, reso, vr);

        std::unordered_map<int, Node> open;
        std::unordered_map<int, Node> closed;
        open[calc_index(ngoal)] = ngoal;

        Eigen::MatrixXd motion = get_motion_model();
        int nmotion = motion.rows();
        std::priority_queue<CostNode, std::vector<CostNode>, CompareNode> pq;
        CostNode cost_goal(ngoal.cost, calc_index(ngoal));
        pq.push(cost_goal);
        int counter = 0;
        while( true ){
            if(open.empty()){
                std::cout << "finish search " << std::endl;
                break;
            }
            CostNode next_node = pq.top();
            pq.pop();
            int c_id = next_node.id;
            Node current = open[c_id];
            open.erase(c_id);
            closed[c_id] = current;
            for(int i=0; i<nmotion; i++){
                Node node(Eigen::Vector2i(current.pose.x() + motion.row(i)[0], current.pose.y() + motion.row(i)[1])
                                , current.cost + motion.row(i)[2], c_id);
                
                if(!verify_node(node)){
                    continue;
                }
                int node_ind = calc_index(node);
                if (closed.find(node_ind) != closed.end()) {
                    continue;
                }
                if (open.find(node_ind) != open.end()) {
                    if(open[node_ind].cost > node.cost){
                        open[node_ind].cost = node.cost;
                        open[node_ind].pind = c_id;
                    }
                }else{
                    open[node_ind] = node;
                    CostNode cost_node(node.cost, calc_index(node));
                    pq.push(cost_node);
                }
            }
        }

        std::cout<< "  grid search length(closed) " << closed.size() << std::endl;
        calc_policy_map(closed, pmap);
    }


    Node GridAStar::search_min_cost_node(const std::unordered_map<int, Node*> open, const Node ngoal){
        double mcost = std::numeric_limits<double>::infinity();
        Node mnode = ngoal;
        for (auto it = open.begin(); it != open.end(); ++it) {
            double cost = it->second->cost + h(it->second->pose.x() - ngoal.pose.x(), it->second->pose.y() - ngoal.pose.y());
            if(mcost > cost){
                mnode = *(it->second);
                mcost = cost;
            }
        }
        return mnode;
    }

    bool GridAStar::verify_node(Node node) const {
        if ((node.pose.x() - minx_) >= xwidth_){
            return false;
        } else if (( node.pose.x() - minx_) < 0 ){
            return false;
        } 
        if ((node.pose.y() - miny_) >= ywidth_){
            return false;
        } else if (( node.pose.y() - miny_) < 0 ){
            return false;
        }       
        if (obs_map(node.pose.x()-minx_, node.pose.y()-miny_) > 0){
            return false;
        }
        return true;                    
    }

    void GridAStar::calc_policy_map(std::unordered_map<int, Node>& closed
                    ,  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& pmap) const {
        pmap.resize(xwidth_, ywidth_);
        pmap.setConstant(std::numeric_limits<double>::infinity());
        for (auto it = closed.begin(); it != closed.end(); ++it) {
            pmap(it->second.pose.x()-minx_, it->second.pose.y()-miny_) = it->second.cost;
        }
    }

    void GridAStar::calc_astar_path(Eigen::Vector2d s, Eigen::Vector2d g
                        , Eigen::MatrixXd obses, double reso, double vr){
        
        Eigen::Vector2i goal_pose (g.x()/reso, g.y()/reso);
        Eigen::Vector2i start_pose (s.x()/reso, s.y()/reso);
        Node ngoal(goal_pose, 0.0, -1);  
        Node nstart(start_pose, 0.0, -1);  

        for(int i=0; i<obses.rows(); i++){
            obses.row(i) = obses.row(i)/reso;
        }   
        calc_obstacle_map(obses, reso, vr);

        std::unordered_map<int, Node> open;
        std::unordered_map<int, Node> closed;

        open[calc_index(nstart)] = nstart;

        Eigen::MatrixXd motion = get_motion_model();
        int nmotion = motion.rows();
        std::priority_queue<CostNode, std::vector<CostNode>, CompareNode> pq;
        CostNode cost_start(calc_cost(nstart, ngoal), calc_index(nstart));
        pq.push(cost_start);
        while( true ){
            if(open.empty()){
                std::cerr << " Error: No open set " << std::endl;
                break;
            }

            if(pq.empty()){
                std::cerr << " Can not find a path " << std::endl;
                break;
            }
            CostNode next_node = pq.top();
            pq.pop();
            int c_id = next_node.id;
            Node current = open[c_id];
            if( (current.pose.x() == ngoal.pose.x())&& (current.pose.y() == ngoal.pose.y())){
                std::cout << " Goal!! " << std::endl;
                closed[c_id] = current;
                break;
            }
            open.erase(c_id);
            closed[c_id] = current;
            for(int i=0; i<nmotion; i++){
                Node node(Eigen::Vector2i(current.pose.x() + motion.row(i)[0], current.pose.y() + motion.row(i)[1])
                                , current.cost + motion.row(i)[2], c_id);
                if(!verify_node(node)){
                    continue;
                }
                int node_ind = calc_index(node);
                if (closed.find(node_ind) != closed.end()) {
                    continue;
                }
                if (open.find(node_ind) != open.end()) {
                    if(open[node_ind].cost > node.cost){
                        open[node_ind].cost = node.cost;
                        open[node_ind].pind = c_id;
                    }
                }else{
                    open[node_ind] = node;
                    CostNode cost_node(calc_cost(node, ngoal), calc_index(node));
                    pq.push(cost_node);
                }
            }
        }
        get_final_path(closed, ngoal, nstart, reso);
    }

    void GridAStar::get_final_path(std::unordered_map<int, Node>& closed, const Node ngoal
                            , const Node nstart, const double reso){

        std::vector<Eigen::Vector2i> path;
        Eigen::Vector2i r = ngoal.pose;
        path.push_back(r);
        int nid = calc_index(ngoal);
        if(closed.empty()){
            std::cerr << " Can not find a path " << std::endl;
            return;
        }
        while(true){
            if (closed.find(nid) == closed.end()) {
                std::cerr << " Can not find the full path " << std::endl;
                break;
            }
            Node next = closed[nid];
            path.push_back(next.pose);
            nid = next.pind;
            if((path.back().x() == nstart.pose.x()) & (path.back().y() == nstart.pose.y())){
                std::cout<< "find grid astar path " << std::endl;
                break;
            }
        }
        
        final_path_.clear();
        for (auto it = path.rbegin(); it != path.rend(); ++it) {
            Eigen::Vector2d next_pose(it->x()*reso, it->y()*reso);
            final_path_.push_back(next_pose);
        }
    }
}

// int main(int argc, char *argv[]){
    
//     Eigen::Vector2d s(10.0, 10.0);
//     Eigen::Vector2d g(50.0, 50.0);

//     Eigen::MatrixXd obss(320, 2);
//     int obs_index = 0;
//     for(int i=0; i< 60; i++){
//         obss.row(obs_index) << i*1.0, 0.0;
//         obs_index++;
//     }
//     for(int i=0; i< 60; i++){
//         obss.row(obs_index) << 60.0, i*1.0;
//         obs_index++;
//     }
//     for(int i=0; i< 60; i++){
//         obss.row(obs_index) << i*1.0, 60.0;
//         obs_index++;
//     }
//     for(int i=0; i< 60; i++){
//         obss.row(obs_index) << 0.0, i*1.0;
//         obs_index++;
//     }
//     for(int i=0; i< 40; i++){
//         obss.row(obs_index) << 20.0, i*1.0;
//         obs_index++;
//     }
//     for(int i=0; i< 40; i++){
//         obss.row(obs_index) << 40.0, 60.0-i*1.0;
//         obs_index++;
//     }

//     double VEHICLE_RADIUS = 5.0;
//     double GRID_RESOLUTION = 1.0;

//     grid_search::GridAStar grid_a_star;
//     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> pmap;
//     // grid_a_star.calc_dist_policy(s, g, obss, GRID_RESOLUTION, VEHICLE_RADIUS, pmap);
//     grid_a_star.calc_astar_path(s, g, obss, GRID_RESOLUTION, VEHICLE_RADIUS);
//     std::vector<Eigen::Vector2d> path_est = grid_a_star.getPath();
//     std::vector<double> path_x, path_y;
//     for(int i=0; i< path_est.size(); i++){
//         path_x.push_back(path_est[i][0]);
//         path_y.push_back(path_est[i][1]);
//     }

//     plt::figure();
//     std::vector<double> ox, oy;
//     for(int i=0; i< obss.rows(); i++){
//         ox.push_back(obss.row(i)[0]);
//         oy.push_back(obss.row(i)[1]);
//     }
//     plt::plot(ox, oy, ".r");
//     plt::plot(path_x, path_y, ".y");
//     plt::plot({s.x()}, {s.y()}, "bo"); // Start point
//     plt::plot({g.x()}, {g.y()}, "go");    // End point
//     plt::show();

// }
