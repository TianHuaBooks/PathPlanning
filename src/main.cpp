#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

    //std::cout << "maps_s.size:" << maps_s.size() << " maps_x:" << maps_x.size() << std::endl;
	while((prev_wp < (int)(maps_s.size()-1)) &&s > maps_s[prev_wp+1])
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

#define SPEED_UNIT 0.224        // 0.5 m/s
#define METER_PER_SECOND 2.24   // 5 m/s
#define MAX_SPEED 49.3          // max speed 50 m/h
#define INIT_SPEED 4            // speed of cold start
#define DISTANCE_UNIT 30        // distance unit 30 meters
#define TOTAL_POINTS 50         // total points projects
#define SAMPLE_TIME 0.02        // 0.02 second (20 ms)
#define LANE_WIDTH 4            // lane width = 4 meters
#define HALF_LANE_WIDTH 2
#define SAFE_DIST 30            // safety distance 30 meters
#define RIGHTEST_LANE 2         // rightest lane index, couting from double-yellow lines as 0
// Fusion data index
#define IDX_VX 3
#define IDX_VY 4
#define IDX_S  5
#define IDX_D  6

class Planner {
public:
    Planner(vector<double>& previous_path_x, vector<double>& previous_path_y, double car_x, double car_y, double car_s, double car_yaw)
    : m_previous_path_x(previous_path_x), m_previous_path_y(previous_path_y),
      m_car_x(car_x), m_car_y(car_y), m_car_s(car_s), m_car_yaw(deg2rad(car_yaw))
    {}
    
    ~Planner() {}
    
    // A test function with frenet coord to get next points
    void get_frenet_points(vector<double>& map_waypoints_x, vector<double>& map_waypoints_y,
                         vector<double>& map_waypoints_s,
                         vector<double>& next_x_vals, vector<double>& next_y_vals){
        double dist_inc = 0.5;
        double d = LANE_WIDTH + HALF_LANE_WIDTH;
        for (int i = 0; i < TOTAL_POINTS; i++) {
            double s = m_car_s + i * dist_inc;
            vector<double> xy = getXY(s, d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            if (i % 10 == 0)
                std::cout << i << ", xy[0]:(" << xy[0] << "," << xy[1] << ")\n";
            next_x_vals.push_back(xy[0]);
            next_y_vals.push_back(xy[1]);
        }
    }
    
    // get next points
    void get_next_points(vector<double>& map_waypoints_x, vector<double>& map_waypoints_y,
         vector<double>& map_waypoints_s, double end_path_s, vector<vector<double>>& sensor_fusion,
         vector<double>& next_x_vals, vector<double>& next_y_vals){

        // array of pts for way pts projection
        vector<double> x_vals, y_vals;
        // reference point and yaw angle
        double x_ref = m_car_x;
        double y_ref = m_car_y;
        double yaw_ref = m_car_yaw;
        
        // first, check distance and speed of the car ahead
        bool too_close = check_car_ahead(m_car_s, sensor_fusion, (int) m_previous_path_x.size());

        // adjust speed
        if (too_close)
            m_ref_vel -= SPEED_UNIT;
        else if (m_ref_vel < MAX_SPEED)
            m_ref_vel += SPEED_UNIT;
        
        // add previous pts into next_x_vals and next_y_vals
        //       and (x_vals, y_vals)
        add_prev_pts(x_vals, y_vals, x_ref, y_ref, yaw_ref, next_x_vals, next_y_vals);
        
        // project three points far away based on waypoints
        add_way_points(x_vals, y_vals, map_waypoints_x, map_waypoints_y, map_waypoints_s);
    
        // transform (x_vals, y_vals) with translation and rotation
        transform_to_car_coord(x_vals, y_vals, x_ref, y_ref, yaw_ref);
        
        // create a spline from (x_vals, y_vals)
        m_spline.set_points(x_vals, y_vals);
        
        // add additional pts to make total (50)
        int addl_pts = TOTAL_POINTS - (int) m_previous_path_x.size();
        if (addl_pts > 0)
            add_addl_points_to_total(next_x_vals, next_y_vals, x_ref, y_ref, yaw_ref, addl_pts);
    }
    
private:
    // utility to add previous pts
    void add_prev_pts(vector<double>& x_vals, vector<double>& y_vals, double& x_ref, double& y_ref, double& yaw_ref, vector<double>& next_x_vals, vector<double>& next_y_vals) {
        if (m_previous_path_x.size() < 2) {
            double prev_x = m_car_x - cos(yaw_ref);
            double prev_y = m_car_y - sin(yaw_ref);
            x_vals.push_back(prev_x);
            y_vals.push_back(prev_y);
            if (m_car_x > prev_x) {
                // spline doesn't like dup or out of order
                x_vals.push_back(m_car_x);
                y_vals.push_back(m_car_y);
            }
        } else {
            // add the last two pts as starting pt
            int last = (int) m_previous_path_x.size() - 1;
            x_ref = m_previous_path_x[last];
            y_ref = m_previous_path_y[last];
            int last_prev = last - 1;
            if (last_prev >= 0) {
                if (x_ref != m_previous_path_x[last_prev]) {
                    // spline doesn't like dup
                    x_vals.push_back(m_previous_path_x[last_prev]);
                    y_vals.push_back(m_previous_path_y[last_prev]);
                    yaw_ref = atan2(y_ref - m_previous_path_y[last_prev], x_ref - m_previous_path_x[last_prev]);
                }
            }
            x_vals.push_back(x_ref);
            y_vals.push_back(y_ref);
        }
        // add all previuos points to next points
        for (auto i = 0; i < m_previous_path_x.size(); i++) {
            next_x_vals.push_back(m_previous_path_x[i]);
            next_y_vals.push_back(m_previous_path_y[i]);
        }
    }
    
    // utility to check safty regards to distance from the car ahead
    // return true if too close
    bool  check_car_ahead(double end_path_s, vector<vector<double>>& sensor_fusion, int sample_count) {
        bool too_close = false;
        double my_lane = m_lane * LANE_WIDTH;
        for (vector<double> data : sensor_fusion) {
            if ((data[IDX_D] >= my_lane) && (data[IDX_D] <= (my_lane + LANE_WIDTH))) {
                if (data[IDX_S] > end_path_s) {
                    double slow_vel = sqrt(data[IDX_VX] * data[IDX_VX] + data[IDX_VY] * data[IDX_VY]);
                    double dist = data[IDX_S] - end_path_s - SAMPLE_TIME * sample_count * slow_vel;
                    if (dist < SAFE_DIST) {
                        if (change_lane(end_path_s, sensor_fusion, sample_count, slow_vel))
                            std::cout << "change lane\n";
                        else
                            std::cout << "front car too close, dist:" << dist << std::endl;
                        //m_ref_vel -= SPEED_UNIT;
                        too_close = true;
                        break;
                    }
                }
            }
        }
        
        return too_close;
    }

    // utility to change lane if feasible
    bool  change_lane(double end_path_s, vector<vector<double>>& sensor_fusion, int sample_count, double slow_car_vel)
    {
        if (m_lane > 0) {
            // check left lane first
            if (is_lane_safe(m_lane - 1, end_path_s, sensor_fusion, sample_count, slow_car_vel))
            {
                m_lane--;
                return true;
            }
        }
        if (m_lane < RIGHTEST_LANE) {
            // check right lane
            if (is_lane_safe(m_lane + 1, end_path_s, sensor_fusion, sample_count, slow_car_vel)) {
                m_lane++;
                return true;
            }
        }
        return false;
    }
    
    bool  is_lane_safe(int idx, double end_path_s, vector<vector<double>>& sensor_fusion, int sample_count, double show_car_vel) {
        bool safe = true;
        double lane = idx * LANE_WIDTH;
        for (vector<double> data : sensor_fusion) {
            if ((data[IDX_D] >= lane) && (data[IDX_D] <= (lane + LANE_WIDTH))) {
                double vel = sqrt(data[IDX_VX]*data[IDX_VX] + data[IDX_VY]*data[IDX_VY]);
                double dist_to_be = SAMPLE_TIME * sample_count * vel;
                if (data[IDX_S] > end_path_s) {
                    double dist = data[IDX_S] - end_path_s - dist_to_be;
                    if (dist < SAFE_DIST) {
                        std::cout << "Lane:" << idx << " not safe to change dist:" << dist << std::endl;
                        safe = false;
                        break;
                    } else if ((dist < SAFE_DIST+SAFE_DIST) && (vel < show_car_vel)) {
                        std::cout << "Lane:" << idx << " not safe to change vel:" << vel << " slow vel:" << show_car_vel << std::endl;
                        safe = false;
                    }
                } else {
                    double dist = end_path_s - dist_to_be - data[IDX_S];
                    if (dist < SAFE_DIST) {
                        std::cout << "Lane:" << idx << " not safe to change (behind) dist:" << dist << std::endl;
                        safe = false;
                        break;
                    }
                }
            }
        }
        return safe;
    }
    
    // utility to add waypoints
    void add_way_points(vector<double>& x_vals, vector<double>& y_vals,
                        vector<double>& map_waypoints_x, vector<double>& map_waypoints_y,
                        vector<double>& map_waypoints_s)
    {
        // project points far away based on waypoints
        double d = LANE_WIDTH * m_lane + HALF_LANE_WIDTH; // lane width is 4 meters
        double s = m_car_s + DISTANCE_UNIT;
        for (int i = 0; i < 3; i++) {
            vector<double> next_wp = getXY(s, d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
            x_vals.push_back(next_wp[0]);
            y_vals.push_back(next_wp[1]);
            s += DISTANCE_UNIT;
        }
    }
    
    // utility to add additional pts to make it to 50
    void add_addl_points_to_total(vector<double>& next_x_vals, vector<double>& next_y_vals, double x_ref, double y_ref, double yaw_ref, int count) {
        double target_x = DISTANCE_UNIT;
        double target_y = m_spline(target_x);
        double target_dist = sqrt(target_x*target_x + target_y*target_y);
        double N = target_dist / (SAMPLE_TIME * m_ref_vel / METER_PER_SECOND);
        double unit_x = target_x / N;
        double x_add_on = 0;
        
        for (auto i = 0; i < count; i++) {
            double x = x_add_on + unit_x;
            double y = m_spline(x);
            x_add_on = x;
            // rotation transformation and translation back to global coord
            double x_pt = x_ref + x * cos(yaw_ref) - y * sin(yaw_ref);
            double y_pt = y_ref + x * sin(yaw_ref) + y * cos(yaw_ref);
            next_x_vals.push_back(x_pt);
            next_y_vals.push_back(y_pt);
        }
    }
    
    // utility to perform translation and rotaion
    void transform_to_car_coord(vector<double>& x_vals, vector<double>& y_vals, double x_org, double y_org, double yaw) {
        for (auto i = 0; i < x_vals.size(); i++) {
            double shift_x = x_vals[i] - x_org;
            double shift_y = y_vals[i] - y_org;
            double angle = 0 - yaw;
            x_vals[i] = shift_x * cos(angle) - shift_y * sin(angle);
            y_vals[i] = shift_x * sin(angle) + shift_y * cos(angle);
            //std::cout << "trax :(" << x_vals[i]  << "," << y_vals[i] << ")\n";
        }
    }
    
    // internal states
    static int m_lane; // lane number starting from double yellow lines as zero
    static double m_ref_vel; // reference velocity miles/hr
    
    // states from simulator
    vector<double>& m_previous_path_x;
    vector<double>& m_previous_path_y;
    double m_car_x;
    double m_car_y;
    double m_car_s;
    double m_car_yaw;
    
    // spline
    tk::spline m_spline;
};

int Planner::m_lane = 1; // middle of the 3 lanes
double Planner::m_ref_vel = INIT_SPEED;

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	vector<double> previous_path_x = j[1]["previous_path_x"];
          	vector<double> previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	vector<vector<double>> sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;
            
            double end_s = previous_path_x.size() > 0 ? end_path_s : car_s;

            // Use planner to figure next points
            Planner myPlanner(previous_path_x, previous_path_y, car_x, car_y, end_s, car_yaw);
            myPlanner.get_next_points(map_waypoints_x, map_waypoints_y, map_waypoints_s,
                                      end_path_s, sensor_fusion,
                                      next_x_vals, next_y_vals);

            // Send (next_x_vals, next_y_vals) to simulator
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
