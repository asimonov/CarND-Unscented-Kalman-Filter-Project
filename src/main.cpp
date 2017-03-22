
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include "Eigen/Dense"
#include "ukf.h"
#include "ground_truth_package.h"
#include "measurement_package.h"
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

void check_arguments(int argc, char* argv[]) {
  string usage_instructions = "Usage instructions: ";
  usage_instructions += argv[0];
  usage_instructions += " path/to/input.txt output.txt";

  bool has_valid_args = false;

  // make sure the user has provided input and output files
  if (argc == 1) {
    cerr << usage_instructions << endl;
  } else if (argc == 2) {
    cerr << "Please include an output file.\n" << usage_instructions << endl;
  } else if (argc == 3) {
    has_valid_args = true;
  } else if (argc > 3) {
    cerr << "Too many arguments.\n" << usage_instructions << endl;
  }

  if (!has_valid_args) {
    exit(EXIT_FAILURE);
  }
}

void check_files(ifstream& in_file, string& in_name,
                 ofstream& out_file, string& out_name) {
  if (!in_file.is_open()) {
    cerr << "Cannot open input file: " << in_name << endl;
    exit(EXIT_FAILURE);
  }

  if (!out_file.is_open()) {
    cerr << "Cannot open output file: " << out_name << endl;
    exit(EXIT_FAILURE);
  }
}



int main(int argc, char* argv[]) {

  check_arguments(argc, argv);

  string in_file_name_ = argv[1];
  ifstream in_file_(in_file_name_.c_str(), ifstream::in);

  string out_file_name_ = argv[2];
  ofstream out_file_(out_file_name_.c_str(), ofstream::out);

  check_files(in_file_, in_file_name_, out_file_, out_file_name_);

  /**********************************************
   *  Set Measurements                          *
   **********************************************/

  vector<MeasurementPackage> measurement_pack_list;
  vector<GroundTruthPackage> gt_pack_list;
  string line;

  // prep the measurement packages (each line represents a measurement at a
  // timestamp)
  while (getline(in_file_, line)) {
    string sensor_type;
    MeasurementPackage meas_package;
    GroundTruthPackage gt_package;
    istringstream iss(line);
    long timestamp;

    // reads first element from the current line
    iss >> sensor_type;

    if (sensor_type.compare("L") == 0) {
      // laser measurement

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::LASER;
      meas_package.raw_measurements_ = VectorXd(2);
      float px;
      float py;
      iss >> px;
      iss >> py;
      meas_package.raw_measurements_ << px, py;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    } else if (sensor_type.compare("R") == 0) {
      // radar measurement

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::RADAR;
      meas_package.raw_measurements_ = VectorXd(3);
      float ro;
      float theta;
      float ro_dot;
      iss >> ro;
      iss >> theta;
      iss >> ro_dot;
      meas_package.raw_measurements_ << ro, theta, ro_dot;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    }

    // read ground truth data to compare later
    float x_gt;
    float y_gt;
    float vx_gt;
    float vy_gt;
    iss >> x_gt;
    iss >> y_gt;
    iss >> vx_gt;
    iss >> vy_gt;
    gt_package.gt_values_ = VectorXd(4);
    gt_package.gt_values_ << x_gt, y_gt, vx_gt, vy_gt;
    gt_pack_list.push_back(gt_package);
  }

  // Create a UKF instance
  UKF ukf;
  // configure what we want to process
  ukf.use_laser_ = true;
  ukf.use_radar_ = true;

  // used to compute the RMSE later
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;

  //for (int i=0; i<100; i++)
  //{
  //  float a = (i-50)*M_PI/10;
  //  cout << a << " " << ukf.NormaliseAngle(a) << endl;
  //}

  size_t number_of_measurements = measurement_pack_list.size();
//  size_t number_of_measurements = 25;


  // start filtering from the second frame (the speed is unknown in the first frame)
  for (size_t k = 0; k < number_of_measurements; ++k) {

    // Call the UKF-based fusion
    cout << endl;
    cout << endl;
    cout << endl;
    cout << "processing measurement " << (k+1) << endl;
    ukf.ProcessMeasurement(measurement_pack_list[k]);

    // output into out file

    // time, secs
    char ch = out_file_name_[out_file_name_.length()-1-4];
    long inittime = (ch=='1' ? 1477010443399637 : 1477010443349642);
    long tdiff = measurement_pack_list[k].timestamp_ - inittime;
    float t = float(tdiff)/1000000.;
    out_file_ << t << "\t";

    // output the ground truth px,py,vx,vy
    out_file_ << gt_pack_list[k].gt_values_(0) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(1) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(2) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(3) << "\t";

    // output the estimation
    out_file_ << ukf.x_(0) << "\t"; // pos1 - est
    out_file_ << ukf.x_(1) << "\t"; // pos2 - est
    out_file_ << ukf.x_(2)*cos(ukf.x_(3)) << "\t"; // vx - est
    out_file_ << ukf.x_(2)*sin(ukf.x_(3)) << "\t"; // vy - est

    // output the measurements
    if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER) {
      // output the estimation

      // p1 - meas
      out_file_ << measurement_pack_list[k].raw_measurements_(0) << "\t";
      // p2 - meas
      out_file_ << measurement_pack_list[k].raw_measurements_(1) << "\t";
      // NIS
      out_file_ << ukf.NIS_laser_ << "\t";

    } else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR) {
      // output the estimation in the cartesian coordinates
      float ro = measurement_pack_list[k].raw_measurements_(0);
      float phi = measurement_pack_list[k].raw_measurements_(1);
      out_file_ << ro * cos(phi) << "\t"; // p1_meas
      out_file_ << ro * sin(phi) << "\t"; // p2_meas
      // NIS
      out_file_ << ukf.NIS_radar_ << "\t";
    }

    // we do not have ground truth for the next state variables, but output their estimates anyway
    out_file_ << ukf.x_(2) << "\t"; // vel_abs -est
    out_file_ << ukf.x_(3) << "\t"; // yaw_angle -est
    out_file_ << ukf.x_(4) << "\t"; // yaw_rate -est

    // output estimates of sigmax, sigmay, correlationxy, anglexy
    MatrixXd v = ukf.P_.topLeftCorner(2, 2); // positions covariance. posterior
    MatrixXd v_prior = ukf.P_prior_.topLeftCorner(2, 2); // positions covariance. prior, before measurement

    // posterior
    // get eigenvalues for ellipse width/height and angle of rotation
    // based on http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
    Eigen::EigenSolver<MatrixXd> es(v);
    VectorXd eval = es.eigenvalues().real();
    cout << "Eigenvalues of v:" << endl << eval << endl;
    MatrixXd evec = es.eigenvectors().real();
    cout << "Eigenvectors of v:" << endl << evec << endl;
    float ellipse_width = 0.;
    float ellipse_height = 0.;
    float ellipse_angle = 0.;
    if (eval(0) >= eval(1))
    {
      // first eigenvalue is largest
      ellipse_width = eval(0);
      ellipse_height = eval(1);
      ellipse_angle = atan2(evec(0,1), evec(0,0));
    }
    else
    {
      // second eigenvalue is largest
      ellipse_width = eval(1);
      ellipse_height = eval(0);
      ellipse_angle = atan2(evec(1,1), evec(1,0));
    }

    // prior
    // get eigenvalues for ellipse width/height and angle of rotation
    // based on http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
    Eigen::EigenSolver<MatrixXd> es2(v_prior);
    VectorXd eval2 = es2.eigenvalues().real();
    cout << "Eigenvalues of v_prior:" << endl << eval2 << endl;
    MatrixXd evec2 = es2.eigenvectors().real();
    cout << "Eigenvectors of v_prior:" << endl << evec2 << endl;
    float ellipse_width2 = 0.;
    float ellipse_height2 = 0.;
    float ellipse_angle2 = 0.;
    if (eval2(0) >= eval2(1))
    {
      // first eigenvalue is largest
      ellipse_width2 = eval2(0);
      ellipse_height2 = eval2(1);
      ellipse_angle2 = atan2(evec2(0,1), evec2(0,0));
    }
    else
    {
      // second eigenvalue is largest
      ellipse_width2 = eval2(1);
      ellipse_height2 = eval2(0);
      ellipse_angle2 = atan2(evec2(1,1), evec2(1,0));
    }

    // get stdev and correlation in x,y coordinates. stupid!
//    MatrixXd d = MatrixXd(2,2);
//    d.fill(0.0);
//    d(0,0) = v(0,0);
//    d(1,1) = v(1,1);
//    MatrixXd sqrtd = d.cwiseSqrt();
//    MatrixXd dsinv = sqrtd.inverse();
//    MatrixXd r = dsinv * v * dsinv;
//    float sig1 = sqrtd(0,0);
//    float sig2 = sqrtd(1,1);
//    float corr = r(0,1);
//    float angle = acos(corr);
    //out_file_ << sig1 << "\t";
    //out_file_ << sig2 << "\t";
    //out_file_ << corr << "\t";
    //out_file_ << angle << "\t";

    // posterior distribution
    out_file_ << ellipse_width << "\t";
    out_file_ << ellipse_height << "\t";
    out_file_ << ellipse_angle << "\t";
    // prior distribution
    out_file_ << ellipse_width2 << "\t";
    out_file_ << ellipse_height2 << "\t";
    out_file_ << ellipse_angle2 << "\t";

    out_file_ << "\n";

    // save gt and estimates for px, py, vx, vy for RMSE calculation later
    VectorXd x = VectorXd(4);
    x(0) = ukf.x_(0);
    x(1) = ukf.x_(1);
    // calculate vx and vy from v and psi
    x(2) = ukf.x_(2)*cos(ukf.x_(3));
    x(3) = ukf.x_(2)*sin(ukf.x_(3));
    estimations.push_back(x);
    ground_truth.push_back(gt_pack_list[k].gt_values_);
  }

  // close files
  if (out_file_.is_open()) {
    out_file_.close();
  }

  if (in_file_.is_open()) {
    in_file_.close();
  }

  // compute the accuracy (RMSE)
  Tools tools;
  cout << "Accuracy - RMSE:" << endl << tools.CalculateRMSE(estimations, ground_truth) << endl;

  cout << "Done!" << endl;
  return 0;
}
