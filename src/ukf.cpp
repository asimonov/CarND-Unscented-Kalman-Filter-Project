#include <iostream>
#include "ukf.h"

using namespace std;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // state of the filter. not initialised yet
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // previous timestamp
  time_us_ = 0;

  // UKF constants
  // State dimension
  n_x_ = 5;
  // Augmented state dimension
  n_aug_ = 7;
  // Sigma point spreading parameter
  lambda_ = 3.0 - n_x_;
  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i=1; i<2*n_aug_+1; i++)
  {
    weights_(i) = 0.5/(lambda_+n_aug_);
  }

  // initial state vector
  x_ = VectorXd(n_x_);
  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);


  // Process noise
  // standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.0; // assuming max acceleration 6m/s2 and take half of that
  // standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/2.;

  // Laser measurement noise
  // standard deviation position1 in m
  std_laspx_ = 0.3;
  // standard deviation position2 in m
  std_laspy_ = 0.3;

  // Radar measurement noise
  // standard deviation radius in m
  std_radr_ = 0.3;
  // standard deviation angle in rad
  std_radphi_ = M_PI/32.; // 5.6 degrees
  // standard deviation radius change in m/s
  std_radrd_ = 0.3;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "UKF first initialization" << endl;

    // Initialize the state x_ with the first measurement.
    float phi; // will calculate from both types of measurement
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates
      float ro = meas_package.raw_measurements_(0);
      phi = meas_package.raw_measurements_(1);
      float rodot = meas_package.raw_measurements_(2);
      x_ <<   ro * cos(phi),
              ro * sin(phi),
              rodot, // assume psi==phi, i.e. object is moving away along the axis from zero to its position
              phi, // assume psi==phi. why not. any direction is ok here.
              0.0; // assume the turn angle is constant
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      float px = meas_package.raw_measurements_(0);
      float py = meas_package.raw_measurements_(1);
      if (abs(px)>1e-5)
        phi = atan2(py, px);
      else
        phi = M_PI/2.; // assume object is straigh ahead
      x_ << meas_package.raw_measurements_(0),
            meas_package.raw_measurements_(1),
            0.0, // assume radial velosity is zero
            phi, // assume object is pointing along the axis from zero to its position
            0.0; // assume the turn angle is constant
    }

    // initialize state covariance matrix
    P_ << std_radr_*std_radr_, 0.,                  0.,            0.,            0., // use radar measurement uncertainty as more imprecise one
            0.,                std_radr_*std_radr_, 0.,            0.,            0., // use radar measurement uncertainty as more imprecise one
            0.,                0.,                  std_a_*std_a_, 0.,            0., // assume velocity uncertainty is on the order of acceleration noise
            0.,                0.,                  0.,            M_PI*M_PI/16., 0., // assume 45 degrees
            0.,                0.,                  0.,            0.,            M_PI*M_PI/16.; // assume 45 degrees

    time_us_ = meas_package.timestamp_;

    is_initialized_ = true;

    // print the output
    cout << "x_ = " << x_ << endl;
    cout << "P_ = " << P_ << endl;

    // done initializing, no need to predict or update
    return;
  }

  //compute the time elapsed between the current and previous measurements
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;     //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;
  cout << "dt = " << dt << endl;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  if (dt > 0.00001) {
    // predict only if processing non-simultaneous measurements
    Prediction(dt); // Kalman Filter prediction step

    cout << "predict done" << endl;
    cout << "x_ = " << x_ << endl;
    cout << "P_ = " << P_ << endl;
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar update
    UpdateRadar(meas_package);
    cout << "radar update done" << endl;
  } else {
    // Laser update
    UpdateLidar(meas_package);
    cout << "laser update done" << endl;
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
//  /**
//     * Update the state transition matrix F according to the new elapsed time.
//   */
//  ekf_.F_ << 1., 0., dt, 0.,
//          0., 1., 0., dt,
//          0., 0., 1., 0.,
//          0., 0., 0., 1.;
//  /**
//     * Update the process noise covariance matrix.
//   */
//  float dtsq = dt * dt;
//  ekf_.Q_ << (dtsq / 4.) * noise_ax_, 0., (dt / 2.) * noise_ax_, 0.,
//          0., (dtsq / 4.) * noise_ay_, 0., (dt / 2.) * noise_ay_,
//          (dt / 2.) * noise_ax_, 0., noise_ax_, 0.,
//          0., (dt / 2.) * noise_ay_, 0., noise_ay_;
//  ekf_.Q_ *= dtsq;
//  /**
//    * predict the state. standard KF equations. assume U==0
//  */
//  // KF Prediction step
//  x_ = F_*x_;
//  P_ = F_*P_*F_.transpose() + Q_;
//
//
//  /**
//    * update the state by using Kalman Filter equations
//  */
//  // KF Measurement update step
//  VectorXd y = z - H_*x_;
//  MatrixXd S = H_*P_*H_.transpose() + R_;
//  MatrixXd K = P_*H_.transpose()*S.inverse();
//  // new state
//  x_ = x_ + K*y;
//  MatrixXd I = MatrixXd::Identity(P_.diagonalSize(), P_.diagonalSize());
//  P_ = (I - K*H_)*P_;
//  /**
//    * update the state by using Extended Kalman Filter equations, but non-linear radar measurement h(x)
//  */
//  float rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
//  float phi = 0.0;
//  float rhorate = 0.0;
//  if (std::abs(x_(0))>1e-6 && std::abs(x_(1))>1e-6) {
//    phi = atan2(x_(1), x_(0));
//    rhorate = (x_(0) * x_(2) + x_(1) * x_(3)) / rho;
//  }
//  VectorXd hx(3);
//  hx << rho, phi, rhorate;

//  // EKF Measurement update step
//  VectorXd y = z - hx;
//  if (y(1) > M_PI)
//    y(1) = y(1) - M_PI;
//  if (y(1) < -M_PI)
//    y(1) = y(1) + M_PI;
//  // rest is same as Update() function
//  MatrixXd S = H_*P_*H_.transpose() + R_;
//  MatrixXd K = P_*H_.transpose()*S.inverse();
//  // new state
//  x_ = x_ + K*y;
//  MatrixXd I = MatrixXd::Identity(P_.diagonalSize(), P_.diagonalSize());
//  P_ = (I - K*H_)*P_;


  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
//  //set measurement covariance
//  ekf_.R_ = R_laser_;
//  //set measurement matrix
//  ekf_.H_ = H_laser_;
//
//  // Kalman Filter measurement update
//  ekf_.Update(meas_package.raw_measurements_);
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
//  //set measurement covariance
//  ekf_.R_ = R_radar_;
//  //set measurement matrix to be Jacobian of non-linear h(x)
//  ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
//
//  // Extended Kalman Filter measurement update
//  ekf_.UpdateEKF(meas_package.raw_measurements_);
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
