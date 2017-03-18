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

  // calculate augmented sigma points vectors
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  AugmentedSigmaPoints(&Xsig_aug);
  cout << "Xsig_aug = " << Xsig_aug << endl;

  // predict sigma points
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  SigmaPointPrediction(Xsig_aug, delta_t, &Xsig_pred);
  cout << "Xsig_pred = " << Xsig_pred << endl;

  // predict mean and covariance
  PredictMeanAndCovariance(Xsig_pred, &x_, &P_);
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













void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0.;
  x_aug(6) = 0.;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  //write result
  *Xsig_out = Xsig_aug;
}



void UKF::SigmaPointPrediction(MatrixXd &Xsig_aug, double delta_t, MatrixXd* Xsig_out) {

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  float dtsq = delta_t*delta_t;

  //predict sigma points
  for (int i=0; i<2*n_aug_+1; i++)
  {
    VectorXd x_k = VectorXd(n_aug_);
    x_k = Xsig_aug.col(i);
    float px = x_k(0);
    float py = x_k(1);
    float v = x_k(2);
    float psi = x_k(3);
    float psidot = x_k(4);
    float nu_a = x_k(5);
    float nu_psidot = x_k(6);

    VectorXd x_k1 = VectorXd(n_x_);
    //avoid division by zero
    if (abs(psidot)>1e-5)
    {
      x_k1(0) = px + v/psidot*(sin(psi+psidot*delta_t)-sin(psi)) + 0.5*dtsq*cos(psi)*nu_a;
      x_k1(1) = py + v/psidot*(-cos(psi+psidot*delta_t)+cos(psi)) + 0.5*dtsq*sin(psi)*nu_a;
    }
    else
    {
      // no change in yaw, going on straight line
      x_k1(0) = px + v*(-sin(psi)) + 0.5*dtsq*cos(psi)*nu_a;
      x_k1(1) = py + v*(cos(psi)) + 0.5*dtsq*sin(psi)*nu_a;
    }
    x_k1(2) = v + 0 + delta_t*nu_a;
    x_k1(3) = psi + psidot*delta_t + 0.5*dtsq*nu_psidot;
    x_k1(4) = psidot + 0 + delta_t*nu_psidot;

    Xsig_pred.col(i) = x_k1;
  }

  //write result
  *Xsig_out = Xsig_pred;
}


void UKF::PredictMeanAndCovariance(MatrixXd &Xsig_pred, VectorXd* x_out, MatrixXd* P_out) {

  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predict state mean
  x = weights_(0)*Xsig_pred.col(0);
  for (int i=1; i<2*n_aug_+1; i++)
  {
    x += weights_(i)*Xsig_pred.col(i);
  }

  //predict state covariance matrix
  P.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++)
  {
    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P += weights_(i)*(x_diff)*(x_diff.transpose());
  }

  //write result
  *x_out = x;
  *P_out = P;
}


