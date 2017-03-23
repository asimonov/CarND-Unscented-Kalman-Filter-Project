#include <iostream>
// for sort
#include <algorithm>
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
  lambda_ = 3.0 - n_aug_;
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
  P_prior_ = MatrixXd(n_x_, n_x_);
  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);


  // Process noise
  // standard deviation longitudinal acceleration in m/s^2
  //std_a_ = 3.0; // assuming max acceleration 6m/s2 and take half of that
  std_a_ = 0.2; // from lessons
  // standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = M_PI/2.;
  //std_yawdd_ = M_PI/20.; about .15
  std_yawdd_ = .2; // from lectures

  // Laser measurement noise
  // standard deviation position1 in m
  std_laspx_ = 0.15;
  // standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise
  // standard deviation radius in m
  std_radr_ = 0.3;
  // standard deviation angle in rad
  //std_radphi_ = M_PI/32.; // 5.6 degrees
  std_radphi_ = 0.0175; // from lectures. ~2.7 degrees
  // standard deviation radius change in m/s
  std_radrd_ = 0.1; // from lectures

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
    double phi; // will calculate from both types of measurement
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates
      double ro = meas_package.raw_measurements_(0);
      phi = meas_package.raw_measurements_(1);
      double rodot = meas_package.raw_measurements_(2);
      x_ <<   ro * cos(phi),
              ro * sin(phi),
              0.0,//rodot, // assume psi==phi, i.e. object is moving away along the axis from zero to its position
              0.0,//NormaliseAngle(phi), // assume psi==phi. why not. any direction is ok here.
              0.0; // assume the turn angle is constant
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      double px = meas_package.raw_measurements_(0);
      double py = meas_package.raw_measurements_(1);
      if (abs(px)>1e-5)
        phi = NormaliseAngle(atan2(py, px));
      else
        phi = M_PI/2.; // assume object is straight ahead
      x_ << meas_package.raw_measurements_(0),
            meas_package.raw_measurements_(1),
            0.0, // assume radial velosity is zero
            0.0,//phi, // assume object is pointing along the axis from zero to its position
            0.0; // assume the turn angle is constant
    }

    // initialize state covariance matrix
    P_ <<   0.01,                0.,              0.,                  0.,            0.,
            0.,                0.01,                  0.,            0.,            0.,
            0.,                0.,                  0.01,            0.,            0.,
            0.,                0.,                  0.,            1.0,            0.,
            0.,                0.,                  0.,            0.,            0.1;
//    P_ << std_laspx_*std_laspx_, 0.,                  0.,            0.,            0., // use radar measurement uncertainty as more imprecise one
//            0.,                std_laspy_*std_laspy_, 0.,            0.,            0., // use radar measurement uncertainty as more imprecise one
//            0.,                0.,                  0.1,            0.,            0., // assume velocity uncertainty is on the order of acceleration noise
//            0.,                0.,                  0.,            0.1,            0., // M_PI*M_PI/16. assume 45 degrees
//            0.,                0.,                  0.,            0.,            0.1; // M_PI*M_PI/16 assume 45 degrees
    P_prior_ = P_;

    time_us_ = meas_package.timestamp_;

    is_initialized_ = true;

    // print the output
    cout << "x_ = " << endl << x_ << endl;
    cout << "P_ = " << endl << P_ << endl;

    // done initializing, no need to predict or update
    return;
  }

  // skip measurement if not supposed to use this type of sensor
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && !use_radar_) {
    cout << "skipping RADAR update" << endl;
    return;
  }
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && !use_laser_) {
    cout << "skipping LIDAR update" << endl;
    return;
  }

  //compute the time elapsed between the current and previous measurements
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;     //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;
  cout << "dt = " << dt << endl;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

//  if (dt > 0.00001) {
    // predict only if processing non-simultaneous measurements
    cout << "predicting process change" << endl;

    // predict in smaller steps if measurements come less frequent than 10Hz
    // it results in better RMSE
    while (dt > 0.1)
    {
      Prediction(0.1); // Kalman Filter prediction step
      dt -= 0.1;
    }
    Prediction(dt); // Kalman Filter prediction step

    cout << "predict done" << endl;
    cout << "x_ = " << endl << x_ << endl;
    cout << "P_ = " << endl << P_ << endl;
//  }
  P_prior_ = P_;

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar update
    cout << "updating radar measurement" << endl;
    UpdateRadar(meas_package);
    cout << "radar update done" << endl;
  } else {
    // Laser update
    cout << "updating laser measurement" << endl;
    UpdateLidar(meas_package);
    cout << "laser update done" << endl;
  }

  // print the output
  cout << "x_ = " << endl << x_ << endl;
  cout << "P_ = " << endl <<  P_ << endl;
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
  cout << "Xsig_aug = " << endl << Xsig_aug << endl;

  // predict sigma points
  SigmaPointPrediction(Xsig_aug, delta_t, &Xsig_pred_); // use our instance variable to later reference it in measurement update
  cout << "Xsig_pred = " << endl << Xsig_pred_ << endl;

  // predict mean and covariance
  PredictMeanAndCovariance(&x_, &P_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  // lidar measurement dimensionality: px, py
  int n_z = 2;

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z,n_z);

  // calculate Zsig, z_pred and S
  PredictLidarMeasurement(n_z, &Zsig, &z_pred, &S);
  cout << "Zsig lidar = " << endl << Zsig << endl;
  cout << "z_pred lidar = " << endl << z_pred << endl;
  cout << "S lidar = " << endl << S << endl;

  // update state and covariance using UKF equations, predicted measurement and actual measurement
  VectorXd z = meas_package.raw_measurements_;
  cout << "z lidar = " << endl << z << endl;
  UpdateState(n_z, Zsig, z_pred, S, z);

  // calculate NIS
  VectorXd z_diff = z - z_pred;
  NIS_laser_ = (z_diff.transpose()) * (S.inverse()) * z_diff;
  cout << "NIS lidar = " << NIS_laser_ << endl;

}




/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  // radar measurement dimensionality: rho, phi, rhodot
  int n_z = 3;

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z,n_z);

  // calculate Zsig, z_pred and S.
  // reuse Xsig_pred_ from prediction step. also no augmentation needed because measurement noise has additive effect
  PredictRadarMeasurement(n_z, &Zsig, &z_pred, &S);
  cout << "Zsig radar = " << endl << Zsig << endl;
  cout << "z_pred radar = " << endl << z_pred << endl;
  cout << "S radar = " << endl << S << endl;

  // update state and covariance using UKF equations, predicted measurement and actual measurement
  VectorXd z = meas_package.raw_measurements_;
  z(1) = NormaliseAngle(z(1));
  cout << "z radar = " << endl << z << endl;
  UpdateState(n_z, Zsig, z_pred, S, z);

  // calculate NIS
  VectorXd z_diff = z - z_pred;
  z_diff(1) = NormaliseAngle(z_diff(1));
  NIS_radar_ = (z_diff.transpose()) * (S.inverse()) * z_diff;
  cout << "NIS radar = " << NIS_radar_ << endl;
}



void UKF::UpdateState(int n_z, MatrixXd &Zsig, VectorXd &z_pred, MatrixXd &S, VectorXd &z) {

  //matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++)
  {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = NormaliseAngle(x_diff(3));

    // measurement difference
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    if (n_z == 3) // awkward check if its radar measurement
      z_diff(1) = NormaliseAngle(z_diff(1));

    Tc += weights_(i)*(x_diff)*(z_diff.transpose());
  }

  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  K = Tc*(S.inverse());

  //residual
  VectorXd z_diff = z - z_pred;
  if (n_z == 3) // awkward check if its radar measurement
    z_diff(1) = NormaliseAngle(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K*z_diff;
  x_(3) = NormaliseAngle(x_(3));
  P_ = P_ - K*S*(K.transpose());

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
  //MatrixXd L = P_aug.llt().matrixL();

  // Take matrix square root
  // 1. compute the Cholesky decomposition of P_aug
  Eigen::LLT<MatrixXd> lltOfPaug(P_aug);
  if (lltOfPaug.info() == Eigen::NumericalIssue) {
    // if decomposition fails, we have numerical issues
    std::cout << "LLT failed!" << std::endl;
    Eigen::EigenSolver<MatrixXd> es(P_aug, false);
    cout << "Eigenvalues of P_aug:" << endl << es.eigenvalues() << endl;

    Eigen::EigenSolver<MatrixXd> es2(P_, false);
    cout << "Eigenvalues of P_:" << endl << es2.eigenvalues() << endl;

    throw std::range_error("LLT failed");
  }
  // 2. get the lower triangle
  MatrixXd L = lltOfPaug.matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug(3,i+1) = NormaliseAngle(Xsig_aug(3,i+1)); // normalise psi
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug(3,i+1+n_aug_) = NormaliseAngle(Xsig_aug(3,i+1+n_aug_)); // normalise psi
  }

  //write result
  *Xsig_out = Xsig_aug;
}



void UKF::SigmaPointPrediction(MatrixXd &Xsig_aug, double delta_t, MatrixXd* Xsig_out) {

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  double dtsq = delta_t*delta_t;

  //predict sigma points
  for (int i=0; i<2*n_aug_+1; i++)
  {
    VectorXd x_k = VectorXd(n_aug_);
    x_k = Xsig_aug.col(i);
    double px = x_k(0);
    double py = x_k(1);
    double v = x_k(2);
    double psi = x_k(3);
    double psidot = x_k(4);
    double nu_a = x_k(5);
    double nu_psidot = x_k(6);

    VectorXd x_k1 = VectorXd(n_x_);
    //avoid division by zero
    if (abs(psidot)>1e-3)
    {
      x_k1(0) = px + v/psidot * ( sin(psi+psidot*delta_t)  - sin(psi) ) ;
      x_k1(1) = py + v/psidot * ( -cos(psi+psidot*delta_t) + cos(psi) );
    }
    else
    {
      // no change in yaw, going on straight line
      //x_k1(0) = px + v*(-sin(psi)) + 0.5*dtsq*cos(psi)*nu_a; // my original implementation
      //x_k1(1) = py + v*(cos(psi)) + 0.5*dtsq*sin(psi)*nu_a;
      x_k1(0) = px + v * delta_t * cos(psi); // from lectures
      x_k1(1) = py + v * delta_t * sin(psi);
    }
    // add noise
    x_k1(0) = x_k1(0) + 0.5*dtsq*cos(psi)*nu_a;
    x_k1(1) = x_k1(1) + 0.5*dtsq*sin(psi)*nu_a;

    x_k1(2) = v + 0 + delta_t*nu_a;
    //x_k1(3) = NormaliseAngle(psi + psidot*delta_t + 0.5*dtsq*nu_psidot);
    //x_k1(4) = NormaliseAngle(psidot + 0 + delta_t*nu_psidot);
    x_k1(3) = psi + psidot*delta_t + 0.5*dtsq*nu_psidot;
    x_k1(4) = psidot + 0 + delta_t*nu_psidot;

    Xsig_pred.col(i) = x_k1;
  }

  //write result
  *Xsig_out = Xsig_pred;
}


void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {

  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //predict state mean
  x.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++)
  {
    x += weights_(i)*Xsig_pred_.col(i);
  }
  //x(3) = NormaliseAngle(x(3));

  //predict state covariance matrix
  P.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++)
  {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    x_diff(3) = NormaliseAngle(x_diff(3));

    P += weights_(i)*(x_diff)*(x_diff.transpose());
  }

  //write result
  *x_out = x;
  *P_out = P;
}







void UKF::PredictMeasurement(int n_z, MatrixXd &Zsig, MatrixXd &R, VectorXd* z_out, MatrixXd* S_out) {
  //calculate mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++)
  {
    z_pred += weights_(i)*Zsig.col(i);
  }

  //calculate measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++)
  {
    // state difference
    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (n_z == 3) // awkward check for radar measurement space
      z_diff(1) = NormaliseAngle(z_diff(1));

    S += weights_(i)*(z_diff)*(z_diff.transpose());
  }
  // add measurement noise
  S = S + R;

  //write result
  *z_out = z_pred;
  *S_out = S;
}


void UKF::PredictLidarMeasurement(int n_z, MatrixXd* Zsig_out, VectorXd* z_out, MatrixXd* S_out) {
  //transform sigma points into LIDAR measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  for (int i=0; i<2*n_aug_+1; i++)
  {
    VectorXd x = Xsig_pred_.col(i);
    double px = x(0);
    double py = x(1);

    // calculate rho
    Zsig(0,i) = px;
    Zsig(1,i) = py;
  }

  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0.,
          0.,                    std_laspy_*std_laspy_;

  PredictMeasurement(n_z, Zsig, R, z_out, S_out);
  *Zsig_out = Zsig;
}


void UKF::PredictRadarMeasurement(int n_z, MatrixXd* Zsig_out, VectorXd* z_out, MatrixXd* S_out) {
  //transform sigma points into RADAR measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  for (int i=0; i<2*n_aug_+1; i++)
  {
    VectorXd x = Xsig_pred_.col(i);
    double px = x(0);
    double py = x(1);
    double v = x(2);
    //double psi = NormaliseAngle(x(3));
    double psi = x(3);

    // calculate rho
    Zsig(0,i) = sqrt(px*px+py*py);

    // calculate phi
    if (abs(px)>1e-5)
      Zsig(1,i) = atan2(py, px);//NormaliseAngle(atan2(py, px));
    else
      Zsig(1,i) = M_PI/2.; // assume object is straight ahead

    // calculate rhodot
    if (abs(Zsig(0,i)) < 1e-5)
      Zsig(2,i) = 0.0; // px and py are zero. rho is zero. cannot have a reasonable change in rho, so let's set it to zero.
    else
      Zsig(2,i) = v * (px*cos(psi) + py*sin(psi)) / Zsig(0,i);
  }

  // RADAR measurement noise matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0.,                      0.,
          0.,                  std_radphi_*std_radphi_, 0.,
          0.,                  0.,                      std_radrd_*std_radrd_;

  PredictMeasurement(n_z, Zsig, R, z_out, S_out);// z_out and S_out will be assigned here

  *Zsig_out = Zsig;
}






double UKF::NormaliseAngle(double angle)
{
  double result = angle;

  // working version from lectures
  if (result>M_PI)
  {
    while (result > M_PI) result -= 2.*M_PI;
  }
  else if (result < -M_PI)
  {
    while (result < -M_PI) result += 2.*M_PI;
  }

  return result;
}

// wolfgang
//static double SNormalizeAngle(double phi)
//{
//  const double Max = M_PI;
//  const double Min = -M_PI;
//
//  return phi < Min
//         ? Max + std::fmod(phi - Min, Max - Min)
//         : std::fmod(phi - Min, Max - Min) + Min;
//}

