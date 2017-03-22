#ifndef UKF_H
#define UKF_H
#include "Eigen/Dense"
#include "measurement_package.h"
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;
  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;
  ///* state covariance matrix
  MatrixXd P_;
  MatrixXd P_prior_;// state after predict, for debug purposes
  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long time_us_;

  ///* Process noise
  ///* standard deviation longitudinal acceleration in m/s^2
  double std_a_;
  ///* standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise
  ///* standard deviation position1 in m
  double std_laspx_;
  ///* standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise
  ///* standard deviation radius in m
  double std_radr_;
  ///* standard deviation angle in rad
  double std_radphi_;
  ///* standard deviation radius change in m/s
  double std_radrd_;

  ///* Weights of sigma points
  VectorXd weights_;
  ///* State dimension
  int n_x_;
  ///* Augmented state dimension
  int n_aug_;
  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  double NIS_radar_;
  ///* the current NIS for laser
  double NIS_laser_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  /**
   * Creates augmented sigma points matrix
   * @param Xsig_out pointer to sigma points matrix where the result is written
   */
  void AugmentedSigmaPoints(MatrixXd* Xsig_out);

  /**
   * Transforms augmented sigma points using process equations
   * @param Xsig_aug augmented sigma points matrix by reference
   * @param delta_t time in seconds
   * @param Xsig_out pointer to predicted sigma points matrix where the result is written
   */
  void SigmaPointPrediction(MatrixXd &Xsig_aug, double delta_t, MatrixXd* Xsig_out);

  /**
   * Predict mean and covariance from predicted sigma points
   * @param x_out pointer to write predicted mean state to
   * @param P_out pointer to write predicted state covariance to
   */
  void PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out);

  /**
   * Predict RADAR measurement mean z_pred and covariance S from predicted sigma points
   * @param n_z   dimensions of RADAR measurements
   * @param Zsig_out pointer to write predicted sigma points in measurement space to
   * @param z_out pointer to write predicted measurement mean to
   * @param S_out pointer to write predicted measurement covariance to
   */
  void PredictRadarMeasurement(int n_z, MatrixXd* Zsig_out, VectorXd* z_out, MatrixXd* S_out);
  /**
   * Predict LIDAR measurement mean z_pred and covariance S from predicted sigma points
   * @param n_z   dimensions of LIDAR measurements
   * @param Zsig_out pointer to write predicted sigma points in measurement space to
   * @param z_out pointer to write predicted measurement mean to
   * @param S_out pointer to write predicted measurement covariance to
   */
  void PredictLidarMeasurement(int n_z, MatrixXd* Zsig_out, VectorXd* z_out, MatrixXd* S_out);

  /**
   * Common part of predicting LIDAR and RADAR measurement mean z_pred and covariance S from predicted sigma points
   * @param n_z   dimensions of LIDAR/RADAR measurements
   * @param Zsig_out pointer to write predicted sigma points in measurement space to
   * @param z_out pointer to write predicted measurement mean to
   * @param S_out pointer to write predicted measurement covariance to
   */
  void PredictMeasurement(int n_z, MatrixXd &Zsig, MatrixXd &R, VectorXd* z_out, MatrixXd* S_out);

  /**
   * Update state and covariance from measurement.
   * @param n_z    dimensions of measurement space
   * @param Zsig   matrix of sigma points in measurement space
   * @param z_pred predicted measurment from sigma points
   * @param S      measurement covariance
   * @param z      actual measurment
   */
  void UpdateState(int n_z, MatrixXd &Zsig, VectorXd &z_pred, MatrixXd &S, VectorXd &z);

  /**
   * Normalize angle to between -PI..PI
   * @param angle  input angle to be normalised
   */
  double NormaliseAngle(double angle);
//  inline double NormaliseAngle( double angle )
//  {
//    //return angle;
//
//    //double twoPi = 2.0 * M_PI;
//    //return angle - twoPi * floor( angle / twoPi );
//
//    return angle - int(angle / M_PI) * M_PI;
//  }
};

#endif /* UKF_H */
