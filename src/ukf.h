#ifndef UKF_H
#define UKF_H

#include <iostream>
#include "Eigen/Dense"
#include "measurement_package.h"
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  UKF();
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

  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // LiDAR measurement matrix
  Eigen::MatrixXd H_;

  // LiDAR measurement covariance
  Eigen::MatrixXd R_;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

  // time when the state is true, in us
  long long time_us_;
  double delta_t_;
  size_t previous_timestamp_;

  //Two vectors to store Sensor data
  VectorXd radarmeasure;
  VectorXd lidarmeasure;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  // Weights of sigma points
  Eigen::VectorXd weights_;

  // State dimension
  int n_x_;
  int n_z_;

  // Augmented state dimension
  int n_aug_;
  // Define Column lenght to simplify the equations
  int col_len;

  // Sigma point spreading parameter
  double lambda_;

  //Created these functions to access privately
private:
  void SigmaPointsGeneration(Eigen::MatrixXd* Xsig_out);
  void SigmaPointsAugmentation(Eigen::MatrixXd* Xsig_out);
  void SigmaPointsPrediction(const Eigen::MatrixXd &Xsig_aug, const double &delta_t, Eigen::MatrixXd* Xsig_out);
  void MeanCovariancePrediction(const Eigen::MatrixXd &Xsig_pred, Eigen::VectorXd* x_pred, Eigen::MatrixXd* P_pred);
  void RadarMeasurementsPrediction(const Eigen::MatrixXd &Xsig_pred,Eigen::MatrixXd* zsig_out, Eigen::VectorXd* z_out, Eigen::MatrixXd* S_out);
  void UpdateRadarState(const Eigen::MatrixXd &Xsig_pred, const Eigen::MatrixXd &Zsig, const Eigen::VectorXd &z_pred, const Eigen::MatrixXd &S);

};

#endif  // UKF_H
