#include "ukf.h"
#include "Eigen/Dense"

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  // initial Lidar related Covariance matrix and Linear functoin H
  P_ = MatrixXd(5, 5);
  P_ << 0.5, 0, 0, 0, 0,
        0, 0.5, 0, 0, 0,
        0, 0, 0.5, 0, 0,
        0, 0, 0, 0.2, 0,
        0, 0, 0, 0, 0.2;

  H_ = MatrixXd(2, 5);
  H_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;


  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;       //I have chosen half of the accelaration 6 m/s^2

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.65; 
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15; 

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03; 

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;  

  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  is_initialized_ = false;
  previous_timestamp_ = 0;

  n_x_ = 5;
  n_z_ = 3;
  n_aug_ = 7;
  col_len = 2*n_aug_ +1;

  lambda_ = 3 - n_x_;
 
  R_ = MatrixXd(2, 2);
  R_ << std_laspx_*std_laspx_,               0,
        0,               std_laspy_*std_laspy_;    


  radarmeasure = VectorXd(3);
  lidarmeasure = VectorXd(2);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) 
{
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

  if(!is_initialized_) {
      //cout << "UKF Initialization " << endl;
      is_initialized_ = true;

      if(meas_package.sensor_type_ == MeasurementPackage::RADAR) 
      {
        //cout << "Initialization from Radar " << endl;

        double rho = meas_package.raw_measurements_[0];
        double psi = meas_package.raw_measurements_[1];
        double rhoDot = meas_package.raw_measurements_[2];

        x_ << rho*cos(psi),
              rho*sin(psi),
              0.4*rhoDot,
              psi,
              0;
                      
        previous_timestamp_ = meas_package.timestamp_;
    
      } 
      else if(meas_package.sensor_type_ == MeasurementPackage::LASER) 
      {
        //cout << "Initialization from Lidar " << endl;
          
        x_ << meas_package.raw_measurements_[0],
              meas_package.raw_measurements_[1],
              0,
              0,
              0;
          
        previous_timestamp_ = meas_package.timestamp_;
      }
    }

    delta_t_ = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = meas_package.timestamp_;

    Prediction(delta_t_);
   
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
       // RADAR measurement update
      UpdateRadar(meas_package);
    } else {  
      // LIDAR measurement update
      UpdateLidar(meas_package);
    }

}

void UKF::Prediction(double delta_t) 
{
  /**
   * Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  // generate sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_,col_len);
	SigmaPointsAugmentation(&Xsig_aug);

  // predict sigma points
	MatrixXd Xsig_pred = MatrixXd(n_x_, col_len);
  SigmaPointsPrediction(Xsig_aug, delta_t, &Xsig_pred);

  // predict mean and covariance of state x_
	VectorXd x_pred = VectorXd(n_x_);
  MatrixXd P_pred = MatrixXd(n_x_, n_x_);
  MeanCovariancePrediction(Xsig_pred, &x_pred, &P_pred);

  // update x_ and P_
  x_ = x_pred;
  P_ = P_pred;
}
/* 
//This part is not used in the project But written here to understand the flow

void UKF::SigmaPointsGeneration(MatrixXd* Xsig_out) 
{

  // create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, col_len);

  // create square root of P
  MatrixXd P_G = P_.llt().matrixL();

  // set first column of sigma point matrix
  Xsig.col(0) = x_;

  // set remaining sigma points
  for(int i = 0; i < n_x_; i++) {
    Xsig.col(i+1) = x_ + sqrt(lambda_ + n_x_) * P_G.col(i);
    Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_ + n_x_) * P_G.col(i);
  }

  cout << "Xsig = " << endl << Xsig << endl;
  *Xsig_out = Xsig;
}
*/


void UKF::SigmaPointsAugmentation(MatrixXd* Xsig_out) {

  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, col_len);

  // create augmented mean state
  
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_ +1) = 0;

  // create augmented covariance matrix
  double std_a_sqr = std_a_ * std_a_;
  double std_yawdd_sqr = std_yawdd_ * std_yawdd_;
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(n_x_,n_x_) = std_a_sqr;
  P_aug(n_x_+1,n_x_ +1) = std_yawdd_sqr;

  // create square root matrix
  MatrixXd P_A = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug_; i++) 
  {
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+ n_aug_) * P_A.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * P_A.col(i);
  }

  cout << "Xsig_aug = " << endl << Xsig_aug << endl;
  *Xsig_out = Xsig_aug;
}


void UKF::SigmaPointsPrediction(const MatrixXd &Xsig_aug, const double &delta_t, MatrixXd* Xsig_out) 
{
 
  // create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, col_len);

  // predict sigma points
  for(int i = 0; i < col_len; i++) {
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    double p_xp;
    double p_yp;
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;


    // avoid division by zero
    if(yawd != 0.0) 
    {
        p_xp = p_x + v/yawd * (sin(yaw_p) - sin(yaw));
        p_yp = p_y + v/yawd * (cos(yaw) - cos(yaw_p));
    } else 
    {
        p_xp = p_x + v*delta_t * cos(yaw);
        p_yp = p_y + v*delta_t * sin(yaw);
    }

    // add noise
    p_xp = p_xp + 0.5*nu_a*delta_t*delta_t*cos(yaw);
    p_yp = p_yp + 0.5*nu_a*delta_t*delta_t*sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma point into right column
    Xsig_pred(0,i) = p_xp;
    Xsig_pred(1,i) = p_yp;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  cout << "Xsig_pred = " << endl << Xsig_pred << endl;
  *Xsig_out = Xsig_pred;
}


void UKF::MeanCovariancePrediction(const MatrixXd &Xsig_pred, VectorXd* x_out, MatrixXd* P_out) 
{
  
  // create vector for weights
  VectorXd weights = VectorXd(col_len);

  // create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  // create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  // set weights
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights(0) = weight_0;
  for(int i = 1; i < col_len; i++) 
  {
    double weight = 0.5/(n_aug_+lambda_);
    weights(i) = weight;
  }

  // predicted state mean
  x.fill(0.0);
  for(int i = 0; i < col_len; i++) 
  {
    x = x + weights(i)*Xsig_pred.col(i);
  }

  // predicted state covariance matrix
  P.fill(0.0);
  for(int i = 0; i < col_len; i++) 
  {
    VectorXd x_diff = Xsig_pred.col(i) - x;
    // angle normalization

    while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
   
    P = P + weights(i) * x_diff * x_diff.transpose();
  }

  cout << "Predicted state" << endl << x << endl;
  cout << "Predicted covariance matrix" << endl << P << endl;

  *x_out = x;
  *P_out = P;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  //Measuremtns from radar Sensor
  radarmeasure <<     meas_package.raw_measurements_[0],
                      meas_package.raw_measurements_[1],
                      meas_package.raw_measurements_[2];

  // generate sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_,col_len);
	SigmaPointsAugmentation(&Xsig_aug);

  // predict sigma points
	MatrixXd Xsig_pred = MatrixXd(n_x_, col_len);
  SigmaPointsPrediction(Xsig_aug, delta_t_, &Xsig_pred);

  MatrixXd Zsig = MatrixXd(n_z_, col_len);
	VectorXd z_out = VectorXd(n_z_);
  MatrixXd S_out = MatrixXd(n_z_, n_z_);
  
  RadarMeasurementsPrediction(Xsig_pred, &Zsig, &z_out, &S_out);
  UpdateRadarState(Xsig_pred, Zsig, z_out, S_out);

}

void UKF::RadarMeasurementsPrediction(const MatrixXd &Xsig_pred, MatrixXd* zsig_out ,VectorXd* z_out, MatrixXd* S_out) 
{
  
  // set measurement dimension, radar can measure r, phi, and r_dot
  
  // set vector for weights
  VectorXd weights = VectorXd(col_len);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights(0) = weight_0;

  for (int i=1; i<col_len; ++i) 
  {  
    double weight = 0.5 / (lambda_ + n_aug_);
    weights(i) = weight;
  }

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, col_len);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_, n_z_);

  // transform sigma points into measurement space
  for(int i = 0; i < col_len; i++) {
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);

    // mesurement model
    double pxx_yy = p_x*p_x + p_y*p_y;
    Zsig(0,i) = sqrt(pxx_yy);
    Zsig(1,i) = atan2(p_y, p_x);
    Zsig(2,i) = (p_x*cos(yaw)*v + p_y*sin(yaw)*v) / sqrt(pxx_yy);  
  }

  // mean predicted measurement
  z_pred.fill(0.0);
  for(int i = 0; i < col_len; i++) 
  {
    z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  // innovation covariance matrix S
  S.fill(0.0);
  for(int i = 0; i < col_len; i++) 
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;
    
    S = S + weights(i)*z_diff*z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_, n_z_);
  R <<  std_radr_*std_radr_,     0,                 0,
        0,      std_radphi_*std_radphi_,            0,
        0,             0,       std_radrd_*std_radrd_;
  
  S = S+R;

  cout << "z_pred: " << endl << z_pred << endl;
  cout << "S: " << endl << S << endl;

  *zsig_out = Zsig;
  *z_out = z_pred;
  *S_out = S;
}

void UKF::UpdateRadarState(const MatrixXd &Xsig_pred, const MatrixXd &Zsig, const VectorXd &z_pred, const MatrixXd &S) 
{

  // set measurement dimension, radar can measure r, phi, and r_dot
  
  // set vector for weights
  VectorXd weights = VectorXd(2*n_aug_ + 1);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights(0) = weight_0;

  for (int i=1; i < col_len; i++) 
  {
    double weight = 0.5 / (lambda_ + n_aug_);  
    weights(i) = weight;
  }

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < col_len; i++) 
  { 
    VectorXd z_diff = Zsig.col(i) - z_pred;

    while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;
   
    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;

    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = radarmeasure - z_pred;

  // angle normalization

  while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;
  
  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  cout << "Updated state x: " << endl << x_ << endl;
  cout << "Updated state covariance P: " << endl << P_ << endl;

  double NIS = z_diff.transpose( ) * S.inverse() * z_diff;
  cout <<"RADAR NIS with 3 dimensions: "<<NIS<< endl;
  cout<<NIS<<" : 95% is less than 7.185"<<endl;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) 
{
  /**
   * Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  lidarmeasure <<      meas_package.raw_measurements_[0],
                       meas_package.raw_measurements_[1];

  
  VectorXd z_pred = H_ * x_;
  VectorXd y = lidarmeasure - z_pred;

  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  x_ = x_ + K * y;

  while (x_(3) > M_PI) x_(3) -= 2.*M_PI;
  while (x_(3) < -M_PI) x_(3) += 2.*M_PI;
   
  size_t x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);

  P_ = (I - K * H_) * P_;

  double NIS = y.transpose( ) * S.inverse() * y;
  cout <<"LIDAR NIS with 2 dimensions: "<<NIS<< endl;
  cout<<NIS<<" :95% is less than 5.991"<<endl;
}

