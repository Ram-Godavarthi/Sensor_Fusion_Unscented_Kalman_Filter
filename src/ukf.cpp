#include "ukf.h"
#include "Eigen/Dense"
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

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

  // initial covariance matrix
  //P_ = MatrixXd::Identity(5, 5);
  
  P_ = MatrixXd(5,5);

  
  P_ <<   0.5, 0,  0,   0,    0,
          0,  0.5, 0,   0,    0,
          0,   0, 0.5,  0,    0,
          0,   0,  0,  0.2,   0,
          0,   0,  0,   0,    0.2;
  

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5; //initially 30, but chnaged to 3.0 by considering half of max acceleration

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.82;
  
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

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  //To store measurement values from sensors
  radarMeasuredInput_ = VectorXd(3); //rho, psi, rhoDot
  lidarMeasuredInput_ = VectorXd(2); //px, py

  is_initialized_ = false;

}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) 
{
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (!is_initialized_)
  {
    is_initialized_ = true;
    if(meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER)
    {
      //x_.fill(0.0);
      x_ << meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1],
            0,
            0,
            0;

      time_us_ = meas_package.timestamp_;
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR)
    {
      double rho = meas_package.raw_measurements_[0];
      double psi = meas_package.raw_measurements_[1];
      double rhoDot = meas_package.raw_measurements_[2];
      //x_.fill(0.0);
      x_ << rho * cos(psi),
            rho * sin(psi),
            0.4 * rhoDot,
            psi,
            0;
      time_us_ = meas_package.timestamp_;
    }

  }

  double delta_t = (meas_package.timestamp_ - time_us_ ) /1000000.0 ;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  if(meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR)
  {
    UpdateRadar(meas_package);
  }
  else
  {
    UpdateLidar(meas_package);
  }
  

}


void UKF::GenerateSigmaPoints(Eigen::MatrixXd* Xsig_out) {
  
  //int n_x = 5;

  // create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  // create square root of P
  MatrixXd A = P_.llt().matrixL();

  // set first column of sigma point matrix
  Xsig.col(0) = x_;

  // set remaining sigma points
  for(int i = 0; i < n_x_; i++) {
    Xsig.col(i+1) = x_ + sqrt(lambda_ + n_x_) * A.col(i);
    Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_ + n_x_) * A.col(i);
  }

  std::cout << "Xsig = " << std::endl << Xsig << std::endl;
  *Xsig_out = Xsig;
}

void UKF::AugmentedSigmaPoints (MatrixXd* Xsig_out)
{
  //Create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  //Create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  //Create sigma point matrix
  MatrixXd Xsig_aug  =MatrixXd(n_aug_ , 2 * n_aug_ +1);

  // This is the code that I have implemented in my classroom exercise
   
   //Create augmented mean state
   x_aug.head(5) = x_;
   x_aug(5) = 0;
   x_aug(6) = 0;
   /*
   x_aug << x_(0),
            x_(1),
            x_(2),
            x_(3),
            x_(4),
             0,
             0;
   */
   //Create augmented covariance matrix
   double sigma_square = std_a_ * std_a_ ;
   double yawdd_square = std_yawdd_ * std_yawdd_ ;
   P_aug.fill(0.0);
   P_aug.topLeftCorner(5,5) = P_;
   P_aug(5,5) =  sigma_square;
   P_aug(6,6) =  yawdd_square;

   //Create square root matrix
   MatrixXd P_sqrt = P_aug.llt().matrixL();

   //Create augmented sigma points
    Xsig_aug.col(0) = x_aug;

    for (int i = 0 ; i < n_aug_ ; i++)
    {
      Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * P_sqrt.col(i);
      Xsig_aug.col(i+1+n_aug_) = x_aug + sqrt(lambda_ + n_aug_) * P_sqrt.col(i);
    }

    std::cout <<"Xsig_aug = " <<std::endl << Xsig_aug <<std::endl;
    *Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(const MatrixXd &Xsig_aug, MatrixXd* Xsig_out, const double &delta_t)
{
  // create matrix with predicted sigma points as columns
  //MatrixXd Xsig_out = 
  MatrixXd Xsig_pred  = MatrixXd(n_x_, 2*n_aug_+1);

  // predict sigma points
  for (int i = 0 ; i < 2 * n_aug_ + 1; i ++)
  {
    // extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);
    
    
    // predicted state values
    double p_xp, p_yp;
    
    // avoid division by zero
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;
    if (yawd > 0.001)
    {
        p_xp = p_x + v/yawd * (sin(yaw_p) - sin(yaw));
        p_yp = p_y + v/yawd * (cos(yaw) - cos(yaw_p));
    } 
    else
    {
        p_xp = p_x + v*delta_t*cos(yaw);
        p_yp = p_y + v*delta_t*sin(yaw);
    }
    // add noise
    p_xp = p_xp + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    p_yp = p_yp + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;
    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma points into right column
  
    Xsig_pred(0,i) = p_xp;
    Xsig_pred(1,i) = p_yp;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }
  std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;
  // write result
  *Xsig_out = Xsig_pred;
  
}

void UKF::PredictMeanAndCovariance(const MatrixXd &Xsig_pred, VectorXd* x_pred, MatrixXd* P_pred)
{
  // create vector for weights
  int col_len = 2* n_aug_+1;
  VectorXd weights = VectorXd(col_len);
  
  // create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  // create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  // set weights
  double w_0 = lambda_/(lambda_+n_aug_);
  weights(0) = w_0;
  for (int i=1; i<col_len; ++i) 
  {  // 2n+1 weights
    double w = 0.5/(n_aug_+lambda_);
    weights(i) = w;
  }

  // predict state mean
  x.fill(0.0);
  for (int i = 0; i < col_len; ++i) 
  {  // iterate over sigma points
    x = x + weights(i) * Xsig_pred.col(i);
  }

  // predict state covariance matrix
  P.fill(0.0);
  
  for (int i = 0; i < col_len; ++i) 
  {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights(i) * x_diff * x_diff.transpose() ;
  }
  // print result
  std::cout << "Predicted state" << std::endl << x << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl << P << std::endl;
  // write result
  *x_pred = x;
  *P_pred = P;
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  //Augmented Sigma Points Generation
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_ + 1);
  AugmentedSigmaPoints (&Xsig_aug);

  //Sigma Points Prediction
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2*n_aug_+1);
  SigmaPointPrediction (Xsig_aug, &Xsig_pred, delta_t);

  /*
  //Mean State & State Covariance Matrix Predition
  VectorXd x_pred = VectorXd(n_x_);
  MatrixXd P_pred = MatrixXd(n_x_, n_x_);
  PredictMeanAndCovariance (Xsig_pred, &x_pred, &P_pred);
  x_ = x_pred;
  P_ = P_pred;
  */

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  lidarMeasuredInput_ <<  meas_package.raw_measurements_[0],
                          meas_package.raw_measurements_[1];
  int lm_size = lidarMeasuredInput_.size();
  MatrixXd H = MatrixXd(lm_size, n_x_);
  H << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0;

  MatrixXd R = MatrixXd(lm_size,lm_size);
  double lpx_sqr = std_laspx_ * std_laspx_;
  double lpy_sqr = std_laspy_ *std_laspy_;

  R <<    lpx_sqr,    0,
            0,      lpy_sqr;

  VectorXd z_pred = H * x_;
  VectorXd y = lidarMeasuredInput_ - z_pred;
  MatrixXd S = H * P_ * H.transpose() + R;
  MatrixXd K = P_ * H.transpose() * S.inverse();

  x_ = x_ + K * y;

  while (x_(3)> M_PI) x_(3)-=2.*M_PI;
  while (x_(3)<-M_PI) x_(3)+=2.*M_PI;
  
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  P_ = (I - K * H) * P_;
  
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  const int col_len = 2* n_aug_+1;
  MatrixXd Xsig_aug =  MatrixXd(n_aug_, col_len);
  AugmentedSigmaPoints(&Xsig_aug);

  MatrixXd Xsig_pred = MatrixXd(n_x_, col_len);
  SigmaPointPrediction(Xsig_aug, &Xsig_pred, delta_t);

  //set measurement dimension, radar can measure rho, psi, rhoDot
  const int n_z = 3;
  VectorXd z_out = VectorXd(n_z);
  MatrixXd S_out = MatrixXd(n_z,n_z);
  MatrixXd Zsig = MatrixXd(n_z, col_len);

  PredictRadarMeasurement(n_z,Xsig_pred, &Zsig, &z_out, &S_out);
  UpdateState(Xsig_pred, Zsig, z_out, S_out);

}
void UKF::PredictRadarMeasurement(const int n_z, const MatrixXd &Xsig_pred, MatrixXd* Zsig_out, VectorXd* z_out, MatrixXd* S_out)
{
  // set vector for weights
  const int col_len = 2*n_aug_+1;
  VectorXd weights = VectorXd(col_len);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;

  for (int i=1; i<col_len; ++i) 
  {  
    double weight = 0.5/(lambda_+n_aug_);
    weights(i) = weight;
  }
  
  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, col_len);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  // transform sigma points into measurement space
  for (int i = 0; i < col_len; ++i) {  // 2n+1 simga points
    // extract values for better readability
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v  = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);

    // measurement model
    double p_xx_yy = p_x *p_x + p_y * p_y;
    Zsig(0,i) = sqrt(p_xx_yy);                       // r
    Zsig(1,i) = atan2(p_y,p_x);                                // phi
    Zsig(2,i) = (p_x*cos(yaw)*v + p_y*sin(yaw)*v) / sqrt(p_xx_yy);   // r_dot
  }

  // mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < col_len; ++i) {
    z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  // innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < col_len; ++i) 
  {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);

  R <<  std_radr_*std_radr_,              0,                 0,
                   0,        std_radphi_*std_radphi_,        0,
                   0,                     0,           std_radrd_*std_radrd_;

  S = S + R;

  // print result
  std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  std::cout << "S: " << std::endl << S << std::endl;

  // write result
  *Zsig_out = Zsig;
  *z_out = z_pred;
  *S_out = S;
}

void UKF::UpdateState(MatrixXd &Xsig_pred, const MatrixXd &Zsig, const VectorXd &z_pred, const MatrixXd &S)
{
  int n_z = 3;
  const int col_len = 2*n_aug_+1;
  VectorXd weights = VectorXd(col_len);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;

  for (int i=1; i< col_len; ++i) 
  {
    double weight = 0.5/(lambda_+n_aug_);
    weights(i) = weight;
  }
   // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0 ; i< col_len; i++)
  {
      VectorXd x_diff = Xsig_pred.col(i) - x_;
      VectorXd z_diff = Zsig.col(i) - z_pred;
      
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
      
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
      
      Tc = Tc + weights(i) *x_diff *z_diff.transpose();
  }

  // calculate Kalman gain K;
   MatrixXd K = Tc * S.inverse();

  // update state mean and covariance matrix
  //residual
  VectorXd z_diff = radarMeasuredInput_ - z_pred;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  // print result
  std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

}