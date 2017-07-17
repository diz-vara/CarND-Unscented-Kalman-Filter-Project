#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd::Zero(5);

  // initial covariance matrix
  P_ = MatrixXd::Zero(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.7;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.7;

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


  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = n_x_ + 2;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  sqrt_lambda_plus_nx_ = sqrt(lambda_ + n_aug_);

  sigma_points_nr_ = 2 * n_aug_ + 1;

  ///* Weights of sigma points

  weights_ = VectorXd(sigma_points_nr_);

  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < sigma_points_nr_; ++i) {
	  weights_(i) = 0.5 / (lambda_ + n_aug_);
  }

  x_ = VectorXd::Zero(n_x_);
  Xsig_pred_ = MatrixXd(n_x_, sigma_points_nr_);

  timestamp_ = 0;

  R_r_ = MatrixXd(3, 3);
  R_r_ << std_radrd_ * std_radrd_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0, std_radrd_ * std_radrd_;

  R_l_ = MatrixXd(2, 2);
  R_l_ << std_laspx_* std_laspx_, 0,
    0, std_laspy_ * std_laspy_;


}

UKF::~UKF() {}

void UKF::Angle_norm(double &angle) const 
{
  //std::cout << "angle: " << angle << " -> ";
	while (angle > M_PI) angle -= 2.*M_PI;
	while (angle < -M_PI) angle += 2.*M_PI;
  //std::cout << angle << std::endl;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:
  -- initialize
  -- calc deltaT
  -- call Prediction
  */
  static int cnt(0);

  std::cout << "Cnt:" << cnt++ << std::endl;

  if (timestamp_ != 0) {
    Prediction( static_cast<double>(meas_package.timestamp_ - timestamp_)/1.e6);
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
      UpdateLidar(meas_package);
    else
      UpdateRadar(meas_package);
  }
  else {
    x_.fill(0);
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }
    else {
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rho_d = meas_package.raw_measurements_(2);

      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
      x_(2) = rho_d;
    }

  }
  timestamp_ = meas_package.timestamp_;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  if (delta_t == 0)
    return;

  VectorXd Nu = VectorXd(n_x_);
  VectorXd x_aug = VectorXd::Zero(n_aug_);
  x_aug.head(n_x_) = x_;

  //create augmented covariance matrix
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);

  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_*std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;


  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  double half_dt_square = delta_t * delta_t * 0.5;
  VectorXd x_a = VectorXd::Zero(n_aug_);;
  VectorXd B = VectorXd::Zero(n_x_);

  //Sigma points prediction 
  for (int col = 0; col < sigma_points_nr_; ++col) {
	  //for col = 0							Xsig_aug.col(0) = x_aug;
	  //for col = 1..n_aug + 1				Xsig_aug.col(col) = x_aug + (sqrt_lambda_plus_nx * A.col(col-1));
	  //for col = n_aug+1 .. 2*n_aug + 1	Xsig_aug.col(col) = x_aug - (sqrt_lambda_plus_nx * A.col(col-n_aug-1));
	  if (col == 0)
		  x_a = x_aug;
	  else if (col <= n_aug_)
		  x_a = x_aug + (sqrt_lambda_plus_nx_ * L.col(col - 1));
	  else
		  x_a = x_aug - (sqrt_lambda_plus_nx_ * L.col(col - n_aug_ - 1));

    //std::cout << col << " : " << x_a << ", ";

	  double yaw = x_a(3);
	  double v = x_a(2);
	  double yaw_dot = x_a(4);

	  double nu_a = x_a(5);
	  double nu_yaw = x_a(6);


	  Nu <<
		  half_dt_square * cos(yaw) *nu_a,
		  half_dt_square * sin(yaw) *nu_a,
		  delta_t * nu_a,
		  half_dt_square * nu_yaw,
		  delta_t * nu_yaw;

	  if (yaw_dot != 0) {
		  double v_by_yaw_dot = v / yaw_dot;
		  double new_yaw = yaw + yaw_dot * delta_t;
		  B <<
			  v_by_yaw_dot * (sin(new_yaw) - sin(yaw)),
			  v_by_yaw_dot * (cos(yaw) - cos(new_yaw)),
			  0,
			  yaw_dot*delta_t,
			  0;
	  }
	  else {
		  //yaw_dot = 0;
		  B <<
			  v * cos(yaw) * delta_t,
			  v * sin(yaw) * delta_t,
			  0, 0, 0;
	  }
	  Xsig_pred_.col(col) = x_a.head(n_x_) + B + Nu;
    //std::cout << "S:" << Xsig_pred_.col(col) << std::endl;
  } //end of sigma point prediction loop


  //predict mean
  x_.fill(0);
  for (int col = 0; col < sigma_points_nr_; ++col) {  //2n+1 weights
	  x_ = x_ + weights_(col) * Xsig_pred_.col(col);
  }

  //predict covariance
  for (int col = 0; col < sigma_points_nr_; ++col) {  //2n+1 weights
	   // state difference
	  VectorXd x_diff = Xsig_pred_.col(col) - x_;
	  //angle normalization
	  Angle_norm(x_diff(3));
	  P_ = P_ + weights_(col) * x_diff * x_diff.transpose();
  }

  std::cout << "P " << x_(0) << ", " << x_(1) << std::endl;

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int n_z = 2;
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff;


  VectorXd z_pred = VectorXd::Zero(n_z);
  VectorXd s_point;
  for (int point = 0; point < sigma_points_nr_; ++point) {
    s_point = Xsig_pred_.col(point).head(n_z);
    z_pred = z_pred + weights_(point) * s_point;
  }
  std::cout << "Zl_pred: " << z_pred(0) << ", " << z_pred(1) << ", " << std::endl;


  //calc covariance matirx 
  MatrixXd S = MatrixXd::Zero(n_z, n_z);
  for (int point = 0; point < sigma_points_nr_; ++point) {
    z_diff = Xsig_pred_.col(point).head(n_z) - z_pred;

    S = S + weights_(point) * z_diff * z_diff.transpose();
  }



  //add measurement noise covariance
  S = S + R_l_;

  //calculate cross correlation matrix
  for (int point = 0; point < sigma_points_nr_; ++point) {
    VectorXd x_diff = Xsig_pred_.col(point) - x_;

    z_diff = Xsig_pred_.col(point).head(n_z) - z_pred;

    Tc = Tc + weights_(point) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  z_diff = z - z_pred;
  Angle_norm(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  std::cout << "Ul " << x_(0) << ", " << x_(1) << std::endl;
  P_ = P_ - K * S * K.transpose();
}

/**
* Transofrmation from process space to the Radar measurement space
* @param in - process state vector
* @param out - radar measurement state vector
*/
VectorXd UKF::H(VectorXd in) const
{
  double p_x = in(0);
  double p_y = in(1);
  double v = in(2);
  double yaw = in(3);

  double v1 = cos(yaw)*v;
  double v2 = sin(yaw)*v;

  VectorXd out = VectorXd::Zero(3);

  double norm2 = sqrt(p_x*p_x + p_y*p_y);
  out(0) = norm2;                           //r
  out(1) = atan2(p_y, p_x);                 //phi
  if (norm2 < 1e-9)
    out(2) = 0;
  else
    out(2) = (p_x*v1 + p_y * v2) / norm2;   //r_dot
  return out;
}
/*
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
	int n_z = 3;
	MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
	VectorXd z = meas_package.raw_measurements_;
	VectorXd z_diff;
  std::cout << "zr     : " << z(0) << ", " << z(1) << ", " << z(2) << std::endl;
  //std::cout << H(x_) << std::endl;


  //sigma-points in radar measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, sigma_points_nr_);
  //mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  VectorXd s_point;
  for (int point = 0; point < sigma_points_nr_; ++point) {
    s_point = H(Xsig_pred_.col(point));
    //std::cout << "p " << point << ", w:" << weights_(point) << ", h(): " << s_point(0) << ", " << s_point(1) << ", " << s_point(2) << std::endl;
    //std::cout << Xsig_pred_.col(point) << std::endl;

    Zsig.col(point) = s_point;
    z_pred = z_pred + weights_(point) * s_point;
  }
  std::cout << "Zr_pred: " << z_pred(0) << ", " << z_pred(1) << ", " << z_pred(2) << std::endl;


  //calc covariance matirx 
  MatrixXd S = MatrixXd::Zero(n_z, n_z);
  for (int point = 0; point < sigma_points_nr_; ++point) {
    z_diff = Zsig.col(point) - z_pred;
    Angle_norm(z_diff(1));

    S = S + weights_(point) * z_diff * z_diff.transpose();
  }



  //add measurement noise covariance
  S = S + R_r_;

  //calculate cross correlation matrix
	for (int point = 0; point < sigma_points_nr_; ++point) {
		VectorXd x_diff = Xsig_pred_.col(point) - x_;
		Angle_norm(x_diff(3));

		z_diff = Zsig.col(point) - z_pred;
		Angle_norm(z_diff(1));

		Tc = Tc + weights_(point) * x_diff * z_diff.transpose();
	}

  //calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();

  z_diff = z - z_pred;
  Angle_norm(z_diff(1));

	//update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  std::cout << "Ur " << x_(0) << ", " << x_(1) << std::endl;
	P_ = P_ - K * S * K.transpose();
}
