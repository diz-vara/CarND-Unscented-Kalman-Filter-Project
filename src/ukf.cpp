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
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
 

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
  Xsig_pred_ = MatrixXd(n_aug_, sigma_points_nr_);

}

UKF::~UKF() {}

void UKF::Angle_norm(double &angle) const 
{
	while (angle > M_PI) angle -= 2.*M_PI;
	while (angle < -M_PI) angle += 2.*M_PI;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
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

  VectorXd Nu = VectorXd(n_x_);
  VectorXd A = VectorXd(n_x_);
  VectorXd x_aug = VectorXd::Zero(n_aug_);
  x_aug.head(n_x_) = x_;

  //create augmented covariance matrix
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);

  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_*std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_*std_yawdd_;


  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();



  double half_dt_square = delta_t * delta_t * 0.5;
  VectorXd x_a = VectorXd::Zero(n_aug_);;

  //Sigma points prediction 
  for (int col = 0; col < sigma_points_nr_; ++col) {
	  //for col = 0							Xsig_aug.col(0) = x_aug;
	  //for col = 1..n_aug + 1				Xsig_aug.col(col) = x_aug + (sqrt_lambda_plus_nx * A.col(col-1));
	  //for col = n_aug+1 .. 2*n_aug + 1	Xsig_aug.col(col) = x_aug - (sqrt_lambda_plus_nx * A.col(col-n_aug-1));
	  if (col == 0)
		  x_a = x_aug;
	  else if (col <= n_aug_)
		  x_a = x_aug + (sqrt_lambda_plus_nx_ * A.col(col - 1));
	  else
		  x_a = x_aug - (sqrt_lambda_plus_nx_ * A.col(col - n_aug_ - 1));


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
		  A <<
			  v_by_yaw_dot * (sin(new_yaw) - sin(yaw)),
			  v_by_yaw_dot * (cos(yaw) - cos(new_yaw)),
			  0,
			  yaw_dot*delta_t,
			  0;
	  }
	  else {
		  //yaw_dot = 0;
		  A <<
			  v * cos(yaw) * delta_t,
			  v * sin(yaw) * delta_t,
			  0, 0, 0;
	  }
	  Xsig_pred_.col(col) = x.head(n_x_) + A + Nu;
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
}

/**
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

	//update state mean and covariance matrix
	z_diff = z - z_pred;
	Angle_norm(z_diff(1));

	x_ = x_ + K * z_diff;
	P_ = P_ - K * S * K.transpose();

}
