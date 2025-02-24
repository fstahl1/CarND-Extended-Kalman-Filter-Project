#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // initialize state transition matrix
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1, 
             0, 0, 1, 0,
             0, 0, 0, 1;

  // set measurement noises
  noise_ax = 9;
  noise_ay = 9;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */
    
    previous_timestamp_ = measurement_pack.timestamp_;

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];

      // convert polar to cartesian coordinates
      ekf_.x_[0] = rho * cos(phi);
      ekf_.x_[1] = rho * sin(phi);
      ekf_.x_[2] = rho_dot * cos(phi);
      ekf_.x_[3] = rho_dot * sin(phi);
    
      /* initialize state covariance matrix: position is certain since based on measurement value, 
      velocity is not that certain, since radar only measures the velocity projected on the line between vehicle and object */
      ekf_.P_ = MatrixXd(4, 4);
      ekf_.P_ << 1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 10, 0,
                 0, 0, 0, 10;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      ekf_.x_[0] = measurement_pack.raw_measurements_[0];
      ekf_.x_[1] = measurement_pack.raw_measurements_[1];

      /* initialize state covariance matrix: position is certain since based on measurement value, 
      velocity is not certain at all, since not directly measured by laser */
      ekf_.P_ = MatrixXd(4, 4);
      ekf_.P_ << 1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 10000, 0,
                 0, 0, 0, 10000;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  // current timestep in microseconds
  dt = measurement_pack.timestamp_ - previous_timestamp_ ;
  // convert to seconds
  dt = dt / 1000000.0;

  // set new previous timestamp
  previous_timestamp_ = measurement_pack.timestamp_;

  // update delta_t in state transition matrix
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  dt2 = dt*dt;
  dt3 = dt2*dt;
  dt4 = dt3*dt;
  ekf_.Q_ = MatrixXd(4,4);
  // update process covariance matrix
  ekf_.Q_ << dt4/4*noise_ax, 0             , dt3/2*noise_ax, 0             ,
             0             , dt4/4*noise_ay, 0             , dt3/2*noise_ay,
             dt3/2*noise_ax, 0             , dt2*noise_ax  , 0             ,
             0             , dt3/2*noise_ay, 0             , dt2*noise_ay;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } 
  else {
    // TODO: Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
