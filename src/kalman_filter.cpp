#include "kalman_filter.h"
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

const float pi2 = 2 * M_PI;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd y = z - (H_ * x_);
  MatrixXd PHt = P_ * H_.transpose();
  MatrixXd S = (H_ * PHt) + R_;
  MatrixXd K = PHt * S.inverse();
  //new estimate / state update
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}


VectorXd Cartesian2Polar(const VectorXd &x_) {

  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];

  float rho, phi, rho_dot;
  
  rho = sqrt(px*px + py*py);
  // avoid dividing by zero
  if (rho < 0.00001) {
    rho = 0.00001;
  }

  phi = atan2(py, px);
  rho_dot = (px*vx + py*vy) / rho;

  VectorXd x_pred(3);
  x_pred << rho, phi, rho_dot;

  return x_pred;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

  VectorXd y = z - Cartesian2Polar(x_);

  // correct error angle phi
  while(y(1) > M_PI){
    y(1) -= pi2;
  }
  while(y(1) < -M_PI){
    y(1) += pi2;
  }

  MatrixXd PHt = P_ * H_.transpose();
  MatrixXd S = (H_ * PHt) + R_;
  MatrixXd K = PHt * S.inverse();
  //new estimate / state update
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * H_) * P_;
}

