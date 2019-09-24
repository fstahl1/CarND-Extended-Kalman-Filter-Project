#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */

  VectorXd rmse_sum(4);
  rmse_sum << 0, 0, 0, 0;
  for (int i = 0; i < estimations.size(); ++i)
  {
    VectorXd residual = estimations[i] - ground_truth[i];
    VectorXd sqr = residual.array() * residual.array();
    rmse_sum += sqr;
  }

  VectorXd rmse_mean = rmse_sum / estimations.size();
  VectorXd rmse = rmse_mean.array().sqrt();

  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float c1 = px*px+py*py;
  // avoid division by zero
  if (fabs(c1) < 0.0001) {
    c1 = 0.0001;
  }
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  MatrixXd jacobian_matrix(3, 4);
  jacobian_matrix << (px/c2)              , (py/c2)              ,     0,     0,
                     -(py/c1)             , (px/c1)              ,     0,     0,
                     py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
  return jacobian_matrix;
}
