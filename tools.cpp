#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if (estimations.size() == 0 ||
    estimations.size() != ground_truth.size()) {
    std::cout << "The estimations data is incorrect!" << std::endl;
    return rmse;
  }

  for (int i = 0; i < estimations.size(); i++) {
    VectorXd residuals;
    residuals = estimations[i] - ground_truth[i];
    residuals = residuals.array() * residuals.array();
    rmse += residuals;
  }

  rmse = rmse / estimations.size();

  rmse = rmse.array().sqrt();

  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */

  MatrixXd Hj(3, 4);
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  const double eps = 1e-3;
  const double c0 = std::max(eps, px * px + py * py);
  const double c1 = sqrt(c0);
  const double c2 = c0 * c1;

  if (fabs(px) < 1e-4 || fabs(py) < 1e-4) {
    cout << "CalculateJacobian Error, division by zero!" << endl;

    return Hj;
  }

  Hj << (px / c1), (py / c1), 0, 0,
        -(py / c0), (px / c0), 0, 0,
        py * (vx * py - vy * px) / c2, px * (vy * px - vx * py) / c2, px / c1, py / c1;

  return Hj;
}
