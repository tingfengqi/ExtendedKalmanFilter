#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
  TODO:
    * predict the state
  */

  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  UpdateCommon(y);
//  MatrixXd Ht = H_.transpose();
//  MatrixXd S = H_ * P_ * Ht + R_;
//  MatrixXd Si = S.inverse();
//  MatrixXd K = P_ * Ht * Si;

//  x_ = x_ + K * y;
//  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
//  P_ = (I - K * H_) * P_;
//  P_ -= K * H_ * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  if (px == 0.0 && py == 0.0) {
    return;
  }

  const double eps = 1e-3;
  const double rho = sqrt(px * px + py * py);
  VectorXd z_pred(3);
  z_pred << std::max(eps, rho), std::atan2(py, px), (px * vx + py * vy) / (std::max(eps, rho));

  VectorXd y = z - z_pred;
  UpdateCommon(y);
//  MatrixXd Ht = H_.transpose();
//  MatrixXd S = H_ * P_ * Ht + R_;
//  MatrixXd Si = S.inverse();
//  MatrixXd K = P_ * Ht * Si;

//  x_ = x_ + K * y;
//  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
//  P_ = (I - K * H_) * P_;
//  P_ -= K * H_ * P_;
}

void KalmanFilter::UpdateCommon(const VectorXd &y) {

  const MatrixXd PHt = P_ * H_.transpose();
  const MatrixXd S = H_ * PHt + R_;
  const MatrixXd K = PHt * S.inverse();

  x_ = x_ + K * y;
  P_ -= K * H_ * P_;
}
