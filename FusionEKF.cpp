#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
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
  R_laser_ << 0.04, 0,
        0, 0.015;

//  R_laser_ << 0.0225, 0,
//        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.20, 0, 0,
        0, 0.001, 0,
        0, 0, 0.10;

  //measurement covariance matrix - radar
//  R_radar_ << 0.09, 0, 0,
//        0, 0.0009, 0,
//        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;

  noise_ax = 16.0;
  noise_ay = 16.0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      cout << "Radar data!" << endl;

      float ro = measurement_pack.raw_measurements_[0];
      float theta = measurement_pack.raw_measurements_[1];
      float ro_dot = measurement_pack.raw_measurements_[2];

      ekf_.x_ << ro * cos(theta), ro * sin(theta), ro_dot * cos(theta), ro_dot * sin(theta);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  const double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  const double dt_2 = dt * dt;
  const double dt_3 = dt_2 * dt;
  const double dt_4 = dt_3 * dt;

  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
             0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
             dt_3 / 2 * noise_ax, 0, dt_2 * noise_ax, 0,
             0, dt_3 / 2 * noise_ay, 0, dt_2 * noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;

    std::ofstream Radar("/home/xiefeng/data/radar.txt", std::ios::app);
    Radar << "ro: " << std::fixed << measurement_pack.raw_measurements_[0] << "\t"
          << "theta: " << std::fixed << measurement_pack.raw_measurements_[1] << "\t"
          << "ro_dot: " << std::fixed << measurement_pack.raw_measurements_[2] << endl;
    Radar << std::endl;

    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

    std::ofstream FilterRadar("/home/xiefeng/data/FilterRadar.txt", std::ios::app);
    FilterRadar << std::fixed << ekf_.x_ << "\t" << endl;
    FilterRadar << endl;
  } else {
    // Laser updates
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;

    std::ofstream Laser("/home/xiefeng/data/laser.txt", std::ios::app);
    Laser << "px: " << std::fixed << measurement_pack.raw_measurements_[0] << "\t"
          << "py: " << std::fixed << measurement_pack.raw_measurements_[1] << endl;
    Laser << endl;

    ekf_.Update(measurement_pack.raw_measurements_);

    std::ofstream FilterLaser("/home/xiefeng/data/FilterLaser.txt", std::ios::app);
    FilterLaser << std::fixed << ekf_.x_ << "\t" << endl;
    FilterLaser << endl;
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
