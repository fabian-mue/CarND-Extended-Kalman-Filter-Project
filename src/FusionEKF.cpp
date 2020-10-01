#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

#define FAC_MICROSEC_IN_SEC 1000000.0F

/**
 * Constructor.
 */
FusionEKF::FusionEKF()
{
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

  // Measurement matrix - laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // State covariance matrix
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;

  // Init state transition matrix
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;

  // Init noise covariance matrix
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << 1, 0, 1, 0,
         0, 1, 0, 1,
         1, 0, 1, 0,
         0, 1, 0, 1;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
  /**
   * Initialization
   */
  if (!is_initialized_)
  {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    double px_d;
    double py_d;
    double vx_d;
    double vy_d;

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      double rho_d = measurement_pack.raw_measurements_[0U];
      double phi_d = measurement_pack.raw_measurements_[1U];
      double rhodot_d = measurement_pack.raw_measurements_[2U];

      px_d = rho_d * cos(phi_d);
      py_d = rho_d * sin(phi_d);
      vx_d = rhodot_d * cos(phi_d);
      vy_d = rhodot_d * sin(phi_d);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      // TODO: Initialize state.
        px_d = measurement_pack.raw_measurements_[0U];
        py_d = measurement_pack.raw_measurements_[1U];
        vx_d = 0.0F;
        vy_d = 0.0F;
    }
    else
    {
        cout << "Sensor unknown" << endl;

        px_d = 1.0F;
        py_d = 1.0F;
        vx_d = 1.0F;
        vy_d = 1.0F;
    }

    ekf_.x_ << px_d, py_d, vx_d, vy_d;

    // Save the timestamp of last measurement
    previous_timestamp_ = measurement_pack.timestamp_;

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

  // Calculate time delta and update last timestamp
  double dt_d = (measurement_pack.timestamp_ - previous_timestamp_) / FAC_MICROSEC_IN_SEC;
  previous_timestamp_ = measurement_pack.timestamp_;

  // Update the state transition matrix;
  ekf_.F_(0, 2) = dt_d;
  ekf_.F_(1, 3) = dt_d;

  // Noise values from the task
  double noise_ax = 9.0;
  double noise_ay = 9.0;

  // Update noise covariance matrix

  double dt2_d = dt_d * dt_d;
  double dt3_d = dt2_d * dt_d;
  double dt4_d = dt3_d * dt_d;

  ekf_.Q_ <<  dt4_d/4*noise_ax, 0, dt3_d/2*noise_ax, 0,
              0, dt4_d/4*noise_ay, 0, dt3_d/2*noise_ay,
              dt3_d/2*noise_ax, 0, dt2_d*noise_ax, 0,
              0, dt3_d/2*noise_ay, 0, dt2_d*noise_ay;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // TODO: Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  else
  {
    // TODO: Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
