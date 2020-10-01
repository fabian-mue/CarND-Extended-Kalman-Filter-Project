#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

#define MIN_ABS_NON_ZERO (0.0001F)

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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

   MatrixXd Ft = F_.transpose();
   x_ = F_ * x_;
   P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z)
{
  /**
   * TODO: update the state by using Kalman Filter equations
   */

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;

  Estimate(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

  double px_d = x_[0U];
  double py_d = x_[1U];
  double vx_d = x_[2U];
  double vy_d = x_[3U];

  double rho_d = sqrt(px_d * px_d + py_d * py_d);
  double phi_d = atan2(py_d, px_d);
  double rhodot_d;

  if (rho_d < MIN_ABS_NON_ZERO)
  {
      rhodot_d = MIN_ABS_NON_ZERO;
  }
  else
  {
      rhodot_d = (px_d * vx_d + py_d * vy_d) / rho_d;
  }

  VectorXd hx(3);
  hx << rho_d, phi_d, rhodot_d;

  VectorXd y = z - hx;
  bool checkPhi_b = true;
  // Check and maybe correct range of phi of vector y

  while(checkPhi_b)
  {
      if (y(1U) > M_PI)
      {
          y(1U) -= 2*M_PI;
      }
      else if (y(1U) < -M_PI)
      {
          y(1U) += 2*M_PI;
      }
      else
      {
          checkPhi_b = false;
      }
  }
  Estimate(y);
}

void KalmanFilter::Estimate(const VectorXd &y)
{
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
