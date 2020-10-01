#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

#define MIN_ABS_NON_ZERO (0.0001F)

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // Check correct estimation vector size
  if (estimations.size() != ground_truth.size() || estimations.size() == 0)
  {
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }

  // accumulate squared differences between ground truth and estimations called delta
  for (unsigned int i=0; i < estimations.size(); ++i)
  {
    VectorXd delta = estimations[i] - ground_truth[i];
    delta = delta.array()*delta.array();
    rmse += delta;
  }

  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */

  MatrixXd Hj(3,4);
  // Get variables for Jacobi matrix calculation
  double px_d = x_state(0U);
  double py_d = x_state(1U);
  double vx_d = x_state(2U);
  double vy_d = x_state(3U);

  double c1_d = px_d * px_d + py_d * py_d;
  double c2_d = sqrt(c1_d);
  double c3_d = c1_d * c2_d;

  // Check for division by zero
  if (fabs(c1_d) < MIN_ABS_NON_ZERO)
  {
    c1_d = MIN_ABS_NON_ZERO;
  }

  if (fabs(c2_d) < MIN_ABS_NON_ZERO)
  {
    c2_d = MIN_ABS_NON_ZERO;
  }

  if (fabs(c3_d) < MIN_ABS_NON_ZERO)
  {
    c3_d = MIN_ABS_NON_ZERO;
  }

  // Generate the Jacobi matrix
  Hj << (px_d/c2_d), (py_d/c2_d), 0.0, 0.0,
       -(py_d/c1_d), (px_d/c1_d), 0.0, 0.0,
         py_d*(vx_d*py_d - vy_d*px_d)/c3_d, px_d*(px_d*vy_d - py_d*vx_d)/c3_d, px_d/c2_d, py_d/c2_d;

  return Hj;
}
