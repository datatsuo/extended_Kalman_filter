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

  /*
  This function computes the root mean square errors for px, py, vx, vy
  */

  VectorXd rmse(4);  /* 4 corresponds to position (px, py) and velocity (vx, vy)*/
  rmse << 0, 0, 0, 0;

  /* in case the data is not correct, display an error message */
  if(estimations.size() != ground_truth.size() || estimations.size() ==0){
    cout << "invalid estimation or ground truth data" << endl;
    return rmse;
  }

  /* computation root mean square error */
  for(int i =0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array()*residual.array();
    rmse += residual;
  }
  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  /**
  TODO:
    * Calculate a Jacobian here.
  */

  /*
  This function computes Jacobian for (rho, phi, vrho) w.r.t  (px, py, vx, vy)
  */

  MatrixXd Hj(3, 4); /* 3 corresponds to (rho, phi, vrho), 4 to (px, py, vx, vy)*/

  /* load the position and velocity data in Cartesian coordinate */
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  /* useful quantities for computing Jacobian: (px^2+py^2)^n with n = 1, 1/2, 3/2 */
  float c1 = px*px + py*py;
  float c2 = sqrt(c1);
  float c3 = c1*c2;
  float c4 = vx*py - vy*px;

  /* in case px^2+py^2 is extremely small, display an error message */
  if(fabs(c1) < 0.001){
    cout << "Error in CalculateJacobian(): division by zero" << endl;
    return Hj;
  }

  /* computation of Jacobian */
  Hj << px/c2, py/c2, 0, 0,
        -py/c1, px/c1, 0, 0,
        py*c4/c3, -px*c4/c3, px/c2, py/c2;
  return Hj;

}
