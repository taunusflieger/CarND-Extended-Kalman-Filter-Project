#include <iostream>
#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

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



VectorXd KalmanFilter::ConvertToPolarCoords(const VectorXd& x_state) {

  VectorXd h = VectorXd(3); // h(x_)

  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float psq = px*px+py*py;

  //check division by zero
  if (psq < 0.000001){
    cout << "ConvertToPolarCoords () - Error - Division by Zero" << endl;
    return h;
  }
  if (fabs(px) < 0.000001) {
    cout << "ConvertToPolarCoords () - Error - Division by Zero" << endl;
    return h;
  }

  float c2 = sqrtf(psq);
  float c3 = atan2f(py, px);
  float c4 = (px * vx + py * vy) / c2;

  h << c2, c3, c4;

  return h;
}


void KalmanFilter::Predict() {
  /**
   * KF Prediction step
   */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * update the state by using Kalman Filter equations
   */

  //VectorXd z = measurements[n];
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K =  P_* H_.transpose()  * S.inverse();

  // new state estimate
  x_ = x_ + K * y;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}


 void KalmanFilter::UpdateEKF(const VectorXd &z) {
  
  VectorXd y = z - ConvertToPolarCoords(x_);

  // Normalize the result to be between -pi and pi
  y[1] -= (2 * M_PI) * floor((y[1] + M_PI) / (2 * M_PI));
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();

  // new state estimate
  x_ = x_ + K * y;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}


