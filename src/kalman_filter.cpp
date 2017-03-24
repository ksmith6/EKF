#include "kalman_filter.h"
#include <iostream>
using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  // R_ = R_in; - deprecated
  Q_ = Q_in;

  // TODO - Re-do Init() to accept R_laser_ and R_radar_.
  R_Laser_ = MatrixXd(2,2);
  R_Radar_ = MatrixXd(3,3);

}

void KalmanFilter::Predict() {
  /**
    * predict the state
  */
  
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  
  P_ = F_ * P_ * Ft + Q_;

  cout << "x_ + = " << endl;
  cout << x_ << endl;
  cout << "P_ + = " << endl;
  cout << P_ << endl;

}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    This is the LASER measurement processing.
    * update the state by using Kalman Filter equations
  */
  cout << "Mapping measurement to state" << endl;
  VectorXd z_pred = H_ * x_;

  cout << "Computing residual" << endl;
  cout << "z.size() = " << z.size() << endl;
  cout << "z_pred.size() = " << z_pred.size() << endl;
  VectorXd y = z - z_pred;
  cout << "y = " << endl;
  cout << y << endl;

  cout << "Transposing H_" << endl;
  MatrixXd Ht = H_.transpose();
  cout << "Ht = " << endl;
  cout << Ht << endl;

  cout << "Computing S" << endl;
  MatrixXd S = H_ * P_ * Ht + R_Laser_;

  cout << "H_ = " << endl;
  cout << H_ << endl;

  cout << "P_ = " << endl;
  cout << P_ << endl;

  cout << "R_Laser_ = " << endl;
  cout << R_Laser_ << endl;

  cout << "S = " << endl;
  cout << S << endl;

  cout << "Inverting S" << endl;
  MatrixXd Si = S.inverse();
  cout << "S^-1 = " << endl;
  cout << Si << endl;

  cout << "Computing P * H'" << endl;
  MatrixXd PHt = P_ * Ht;

  cout << "Computing Kalman Gain K" << endl;
  MatrixXd K = PHt * Si;
  cout << "K = " << endl;
  cout << K << endl;

  //new estimate
  cout << "Updating state estimate " << endl;
  x_ = x_ + (K * y);
  long x_size = x_.size();

  MatrixXd I = MatrixXd::Identity(x_size, x_size);

  cout << "Updating covariance" << endl;
  cout << "P- = " << endl;
  cout << P_ << endl;
  P_ = (I - K * H_) * P_;

  cout << "P+ = " << endl;
  cout << P_ << endl;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    Update the state by using Extended Kalman Filter equations
  */

  // TODO: 
  // tools = Tools();
  cout << "Measurement z = " << endl;
  cout << z << endl;

  cout << "x_ = " << endl;
  cout << x_ << endl;

  MatrixXd Hj = tools.CalculateJacobian(x_); 
  cout << "Jacobian Hj = " << endl;
  cout << Hj << endl;

  VectorXd z_pred = Hj * x_;

  cout << "z_pred= " << endl;
  cout << z_pred << endl;

  cout << "Computing residual y" << endl;
  VectorXd y = z - z_pred;

  cout << "Transposing Hj" << endl;
  MatrixXd Ht = Hj.transpose();

  cout << "Computing S" << endl;
  MatrixXd S = Hj * P_ * Ht + R_Radar_;

  cout << "Inverting S" << endl;
  MatrixXd Si = S.inverse();

  cout << "P * Ht" << endl;
  MatrixXd PHt = P_ * Ht;

  cout << "Kalman Gain K" << endl;
  MatrixXd K = PHt * Si;

  //new estimate
  cout << "Updating state" << endl;
  x_ = x_ + (K * y);

  cout << "Updating covariance" << endl;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;

}

