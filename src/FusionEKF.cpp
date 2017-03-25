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
  R_laser_ << 0,0,0,0;
  R_radar_ = MatrixXd(3, 3);

  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  MatrixXd H_in = MatrixXd(2,4);
  H_in << 1, 0, 0, 0,
          0, 1, 0, 0;

  VectorXd x_in = VectorXd(4);
  x_in << 0.1,0.1,0.1,0.1;
  MatrixXd P_in = MatrixXd(4,4);
  P_in << 0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0;
  MatrixXd F_in = MatrixXd(4,4);
  F_in << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1;
  MatrixXd Q_in = MatrixXd(4,4);
  Q_in << 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0;
  ekf_.Init(x_in, P_in, F_in, H_in, R_laser_, Q_in);

  // Initialize the Measurement Noise Covariance (R) for the RADAR measurements.
  ekf_.R_Radar_ = MatrixXd(3,3);
  ekf_.R_Radar_(0,0) = 0.09;
  ekf_.R_Radar_(1,1) = 0.0009;
  ekf_.R_Radar_(2,2) = 0.09;

  // Initialize the Measurement Noise Covariance (R) for the LASER measurements
  const float LASER_VAR = 0.0225;
  ekf_.R_Laser_ = MatrixXd(2,2);
  ekf_.R_Laser_ <<  LASER_VAR, 0, 0, LASER_VAR;

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
    cout << "=========== FILTER INITIALIZATION ============ " << endl;

    cout << "Initializing filter with";
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float range = measurement_pack.raw_measurements_[0];
      float bearing = measurement_pack.raw_measurements_[1];
      float rangeRate = measurement_pack.raw_measurements_[2]; 

      // TODO - Convert polar measuremnts to Cartesian coordinates
      float px = range * cos(bearing);
      float py = range * sin(bearing);
      float vx = rangeRate * cos(bearing);
      float vy = rangeRate * sin(bearing);
      ekf_.x_ << px, py, vx, vy;
      cout << " Radar measurement." << endl;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize states.
      NOTE: Initial velocity states initialized with the assumption that the object is stationary.
      */
      float px = measurement_pack.raw_measurements_[0];
      float py = measurement_pack.raw_measurements_[1];
      float vx = 0.0;   // Assumes object is stationary
      float vy = 0.0;   // Assumes object is stationary

      // Initialize the states with the first laser measurement.
      ekf_.x_ << px, py, vx, vy;
      cout << " Laser measurement." << endl;

      
    }
    // Initialize the covariance matrix
    float pVar = 1000;  // Position Variance
    float vVar = 100; // Velocity Variance

    // TODO: Change pVar to measurement variance for position & velocity.
    ekf_.P_ << pVar, 0, 0, 0,
             0, pVar, 0, 0, 
             0, 0, vVar, 0,
             0, 0, 0, vVar; 


     

    // Save off the initial measurement time-step as the "previous"
    previous_timestamp_ = measurement_pack.timestamp_;


    // done initializing, no need to predict or update
    is_initialized_ = true;
    
  } else {

    /*
    NOT the first time-step (GENERAL CASE).  
    In this case, perform a state update from the previous state estimate 
    (and its associated time) to the current state, and then perform a measurement update.
    */
    /*****************************************************************************
     *  Prediction
      * Update the state transition matrix F according to the new elapsed time.
        - Time is measured in seconds.
      * Update the process noise covariance matrix.
     ****************************************************************************/
    if (ekf_.Debug) {
      cout << "Filter Already Initialized" << endl;
      cout << "     <------------- PREDICTION ----------------->" << endl;  
    }
    
    // Compute elapsed time from previous state estimate.
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;  //dt - expressed in seconds
    previous_timestamp_ = measurement_pack.timestamp_;

    /* Update the time-dependent elements of the State Transition Matrix 
     (F) to reflect the updated elapsed time. 
     */
    if (ekf_.Debug) {
      cout << "dt = " << dt << endl;
    }
    ekf_.F_(0,2) = dt;
    ekf_.F_(1,3) = dt;

    // Cache some shorthand variables to avoid repeated multiplication.
    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;

    // Cache shorthand variables from EKF tuning parameters to improve readability of updates to process noise matrix.
    float qx = ekf_.noise_ax;
    float qy = ekf_.noise_ay;

    // Update the process noise covariance matrix, Q.
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ <<  dt_4/4*qx,  0,          dt_3/2*qx,  0,
                0,          dt_4/4*qy,  0,          dt_3/2*qy,
                dt_3/2*qx,  0,          dt_2*qx,    0,
                0,          dt_3/2*qy,  0,           dt_2*qy;

    
    // Propagate the state to the current time.
    ekf_.Predict();

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    /**
       * Use the sensor type to perform the update step.
       * Update the state and covariance matrices.
     */
    if (ekf_.Debug) {
      cout << "    <------------- MEASUREMENT UPDATE ----------------->" << endl;
      cout << "Meas Update with ";
    }
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Radar updates - process with EKF logic due to nonlinearity of measurements.
      if (ekf_.Debug) {
        cout << "Radar" << endl;
      }
      // Check for zero measurements
      bool validMeas = measurement_pack.raw_measurements_[0] != 0.0 && measurement_pack.raw_measurements_[1] != 0.0 && measurement_pack.raw_measurements_[2] != 0.0;
      if (validMeas) {
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);  
      } else if (ekf_.Debug) {
        cout << "Invalid Radar measurement." << endl;
      }
      
    } else {
      // Laser updates - linear updates
      if (ekf_.Debug) {
        cout << "Laser";
      }
      bool validMeas = measurement_pack.raw_measurements_[0] != 0.0 && measurement_pack.raw_measurements_[1] != 0.0; 
      if (validMeas) {
        ekf_.Update(measurement_pack.raw_measurements_);
      } else if (ekf_.Debug) {
        cout << "Invalid Laser measurement." << endl;
      }
    }
    if (ekf_.Debug) {
      cout << "...complete" << endl;
    }
  }

  

  // print the output
  if (ekf_.Debug) {
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
  }
}
