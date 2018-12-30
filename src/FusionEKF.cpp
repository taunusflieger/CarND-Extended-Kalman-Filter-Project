#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
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


  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

  // INFO: MZILL
  ekf_.F_ = MatrixXd(4, 4);

	// Process noise covariance matrix
	ekf_.Q_ = MatrixXd(4, 4);
	//create a 4D state vector, we don't know yet the values of the x state
	ekf_.x_ = VectorXd(4);
	
	//state covariance matrix P
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
						 0, 1, 0, 0,
						 0, 0, 1, 0,
						 0, 0, 0, 1;
	
	//measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
              
  //projection  matrix laser
	H_laser_ << 1, 0, 0, 0,
			  		  0, 1, 0, 0;

	//set the acceleration noise components
	noise_ax = 10;
	noise_ay = 10;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * Initialize the state ekf_.x_ with the first measurement.
     * Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF initialize: " << endl;
    //ekf_.x_ = VectorXd(4);
    //ekf_.x_ << 1, 1, 1, 1;

    // initialize state convariance matrix
    Eigen::MatrixXd P_init(4,4);
    P_init << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;

    // initialize state transition matrix
    Eigen::MatrixXd F_init(4,4);
    F_init << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

    // initialize process covariance matrix
    Eigen::MatrixXd Q_init(4,4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      cout << "Radar - step 1" << endl;

      // Convert radar from polar to cartesian coordinates and initialize state.
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
    
      Eigen::VectorXd x_init(4);
      x_init << rho * cos(phi), rho * sin(phi), 0, 0;

      ekf_.Init(x_init, P_init, F_init, Hj_, R_radar_, Q_init);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      cout << "Lidar - step 1" << endl;

      // Initialize state
      Eigen::VectorXd x_init(4);
      x_init << measurement_pack.raw_measurements_[0], 
                measurement_pack.raw_measurements_[1], 
                0,
                0;

      ekf_.Init(x_init, P_init, F_init, H_laser_, R_laser_, Q_init);
    }

    // save first timestamp
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }  // not initialized

  /*///////////////////////////////////////////////////////////////////////////////
   * Prediction
   *///////////////////////////////////////////////////////////////////////////////

  /**
   * Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * pdate the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

 // update state transition matrix
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  // update process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
            0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
            dt_3 / 2 * noise_ax, 0, dt_2 * noise_ax, 0,
            0, dt_3 / 2 * noise_ay, 0, dt_2 * noise_ay;

  ekf_.Predict();

  /*///////////////////////////////////////////////////////////////////////////////
   * Update
   *///////////////////////////////////////////////////////////////////////////////


  /**
   * 
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

    // Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;

    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;

     ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  // << "x_ = " << ekf_.x_ << endl;
  // cout << "P_ = " << ekf_.P_ << endl;
}
