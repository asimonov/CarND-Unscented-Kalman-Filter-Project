#include <iostream>
#include "tools.h"

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                       const vector<VectorXd> &ground_truth){

    VectorXd rmse(4);
    rmse << 0,0,0,0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if(estimations.size() != ground_truth.size()
       || estimations.size() == 0){
        cout << "Invalid estimation or ground_truth data" << endl;
        return rmse;
    }

    //accumulate squared residuals
    for(unsigned int i=0; i < estimations.size(); ++i){

        VectorXd residual = estimations[i] - ground_truth[i];

        //coefficient-wise multiplication
        residual = residual.array()*residual.array();
        rmse += residual;
    }

    //calculate the mean
    rmse = rmse/estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}


MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate Jacobian of our h(x).
  */
    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    //check division by zero
    float sumsq = px*px + py*py;
    if (sumsq <= 1e-6)
    {
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        // for EKF if our position is close to (0,0) let's use Hj of zeros so that
        // we regard the measurement as error and keep the uncertainty matrix as it was
        Hj << 0.,                       0.,                       0.,         0.,
              0.,                       0.,                       0.,         0.,
              0.,                       0.,                       0.,         0.;
    }
    //compute the Jacobian matrix
    else
    {
        float sqsumsq = sqrt(sumsq);
        Hj << px/sqsumsq,                       py/sqsumsq,                       0.,         0.,
              -py/sumsq,                        px/sumsq,                         0.,         0.,
              py*(vx*py-vy*px)/(sumsq*sqsumsq), px*(vy*px-vx*py)/(sumsq*sqsumsq), px/sqsumsq, py/sqsumsq;
    }

    return Hj;
}
