#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

VectorXd Tools::CalculateRMSE(const std::vector<VectorXd> &estimations,
                              const std::vector<VectorXd> &ground_truth) {

  VectorXd rmse = VectorXd::Zero(4);
  
  if (estimations.size() == 0) {
    std::cerr << "Zero vector length" << std::endl;
    return rmse;
  }

  int len = estimations.size();
  if (len != ground_truth.size()) {
    std::cerr << "Vector lengths are not equal" << std::endl;
    return rmse;

  }

  //accumulate squared residuals
  for (int i = 0; i < len; ++i) {
    VectorXd diff = estimations[i] - ground_truth[i];
    diff = diff.array() * diff.array();
    rmse = rmse + diff;
  }

  //calculate the mean
  rmse = rmse / len;

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;

  
  return rmse;
}