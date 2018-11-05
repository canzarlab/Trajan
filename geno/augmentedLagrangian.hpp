// Copyright Jens K. Mueller
// Friedrich-Schiller-University Jena
//

#pragma once
#ifndef AUGMENTED_LAGRANGIAN_H
#define AUGMENTED_LAGRANGIAN_H

#include "lbfgsb.hpp"

// using Scalar = double;
// using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
// using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

class AugmentedLagrangian
{
public:
  AugmentedLagrangian(GenoNLP& genoNLP, size_t correctionPairs);

  SolverStatus solve();

  // Set x for warm start.
  //  const Vector& x(const Vector& x);
  const Scalar* x() const;
  const Scalar* y() const;
  Scalar f() const;

  // Set parameters
  bool setParameter(std::string parameter, Scalar value);

private:
  // disable default constructor
  AugmentedLagrangian();
  // disable copy constructor
  AugmentedLagrangian(const AugmentedLagrangian& other);
  
  GenoNLP& _genoNLP;
  double _constraintsTol;
  double _c;
  size_t _correctionPairs;
  Scalar _f;
  Vector _x;
  Vector _y;
  size_t _maxIter;

  Scalar _tol;
  Scalar _tolFun;
  bool _verbose;
};

#endif // AUGMENTED_LAGRANGIAN_H
