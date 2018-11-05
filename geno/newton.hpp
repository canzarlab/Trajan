// Copyright Soeren Laue, Jens K. Mueller, Lars Kuehne
// Friedrich-Schiller-University Jena

#ifndef NEWTON_HPP_
#define NEWTON_HPP_

#include <Eigen/Dense>

#include "genoNLP.hpp"
#include "lbfgsb.hpp"

typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

//enum class SolverStatus {SOLVED = 0, SUBOPTIMAL, UNBOUNDED, INFEASIBLE, NUM_ERROR};


class Newton {
  private:
    GenoNLP& _genoNLP;
    Scalar _f;
    Vector _x;
    Vector _lb;
    Vector _ub;
    Vector _g;
    Matrix _h; // hessian matrix
    size_t _funEval; 
    bool _verbose;
    size_t _maxIter;
    Scalar _tol;
    Scalar _tolFun;
    Vector _s;
    Vector _y;
    Scalar _fNew;
    Vector _xNew;
    Vector _gNew;
    bool _constrained;
    bool _boxed;
  public:
    Newton(GenoNLP& genoNLP):
      _genoNLP(genoNLP), 
      _f(std::numeric_limits<Scalar>::infinity()),
      _x(Vector(genoNLP.getN())),
      _lb(Vector(genoNLP.getN())),
      _ub(Vector(genoNLP.getN())),
      _g(Vector(genoNLP.getN())),
      _h(Matrix(genoNLP.getN(), genoNLP.getN())),
      _funEval(0),
      _verbose(true),
      _maxIter(10000),
      _tol(1e-5),
      //_tolFun(1e+7 * std::numeric_limits<Scalar>::epsilon()),
      _tolFun(0 * std::numeric_limits<Scalar>::epsilon()),
      _s(Vector(genoNLP.getN())), 
      _y(Vector(genoNLP.getN())),
      _fNew(std::numeric_limits<Scalar>::infinity()),
      _xNew(Vector(genoNLP.getN())),
      _gNew(Vector(genoNLP.getN())),
      _constrained(false),
      _boxed(true)
      {     
        _genoNLP.getBounds(_lb.data(), _ub.data());
        _genoNLP.getStartingPoint(_x.data());   
      }

    SolverStatus solve();
     // Set x for warm start.
    const Scalar* x() const;
    const Scalar* x(const Scalar* x);
    Scalar f() const;
    Scalar helper();
    bool initFeasible();
    // lineSearch should not be a member of a class (rather friend?)
    Scalar lineSearchMT(const Vector& deltaX, size_t iter);
    bool setParameter(std::string parameter, Scalar value);    
    
};

#endif // NEWTON_HPP_

