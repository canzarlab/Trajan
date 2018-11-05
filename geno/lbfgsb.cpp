#include "lbfgsb.hpp"
#include "lineSearch.hpp"

#include <cassert>
#include <cstdio>
#include <cmath>

#include <limits>
#include <vector>
#include <algorithm>
#include <utility>

#include <iostream>

#include <Eigen/LU>

#ifdef MATLAB_MEX_FILE
#include <mex.h>
#endif

void LBFGSB::boundInit() {
  _genoNLP.getBounds(_lb.data(), _ub.data());
  Scalar const INF = std::numeric_limits<Scalar>::infinity();
  for (int i = 0; !_constrained && i < _lb.size(); ++i)
    _constrained = (_lb[i] != -INF || _ub[i] != INF);
  for (int i = 0; _boxed && i < _x.size(); ++i)
    _boxed = (_lb[i] != -INF && _ub[i] != INF);
}

void LBFGSB::updateB() {
  Scalar ys = _y.dot(_s);
  Scalar yy = _y.dot(_y);
  Scalar const eps = std::numeric_limits<Scalar>::epsilon();
  if (ys > (eps * yy)) {
    //		  _theta = _y.dot(_y) / ys / 2;  // TODO redo
    _theta = _y.dot(_y) / ys;
    Vector yws = _y;
    Vector sws = _s;
    yws = yws.cwiseProduct(_workingSet);
    sws = sws.cwiseProduct(_workingSet);

    if (_S.cols() < _m) {
      _S.conservativeResize(_S.rows(), _S.cols() + 1);
      _Y.conservativeResize(_Y.rows(), _Y.cols() + 1);
      _SY.conservativeResize(_SY.rows() + 1, _SY.cols() + 1);
      _SS.conservativeResize(_SS.rows() + 1, _SS.cols() + 1);
      _SYws.conservativeResize(_SYws.rows() + 1, _SYws.cols() + 1);
      _SSws.conservativeResize(_SSws.rows() + 1, _SSws.cols() + 1);
      _YYws.conservativeResize(_YYws.rows() + 1, _YYws.cols() + 1);
    }  
    _S.col(_index) = _s;
    _Y.col(_index) = _y;
    
    _SY.col(_index) = _S.transpose()*_y;
    _SS.col(_index) = _S.transpose()*_s;
    _SYws.col(_index) = _S.transpose() * yws;
    _SSws.col(_index) = _S.transpose() * sws;
    _SYws.row(_index) = _Y.transpose() * sws;
    _SY.row(_index) = _Y.transpose()*_s;
    _YYws.col(_index) = _Y.transpose() * yws;
    _SS.row(_index) = _SS.col(_index).eval();
    _SSws.row(_index) = _SSws.col(_index).eval();
    _YYws.row(_index) = _YYws.col(_index).eval();
  
    _index = (_index + 1) % _m;
   } else {
    std::clog << "skipping LBFGS-B update\n";
    return;
  }
  _Q = Matrix::Identity(_S.cols(), _S.cols());
  if (_S.cols() >= _m) {  // rotate by index % m to the right
    std::rotate(_Q.data(), _Q.data() + ((_m - _index) % _m) * _m,
		_Q.data() + _m * _m);
  }
  Matrix SS2 = _Q * _SS * _Q.transpose();
  Matrix SY2 = _Q * _SY * _Q.transpose();
  Matrix D = Matrix(_SY.rows(), _SY.cols());
  D.setZero();
  D.diagonal() = SY2.diagonal();
  Matrix L = SY2.triangularView<Eigen::StrictlyLower>();

  Matrix tmp(D.rows() * 2, 2 * D.cols());
  tmp.topLeftCorner(D.rows(), D.cols()) = _theta * SS2;
  tmp.topRightCorner(D.rows(), D.cols()) = L;
  tmp.bottomLeftCorner(D.rows(), D.cols()) = L.transpose();
  tmp.bottomRightCorner(D.rows(), D.cols()) = -D;

  Eigen::FullPivLU<Matrix> lu(tmp);
  _M = lu.inverse();

  _M.topLeftCorner(D.rows(), D.cols()) *= pow(_theta, 2);
  _M.topRightCorner(D.rows(), D.cols()) *= _theta;
  _M.bottomLeftCorner(D.rows(), D.cols()) *= _theta;

  _M.topLeftCorner(D.rows(), D.cols()) = _Q.transpose() * _M.topLeftCorner(D.rows(), D.cols()) * _Q;
  _M.topRightCorner(D.rows(), D.cols()) = _Q.transpose() * _M.topRightCorner(D.rows(), D.cols()) * _Q;
  _M.bottomLeftCorner(D.rows(), D.cols()) = _Q.transpose() * _M.bottomLeftCorner(D.rows(), D.cols()) * _Q;
  _M.bottomRightCorner(D.rows(), D.cols()) = _Q.transpose() * _M.bottomRightCorner(D.rows(), D.cols()) * _Q;
}

void LBFGSB::choleskyFactorK() {
  Index n = _K.rows() / 2;
  Matrix A1 = -_K.topLeftCorner(n, n);
  Matrix A2 = -_K.bottomLeftCorner(n, n);
  Matrix A3 = _K.bottomRightCorner(n, n);

  Matrix L1 = A1.llt().matrixL();
  Matrix U1 = L1.transpose();
  Matrix B2 = A2;
  U1.triangularView<Eigen::Upper>().solveInPlace<Eigen::OnTheRight>(B2);
  Matrix B3 = A3 + B2*B2.transpose();
  Matrix L3 = B3.llt().matrixL();

  _L = Matrix::Zero(n * 2, 2 * n);
  _L.topLeftCorner(n, n) = L1;
  _L.bottomLeftCorner(n, n) = B2;
  _L.bottomRightCorner(n, n) = L3;

  _E = Matrix::Zero(n * 2, 2 * n);
  _E.topLeftCorner(n, n) = -Matrix::Identity(n, n);
  _E.bottomRightCorner(n, n) = Matrix::Identity(n, n);
}

void LBFGSB::formK() {
  Vector diffWS = _workingSet - _workingSetOld;
  for (Index i = 0; i < _workingSet.size(); ++i) {
    if (diffWS[i]) {
      _SSws += diffWS[i] * _S.row(i).transpose() * _S.row(i);
      _SYws += diffWS[i] * _S.row(i).transpose() * _Y.row(i);
      _YYws += diffWS[i] * _Y.row(i).transpose() * _Y.row(i);
    }
  }

  Matrix SY2 = _Q * _SY * _Q.transpose();
  Matrix D = Matrix(_SY.rows(), _SY.cols());
  D.setZero();
  D.diagonal() = SY2.diagonal();

  Matrix M1 = SY2.triangularView<Eigen::StrictlyLower>();
  Matrix SYws2 = _Q * _SYws * _Q.transpose();
  M1 = M1 - SYws2;
  Matrix YYws2 = _Q * _YYws * _Q.transpose();

  Matrix SS2 = _Q * _SS * _Q.transpose();
  Matrix SSws2 = _Q * _SSws * _Q.transpose();

  _K = Matrix::Zero(D.rows() * 2, 2 * D.cols());

  _K.topLeftCorner(D.rows(), D.cols()) = -D-YYws2/_theta;
  _K.topRightCorner(D.rows(), D.cols()) = M1.transpose();
  _K.bottomLeftCorner(D.rows(), D.cols()) = M1;
  _K.bottomRightCorner(D.rows(), D.cols()) = _theta*(SS2-SSws2);

  choleskyFactorK();
}

void LBFGSB::cauchyPoint() {
  _c = Vector::Zero(2 * _S.cols());

  Vector bnd(_x.size());
  Vector d(_x.size());
  typedef std::pair<Scalar, Index> Pair;
  std::vector<Pair> vector;
  vector.reserve(_x.size());

  for (Index i = 0; i < _g.size(); ++i) {
    Scalar t = std::numeric_limits<Scalar>::infinity();
    if (_g[i] < 0) {
      t = (_x[i] - _ub[i]) / _g[i];
      bnd[i] = _ub[i];
    } else if (_g[i] > 0) {
      t = (_x[i] - _lb[i]) / _g[i];
      bnd[i] = _lb[i];
    }
    if (t == 0) {
      d[i] = 0;
      _workingSet[i] = 0;
    } else {
      d[i] = -_g[i];
      _workingSet[i] = 1;
    }
    //	    vector.emplace_back(t, i);
    vector.push_back(Pair(t, i));
  }

  auto nonzero_begin = std::partition(vector.begin(), vector.end(),
      [](Pair const& p) {return p.first == 0.0;});
  auto heapComp = [](Pair const& p1, Pair const& p2)
  {return p1.first >= p2.first;};
  std::make_heap(nonzero_begin, vector.end(), heapComp);

  // Initialize
  Vector p(2 * _S.cols());
  p.head(_S.cols()) = _S.transpose() * d;
  p.tail(_S.cols()) = _Y.transpose() * d;

  Scalar f1 = -d.dot(d);
  Scalar f2 = -_theta * f1 - p.transpose() * _M * p;

  _xCP = _x;
  Scalar tOld = 0;

  Scalar deltaTMin = 0;
  auto it = nonzero_begin;
  auto heap_end = vector.end();
  for (; it != heap_end; --heap_end) {
    deltaTMin = (0 == f1 && 0 == f2) ? 0 : -f1/f2;  // TODO think about just f2 == 0
    Index b = it->second;
    Scalar tThis = it->first;
    Scalar deltaT = tThis - tOld;

    if (deltaTMin < deltaT) break;
    std::pop_heap(it, heap_end, heapComp);

    assert(0 != _g[b]);
    _xCP(b) = bnd(b);
    Scalar z = _xCP(b) - _x(b);
    _c = _c + deltaT*p;
    Vector tmp(2 * _S.cols());// = _W.row(b);
    tmp.head(_S.cols()) = _S.row(b);
    tmp.tail(_S.cols()) = _Y.row(b);

    f1 = f1 + deltaT*f2 + _g(b)*_g(b) + _theta*_g(b)*z - _g(b)*(tmp.dot(_M*_c));
    // TODO
    // -1 should + 3
    f2 = f2 - 1 * _theta*_g(b)*_g(b) - 2*_g(b)*(tmp.dot(_M*p)) -
      _g(b) * _g(b)*(tmp.dot(_M*tmp));
    p = p + _g(b)*tmp;
    d(b) = 0;
    tOld = tThis;
    _workingSet[b] = 0;
  }

  deltaTMin = std::max(deltaTMin, Scalar(0.0));
  tOld = tOld + deltaTMin;
  _xCP += tOld * d;

  _c = _c + deltaTMin*p;
}

void LBFGSB::subsm() {
  Vector tmp2 = _M * _c;
  Vector tmp3 = _S * tmp2.head(_S.cols()) + _Y * tmp2.tail(_Y.cols());
  Vector r = (tmp3 - (_g + _theta*(_xCP - _x))).cwiseProduct(_workingSet);
  Vector wv(_M.rows());
  wv.head(_Y.cols()) = _Y.transpose() * r;
  wv.tail(_S.cols()) = _S.transpose() * r;
  wv.tail(_S.cols()) *= _theta;

  wv.head(_Y.cols()) = _Q * wv.head(_Y.cols());
  wv.tail(_S.cols()) = _Q * wv.tail(_S.cols());
  Vector wvl = wv;
  _L.triangularView<Eigen::Lower>().solveInPlace(wvl);
  wvl = _E*wvl;
 _L.transpose().triangularView<Eigen::Upper>().solveInPlace(wvl);
  wvl.tail(_S.cols()) *= _theta;
  wvl.head(_Y.cols()) = _Q.transpose() * wvl.head(_Y.cols());
  wvl.tail(_S.cols()) = _Q.transpose() * wvl.tail(_S.cols());


  Vector tmp8 = _Y * wvl.head(_Y.cols()) + _S * wvl.tail(_S.cols());
  Vector d = 1.0/_theta * r + 1.0 / (_theta * _theta) * tmp8;
  d = d.cwiseProduct(_workingSet);


  _xSubspace = _xCP + d;
  // Project _xSubspace back to bounds.
  _xSubspace = _xSubspace.array().max(_lb.array()).min(_ub.array());

  const Scalar machinePrecision = std::numeric_limits<Scalar>::epsilon();
  Scalar dummy = (_xSubspace - _x).dot(_g);
  if (dummy < 0 || d.norm() < machinePrecision)
  {
    return;
  }

  if (_verbose)
    printf("\npositive directional derivative, using backtracking line search\n");

  // Backtrack solution
  Vector alphaLB = (d.array() >= 0).select(std::numeric_limits<Scalar>::infinity(), (_lb - _xCP).cwiseQuotient(d));
  Vector alphaUB = (d.array() <= 0).select(std::numeric_limits<Scalar>::infinity(), (_ub - _xCP).cwiseQuotient(d));

  Index indexLB;
  Scalar minAlphaLB = alphaLB.minCoeff(&indexLB);
  Index indexUB;
  Scalar minAlphaUB = alphaUB.minCoeff(&indexUB);

  Index index2 = 0;
  Scalar alpha = std::min(minAlphaLB, minAlphaUB);
  if (alpha == minAlphaUB)
    index2 = 1;

  if (alpha < 1) {
    _xSubspace = _xCP + d * alpha;
    // Just to eliminate roundoff errors.
    if (0 == index2)  // Did we hit a lower bound?
      _xSubspace(indexLB) = _lb(indexLB);
    else
      _xSubspace(indexUB) = _ub(indexUB);
  } else {
    _xSubspace = _xCP + d;
  }
}

Scalar LBFGSB::lineSearchMT(const Vector& deltaX, size_t iter) {
  Scalar const INF = std::numeric_limits<Scalar>::infinity();
  double stp;
  double a1, a2;
  double big  = 1.0e10;
  double ftol = 1.0e-3,
         gtol = 0.9,
         xtol = 0.1;
  Scalar dnorm = deltaX.norm();
  // Determine the maximum step length.
  Scalar stpmx = big;
  if (_constrained) {
    if (iter == 1) {
      stpmx = 1.0;
    } else {
      for (Index i = 0; i < _x.size(); ++i) {
        a1 = deltaX[i];
	//        if (_nbd[i] != 0) {     
	if (_lb[i] != -INF || _ub[i] != INF) {
	  //      if (a1 < zero && _nbd[i] <= 2) {
	  if (a1 < 0.0 && _lb[i] != -INF) {
	    a2 = _lb(i) - _x(i);
            if (a2 >= 0.0) {
              stpmx = 0.0;
            } else if (a1*stpmx < a2) {
              stpmx = a2/a1;
            }
	    //          } else if (a1 > zero && _nbd[i] >= 2) {
          } else if (a1 > 0.0 && _ub[i] != INF) {
            a2 = _ub(i) - _x(i);
            if (a2 <= 0.0) {
              stpmx = 0.0;
            } else if (a1*stpmx > a2) {
              stpmx = a2/a1;
            }
          }
        }
      }
    }
  }
  if (iter == 1 && !_boxed)
    stp = std::min(1.0/dnorm, stpmx);
  else
    stp = 1.0;
  _gNew = _g;
  _fNew = _f;
  TaskType csave = TaskType::START;
  Scalar gd = _gNew.dot(deltaX);
  if (gd >= 0.0) {
    std::cerr << "ascent direction in projection gd = " << gd << std::endl;
    stp = -1;
    return stp;
    //    std::exit(EXIT_FAILURE);
  }
  int funIter = 0;
  TMPContainer TCON;
  do {
    dcsrch(_fNew, gd, stp, ftol, gtol, xtol, 0.0, stpmx, csave, TCON);
    if (TaskType::ERROR == csave)
      std::exit(EXIT_FAILURE);
    if (csave == TaskType::CONVERGENCE || csave == TaskType::WARNING)
      break;
    _xNew = stp * deltaX + _x;
    _genoNLP.functionValueAndGradient(_xNew.data(), _fNew, _gNew.data());
    ++_funEval;
    gd = _gNew.dot(deltaX);
    ++funIter;
    if (funIter >= 19) {
      stp = -1;  // something went wrong
      break;
    }
  } while(true);
  return stp;
}

void LBFGSB::projectedGradient()
{
  _projectedGradient = _g.array().max((_x - _ub).array());
  _projectedGradient = _projectedGradient.array().min((_x - _lb).array());

  // Compute projected gradient.
  //  for (Index i = 0; i < _x.size(); ++i) {
  //    _projectedGradient[i] = 
  //      ((_g[i] > 0 && _x[i] == _lb[i]) || (_g[i] < 0 && _x[i] == _ub[i])) ?
  //      0 : _g[i];
  //  }
}

Scalar LBFGSB::helper() {
  _s = _xNew - _x;
  _x.swap(_xNew);
  _y = _gNew - _g;
  _g.swap(_gNew);

  Scalar deltaF = (_f - _fNew) / (fabs(_fNew) + 1);  
  // (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} in original implementation
  //  Scalar deltaF = (_f - _fNew) / std::max(std::abs(_f), std::max(std::abs(_fNew), 1.0));
  _f = _fNew;

  projectedGradient();
  return deltaF;
}

void LBFGSB::restart() {
  refresh();
  _f = std::numeric_limits<Scalar>::infinity();
  _funEval = 0;
}

void LBFGSB::refresh() {
  _index = 0;
  _S.resize(_x.size(), 0);
  _Y.resize(_x.size(), 0);
  _theta = 1;
  _SY.resize(0, 0);
  _SS.resize(0, 0);
  _SYws.resize(0, 0);
  _SSws.resize(0, 0);
  _YYws.resize(0, 0);
  _M.resize(0, 0);
  _Q.resize(0, 0);
  _K.resize(0, 0);
  _L.resize(0, 0);
  _E.resize(0, 0);
  _workingSet.setOnes();
  _fNew = (std::numeric_limits<Scalar>::infinity());
}

bool LBFGSB::initFeasible() {
  for (Index i = 0; i < _x.size(); ++i) {
    if (_lb[i] > _ub[i])
      return false;    // infeasible
    if (_lb[i] > _x[i])
      _x[i] = _lb[i];
    if (_ub[i] < _x[i])
      _x[i] = _ub[i];
  }
  return true;
}

SolverStatus LBFGSB::solve() {
  Scalar const INF = std::numeric_limits<Scalar>::infinity();
  if (!initFeasible()) {
    std::cout << "Prolem is not feasible." << std::endl;
    return INFEASIBLE;
  }

  _funEval = 1;
  _genoNLP.functionValueAndGradient(_x.data(), _f, _g.data());
  //    std::cout << "x " << _x.transpose() << std::endl;
  //    std::cout << "g " << _g.transpose() << std::endl;
  // check for NANs and INFs
  double sum = _f + _g.sum();
  if ((-INF < sum) && (sum < INF)) {
  } else {
    std::cout << "encountered INF or NAN in function value or gradient" << std::endl;
    return NUM_ERROR;      
  }
  projectedGradient();
  if (_projectedGradient.lpNorm<Eigen::Infinity>() < _tol) {
    std::clog << "already optimal" << std::endl;
    return SOLVED;
  }

  if (_verbose)
    printf("\n%8s %8s %14s %14s %14s\n",
        "Iteration", "FunEval", "Step Length", "FunValue", "Proj.Grad.");
  size_t k = 1;
  for (;; ++k) {
    // Backup the old state
    _workingSetOld = _workingSet;
    if (_constrained || 1 == k) {
      cauchyPoint();
    } else {
      _xCP = _x;
      _c = Vector::Zero(_S.cols() * 2);
    }

    formK();
    subsm();
    Vector deltaX = _xSubspace - _x;
    Scalar t = lineSearchMT(deltaX, k);
    if (-1 == t) {
      if (0 != _S.cols()) {
        refresh();
        std::clog << "refresh() called\n";
        --k;
        --_funEval;  // TODO remove
        continue;
      }
      //        std::cerr << "error in line search - abnormal termination" << std::endl;
      //        std::exit(EXIT_FAILURE);
      return SUBOPTIMAL;
    }
    Scalar deltaF = helper();

    if (_verbose) {
      printf(" %8lu %8lu %14.5g %14.5E %14.5E\n",
          k, _funEval, t, _fNew, _projectedGradient.lpNorm<Eigen::Infinity>());
    }

    if (_projectedGradient.lpNorm<Eigen::Infinity>() < _tol) {
      if (_verbose) printf("Directional derivative below tol\n");
      return SOLVED;
    }

    if (deltaF < _tolFun) {
      if (_verbose) printf("Change in function value below tolFun\n");
      return SOLVED;
    }

    if (k >= _maxIter) {
      if (_verbose) printf("Maximum iteration count reached\n");
      return SUBOPTIMAL;
    }

    updateB();
  }
}


LBFGSB::LBFGSB(GenoNLP& genoNLP, Index m, bool verbose)
  : _genoNLP(genoNLP),
  _f(std::numeric_limits<Scalar>::infinity()),
  _x(Vector(genoNLP.getN())),
  _lb(Vector(genoNLP.getN())),
  _ub(Vector(genoNLP.getN())),
  _g(Vector(genoNLP.getN())),
  _funEval(0),
  _verbose(verbose),
  _maxIter(50000),
    //    _maxIter(100),
  _tol(1e-8),
  //_tolFun(1E7 * std::numeric_limits<Scalar>::epsilon()),
    //  _tolFun(1E5 * std::numeric_limits<Scalar>::epsilon()),
  _tolFun(1E2 * std::numeric_limits<Scalar>::epsilon()),
  _m(m),
  _index(0),
  _s(Vector(genoNLP.getN())),
  _y(Vector(genoNLP.getN())),
  _S(Matrix(genoNLP.getN(), 0)),
  _Y(Matrix(genoNLP.getN(), 0)),
  _theta(1),
  _SY(Matrix(0, 0)),
  _SS(Matrix(0, 0)),
  _M(Matrix(0, 0)),
  _Q(Matrix(0, 0)),
  _K(Matrix(0, 0)),
  _L(Matrix(0, 0)),
  _E(Matrix(0, 0)),
  _SYws(Matrix(0, 0)),
  _SSws(Matrix(0, 0)),
  _YYws(Matrix(0, 0)),
  _xCP(Vector(genoNLP.getN())),
  _workingSet(Vector(genoNLP.getN())),
  _workingSetOld(Vector(genoNLP.getN())),
  _c(Vector(genoNLP.getN())),
  _xSubspace(Vector(genoNLP.getN())),
  _projectedGradient(Vector(genoNLP.getN())),
  _fNew(std::numeric_limits<Scalar>::infinity()),
  _xNew(Vector(genoNLP.getN())),
  _gNew(Vector(genoNLP.getN())),
  _constrained(false),
  _boxed(true)
{
  _workingSet.setOnes();
  boundInit();
  _genoNLP.getStartingPoint(_x.data());
}

const Scalar* LBFGSB::x(const Scalar* x) {
  _x = Vector::Map(x, _x.size());
  return _x.data();
}

const Scalar* LBFGSB::x() const {
  return _x.data();
}

Scalar LBFGSB::f() const {
  return _f;
}

bool LBFGSB::setParameter(std::string parameter, Scalar value)
{
  // TODO
  // lower case parameter

  if (parameter == "pgtol")
  {
    _tol = value;
    return true;
  }

  if (parameter == "factr")
  {
    _tolFun = value * std::numeric_limits<Scalar>::epsilon();
    return true;
  }

  return false;
}

