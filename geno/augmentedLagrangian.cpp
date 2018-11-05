#include <iostream>
#include <vector>

#include "lbfgsb.hpp"
#include "augmentedLagrangian.hpp"

static const Scalar INF = std::numeric_limits<Scalar>::infinity();

class AugmentedNLP : public GenoNLP
{
public:
  AugmentedNLP(GenoNLP& genoNLP, Scalar rho, const Vector& y)
    :
    _genoNLP(genoNLP),
    _n(genoNLP.getN()),
    _m(genoNLP.getM()),
    _mEqualities(0),
    _mInequalities(0),
    _constraintValues(Vector(genoNLP.getM())),
//    _jacobian(Matrix(genoNLP.getM(), genoNLP.getN())),
    _rho(rho),
    _y(y),
    _cl(Vector(genoNLP.getM())),
    _cu(Vector(genoNLP.getM())),
    _cType(std::vector<int>(genoNLP.getM())),
    _constraintError(Vector(genoNLP.getM()))
  {
    _genoNLP.getBoundsConstraints(_cl.data(), _cu.data());

    _mEqualities = 0;
    _mInequalities = 0;
    for (Index i = 0; i < _genoNLP.getM(); ++i)
    {
      if (_cl(i) == _cu(i)) {
	_cType[i] = 1;   // equality constraint
	++_mEqualities;
	if (_mInequalities) {
	  std::cout << "Error : Equalitites must precede inequalitites!" 
		    << std::endl;
	  exit(1);
	}
      } else {
	if (_cl(i) != -INF) {
	  _cType[i] = 2;   // has lower bound
	  ++_mInequalities;
	  //	  std::cout << "Error : No lower bounds allowed!"
	  //		    << std::endl;
	  //	  exit(1);
 
	  if ((_cu(i) != INF) && (_cl(i) != _cu(i))) {
	    _cType[i] = 3;
	    //   	    std::cout << "Error : Upper bounds and lower bound are not allowed at the same time!" << std::endl;
		    //	    	    exit(1);
		    //		    		    _cType[i] = 2;
	  }
	}
	if (_cl(i) == -INF) {
	  _cType[i] = 4;      // no lower bound
	  ++_mInequalities;
	}
      }
    }
    //    _s = Vector::Zero(_mInequalities);
    assert(_m = _mEqualities + _mInequalities);
  }

  virtual ~AugmentedNLP() 
  {
  }

  virtual bool getInfo(Index& n, Index& m) 
  {
    n = _n; 
    m = 0;
    return true;
  }

  virtual bool getBounds(Scalar* lb, Scalar* ub)
  {
    _genoNLP.getBounds(lb, ub);
    return true;
  }

  virtual bool getBoundsConstraints(Scalar* cl, Scalar* cu) 
  {
    (void)cl;
    (void)cu;
    return false;
  };

  virtual bool getStartingPoint(Scalar* x)
  {
    _genoNLP.getStartingPoint(x);
    return true;
  }

  virtual bool getStartingPointDual(Scalar* y) 
  {
	(void)y;
    return false;
  };

  void computeConstraintError()
  {
    _constraintError = _constraintValues;

    // for equality constraints
    _constraintError.head(_mEqualities) -= _cl.head(_mEqualities);

    for (size_t i = 0; i < _mInequalities; ++i)
    {
      if (_cType[i+_mEqualities] == 2)
      {
	// lower bound not -INF
	_constraintError(i+_mEqualities) -= _cl(i+_mEqualities);
      }
      if (_cType[i+_mEqualities] == 3)
      {
	// lower bound not -INF and upper not INF
	double err = std::max(_constraintError(i+_mEqualities) - _cu(i + _mEqualities), 0.0);
	err = std::min(_constraintError(i+_mEqualities) - _cl(i+_mEqualities), err);
	_constraintError(i+_mEqualities) = err;
	//	std::cout << "c = " << _constraintError(i+_mEqualities) << std::endl;
	//	std::cout << "e = " << err << std::endl;
	//	std::cout << "lb= " << _cl(i+_mEqualities) << std::endl;
	//	std::cout << "ub= " << _cu(i+_mEqualities) << std::endl;
      }
      if (_cType[i+_mEqualities] == 4)
      {
	// upper bound not INF
	_constraintError(i+_mEqualities) -= _cu(i+_mEqualities);
      }
    }
  }

  virtual bool functionValueAndGradient(const Scalar* x,
					Scalar& functionValue,
					Scalar* gradient)
  {
    _genoNLP.functionValueAndGradient(x, functionValue, gradient);

    _genoNLP.functionValueConstraints(x, _constraintValues.data());

    computeConstraintError();
    
    Vector toAdd(_constraintError);
    // for equality constraints
    toAdd.head(_mEqualities) += _y.head(_mEqualities) / _rho;

    // for inequality constraints
    toAdd.tail(_mInequalities) += _y.tail(_mInequalities) / _rho;
    for (size_t i = 0; i < _mInequalities; ++i) {
      if (_cType[i+_mEqualities] == 2)
        if (toAdd(_mEqualities + i) > 0)
          toAdd(_mEqualities + i) = 0;

      if (_cType[i+_mEqualities] == 4)
	if (toAdd(_mEqualities + i) < 0)
	  toAdd(_mEqualities + i) = 0;
    }
    functionValue += _rho/2 * toAdd.squaredNorm();

    Vector::MapType gradientMapX(gradient, _n);
    Vector gradientX(_n);

    Vector v = _rho * toAdd;
    _genoNLP.gradientConstraintsTimesVector(x, v.data(), gradientX.data());

    gradientMapX += gradientX;
    //    std::cout << "const " << _constraintValues.transpose() << std::endl;
    //    std::cout << "toAdd " << toAdd.transpose() << std::endl;
    return true;
  }


  virtual bool hessian(const Scalar* x,
		       Scalar* hessian)
  {
    (void)x;
    (void)hessian;
    /*
    Matrix::MapType hessianMap(hessian, 
			       _n + _mInequalities, 
			       _n + _mInequalities);

    //    _genoNLP.functionValueAndGradientConstraints(x,
    //                                                 _constraintValues.data(),
    //                                                 _jacobian.data());
    computeConstraintError();
    //    Matrix funHxx = Matrix(_n, _n);
    //    _genoNLP.hessian(x, funHxx.data());
    //    std::cout << x[0] << "  " << x[1] << std::endl;
    Matrix Hxx = Matrix(_n, _n);
    Vector v = _y + _rho * _constraintError;
    //    std::cout << "v= " << v.transpose() << std::endl;
    _genoNLP.lagrangianHessian(x, v.data(), Hxx.data());
    //    std::cout << "Hxx1 =\n" << Hxx << std::endl;
    Matrix JJ = (_jacobian.transpose() * _jacobian);
    Hxx += _rho * (_jacobian.transpose() * _jacobian);
    */
    /*    
    std::cout << "JJ =\n" << JJ << std::endl;
    std::cout << "jacobian =\n" << _jacobian << std::endl;
    std::cout << "Hxx2 =\n" << Hxx << std::endl;
    */
    
    /*
    Matrix Hxs = -_rho * _jacobian.bottomLeftCorner(_mInequalities, _n);
    
    Matrix Hss = _rho * Matrix::Identity(_mInequalities, _mInequalities);

    hessianMap.topLeftCorner(_n, _n) = Hxx; 
    */
    /*
    std::cout << "Hxs = \n" << Hxs << std::endl;
    std::cout << "n = " << _n 
	      << " _mInequalities = " << _mInequalities << std::endl;
    */
    /*
    hessianMap.topRightCorner(_n, _mInequalities) = Hxs.transpose();
    hessianMap.bottomLeftCorner(_mInequalities, _n) = Hxs;
    hessianMap.bottomRightCorner(_mInequalities, _mInequalities) = Hss;
    //    std::cout << "Hessian =\n" << hessianMap << std::endl;
    return true;
    */
    return false;
  }


  virtual bool functionValueAndGradientConstraints(const Scalar* x,
						   Scalar* constraintValues,
						   Scalar* jacobian)
  {
    (void)x;
    (void)constraintValues;
    (void)jacobian;
    return false;
  }

  const Vector& cu()
  {
    return _cu;
  }

  const Vector& cl()
  {
    return _cl;
  }

  const Vector& constraintValues()
  {
    return _constraintValues;
  }

  const Vector& constraintError()
  {
    return _constraintError;
  }

  void setParams(Scalar rho, const Vector& y)
  {
    _rho = rho;
    _y = y;
  }

  size_t mEqualities() {
    return _mEqualities;
  }

  size_t mInequalities() {
    return _mInequalities;
  }

  const std::vector<int>& cType() const {
    return _cType;
  }

private:
  AugmentedNLP(const AugmentedNLP&);
  void operator=(const AugmentedNLP&);
  
  GenoNLP& _genoNLP;
  size_t _n;
  size_t _m;
  size_t _mEqualities;
  size_t _mInequalities;
  Vector _constraintValues;
//  Matrix _jacobian;
  Scalar _rho;
  Vector _y;
  Vector _cl;
  Vector _cu;
  std::vector<int> _cType;
  Vector _constraintError;
};


AugmentedLagrangian::AugmentedLagrangian(GenoNLP& genoNLP, 
					 size_t correctionPairs)
  : _genoNLP(genoNLP),
    _constraintsTol(1e-4),
    _c(1),
    _correctionPairs(correctionPairs),
    _f(0),
    _x(Vector(genoNLP.getN())),
    _y(Vector::Zero(_genoNLP.getM())),
    _maxIter(50),
    _tol(-1),
    _tolFun(-1),
    _verbose(false)
{
}

SolverStatus AugmentedLagrangian::solve() 
{  
  static const double gamma = 2;
  //  static const double rhoMin = 1E-6;
  //  static const double rhoMax = 10;
    
  Scalar rho = 0;
//  Vector y = Vector::Zero(_genoNLP.getM());
  _genoNLP.getStartingPointDual(_y.data());
  //  std::cout << "y= " << y.transpose() << std::endl;
  AugmentedNLP augmentedNLP(_genoNLP, rho, _y);
  size_t mEqualities = augmentedNLP.mEqualities();
  size_t mInequalities = augmentedNLP.mInequalities();
  std::vector<int> cType = augmentedNLP.cType();
  

  /*
  _genoNLP.getStartingPoint(_x.data());
  Vector dummyG(_genoNLP.getN());
  _genoNLP.functionValueAndGradient(_x.data(), _f, dummyG.data());
  double dummyF;
  augmentedNLP.functionValueAndGradient(_x.data(), dummyF, dummyG.data());
  Vector constraintError = augmentedNLP.constraintError();

  for (size_t i = 0; i < mInequalities; ++i)
    if (constraintError(mEqualities + i) < 0)
      constraintError(mEqualities + i) = 0;
  std::cout << "constraint Error = " << constraintError.transpose() << std::endl;
  rho = 2 * (std::abs(_f) + 1) / constraintError.squaredNorm();
  if (rho < rhoMin)
    rho = rhoMin;
  if(rho > rhoMax)
    rho = rhoMax;
  */

  Vector constraintError;
  Vector cl = augmentedNLP.cl();
  Vector cu = augmentedNLP.cu();
  Vector dummyG(_genoNLP.getN());
  rho = 1.0;

  LBFGSB solver(augmentedNLP, _correctionPairs, _verbose);
  //NEWTONB solver(augmentedNLP, _correctionPairs);
  if (_tol != -1) solver.setParameter("pgtol", _tol);
  if (_tolFun != -1) solver.setParameter("factr", _tolFun);
  double oldFactor = INF;
  double tau = 0.5;
  double lambdaMin = -1E20;
  double lambdaMax = 1E20;
  double muMax = 1E20;
  double muMin = -1E20;
  SolverStatus status = SOLVED;
  
  for (size_t iter = 0; iter < _maxIter; ++iter) 
  {
    if (_verbose)
      std::cout << "rho = " << rho << std::endl;

    augmentedNLP.setParams(rho, _y);
    solver.restart();
    status = solver.solve();
    if ((status != SOLVED) && (status != SUBOPTIMAL))
      break;
    constraintError = augmentedNLP.constraintError();
    Vector trueConstraintError(constraintError);
    for (size_t i = 0; i < mInequalities; ++i) {
      if (cType[i+mEqualities] == 2)
        if (trueConstraintError(mEqualities + i) > 0)
          trueConstraintError(mEqualities + i) = 0;

      if (cType[i+mEqualities] == 4)
	if (trueConstraintError(mEqualities + i) < 0)
	  trueConstraintError(mEqualities + i) = 0;
    }
    double constraintErrorNorm = 
      trueConstraintError.lpNorm<Eigen::Infinity>();

    _x = Vector::Map(solver.x(), _x.size());

    //    std::cout << "cl " << cl.transpose() << std::endl;
    //    std::cout << "cu " << cu.transpose() << std::endl;
    //    std::cout << "constraint Values " << augmentedNLP.constraintValues().transpose() << std::endl;
    //    std::cout << "y " << y.transpose() << std::endl;
    if (constraintErrorNorm < _constraintsTol)
      break;

    //    y = y + rho * constraintError;

    _y.head(mEqualities) +=  rho * constraintError.head(mEqualities);

    //safeguard Lagrange multiplier for equality constraints
    for (size_t i = 0; i < mEqualities; ++i) {
      if (_y(i) < lambdaMin)
	_y(i) = lambdaMin;
      if (_y(i) > lambdaMax)
	_y(i) = lambdaMax;
    }


    Vector constraintValues = augmentedNLP.constraintValues();

    // compute new Lagrange multiplier and safeguard it
    for (size_t i = 0; i < mInequalities; ++i) {
      if (cType[i+mEqualities] == 2) {
	_y(i+mEqualities) += rho * constraintError(i+mEqualities);
        if (_y(i+mEqualities) > 0)
          _y(i+mEqualities) = 0;
      }

      if (cType[i+mEqualities] == 3) {
	double yLB = std::min(0.0, _y(i+mEqualities));
	double yUB = std::max(0.0, _y(i+mEqualities));
	yLB += rho * (constraintValues(i+mEqualities) - cl(i));
	yUB += rho * (constraintValues(i+mEqualities) - cu(i));
	yLB = std::min(0.0, yLB);
	yUB = std::max(0.0, yUB);

	_y(i+mEqualities) = yLB + yUB;

      }

      if (cType[i+mEqualities] == 4) {
        _y(i+mEqualities) += rho * constraintError(i+mEqualities);
 	if (_y(i+mEqualities) < 0)
	  _y(i+mEqualities) = 0;
      }

      if (_y(i+mEqualities) < muMin)
        _y(i+mEqualities) = muMin;
      if (_y(i+mEqualities) > muMax)
	_y(i+mEqualities) = muMax;
    }
    
    //    double factor = constraintErrorNorm;

    double factor = 0;
    // need to have the if there since if there are no equalities Eigen will make a mistake
    if (mEqualities)  factor = constraintError.head(mEqualities).lpNorm<Eigen::Infinity>();

    //    std::cout << "factor = " << factor << std::endl;
    for (size_t i = 0; i < mInequalities; ++i) {
      factor = std::max(std::abs(std::max(constraintError(mEqualities + i), -_y(mEqualities + i) / rho)), factor);
    }
    if (factor > tau * oldFactor) {
      rho *= gamma;
    }
    oldFactor = factor;
  }
  _genoNLP.functionValueAndGradient(_x.data(), _f, dummyG.data());
  return status;
}

/*
const Vector&  AugmentedLagrangian::x(const Vector& x)
{
  return _x;
}
*/

const Scalar* AugmentedLagrangian::x() const
{
  return _x.data();
}

const Scalar* AugmentedLagrangian::y() const
{
  return _y.data();
}
  
Scalar AugmentedLagrangian::f() const
{
  return _f;
}

bool AugmentedLagrangian::setParameter(std::string parameter, Scalar value)
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
		_tolFun = value;
		return true;
	}

	if (parameter == "constraintsTol")
	{
		_constraintsTol = value;
		return true;
	}
	if (parameter == "verbose")
	{
	        _verbose = value;
	        return true;
        }
	        

	return false;
}

