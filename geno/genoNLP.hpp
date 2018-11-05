#ifndef __GENONLP_HPP__
#define __GENONLP_HPP__

#include <cstdint>

typedef int_fast64_t Index;
typedef double Scalar;

/** Base class for all NLP's that use standard triplet matrix form
 *  and dense vectors.  This is the standard base class for all
 *  NLP's that use the standard triplet matrix form (as for Harwell
 *  routines) and dense vectors. The class TNLPAdapter then converts
 *  this interface to an interface that can be used directly by
 *  ipopt.
 *
 *  This interface presents the problem form:
 *  
 *     min f(x)
 *
 *     s.t. gL <= g(x) <= gU
 *
 *          xL <=  x   <= xU
 *
 *  In order to specify an equality constraint, set gL_i = gU_i =
 *  rhs.  The value that indicates "infinity" for the bounds
 *  (i.e. the variable or constraint has no lower bound (-infinity)
 *  or upper bound (+infinity)) is set through the option
 *  nlp_lower_bound_inf and nlp_upper_bound_inf.  To indicate that a
 *  variable has no upper or lower bound, set the bound to
 *  -ipopt_inf or +ipopt_inf respectively
 */
class GenoNLP 
{
public:
  GenoNLP() {}
  virtual ~GenoNLP() 
  {
  }

  virtual bool getInfo(Index& n, Index& m) = 0;
  virtual Index getN()
  {
    Index n;
    Index m;
    getInfo(n, m);
    return n;
  }
  virtual Index getM()
  {
    Index n;
    Index m;
    getInfo(n, m);
    return m;
  }

  virtual bool getBounds(Scalar* lb, Scalar* ub) = 0;

  virtual bool getBoundsConstraints(Scalar* cl, Scalar* cu) 
  {
    (void)cl;
    (void)cu;
    return false;
  };

  virtual bool getStartingPoint(Scalar* x) = 0;

  virtual bool getStartingPointDual(Scalar* y) 
  {
    (void)y;
    return false;
  };

  // TODO
  // add bool 
  // if bool is true gradient needs to be computed
  virtual bool functionValueAndGradient(const Scalar* x,
					Scalar& functionValue,
					Scalar* gradient) = 0;

  virtual bool hessian(const Scalar* x,
		       Scalar* hessian)
  {
    (void)x;
    (void)hessian;
    return false;
  }

  virtual bool lagrangianHessian(const Scalar* x, 
				 const Scalar* y,
				 Scalar* lHessian)
  {
    (void)x;
    (void)y;
    (void)lHessian;
    return false;
  }

  

  virtual bool augmentedLagrangianHessian(const Scalar* x, 
					  const Scalar* y,
					  Scalar rho,
					  Scalar* augLHessian)
  {
    (void)x;
    (void)y;
    (void)rho;
    (void)augLHessian;
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

  virtual bool functionValueConstraints(const Scalar* x,
						   Scalar* constraintValues)
  {
    (void)x;
    (void)constraintValues;
    return false;
  }

  virtual bool gradientConstraintsTimesVector(const Scalar* x,
						   const Scalar* y,
						   Scalar* gradient)
  {
    (void)x;
    (void)y;
    (void)gradient;
    return false;
  }
  
private:
  GenoNLP(const GenoNLP&);
  void operator=(const GenoNLP&);
};

#endif
