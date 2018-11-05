#ifndef GENO_H
#define GENO_H

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "geno/genoNLP.hpp"
#include "geno/augmentedLagrangian.hpp"
#include <limits>
using namespace std;

Scalar const INF = numeric_limits<Scalar>::infinity();

// declares a column-major sparse matrix type of double
typedef Eigen::SparseMatrix<Scalar> SpMat;
typedef Eigen::Triplet<Scalar> ET;

class SimpleJRF : public GenoNLP
{
 public:
    SimpleJRF(const SpMat& A,
              const Vector& b,
              const Vector& c,
              Vector& x,
              Vector& y)
        : _A(A), _b(b), _c(c), _x(x), _y(y), _n(A.cols()), _m(A.rows())
     {
     }


    virtual bool getInfo(Index& n, Index& m)
    {
        // number of variables
        n = _n;

        // number of constraints (only real constraints, bound constraints do not count)
        m = _m;
        return true;
    }

    // bounds on the variables
    virtual bool getBounds(Scalar* lb, Scalar* ub)
    {
        Vector::MapType(lb, _n) = Vector::Constant(_n, 0.0);
        Vector::MapType(ub, _n) = Vector::Constant(_n, INF);
        return true;
    }

/*    virtual bool getBoundsConstraints(Scalar* cl, Scalar* cu)
    {
        // we have equality constraints here
        Vector::MapType(cl, _m) = _b;
        Vector::MapType(cu, _m) = Vector::Constant(_m, INF);
        return true;
    };
*/
    virtual bool getStartingPoint(Scalar* x)
    {
//        Vector::MapType(x, _n) = Vector::Zero(_n);
        Vector::MapType(x, _n) = _x;
        return true;
    }

    virtual bool getStartingPointDual(Scalar* y)
    {
//        Vector::MapType(y, _m) = Vector::Zero(_m);
        Vector::MapType(y, _m) = _y;
        return true;
    };

    virtual bool functionValueAndGradient(const Scalar* variablesPtr,
                                          Scalar& functionValue,
                                          Scalar* gradientPtr)
    {
        Vector::ConstMapType x = Vector::ConstMapType(variablesPtr, _n);
        Vector::MapType g = Vector::MapType(gradientPtr, _n);
        functionValue = _c.dot(x);
        g = _c;
        return true;
    }

    virtual bool functionValueConstraints(const Scalar* variablesPtr,
                                          Scalar* constraintValuesPtr)
    {
        Vector::ConstMapType x = Vector::ConstMapType(variablesPtr, _n);
        Vector::MapType constraintValues = Vector::MapType(constraintValuesPtr, _m);
        constraintValues = _A * x;
        return true;
    }

    virtual bool gradientConstraintsTimesVector(const Scalar* variablesPtr,
                                                const Scalar* dualVariablesPtr,
                                                Scalar* gradientPtr)
    {
        Vector::ConstMapType x = Vector::ConstMapType(variablesPtr, _n);
        Vector::ConstMapType y = Vector::ConstMapType(dualVariablesPtr, _m);
        Vector::MapType gradient = Vector::MapType(gradientPtr, _n);
        gradient = _A.transpose() * y;
        return true;
    }

protected:
    const SpMat& _A;
    const Vector& _b;
    const Vector& _c;
    Vector& _x;
    Vector& _y;
    Index _n;
    Index _m;
};

class CoveringJRF : public SimpleJRF
{
 public:
    CoveringJRF(const SpMat& A,
                const Vector& b,
                const Vector& c,
                Vector& x,
                Vector& y): SimpleJRF(A, b, c, x, y) { }
    bool getBoundsConstraints(Scalar* cl, Scalar* cu) override
    {
    // we have equality constraints here
        Vector::MapType(cl, _m) = _b;
        Vector::MapType(cu, _m) = Vector::Constant(_m, INF);
        return true;
    }

};


class PackingJRF : public SimpleJRF
{
 public:
    PackingJRF(const SpMat& A,
                const Vector& b,
                const Vector& c,
                Vector& x,
                Vector& y): SimpleJRF(A, b, c, x, y) { }
    bool getBoundsConstraints(Scalar* cl, Scalar* cu) override
    {
    // we have equality constraints here
        Vector::MapType(cu, _m) = _b;
        Vector::MapType(cl, _m) = Vector::Constant(_m, -INF);
        return true;
    }
};


/*
 min c^tx
     x_i(1-x_i) = 0, forall i=1,...,n
     Ax <= b
     x > 0
*/
class IntegerPackingJRF : public SimpleJRF {
 public:
     IntegerPackingJRF(const SpMat& A,
                       const Vector& b,
                       const Vector& c,
                       Vector& x,
                       Vector& y) : SimpleJRF(A,b,c,x,y){ _m += _n; }

     bool getBoundsConstraints(Scalar* cl, Scalar* cu) override
     {
         Vector::MapType cLower = Vector::MapType(cl, _m);
         Vector::MapType cUpper = Vector::MapType(cu, _m);
         cLower.head(_n) = Vector::Zero(_n);
         cUpper.head(_n) = Vector::Zero(_n);
         cUpper.tail(_A.rows()) = _b;
         cLower.tail(_A.rows()) = Vector::Constant(_A.rows(), -INF);
         return true;
     }
     bool functionValueConstraints(const Scalar* variablesPtr,
                                   Scalar* constraintValuesPtr) override
     {
          Vector::ConstMapType x = Vector::ConstMapType(variablesPtr, _n);
          Vector::MapType constraintValues = Vector::MapType(constraintValuesPtr, _m);
          Vector t = Vector::Ones(_n);
          constraintValues.head(_n) = x.cwiseProduct(x-t);
          constraintValues.tail(_A.rows()) = _A * x;
          return true;
      }
      bool gradientConstraintsTimesVector(const Scalar* variablesPtr,
                                          const Scalar* dualVariablesPtr,
                                          Scalar* gradientPtr) override
      {
           Vector::ConstMapType x = Vector::ConstMapType(variablesPtr, _n);
           Vector::ConstMapType y = Vector::ConstMapType(dualVariablesPtr, _m);
           Vector::MapType gradient = Vector::MapType(gradientPtr, _n);
           gradient = (2*x).cwiseProduct(y.head(_n)) - y.head(_n);
           gradient += _A.transpose() * y.tail(_A.rows());
           return true;
       }

};

class BranchingJRF : public SimpleJRF
{
     public:

    BranchingJRF(const SpMat& A, const Vector& b, const Vector& c, Vector& x, Vector& y) : 
    SimpleJRF(A, b, c, x, y) 
    { 
        lo.conservativeResizeLike(Vector::Zero(_n));
        hi.conservativeResizeLike(Vector::Ones(_n));
    }

    virtual bool getBounds(Scalar* lb, Scalar* ub) override
    {
        Vector::MapType(lb, _n) = lo;
        Vector::MapType(ub, _n) = hi;
        return true;
    }

    bool getBoundsConstraints(Scalar* cl, Scalar* cu) override
    {
        Vector::MapType(cu, _m) = _b;
        Vector::MapType(cl, _m) = Vector::Constant(_m, -INF);
        return true;
    }

    Vector hi, lo;
};






#endif
