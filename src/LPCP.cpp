#include "LPCP.h"
#include "Timer.h"
#include <iostream>

LPCP::LPCP(Graph& t1, Graph& t2, string d, double k, bool dag) : LPInt(t1, t2, d, k, dag)
{
}

bool LPCP::SolveLP()
{
    clog << "nr_rows = " << nr_rows << " and nr_cols = " << nr_cols << endl;
    warm_y.conservativeResizeLike(Vector::Zero(nr_rows)); // resizes y with 0's, but keeping old values intact.

    SpMat A(nr_rows, nr_cols);
    A.setFromTriplets(Triplets.begin(), Triplets.end());
    SpMat A_t = A.transpose();
    Vector b = Vector::Ones(nr_rows);

    //x = Vector::Zero(nr_cols);
    //y = Vector::Zero(nr_rows);

    CoveringJRF simpleJRF(A_t, c, b, warm_y, warm_x);

    Vector c1 = -c;
//        PackingJRF simpleJRF(A, b, c1, x, y);
    AugmentedLagrangian solver(simpleJRF, 15);
    solver.setParameter("verbose", false);
    solver.setParameter("pgtol", 1e-1); // should influence running time a lot
    solver.setParameter("constraintsTol", 1e-3);
    Timer timeGeno;
    timeGeno.start();
    solver.solve();
    timeGeno.stop();

    clog << "f = " << solver.f() << " computed in time: " << timeGeno.secs() << " secs" << endl;

    /* when Packing: x->x, y->y, when Covering: -y -> x, x -> y */
    x = -Vector::ConstMapType(solver.y(), nr_cols);
    y = Vector::ConstMapType(solver.x(), nr_rows);

    warm_x = Vector::ConstMapType(solver.y(), nr_cols);
    warm_y = Vector::ConstMapType(solver.x(), nr_rows);


    Vector t = A_t*y-c;
    int nr_tight_constr =  nr_cols - (t.array() > 0.1).count();
    clog << "Number of tight constraints in the dual: " << nr_tight_constr << endl;
    //idea, truncate matrix A for columns that correspond to non-tight constraints in dual
    SpMat truncA(nr_rows, nr_tight_constr);
    Vector truncc(nr_tight_constr);
    Vector truncx(nr_tight_constr);
    int truncA_col = 0;
    for (int i=0; i<nr_cols; i++)
        if (t(i)> 0.1)
            x(i) = 0.0;
        else {
            for (SpMat::InnerIterator it(A,i); it; ++it)
            {
                truncA.coeffRef(it.row(), truncA_col) = it.value();
            }
            truncc(truncA_col) = - c(i);
            truncx(truncA_col) = x(i);
            truncA_col ++;
        }

    assert(truncA_col == nr_tight_constr);

    clog << "Truncated matrix formed ... resolve" << endl;
    Vector expy = Vector::Zero(truncA.rows() + truncA.cols());
    IntegerPackingJRF simpleJRF1(truncA, b, truncc, truncx, expy);
    //PackingJRF simpleJRF1(truncA, b, truncc, truncx, y);
    AugmentedLagrangian solver1(simpleJRF1, 15);
    solver1.setParameter("verbose", false);
    solver1.setParameter("pgtol", 1e-1); // should influence running time a lot
    solver1.setParameter("constraintsTol", 1e-5);
    Timer timeGeno1;
    timeGeno1.start();
    solver1.solve();
    timeGeno1.stop();
    clog << "trunc f = " << solver1.f() << " computed in time: " << timeGeno1.secs() << " secs" << endl;

    // map the solution back to vector x
    truncx = Vector::ConstMapType(solver1.x(), nr_tight_constr);
    clog <<"MAX INTEGER PACKING VALUE: " <<  truncx.maxCoeff() << endl;
    truncA_col = 0;
    for (int i=0; i<nr_cols; i++)
        if (t(i)<= 0.1){
            x(i) = truncx(truncA_col);
            truncA_col++;
        }

    assert(truncA_col == nr_tight_constr);
    return true;
}
