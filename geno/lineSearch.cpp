#include "lineSearch.hpp"

#include <cmath>
#include <iostream>

void dcstep(double & stx,
            double & fx,
            double & dx,
            double & sty,
            double & fy,
            double & dy,
            double & stp,
            double fp,
            double dp,
            bool & brackt,
            double const stpmin,
            double const stpmax) {
  double const zero  = 0.0,
               p66   = 0.66,
               two   = 2.0,
               three = 3.0;
  double gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta;
  sgnd = dp * (dx / std::fabs(dx));
  if (fp > fx) {
    //    std::cout << "dcstep 1" << std::endl;
     theta = three * (fx - fp) / (stp - stx) + dx + dp;
     s = std::max(std::fabs(theta), std::max(std::fabs(dx), std::fabs(dp)));
     gamma = s * std::sqrt(std::pow(theta / s, 2) - (dx / s) * (dp / s));
     if (stp < stx) gamma = -gamma;
     p = (gamma - dx) + theta;
     q = ((gamma - dx) + gamma) + dp;
     r = p / q;
     stpc = stx + r*(stp - stx);
     stpq = stx + ((dx/((fx - fp)/(stp - stx) + dx))/two) * (stp - stx);
     if (std::fabs(stpc-stx) < std::fabs(stpq-stx)) {
        stpf = stpc;
     } else {
        stpf = stpc + (stpq - stpc)/two;
     }
     brackt = true;
  } else if (sgnd < zero) {
    //    std::cout << "dcstep 2" << std::endl;
     theta = three*(fx - fp)/(stp - stx) + dx + dp;
     s = std::max(std::fabs(theta), std::max(std::fabs(dx),std::fabs(dp)));
     gamma = s*std::sqrt(std::pow(theta/s, 2) - (dx/s)*(dp/s));
     if (stp > stx) gamma = -gamma;
     p = (gamma - dp) + theta;
     q = ((gamma - dp) + gamma) + dx;
     r = p/q;
     stpc = stp + r*(stx - stp);
     stpq = stp + (dp/(dp - dx))*(stx - stp);
     if (std::fabs(stpc-stp) > std::fabs(stpq-stp)) {
        stpf = stpc;
     } else {
        stpf = stpq;
     }
     brackt = true;
  } else if (std::fabs(dp) < std::fabs(dx)) {
    //    std::cout << "dcstep 3 ***************" << std::endl;
     theta = three*(fx - fp)/(stp - stx) + dx + dp;
     s = std::max(std::fabs(theta), std::max(std::fabs(dx),std::fabs(dp)));

     gamma = s*std::sqrt(std::max(zero,std::pow(theta/s, 2)-(dx/s)*(dp/s)));
     //     std::cout << "theta = " << theta << std::endl;
     //     std::cout << "gamma = " << gamma << std::endl;
     if (stp > stx) gamma = -gamma;
     p = (gamma - dp) + theta;
     q = (gamma + (dx - dp)) + gamma;
     r = p/q;
     if (r < zero && gamma != zero) {
        stpc = stp + r*(stx - stp);
     } else if (stp > stx) {
        stpc = stpmax;
     } else {
        stpc = stpmin;
     }
     //     std::cout << "stpc = " << stpc << std::endl;
     
     stpq = stp + (dp/(dp - dx))*(stx - stp);

     //     std::cout << "--------------------------------------------" << std::endl;
     //     std::cout << "stpq = " << stpq << std::endl;
     //     std::cout << "stp  = " << stp << std::endl;
     //     std::cout << "dp   = " << dp << std::endl;
     //     std::cout << "dx   = " << dx << std::endl;
     //     std::cout << "stx  = " << stx << std::endl;
     //     std::cout << "dp/(dp-dx) = " << dp/(dp-dx) << std::endl;
     //     std::cout << "--------------------------------------------" << std::endl;



     if (brackt) {
        if (std::fabs(stpc-stp) < std::fabs(stpq-stp)) {
           stpf = stpc;
        } else {
           stpf = stpq;
        }
        if (stp > stx) {
           stpf = std::min(stp+p66*(sty-stp),stpf);
        } else {
           stpf = std::max(stp+p66*(sty-stp),stpf);
        }
     } else {
        if (std::fabs(stpc-stp) > std::fabs(stpq-stp)) {
           stpf = stpc;
        } else {
           stpf = stpq;
        }
        stpf = std::min(stpmax,stpf);
        stpf = std::max(stpmin,stpf);
     }
  } else { 
    //    std::cout << "dcstep 4" << std::endl;
    if (brackt) {

        theta = three*(fp - fy)/(sty - stp) + dy + dp;
        s = std::max(std::fabs(theta), std::max(std::fabs(dy),std::fabs(dp)));
        gamma = s*std::sqrt(std::pow(theta/s, 2) - (dy/s)*(dp/s));
        if (stp > sty) gamma = -gamma;
        p = (gamma - dp) + theta;
        q = ((gamma - dp) + gamma) + dy;
        r = p/q;
        stpc = stp + r*(sty - stp);
        stpf = stpc;
     } else if (stp > stx) {
        stpf = stpmax;
     } else {
        stpf = stpmin;
     }
  }

  if (fp > fx) {
     sty = stp;
     fy = fp;
     dy = dp;
  } else {
     if (sgnd < zero) {
        sty = stx;
        fy = fx;
        dy = dx;
     }
     stx = stp;
     fx = fp;
     dx = dp;
  }
  stp = stpf;
}

void dcsrch(double & f,
            double & g,
            double & stp,
            double const ftol,
            double const gtol,
            double const xtol,
            double const stpmin,
            double const stpmax,
            TaskType & task, TMPContainer& TCON) {
  double const zero = 0.0,
               p5   = 0.5,
               p66  = 0.66,
               xtrapl = 1.1,
               xtrapu = 4.0;
  //static bool brackt;
  //static int stage;
  //static double finit, fx, fy, ginit,gtest, gx, gy,
  //              stx,sty,stmin,stmax,width,width1;
  double ftest,fm,fxm,fym, gm,gxm,gym;

      if (task == TaskType::START) {

//       Check the input arguments for errors.

         if (stp < stpmin) {
          task = TaskType::ERROR;
          std::cerr << "ERROR: STP < STPMIN" << std::endl;
         }
         if (stp > stpmax) {
          //task = TaskType::ERROR;
          stp = stpmax;
          std::cerr << "ERROR: STP > STPMAX" << std::endl;
         } 
         if (g >= zero) {
          task = TaskType::ERROR;
          std::cerr << "ERROR: INITIAL G >= ZERO" << std::endl;
         }
         if (ftol < zero) {
          task = TaskType::ERROR;
          std::cerr << "ERROR: FTOL < ZERO" << std::endl;
         }
         if (gtol < zero) {
          task = TaskType::ERROR;
          std::cerr << "ERROR: GTOL < ZERO" << std::endl;
         }
         if (xtol < zero) {
          task = TaskType::ERROR;
          std::cerr << "ERROR: XTOL < ZERO" << std::endl;
         }
         if (stpmin < zero) {
          task = TaskType::ERROR;
          std::cerr << "ERROR: STPMIN < ZERO" << std::endl;
         }
         if (stpmax < stpmin) {
          task = TaskType::ERROR;
          std::cerr << "ERROR: STPMAX < STPMIN" << std::endl;
         }

//        Exit if there are errors on input.

         if (task == TaskType::ERROR)
          return;

//        Initialize local variables.

         TCON.brackt = false;
         TCON.stage = 1;
         TCON.finit = f;
         TCON.ginit = g;
         TCON.gtest = ftol*TCON.fy;
         TCON.width = stpmax - stpmin;
         TCON.width1 = TCON.width/p5;

//        The variables stx, fx, gx contain the values of the step, 
//        function, and derivative at the best step. 
//        The variables sty, fy, gy contain the value of the step, 
//        function, and derivative at sty.
//        The variables stp, f, g contain the values of the step, 
//        function, and derivative at stp.

         TCON.stx = zero;
         TCON.fx = TCON.finit;
         TCON.gx = TCON.ginit;
         TCON.sty = zero;
         TCON.fy = TCON.finit;
         TCON.gy = TCON.ginit;
         TCON.stmin = zero;
         TCON.stmax = stp + xtrapu*stp;
         task = TaskType::FG;

         return;

      }

//     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
//     algorithm enters the second stage.

      ftest = TCON.finit + stp*TCON.gtest;
      if (TCON.stage == 1 && f <= ftest && g >= zero) 
        TCON.stage = 2;

//     Test for warnings.

      if (TCON.brackt && (stp <= TCON.stmin || stp >= TCON.stmax)) {
        task = TaskType::WARNING;
//        std::cerr << "ROUNDING ERRORS PREVENT PROGRESS\n";
      }
      if (TCON.brackt && TCON.stmax - TCON.stmin <= xtol*TCON.stmax) {
        task = TaskType::WARNING;
//        std::cerr << "XTOL TEST SATISFIED\n";
      }
      if (stp == stpmax && f <= ftest && g <= TCON.gtest) {
        task = TaskType::WARNING;
//        std::cerr << "STP = STPMAX\n";
      }
      if (stp == stpmin && (f > ftest || g >= TCON.gtest)) {
        task = TaskType::WARNING;
//        std::cerr << "STP = STPMIN\n";
      }
//     Test for convergence.

/*
      printf("f = %15.15g \n", f);
      printf("ftest = %15.15g\n", ftest);
      std::cout << "f " << f << std::endl;
      std::cout << "ftest " << ftest << std::endl;
      std::cout << "fabs(g) " << std::fabs(g) << std::endl;
      std::cout << "gtol*ginit " << gtol * (-ginit) << std::endl;
*/
//      if (f <= ftest && std::fabs(g) <= gtol * (-TCON.ginit))
      double eps = 1E-6;
      if (f <= ftest + eps*(std::fabs(ftest) + 1) && std::fabs(g) <= gtol * (-TCON.ginit)) 
        task = TaskType::CONVERGENCE;

//     Test for termination.

      if (task == TaskType::WARNING || task == TaskType::CONVERGENCE)
        return;

//     A modified function is used to predict the step during the
//     first TCON.stage if a lower function value has been obtained but 
//     the decrease is not sufficient.

      if (TCON.stage == 1 && f <= TCON.fx && f > ftest) {
//        Define the modified function and derivative values.
         fm = f - stp*TCON.gtest;
         fxm = TCON.fx - TCON.stx*TCON.gtest;
         fym = TCON.fy - TCON.sty*TCON.gtest;
         gm = g - TCON.gtest;
         gxm = TCON.gx - TCON.gtest;
         gym = TCON.gy - TCON.gtest;
//        Call dcstep to update stx, sty, and to compute the new step.
//	 std::cout << "called dcstep 1" << std::endl;
         dcstep(TCON.stx,fxm,gxm,TCON.sty,fym,gym,stp,fm,gm,TCON.brackt,TCON.stmin,TCON.stmax);
//        Reset the function and derivative values for f.
         TCON.fx = fxm + TCON.stx*TCON.gtest;
         TCON.fy = fym + TCON.sty*TCON.gtest;
         TCON.gx = gxm + TCON.gtest;
         TCON.gy = gym + TCON.gtest;
      } else {

//       Call dcstep to update stx, sty, and to compute the new step.
//	std::cout << "called dcstep 2" << std::endl;
        dcstep(TCON.stx,TCON.fx,TCON.gx,TCON.sty,TCON.fy,TCON.gy,stp,f,g,TCON.brackt,TCON.stmin,TCON.stmax);
	//	std::cout << "stx = " << stx << std::endl;
	//	std::cout << "fx  = " << fx << std::endl;
	//	std::cout << "gx  = " << gx << std::endl;
	//	std::cout << "sty = " << sty << std::endl;
	//	std::cout << "fy  = " << fy << std::endl;
	//	std::cout << "gy  = " << gy << std::endl;
	//	std::cout << "stp = " << stp << std::endl;
	//	std::cout << "f   = " << f << std::endl;
	//	std::cout << "g   = " << g << std::endl;
	//	std::cout << "stmin = " << stmin << std::endl;
	//	std::cout << "stmax = " << stmax << std::endl;
      }
//     Decide if a bisection step is needed.

      if (TCON.brackt) {
	//	std::cout << "brackt true" << std::endl;
         if (std::fabs(TCON.sty-TCON.stx) >= p66*TCON.width1) stp = TCON.stx + p5*(TCON.sty - TCON.stx);
         TCON.width1 = TCON.width;
         TCON.width = std::fabs(TCON.sty-TCON.stx);
      }

//     Set the minimum and maximum steps allowed for stp.

      if (TCON.brackt) {
         TCON.stmin = std::min(TCON.stx,TCON.sty);
         TCON.stmax = std::max(TCON.stx,TCON.sty);
      } else {
         TCON.stmin = stp + xtrapl*(stp - TCON.stx);
         TCON.stmax = stp + xtrapu*(stp - TCON.stx);
      }
// 
//     Force the step to be within the bounds stpmax and stpmin.
// 
      stp = std::max(stp,stpmin);
      stp = std::min(stp,stpmax);

//     If further progress is not possible, let stp be the best
//     point obtained during the search.

      if ((TCON.brackt && (stp <= TCON.stmin || stp >= TCON.stmax))
          || (TCON.brackt && TCON.stmax-TCON.stmin <= xtol*TCON.stmax)) {
	//	std::cout << "no further progress" << std::endl;
	stp = TCON.stx;
      }

//     Obtain another function and derivative.

      task = TaskType::FG;
}

