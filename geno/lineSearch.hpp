class TMPContainer
{
  public:
  bool brackt;
  int stage;
  double finit, fx, fy, ginit, gtest, gx, gy, stx, sty, stmin, stmax, width, width1;
};

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
            double const stpmax);
            
enum TaskType {START, ERROR, FG, WARNING, CONVERGENCE};

void dcsrch(double & f,
            double & g,
            double & stp,
            double const ftol,
            double const gtol,
            double const xtol,
            double const stpmin,
            double const stpmax,
            TaskType & task, TMPContainer& TCON);
