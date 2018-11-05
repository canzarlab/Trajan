#ifndef LPCP_H
#define LPCP_H

#include "LPInt.h"

class LPCP : public LPInt
{
public:
    LPCP(Graph& t1, Graph& t2, string d, double k, bool dag);

private:
    bool SolveLP() override;
};

#endif
