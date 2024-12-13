#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdlib.h>
#include <stdio.h>
#include <cstring>

#include "shared.h"
#include "exact-method.h"

MethodResoltion methodResoltion = GUROBI;
char sourceFile[] = "instances/inst5.txt";
char targetFile[] = "interactive-optimizers/inst5.lp";
int hubQuantity = 1;
int showLogs = 0;
double alpha = 0.75;
InstanceEntries instanceEntries;
double** costMatrix;
double** costAggMatrix;

int main()
{
    switch (methodResoltion)
    {
    case CPLEX:
        cplexSolver();
        break;
    case GUROBI:
        gurobiSolver();
        break;
    }

    freeCoordinatesInstance();
    freeCostMatriz();
    freeCostAggMatriz();

    return 0;
}