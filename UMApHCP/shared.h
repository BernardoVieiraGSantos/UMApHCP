#ifndef SHARED_H
#define SHERED_H
#define MAX_PATH_LENGTH 100

typedef struct
{
    double x;
    double y;
} Coordinate;

typedef struct
{
    int nodeQuantity;
    Coordinate* coordinates;
} InstanceEntries;

enum MethodResoltion {
    CPLEX,
    GUROBI
};

extern MethodResoltion methodResoltion;
extern char sourceFile[];
extern char targetFile[];
extern int hubQuantity;
extern int showLogs;
extern double alpha;
extern InstanceEntries instanceEntries;
extern double** costMatrix;
extern double** costAggMatrix;

int readInstance();
int freeCoordinatesInstance();
int generateCostMatriz();
int freeCostMatriz();
int generateCostAggMatriz();
int freeCostAggMatriz();

#endif