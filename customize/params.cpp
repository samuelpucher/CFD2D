////////////////////////////////////////////////////////
//
//  Author:     Samuel
//  Date:       6th April 2025
//  Desc:       This file contains all 
//              relevant numerical parameters
//              for the actual simulation
//
////////////////////////////////////////////////////////

#include <cmath>

// Physical constants
double NU = 0.005;
double AMP = 0.05;

// Simulation constants
double LENGTH_X = 4 * M_PI;
double LENGTH_Y = 4 * M_PI;
double T_TOTAL = 100.0;

// Numerical constants
double T_STEP = 0.001;
int GRID_X = 256;
int GRID_Y = 256;

// Helping varibles
double DELTA_X = LENGTH_X/GRID_X;
double DELTA_Y = LENGTH_Y/GRID_Y;
double STEPS_TOTAL = T_TOTAL/T_STEP;