#include <iostream>
#include <vector>
#include <solver.h>
#include <omp.h> 
#include "vtk_writer.h"

int main(){
    const int N = 100;
    const int ibd = 4;
    const int N_total = N + 2*ibd;
    const double x0 = 0.0, x1 = 1.0;
    const double dx = (x1 - x0) / N;
    const double CFL = 0.2;
    const double T_end = 1.0;
    

    std::vector<double> u(N+2*ibd,0.0);

    //Initial condition: sin wave function
    burgers_initialize(u,N,x0,x1,dx,ibd);

    //Time step
    const double dt = CFL * dx;
    const int n_steps = T_end / dt;

    //Solve the equation
    simulate_burgers1d(u,dx,CFL,T_end,ibd);

    return 0;
}