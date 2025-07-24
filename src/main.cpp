#include <iostream>
#include <vector>
#include <solver.h>
#include "vtk_writer.h"

int main(){
    const int N = 100;
    const double x0 = 0.0, x1 = 1.0;
    const double dx = (x1 - x0) / N;
    const double CFL = 0.5;
    const double T_end = 1.0;

    std::vector<double> u(N,0.0);

    //Initial condition: sin wave function
    for(int i = 0; i < N; i++){
        u[i] = std::sin(2.0 * M_PI * (i * dx));
    }

    //Time step
    const double dt = CFL * dx;
    const int n_steps = T_end / dt;

    //Solve the equation
    simulate_burgers1d(u,dx,CFL,T_end);

    return 0;
}