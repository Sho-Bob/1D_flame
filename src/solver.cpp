#include "solver.h"
// #include "limiter.h"
#include "vtk_writer.h"
#include <vector>
#include <cmath>
#include <algorithm>

void simulate_burgers1d(std::vector<double>& u, double dx, double CFL, double T_end){
    const int N = u.size();
    int n_steps = 0;
    std::vector<double> u_new(N,0.0);
    std::vector<double> u_flux(N+1,0.0);
    std::vector<double> uL(N+1,0.0), uR(N+1,0.0);

    double t = 0.0;

    while(t < T_end){

        n_steps++;
        //Compute dt from CFL condition
        double max_speed = 0.0;
        for (int i = 0;i < N;i++){
            max_speed = std::max(max_speed,std::abs(u[i]));
        }
        double dt = CFL * dx / max_speed;

        //Reconstruct left and right states (1st order)
        for (int i=1;i < N;i++){
        // MUSCL for future work
            uL[i] = u[i-1];
            uR[i] = u[i];
        }

        // // Boundary conditions (periodic)
        //     uL[0] = u[N-1];
        //     uR[0] = u[0];
        //     uL[N] = u[N-1];
        //     uR[N] = u[0];

        // Boundary conditions (periodic)
            uL[0] = u[0];
            uR[0] = u[0];
            uL[N] = u[N-1];
            uR[N] = u[N-1];
    
        //Compute fluxes
        for (int i=0; i<N+1; i++){
            double qm = 0.5 * (uL[i] + std::abs(uL[i]));
            double qp = 0.5 * (uR[i] - std::abs(uR[i]));
            u_flux[i] = std::max(0.5 * qm*qm, 0.5 * qp*qp);
        }

        // Update solution
        for (int i=0; i<N; i++){
            u_new[i] = u[i] - dt/dx * (u_flux[i+1]-u_flux[i]);
        }

        // Update solution
        for (int i=0; i<N; i++){
            u[i] = u_new[i];
        }

        //Write solution to VTK file
        if(n_steps % 10 == 0){
            write_vtk(u,dx,"burgers1d_step_" + std::to_string(n_steps) + ".vtk");
        }

        t += dt;
    }
    
}