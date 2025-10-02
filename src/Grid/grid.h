#ifndef GRID_H
#define GRID_H

class Grid {
  public:
    Grid(int nx, int num_boundary_points, double x0, double x1);
    ~Grid();

    double** get_coordinates() { return coordinates; }
    int* get_n_ptr() { return n; }
    int* get_n_tot_ptr() { return n_tot; }
    double* get_Delta_ptr() { return Delta; }
    void get_n(int* n_out) { n_out = this->n; }
    void get_n_tot(int* n_tot_out) { n_tot_out = this->n_tot; }
    int get_num_dim() { return num_dim; }
    int get_num_boundary_points() { return num_boundary_points; }

    int get_icv0(int ifa) { return this->num_boundary_points + ifa - 1; }
    int get_icv1(int ifa) { return this->num_boundary_points + ifa ; }

  private:
    int num_dim; // number of dimensions
    int* n; // number of grid points in each dimension EXCLUDING ghost cells
    int* n_tot; // number of grid points in each dimension INLCUDING ghost cells
    int num_boundary_points; // number of ghost cells, in each direction and on each side
    // --------------------------------------------------
    // 1D illustration of the grid:
    // o---x---o---x---o---x---o---x---o---x---o---x---o
    // |       |       |       |       |       |       |
    // x0     x1      x2      x3      x4      x5     x6
    // --------------------------------------------------
    // where 'o' are the cell centers and 'x' are the cell faces (interfaces)
    // nx is the number of cells (centers); in this case nx = 5, and nbd = 1
    // therefore, x0 and x6 are ghost cells
    // --------------------------------------------------
    // The total number of grid points including ghost cells is:
    // N_total = nx + 2 * num_boundary_points
    // --------------------------------------------------
    double x0, x1; // domain [x0,x1]
    double xmin, xmax; // domain [xmin,xmax] excluding ghost cells
    double dx; // grid spacing
    double* Delta; // grid spacing in each dimension
    double** coordinates; // coordinates of the grid points
};

#endif // GRID_H
