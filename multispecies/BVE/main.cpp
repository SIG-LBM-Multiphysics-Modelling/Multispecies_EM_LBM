// 2-LBM solver with barycentric velocity equilibrium formulation
#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "VTIWriter.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// 1m x 1m physical domain
const double Lx = 1.0, Ly = Lx;
const int NX = 256, NY = Ly / Lx * NX;
const int Q = 9;

// --------------------- Fields ---------------------
Vec2D phi(NX, Vec1D(NY, 0.));
Vec2D Bz(NX, Vec1D(NY, 0.));
Vec2D q_rho(NX, Vec1D(NY, 0.));
Vec2D Jx(NX, Vec1D(NY, 0.));
Vec2D Jy(NX, Vec1D(NY, 0.));
Vec2D Ex(NX, Vec1D(NY, 0.));
Vec2D Ey(NX, Vec1D(NY, 0.));
Vec2D rhs_phi(NX, Vec1D(NY, 0.));
Vec2D rhs_Bz(NX, Vec1D(NY, 0.));

// --------------------- Fluid ---------------------
const double cs2 = 1. / 3.; // lattice speed of sound squared
const int cx[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int cy[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
const double wf[Q] = {4/9., 1/9., 1/9., 1/9., 1/9., 
                        1/36., 1/36., 1/36., 1/36.};

const double dx_phys = Lx / (double)NX;
const double dt_phys = dx_phys / 1.0;  // lattice speed c = 1.0
const double nu_phys = 0.01 / 256.;
const double nu_lb = nu_phys * dt_phys / (dx_phys * dx_phys);

const double tau = 0.5 + nu_lb / cs2;
const double omega = 1. / tau, omega1 = 1. - omega;

// Helper to generate filename inside "VTI_Poisson"
std::string getVTIFileName(int timestep) {
    std::ostringstream oss;
    oss << "VTI_Poisson/poisson_t_" 
        << std::setw(6) << std::setfill('0') << timestep 
        << ".vti";
    return oss.str();
}

void initialise(Vec2D &rho_m, Vec2D &ux_m, Vec2D &uy_m, Vec3D &f1_m, Vec3D &f2_m, 
            Vec2D &rho_p, Vec2D &ux_p, Vec2D &uy_p, Vec3D &f1_p, Vec3D &f2_p, 
            Vec2D &q_rho, Vec2D &Jx, Vec2D &Jy, Vec2D &Ex, Vec2D &Ey, Vec2D &Bz) {

    const double E0 = 1.e-4;

    for (int i = 0; i < NX; i++) {
        double x = double(i) / double(NX);
        for (int j = 0; j < NY; j++) {
            double y = double(j) / double(NY);

            // --- Densities ---
            rho_p[i][j] = 1.0;
            rho_m[i][j] = 1.0;

            // --- Cations velocities ---
            ux_p[i][j] = 0.;
            uy_p[i][j] = 0.;

            // --- Anions velocities ---
            ux_m[i][j] = 0.;
            uy_m[i][j] = 0.;

            // --- Electromagnetic fields ---
            Ex[i][j] = E0 * std::sin(2.0 * M_PI * y);
            Ey[i][j] = 0.0;
            Bz[i][j] = 0.0;

            // --- Charge and current densities ---
            q_rho[i][j] = rho_p[i][j] - rho_m[i][j];
            Jx[i][j] = rho_p[i][j] * ux_p[i][j] - rho_m[i][j] * ux_m[i][j];
            Jy[i][j] = rho_p[i][j] * uy_p[i][j] - rho_m[i][j] * uy_m[i][j];

            // Initialise distribution functions
            for (int k = 0; k < Q; k++) {
                // Anions
                double first_order = ux_m[i][j] * cx[k] + uy_m[i][j] * cy[k];
                f1_m[i][j][k] = f2_m[i][j][k] = wf[k] * rho_m[i][j] * (1. + 3. * first_order +
                    4.5 * first_order * first_order - 1.5 * (ux_m[i][j] * ux_m[i][j] + uy_m[i][j] * uy_m[i][j]));

                // Cations
                first_order = ux_p[i][j] * cx[k] + uy_p[i][j] * cy[k];
                f1_p[i][j][k] = f2_p[i][j][k] = wf[k] * rho_p[i][j] * (1. + 3. * first_order +
                    4.5 * first_order * first_order - 1.5 * (ux_p[i][j] * ux_p[i][j] + uy_p[i][j] * uy_p[i][j]));
            }
        }
    }
}

int algoLB (Vec2D &rho_m, Vec2D &ux_m, Vec2D &uy_m, Vec3D &f1_m, Vec3D &f2_m, 
            Vec2D &rho_p, Vec2D &ux_p, Vec2D &uy_p, Vec3D &f1_p, Vec3D &f2_p, 
            Vec2D &q_rho, Vec2D &Jx, Vec2D &Jy, const Vec2D &Ex, const Vec2D &Ey, const Vec2D &Bz) {

    double Rm, Um, Vm, Rp, Up, Vp;
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {

            Rm = Um = Vm = Rp = Up = Vp = 0.;
            for (int k = 0; k < Q; k++) {
                Rm += f1_m[i][j][k];
                Um += f1_m[i][j][k] * cx[k];
                Vm += f1_m[i][j][k] * cy[k];
                Rp += f1_p[i][j][k];
                Up += f1_p[i][j][k] * cx[k];
                Vp += f1_p[i][j][k] * cy[k];
            }

            double Fx_m = 0., Fy_m = 0.;
            double Fx_p = 0., Fy_p = 0.;
            Fx_m = -(Ex[i][j] + uy_m[i][j] * Bz[i][j]);
            Fy_m = -(Ey[i][j] - ux_m[i][j] * Bz[i][j]);

            Fx_p = (Ex[i][j] + uy_p[i][j] * Bz[i][j]);
            Fy_p = (Ey[i][j] - ux_p[i][j] * Bz[i][j]);

            Um = (Um + 0.5 * Fx_m) / Rm;
            Vm = (Vm + 0.5 * Fy_m) / Rm;
            rho_m[i][j] = Rm;
            ux_m[i][j] = Um;
            uy_m[i][j] = Vm;

            Up = (Up + 0.5 * Fx_p) / Rp;
            Vp = (Vp + 0.5 * Fy_p) / Rp;
            rho_p[i][j] = Rp;
            ux_p[i][j] = Up;
            uy_p[i][j] = Vp;

            double Ux_bary = (Rm * Um + Rp * Up) / (Rm + Rp);
            double Uy_bary = (Rm * Vm + Rp * Vp) / (Rm + Rp);

            Um = Ux_bary;
            Vm = Uy_bary;
            Up = Ux_bary;
            Vp = Uy_bary;

            q_rho[i][j] = rho_p[i][j] - rho_m[i][j];
            Jx[i][j] = rho_p[i][j] * ux_p[i][j] - rho_m[i][j] * ux_m[i][j];
            Jy[i][j] = rho_p[i][j] * uy_p[i][j] - rho_m[i][j] * uy_m[i][j];

            double Um2 = Um * Um,
                   Up2 = Up * Up;
            double Vm2 = Vm * Vm,
                   Vp2 = Vp * Vp;
            double UVm = Um * Vm,
                   UVp = Up * Vp;

            // Central and raw moments
            double k1_m, k2_m, k3_m, k4_m, k5_m, k6_m, k7_m, k8_m;
            double k1_p, k2_p, k3_p, k4_p, k5_p, k6_p, k7_p, k8_p;
            double r1_m, r2_m, r3_m, r4_m, r5_m, r6_m, r7_m, r8_m;
            double r1_p, r2_p, r3_p, r4_p, r5_p, r6_p, r7_p, r8_p;

            // Pre-collision moments transformed from populations
            r4_m = f1_m[i][j][1] - f1_m[i][j][2] + f1_m[i][j][3] - f1_m[i][j][4];
            r5_m = f1_m[i][j][5] - f1_m[i][j][6] + f1_m[i][j][7] - f1_m[i][j][8];
            k4_m = r4_m - Rm * (Um2 - Vm2);
            k5_m = r5_m - Rm * UVm;
            k1_m = 0.5 * Fx_m;
            k2_m = 0.5 * Fy_m;
            k3_m = 2. * Rm * cs2;
            k4_m = omega1 * k4_m;
            k5_m = omega1 * k5_m;
            k6_m = 0.5 * Fy_m * cs2;
            k7_m = 0.5 * Fx_m * cs2;
            k8_m = Rm * cs2 * cs2;
            r1_m = k1_m + Rm * Um;
            r2_m = k2_m + Rm * Vm;
            r3_m = k3_m + 2. * Um * k1_m + 2. * Vm * k2_m + Rm * (Um2 + Vm2);
            r4_m = k4_m + 2. * Um * k1_m - 2. * Vm * k2_m + Rm * (Um2 - Vm2);
            r5_m = k5_m + Um * k2_m + Vm * k1_m + Rm * (UVm);
            r6_m = k6_m + 2. * Um * k5_m + 0.5 * Vm * (k3_m + k4_m) + Um2 * k2_m + 2. * UVm * k1_m + Rm * Um2 * Vm;
            r7_m = k7_m + 0.5 * Um * (k3_m - k4_m) + 2. * Vm * k5_m + Vm2 * k1_m + 2. * UVm * k2_m + Rm * Um * Vm2;
            r8_m = k8_m + 2. * Um * k7_m + 2. * Vm * k6_m + 0.5 * k3_m * (Um2 + Vm2) - 0.5 * k4_m * (Um2 - Vm2) + 
                    Rm * Um2 * Vm2 + 4. * UVm * k5_m + 2. * Um * Vm2 * k1_m + 2. * Um2 * Vm * k2_m;

            f1_m[i][j][0] = Rm - r3_m + r8_m;
            f1_m[i][j][1] = 0.25 * (r3_m + r4_m) + 0.5 * (r1_m - r7_m - r8_m);
            f1_m[i][j][2] = 0.25 * (r3_m - r4_m) + 0.5 * (r2_m - r6_m - r8_m);
            f1_m[i][j][3] = 0.25 * (r3_m + r4_m) + 0.5 * (-r1_m + r7_m - r8_m);
            f1_m[i][j][4] = 0.25 * (r3_m - r4_m) + 0.5 * (-r2_m + r6_m - r8_m);
            f1_m[i][j][5] = 0.25 * (r5_m + r6_m + r7_m + r8_m);
            f1_m[i][j][6] = 0.25 * (-r5_m + r6_m - r7_m + r8_m);
            f1_m[i][j][7] = 0.25 * (r5_m - r6_m - r7_m + r8_m);
            f1_m[i][j][8] = 0.25 * (-r5_m - r6_m + r7_m + r8_m);

            r4_p = f1_p[i][j][1] - f1_p[i][j][2] + f1_p[i][j][3] - f1_p[i][j][4];
            r5_p = f1_p[i][j][5] - f1_p[i][j][6] + f1_p[i][j][7] - f1_p[i][j][8];
            k4_p = r4_p - Rp * (Up2 - Vp2);
            k5_p = r5_p - Rp * UVp;
            k1_p = 0.5 * Fx_p;
            k2_p = 0.5 * Fy_p;
            k3_p = 2. * Rp * cs2;
            k4_p = omega1 * k4_p;
            k5_p = omega1 * k5_p;
            k6_p = 0.5 * Fy_p * cs2;
            k7_p = 0.5 * Fx_p * cs2;
            k8_p = Rp * cs2 * cs2;
            r1_p = k1_p + Rp * Up;
            r2_p = k2_p + Rp * Vp;
            r3_p = k3_p + 2. * Up * k1_p + 2. * Vp * k2_p + Rp * (Up2 + Vp2);
            r4_p = k4_p + 2. * Up * k1_p - 2. * Vp * k2_p + Rp * (Up2 - Vp2);
            r5_p = k5_p + Up * k2_p + Vp * k1_p + Rp * (UVp);
            r6_p = k6_p + 2. * Up * k5_p + 0.5 * Vp * (k3_p + k4_p) + Up2 * k2_p + 2. * UVp * k1_p + Rp * Up2 * Vp;
            r7_p = k7_p + 0.5 * Up * (k3_p - k4_p) + 2. * Vp * k5_p + Vp2 * k1_p + 2. * UVp * k2_p + Rp * Up * Vp2;
            r8_p = k8_p + 2. * Up * k7_p + 2. * Vp * k6_p + 0.5 * k3_p * (Up2 + Vp2) - 0.5 * k4_p * (Up2 - Vp2) + 
                    Rp * Up2 * Vp2 + 4. * UVp * k5_p + 2. * Up * Vp2 * k1_p + 2. * Up2 * Vp * k2_p;
            
            f1_p[i][j][0] = Rp - r3_p + r8_p;
            f1_p[i][j][1] = 0.25 * (r3_p + r4_p) + 0.5 * (r1_p - r7_p - r8_p);
            f1_p[i][j][2] = 0.25 * (r3_p - r4_p) + 0.5 * (r2_p - r6_p - r8_p);
            f1_p[i][j][3] = 0.25 * (r3_p + r4_p) + 0.5 * (-r1_p + r7_p - r8_p);
            f1_p[i][j][4] = 0.25 * (r3_p - r4_p) + 0.5 * (-r2_p + r6_p - r8_p);
            f1_p[i][j][5] = 0.25 * (r5_p + r6_p + r7_p + r8_p);
            f1_p[i][j][6] = 0.25 * (-r5_p + r6_p - r7_p + r8_p);
            f1_p[i][j][7] = 0.25 * (r5_p - r6_p - r7_p + r8_p);
            f1_p[i][j][8] = 0.25 * (-r5_p - r6_p + r7_p + r8_p);

            for (int k = 0; k < Q; k++) {
                int in = (i + cx[k] + NX) % NX;
                int jn = (j + cy[k] + NY) % NY;
                f2_m[in][jn][k] = f1_m[i][j][k];
                f2_p[in][jn][k] = f1_p[i][j][k];
            }

            if (Um > 0.1 || Up > 0.1 || Vm > 0.1 || Vp > 0.1) {
                return 1;
            }
        }
    }

    return 0;
}

void regularisedBC(Vec2D &rho_m, Vec2D &ux_m, Vec2D &uy_m, Vec3D &f1_m, Vec3D &f2_m, 
            Vec2D &rho_p, Vec2D &ux_p, Vec2D &uy_p, Vec3D &f1_p, Vec3D &f2_p) {

    // Top and Bottom walls
    // No-slip
    for (int i = 0; i < NX; i++) {

        double Rm, Um, Vm, Pxxm, Pxym, Pyym, fneq, QPI;
        double Rp, Up, Vp, Pxxp, Pxyp, Pyyp;

        // Bottom (y = 0)
        int j = 0;
        ux_m[i][j] = Um = 0.;
        uy_m[i][j] = Vm = 0.;
        ux_p[i][j] = Up = 0.;
        uy_p[i][j] = Vp = 0.;
        Rm = 1. / (1. - uy_m[i][j]) * (f1_m[i][j][0] + f1_m[i][j][1] + f1_m[i][j][3] +
                2. * (f1_m[i][j][4] + f1_m[i][j][7] + f1_m[i][j][8]));
        Rp = 1. / (1. - uy_p[i][j]) * (f1_p[i][j][0] + f1_p[i][j][1] + f1_p[i][j][3] +
                2. * (f1_p[i][j][4] + f1_p[i][j][7] + f1_p[i][j][8]));

        Pxxm = Pxym = Pyym = Pxxp = Pxyp = Pyyp = 0.;
        Vec1D Feqm(Q, 0.), Feqp(Q, 0.);
        for (int k = 0; k < Q; k++) {

            double A = Um * cx[k] + Vm * cy[k];
            Feqm[k] = wf[k] * Rm * (1. + 3. * A + 4.5 * A * A
                 - 1.5 * (Um * Um + Vm * Vm));
            fneq = f1_m[i][j][k] - Feqm[k];

            Pxxm += fneq * cx[k] * cx[k];
            Pxym += fneq * cx[k] * cy[k];
            Pyym += fneq * cy[k] * cy[k];

            A = Up * cx[k] + Vp * cy[k];
            Feqp[k] = wf[k] * Rp * (1. + 3. * A + 4.5 * A * A
                 - 1.5 * (Up * Up + Vp * Vp));
            fneq = f1_p[i][j][k] - Feqp[k];

            Pxxp += fneq * cx[k] * cx[k];
            Pxyp += fneq * cx[k] * cy[k];
            Pyyp += fneq * cy[k] * cy[k];
        }

        for (int k = 0; k < Q; k++) {
            QPI = Pxxm * (cx[k] * cx[k] - cs2) + 2. * Pxym * cx[k] * cy[k] + Pyym * (cy[k] * cy[k] - cs2);
            f2_m[i][j][k] = Feqm[k] + 4.5 * wf[k] * QPI;

            QPI = Pxxp * (cx[k] * cx[k] - cs2) + 2. * Pxyp * cx[k] * cy[k] + Pyyp * (cy[k] * cy[k] - cs2);
            f2_p[i][j][k] = Feqp[k] + 4.5 * wf[k] * QPI;
        }

        // Top (y = NY - 1)
        j = NY - 1;
        ux_m[i][j] = Um = 0.;
        uy_m[i][j] = Vm = 0.;
        ux_p[i][j] = Up = 0.;
        uy_p[i][j] = Vp = 0.;
        Rm = 1. / (1. + uy_m[i][j]) * (f1_m[i][j][0] + f1_m[i][j][1] + f1_m[i][j][3] +
                2. * (f1_m[i][j][2] + f1_m[i][j][5] + f1_m[i][j][6]));
        Rp = 1. / (1. + uy_p[i][j]) * (f1_p[i][j][0] + f1_p[i][j][1] + f1_p[i][j][3] +
                2. * (f1_p[i][j][2] + f1_p[i][j][5] + f1_p[i][j][6]));

        Pxxm = Pxym = Pyym = Pxxp = Pxyp = Pyyp = 0.;
        for (int k = 0; k < Q; k++) {

            double A = Um * cx[k] + Vm * cy[k];
            Feqm[k] = wf[k] * Rm * (1. + 3. * A + 4.5 * A * A
                 - 1.5 * (Um * Um + Vm * Vm));
            fneq = f1_m[i][j][k] - Feqm[k];

            Pxxm += fneq * cx[k] * cx[k];
            Pxym += fneq * cx[k] * cy[k];
            Pyym += fneq * cy[k] * cy[k];

            A = Up * cx[k] + Vp * cy[k];
            Feqp[k] = wf[k] * Rp * (1. + 3. * A + 4.5 * A * A
                 - 1.5 * (Up * Up + Vp * Vp));
            fneq = f1_p[i][j][k] - Feqp[k];

            Pxxp += fneq * cx[k] * cx[k];
            Pxyp += fneq * cx[k] * cy[k];
            Pyyp += fneq * cy[k] * cy[k];
        }

        for (int k = 0; k < Q; k++) {
            QPI = Pxxm * (cx[k] * cx[k] - cs2) + 2. * Pxym * cx[k] * cy[k] + Pyym * (cy[k] * cy[k] - cs2);
            f2_m[i][j][k] = Feqm[k] + 4.5 * wf[k] * QPI;

            QPI = Pxxp * (cx[k] * cx[k] - cs2) + 2. * Pxyp * cx[k] * cy[k] + Pyyp * (cy[k] * cy[k] - cs2);
            f2_p[i][j][k] = Feqp[k] + 4.5 * wf[k] * QPI;
        }
    }

    // Side walls
    // No-slip
    for (int j = 0; j < NY; j++) {

        double Rm, Um, Vm, Pxxm, Pxym, Pyym, fneq, QPI;
        double Rp, Up, Vp, Pxxp, Pxyp, Pyyp;

        // Left (x = 0)
        int i = 0;
        ux_m[i][j] = Um = 0.;
        uy_m[i][j] = Vm = 0.;
        ux_p[i][j] = Up = 0.;
        uy_p[i][j] = Vp = 0.;
        Rm = 1. / (1. - ux_m[i][j]) * (f1_m[i][j][0] + f1_m[i][j][2] + f1_m[i][j][4] +
                2. * (f1_m[i][j][3] + f1_m[i][j][6] + f1_m[i][j][7]));
        Rp = 1. / (1. - ux_p[i][j]) * (f1_p[i][j][0] + f1_p[i][j][2] + f1_p[i][j][4] +
                2. * (f1_p[i][j][3] + f1_p[i][j][6] + f1_p[i][j][7]));

        Pxxm = Pxym = Pyym = Pxxp = Pxyp = Pyyp = 0.;
        Vec1D Feqm(Q, 0.), Feqp(Q, 0.);
        for (int k = 0; k < Q; k++) {

            double A = Um * cx[k] + Vm * cy[k];
            Feqm[k] = wf[k] * Rm * (1. + 3. * A + 4.5 * A * A
                 - 1.5 * (Um * Um + Vm * Vm));
            fneq = f1_m[i][j][k] - Feqm[k];

            Pxxm += fneq * cx[k] * cx[k];
            Pxym += fneq * cx[k] * cy[k];
            Pyym += fneq * cy[k] * cy[k];

            A = Up * cx[k] + Vp * cy[k];
            Feqp[k] = wf[k] * Rp * (1. + 3. * A + 4.5 * A * A
                 - 1.5 * (Up * Up + Vp * Vp));
            fneq = f1_p[i][j][k] - Feqp[k];

            Pxxp += fneq * cx[k] * cx[k];
            Pxyp += fneq * cx[k] * cy[k];
            Pyyp += fneq * cy[k] * cy[k];
        }

        for (int k = 0; k < Q; k++) {
            QPI = Pxxm * (cx[k] * cx[k] - cs2) + 2. * Pxym * cx[k] * cy[k] + Pyym * (cy[k] * cy[k] - cs2);
            f2_m[i][j][k] = Feqm[k] + 4.5 * wf[k] * QPI;

            QPI = Pxxp * (cx[k] * cx[k] - cs2) + 2. * Pxyp * cx[k] * cy[k] + Pyyp * (cy[k] * cy[k] - cs2);
            f2_p[i][j][k] = Feqp[k] + 4.5 * wf[k] * QPI;
        }

        // Right (x = NX - 1)
        i = NX - 1;
        ux_m[i][j] = Um = 0.;
        uy_m[i][j] = Vm = 0.;
        ux_p[i][j] = Up = 0.;
        uy_p[i][j] = Vp = 0.;
        Rm = 1. / (1. + ux_m[i][j]) * (f1_m[i][j][0] + f1_m[i][j][2] + f1_m[i][j][4] +
                2. * (f1_m[i][j][1] + f1_m[i][j][5] + f1_m[i][j][8]));
        Rp = 1. / (1. + ux_p[i][j]) * (f1_p[i][j][0] + f1_p[i][j][2] + f1_p[i][j][4] +
                2. * (f1_p[i][j][1] + f1_p[i][j][5] + f1_p[i][j][8]));

        Pxxm = Pxym = Pyym = Pxxp = Pxyp = Pyyp = 0.;
        for (int k = 0; k < Q; k++) {

            double A = Um * cx[k] + Vm * cy[k];
            Feqm[k] = wf[k] * Rm * (1. + 3. * A + 4.5 * A * A
                 - 1.5 * (Um * Um + Vm * Vm));
            fneq = f1_m[i][j][k] - Feqm[k];

            Pxxm += fneq * cx[k] * cx[k];
            Pxym += fneq * cx[k] * cy[k];
            Pyym += fneq * cy[k] * cy[k];

            A = Up * cx[k] + Vp * cy[k];
            Feqp[k] = wf[k] * Rp * (1. + 3. * A + 4.5 * A * A
                 - 1.5 * (Up * Up + Vp * Vp));
            fneq = f1_p[i][j][k] - Feqp[k];

            Pxxp += fneq * cx[k] * cx[k];
            Pxyp += fneq * cx[k] * cy[k];
            Pyyp += fneq * cy[k] * cy[k];
        }

        for (int k = 0; k < Q; k++) {
            QPI = Pxxm * (cx[k] * cx[k] - cs2) + 2. * Pxym * cx[k] * cy[k] + Pyym * (cy[k] * cy[k] - cs2);
            f2_m[i][j][k] = Feqm[k] + 4.5 * wf[k] * QPI;

            QPI = Pxxp * (cx[k] * cx[k] - cs2) + 2. * Pxyp * cx[k] * cy[k] + Pyyp * (cy[k] * cy[k] - cs2);
            f2_p[i][j][k] = Feqp[k] + 4.5 * wf[k] * QPI;
        }
    }
}

double maxValue(const Vec2D &v) {
    double max_val = 0.0, original = 0.0;
    for (int i = 0; i < v.size(); i++)
        for (int j = 0; j < v[i].size(); j++)
            if (std::abs(v[i][j]) > max_val) {
                max_val = std::abs(v[i][j]);
                original = v[i][j];
            }
    return original;
}

double minValue(const Vec2D &v) {
    double min_val = 0.0, original = 0.0;
    for (int i = 0; i < v.size(); i++)
        for (int j = 0; j < v[i].size(); j++)
            if (std::abs(v[i][j]) < min_val) {
                min_val = std::abs(v[i][j]);
                original = v[i][j];
            }
    return original;
}

int main() {

    system("mkdir VTI_Poisson");

    Vec2D rho_m(NX, Vec1D(NY, 0.));
    Vec2D ux_m(NX, Vec1D(NY, 0.));
    Vec2D uy_m(NX, Vec1D(NY, 0.));
    Vec3D f1_m(NX, Vec2D(NY, Vec1D(Q, 0.)));
    Vec3D f2_m(NX, Vec2D(NY, Vec1D(Q, 0.)));

    Vec2D rho_p(NX, Vec1D(NY, 0.));
    Vec2D ux_p(NX, Vec1D(NY, 0.));
    Vec2D uy_p(NX, Vec1D(NY, 0.));
    Vec3D f1_p(NX, Vec2D(NY, Vec1D(Q, 0.)));
    Vec3D f2_p(NX, Vec2D(NY, Vec1D(Q, 0.)));

    double h = 1.;

    initialise(rho_m, ux_m, uy_m, f1_m, f2_m,
                rho_p, ux_p, uy_p, f1_p, f2_p,
                q_rho, Jx, Jy, Ex, Ey, Bz);

    VTIWriter::write2DVTI(getVTIFileName(0),
                            Ex, Ey, Bz, phi, Jx, Jy, q_rho,
                            rho_m, rho_p, ux_m, uy_m, ux_p, uy_p);

    // --- Build CSV filename ---
    std::string filename = "conservation_log(" + std::to_string(NX) + "x" + std::to_string(NY) + ").csv";
    std::ofstream csvFile(filename);
    csvFile << "timestep,M_total,Px_total,Py_total,"
            << "ux_p_avg,uy_p_avg,ux_m_avg,uy_m_avg\n";

    // --- Log initial totals at t=0 ---
    double M_total = 0.0, Px_total = 0.0, Py_total = 0.0;
    double ux_p_sum = 0.0, uy_p_sum = 0.0;
    double ux_m_sum = 0.0, uy_m_sum = 0.0;

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            M_total += rho_p[i][j] + rho_m[i][j];

            Px_total += rho_p[i][j] * ux_p[i][j]
                    + rho_m[i][j] * ux_m[i][j];

            Py_total += rho_p[i][j] * uy_p[i][j]
                    + rho_m[i][j] * uy_m[i][j];

            ux_p_sum += ux_p[i][j];
            uy_p_sum += uy_p[i][j];
            ux_m_sum += ux_m[i][j];
            uy_m_sum += uy_m[i][j];
        }
    }

    const double invN = 1.0 / (NX * NY);
    double ux_p_avg = ux_p_sum * invN;
    double uy_p_avg = uy_p_sum * invN;
    double ux_m_avg = ux_m_sum * invN;
    double uy_m_avg = uy_m_sum * invN;

    csvFile << 0 << ","
        << std::setprecision(15) << M_total << ","
        << std::setprecision(15) << Px_total << ","
        << std::setprecision(15) << Py_total << ","
        << std::setprecision(15) << ux_p_avg << ","
        << std::setprecision(15) << uy_p_avg << ","
        << std::setprecision(15) << ux_m_avg << ","
        << std::setprecision(15) << uy_m_avg << "\n";
    
    // ------------------------ Time loop ------------------------
    const double max_phys_time = 5.; // in seconds
    const int Nsteps = std::ceil(max_phys_time / dt_phys), n_out = 20;
    int check = 0;

    std::cout << "Total number of time-steps: " << Nsteps << std::endl;

    for (int t = 1; t <= Nsteps; t++) {

        // --- LBM step (collision + streaming) ---
        check = algoLB(rho_m, ux_m, uy_m, f1_m, f2_m, 
                        rho_p, ux_p, uy_p, f1_p, f2_p, 
                        q_rho, Jx, Jy, Ex, Ey, Bz);


        f1_m = f2_m;
        f1_p = f2_p;

        // --- Output VTI files every n_out steps ---
        if (t % n_out == 0) {
            VTIWriter::write2DVTI(getVTIFileName(t),
                                    Ex, Ey, Bz, phi, Jx, Jy, q_rho,
                                    rho_m, rho_p, ux_m, uy_m, ux_p, uy_p);
        }

        // --- Compute totals ---
        M_total = 0.0; Px_total = 0.0; Py_total = 0.0;
        ux_p_sum = 0.0, uy_p_sum = 0.0;
        ux_m_sum = 0.0, uy_m_sum = 0.0;
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                M_total += rho_p[i][j] + rho_m[i][j];

                Px_total += rho_p[i][j] * ux_p[i][j]
                        + rho_m[i][j] * ux_m[i][j];

                Py_total += rho_p[i][j] * uy_p[i][j]
                        + rho_m[i][j] * uy_m[i][j];

                ux_p_sum += ux_p[i][j];
                uy_p_sum += uy_p[i][j];
                ux_m_sum += ux_m[i][j];
                uy_m_sum += uy_m[i][j];
            }
        }

        const double invN = 1.0 / (NX * NY);
        double ux_p_avg = ux_p_sum * invN;
        double uy_p_avg = uy_p_sum * invN;
        double ux_m_avg = ux_m_sum * invN;
        double uy_m_avg = uy_m_sum * invN;

        std::cout << "Time-step: " << t << ", Physical Time: " << t * dt_phys
                << ", Total Mass = " << M_total
                << ", Total Momentum: Px = " << Px_total
                << ", Py = " << Py_total << std::endl;

        // --- Write to CSV ---
        csvFile << t << ","
                << std::setprecision(15) << M_total << ","
                << std::setprecision(15) << Px_total << ","
                << std::setprecision(15) << Py_total << ","
                << std::setprecision(15) << ux_p_avg << ","
                << std::setprecision(15) << uy_p_avg << ","
                << std::setprecision(15) << ux_m_avg << ","
                << std::setprecision(15) << uy_m_avg << "\n";
        
        // --- Instability check ---
        if (check == 1) {
            std::cerr << "Simulation became unstable at step " << t << std::endl;
            double max_um = 0., max_vm = 0.;
            double max_up = 0., max_vp = 0.;

            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    if (std::abs(ux_m[i][j]) > max_um) max_um = std::abs(ux_m[i][j]);
                    if (std::abs(uy_m[i][j]) > max_vm) max_vm = std::abs(uy_m[i][j]);
                    if (std::abs(ux_p[i][j]) > max_up) max_up = std::abs(ux_p[i][j]);
                    if (std::abs(uy_p[i][j]) > max_vp) max_vp = std::abs(uy_p[i][j]);
                }
            }

            std::cerr << "Maximum velocities:\n";
            std::cerr << "  Negative species: ux = " << max_um << ", uy = " << max_vm << "\n";
            std::cerr << "  Positive species: ux = " << max_up << ", uy = " << max_vp << "\n";

            VTIWriter::write2DVTI(getVTIFileName(t),
                                    Ex, Ey, Bz, phi, Jx, Jy, q_rho,
                                    rho_m, rho_p, ux_m, uy_m, ux_p, uy_p);

            break;
        }
    }

    // Close the CSV file after the loop
    csvFile.close();

    return 0;
}

