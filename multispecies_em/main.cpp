// 2D Poisson-LB solver for multispecies fluids with electromagnetic forces
#include "VTIWriter.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct D2Q9 {
    static constexpr int Q = 9;
    const int cx[Q], cy[Q], opp[Q];
    const double wf[Q];

    static constexpr double cs2 = 1/3.;
    static constexpr double LX = 1.0;
    static constexpr int NX = 256;
    static constexpr int NY = NX;
    static constexpr double dx = LX / NX;
    static constexpr double dt = dx;
    static constexpr double nu_phys = 1./6400;
    static constexpr double nu = nu_phys * dt / (dx * dx);
    static constexpr double tau = 0.53;
    static constexpr double omega = 1. / tau, omega1 = 1. - omega;

    static constexpr double q_charge = 0.025;
    static constexpr double collision_freq = 0.0006;

    D2Q9() : cx{0, 1, 0, -1, 0, 1, -1, -1, 1},
             cy{0, 0, 1, 0, -1, 1, 1, -1, -1},
             opp{0, 3, 4, 1, 2, 7, 8, 5, 6},
             wf{4/9., 1/9., 1/9., 1/9., 1/9., 1/36., 1/36., 1/36., 1/36.} {}
};

std::string getVTIFileName(int timestep) {
    std::ostringstream oss;
    oss << "VTI_Poisson/poisson_t"
        << std::setw(6) << std::setfill('0') << timestep
        << ".vti";
    return oss.str();
}

void initialise(Vec2D &rho_p, Vec2D &rho_m, Vec2D &ux_p, Vec2D &uy_p, Vec2D &ux_m, Vec2D &uy_m,
                Vec3D &f1_p, Vec3D &f2_p, Vec3D &f1_m, Vec3D &f2_m, const D2Q9 &lattice,
                Vec2D &phi, Vec2D &Ex, Vec2D &Ey, Vec2D &q_rho) {

    for (int i = 0; i < lattice.NX; i++) {
        double x = double(i) / double(lattice.NX);
        for (int j = 0; j < lattice.NY; j++) {
            double y = double(j) / double(lattice.NY);

            static constexpr double U0 = 0.06;
            rho_p[i][j] = 1. + 0.01 * std::sin(2 * M_PI * x);
            rho_m[i][j] = 1. - 0.01 * std::sin(2 * M_PI * x);
            ux_p[i][j] = 0.0;
            uy_p[i][j] = 0.0;
            ux_m[i][j] = 0.0;
            uy_m[i][j] = 0.0;

            phi[i][j] = 0.0;
            Ex[i][j]  = 0.0;
            Ey[i][j]  = 0.0;
            q_rho[i][j] = lattice.q_charge * (rho_p[i][j] - rho_m[i][j]);

            for (int k = 0; k < lattice.Q; k++) {
                double A_p = ux_p[i][j] * lattice.cx[k] + uy_p[i][j] * lattice.cy[k];
                f1_p[i][j][k] = lattice.wf[k] * rho_p[i][j] * (1. + 3. * A_p + 4.5 * A_p * A_p - 1.5 * (ux_p[i][j] * ux_p[i][j] + uy_p[i][j] * uy_p[i][j]));
                f2_p[i][j][k] = f1_p[i][j][k];

                double A_m = ux_m[i][j] * lattice.cx[k] + uy_m[i][j] * lattice.cy[k];
                f1_m[i][j][k] = lattice.wf[k] * rho_m[i][j] * (1. + 3. * A_m + 4.5 * A_m * A_m - 1.5 * (ux_m[i][j] * ux_m[i][j] + uy_m[i][j] * uy_m[i][j]));
                f2_m[i][j][k] = f1_m[i][j][k];
            }
        }
    }
}

int algoLB(Vec2D &rho_p, Vec2D &rho_m, Vec2D &ux_p, Vec2D &uy_p, Vec2D &ux_m, Vec2D &uy_m,
           Vec3D &f1_p, Vec3D &f2_p, Vec3D &f1_m, Vec3D &f2_m, const D2Q9 &lattice,
           Vec2D &Ex, Vec2D &Ey, Vec2D &q_rho, Vec2D &Bz) {

    for (int i = 0; i < lattice.NX; i++) {
        for (int j = 0; j < lattice.NY; j++) {

            double Rp, Rm, Up, Vp, Um, Vm;
            Rp = 0.0, Rm = 0.0, Up = 0.0, Vp = 0.0, Um = 0.0, Vm = 0.0;
            for (int k = 0; k < lattice.Q; k++) {
                Rp += f1_p[i][j][k];
                Rm += f1_m[i][j][k];
                Up += f1_p[i][j][k] * lattice.cx[k];
                Vp += f1_p[i][j][k] * lattice.cy[k];
                Um += f1_m[i][j][k] * lattice.cx[k];
                Vm += f1_m[i][j][k] * lattice.cy[k];
            }

            rho_p[i][j] = Rp;
            rho_m[i][j] = Rm;
            double UP = Up / Rp;
            double VP = Vp / Rp;
            double UM = Um / Rm;
            double VM = Vm / Rm;

            double reduced_mass  = 2. * Rp * Rm / (Rp + Rm);

            // Electric + drag forces
            double Fx_p = Rp * lattice.q_charge * Ex[i][j] + reduced_mass * lattice.collision_freq * (UM - UP);
            double Fy_p = Rp * lattice.q_charge * Ey[i][j] + reduced_mass * lattice.collision_freq * (VM - VP);
            double Fx_m = Rm * -lattice.q_charge * Ex[i][j] + reduced_mass * lattice.collision_freq * (UP - UM);
            double Fy_m = Rm * -lattice.q_charge * Ey[i][j] + reduced_mass * lattice.collision_freq * (VP - VM);

            // Magnetic Lorentz force
            double Bz_loc = Bz[i][j];
            Fx_p += Rp * lattice.q_charge * VP * Bz_loc;
            Fy_p += Rp * -lattice.q_charge * UP * Bz_loc;
            Fx_m += Rm * -lattice.q_charge * VM * Bz_loc;
            Fy_m += Rm * lattice.q_charge * UM * Bz_loc;

            Up = (Up + 0.5 * Fx_p) / Rp;
            Vp = (Vp + 0.5 * Fy_p) / Rp;
            Um = (Um + 0.5 * Fx_m) / Rm;
            Vm = (Vm + 0.5 * Fy_m) / Rm;
            ux_p[i][j] = Up; uy_p[i][j] = Vp; ux_m[i][j] = Um; uy_m[i][j] = Vm;
            q_rho[i][j] = lattice.q_charge * (Rp - Rm);

            double Um2 = Um * Um, Up2 = Up * Up;
            double Vm2 = Vm * Vm, Vp2 = Vp * Vp;
            double UVm = Um * Vm, UVp = Up * Vp;

            double k1_m, k2_m, k3_m, k4_m, k5_m, k6_m, k7_m, k8_m;
            double k1_p, k2_p, k3_p, k4_p, k5_p, k6_p, k7_p, k8_p;
            double r1_m, r2_m, r3_m, r4_m, r5_m, r6_m, r7_m, r8_m;
            double r1_p, r2_p, r3_p, r4_p, r5_p, r6_p, r7_p, r8_p;

            r4_m = f1_m[i][j][1] - f1_m[i][j][2] + f1_m[i][j][3] - f1_m[i][j][4];
            r5_m = f1_m[i][j][5] - f1_m[i][j][6] + f1_m[i][j][7] - f1_m[i][j][8];
            k4_m = r4_m - Rm * (Um2 - Vm2);
            k5_m = r5_m - Rm * UVm;
            k1_m = 0.5 * Fx_m;
            k2_m = 0.5 * Fy_m;
            k3_m = 2. * Rm * lattice.cs2;
            k4_m = lattice.omega1 * k4_m;
            k5_m = lattice.omega1 * k5_m;
            k6_m = 0.5 * Fy_m * lattice.cs2;
            k7_m = 0.5 * Fx_m * lattice.cs2;
            k8_m = Rm * lattice.cs2 * lattice.cs2;
            r1_m = k1_m + Rm * Um;
            r2_m = k2_m + Rm * Vm;
            r3_m = k3_m + 2. * Um * k1_m + 2. * Vm * k2_m + Rm * (Um2 + Vm2);
            r4_m = k4_m + 2. * Um * k1_m - 2. * Vm * k2_m + Rm * (Um2 - Vm2);
            r5_m = k5_m + Um * k2_m + Vm * k1_m + Rm * UVm;
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
            k3_p = 2. * Rp * lattice.cs2;
            k4_p = lattice.omega1 * k4_p;
            k5_p = lattice.omega1 * k5_p;
            k6_p = 0.5 * Fy_p * lattice.cs2;
            k7_p = 0.5 * Fx_p * lattice.cs2;
            k8_p = Rp * lattice.cs2 * lattice.cs2;
            r1_p = k1_p + Rp * Up;
            r2_p = k2_p + Rp * Vp;
            r3_p = k3_p + 2. * Up * k1_p + 2. * Vp * k2_p + Rp * (Up2 + Vp2);
            r4_p = k4_p + 2. * Up * k1_p - 2. * Vp * k2_p + Rp * (Up2 - Vp2);
            r5_p = k5_p + Up * k2_p + Vp * k1_p + Rp * UVp;
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

            for (int k = 0; k < lattice.Q; k++) {
                int in = (i + lattice.cx[k] + lattice.NX) % lattice.NX;
                int jn = (j + lattice.cy[k] + lattice.NY) % lattice.NY;
                f2_m[in][jn][k] = f1_m[i][j][k];
                f2_p[in][jn][k] = f1_p[i][j][k];
            }

            if (fabs(Um) > 0.1 || fabs(Up) > 0.1 || fabs(Vm) > 0.1 || fabs(Vp) > 0.1) {
                return 1;
            }
        }
    }

    return 0;
}

void dirichlet_boundary(Vec2D &rho_p, Vec2D &rho_m, Vec2D &ux_p, Vec2D &uy_p, Vec2D &ux_m, Vec2D &uy_m,
              Vec3D &f1_p, Vec3D &f2_p, Vec3D &f1_m, Vec3D &f2_m, const D2Q9 &lattice,
              Vec2D &phi, Vec2D &Ex, Vec2D &Ey, Vec2D &q_rho) {

    int i = 0;
    for (int j = 0; j < lattice.NY; j++) {
        f2_p[i][j][1] = f1_p[i][j][3]; f2_p[i][j][5] = f1_p[i][j][7]; f2_p[i][j][8] = f1_p[i][j][6];
        f2_m[i][j][1] = f1_m[i][j][3]; f2_m[i][j][5] = f1_m[i][j][7]; f2_m[i][j][8] = f1_m[i][j][6];
    }
    i = lattice.NX - 1;
    for (int j = 0; j < lattice.NY; j++) {
        f2_p[i][j][3] = f1_p[i][j][1]; f2_p[i][j][7] = f1_p[i][j][5]; f2_p[i][j][6] = f1_p[i][j][8];
        f2_m[i][j][3] = f1_m[i][j][1]; f2_m[i][j][6] = f1_m[i][j][8]; f2_m[i][j][7] = f1_m[i][j][5];
    }
    int j = 0;
    for (int i = 0; i < lattice.NX; i++) {
        f2_p[i][j][2] = f1_p[i][j][4]; f2_p[i][j][5] = f1_p[i][j][7]; f2_p[i][j][6] = f1_p[i][j][8];
        f2_m[i][j][2] = f1_m[i][j][4]; f2_m[i][j][5] = f1_m[i][j][7]; f2_m[i][j][6] = f1_m[i][j][8];
    }
    j = lattice.NY - 1;
    for (int i = 0; i < lattice.NX; i++) {
        f2_p[i][j][4] = f1_p[i][j][2]; f2_p[i][j][7] = f1_p[i][j][5]; f2_p[i][j][8] = f1_p[i][j][6];
        f2_m[i][j][4] = f1_m[i][j][2]; f2_m[i][j][7] = f1_m[i][j][5]; f2_m[i][j][8] = f1_m[i][j][6];
    }
    for (int i = 0; i < lattice.NX; i++)
        for (int j = 0; j < lattice.NY; j++)
            if (i == 0 || i == lattice.NX - 1 || j == 0 || j == lattice.NY - 1)
                ux_p[i][j] = uy_p[i][j] = ux_m[i][j] = uy_m[i][j] = 0.0;
}

// Periodic Poisson CG solver
void Poisson_CG(Vec2D &sol, const Vec2D &rhs_in, int NX, int NY, double tol = 1.e-6, int max_iter = 10000) {

    Vec2D rhs(NX, Vec1D(NY, 0.));
    double mean = 0.;
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++) {
            rhs[i][j] = rhs_in[i][j];
            mean += rhs[i][j];
        }
    mean /= (NX * NY);
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
            rhs[i][j] -= mean;

    auto apply_L = [&](const Vec2D &p, Vec2D &Ap) {
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++) {
                int ip = (i+1)%NX, im = (i-1+NX)%NX;
                int jp = (j+1)%NY, jm = (j-1+NY)%NY;
                Ap[i][j] = p[ip][j] + p[im][j] + p[i][jp] + p[i][jm] - 4.*p[i][j];
            }
    };

    Vec2D Ap(NX, Vec1D(NY, 0.));
    apply_L(sol, Ap);
    Vec2D r(NX, Vec1D(NY, 0.));
    double rr = 0.;
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++) {
            r[i][j] = -rhs[i][j] - Ap[i][j];
            rr += r[i][j] * r[i][j];
        }

    Vec2D p = r;
    for (int iter = 0; iter < max_iter; iter++) {
        apply_L(p, Ap);
        double pAp = 0.;
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++)
                pAp += p[i][j] * Ap[i][j];
        if (std::abs(pAp) < 1e-30) break;
        double alpha = rr / pAp;
        double rr_new = 0.;
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++) {
                sol[i][j] += alpha * p[i][j];
                r[i][j]   -= alpha * Ap[i][j];
                rr_new    += r[i][j] * r[i][j];
            }
        if (std::sqrt(rr_new) < tol) break;
        double beta = rr_new / rr;
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++)
                p[i][j] = r[i][j] + beta * p[i][j];
        rr = rr_new;
    }

    double pmean = 0.;
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
            pmean += sol[i][j];
    pmean /= (NX * NY);
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
            sol[i][j] -= pmean;
}

// Electric Poisson
void Poisson_E(Vec2D &phi, Vec2D &Ex, Vec2D &Ey, const Vec2D &q_rho, const D2Q9 &lattice,
               const Vec2D *q_rho_fixed = nullptr) {

    int NX = lattice.NX, NY = lattice.NY;
    Vec2D rhs(NX, Vec1D(NY, 0.));
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
            rhs[i][j] = q_rho[i][j] + (q_rho_fixed ? (*q_rho_fixed)[i][j] : 0.);

    Poisson_CG(phi, rhs, NX, NY);

    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++) {
            Ex[i][j] = -(phi[(i + 1) % NX][j] - phi[(i - 1 + NX) % NX][j]) / 2.;
            Ey[i][j] = -(phi[i][(j + 1) % NY] - phi[i][(j - 1 + NY) % NY]) / 2.;
        }
}

// Magnetic Poisson
void Poisson_B(Vec2D &Bz, const Vec2D &Jx, const Vec2D &Jy, const D2Q9 &lattice) {

    int NX = lattice.NX, NY = lattice.NY;
    Vec2D curl_J(NX, Vec1D(NY, 0.));
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++) {
            double dJy_dx = (Jy[(i + 1) % NX][j] - Jy[(i - 1 + NX) % NX][j]) / 2.;
            double dJx_dy = (Jx[i][(j + 1) % NY] - Jx[i][(j - 1 + NY) % NY]) / 2.;
            curl_J[i][j] = dJy_dx - dJx_dy;
        }

    Poisson_CG(Bz, curl_J, NX, NY);
}

int main() {

    const D2Q9 lattice;
    system("mkdir VTI_Poisson");

    int NX = lattice.NX, NY = lattice.NY;

    // Fluid fields
    Vec2D rho_p(NX, Vec1D(NY, 0.));
    Vec2D rho_m(NX, Vec1D(NY, 0.));
    Vec2D ux_p(NX, Vec1D(NY, 0.));
    Vec2D uy_p(NX, Vec1D(NY, 0.));
    Vec2D ux_m(NX, Vec1D(NY, 0.));
    Vec2D uy_m(NX, Vec1D(NY, 0.));
    Vec3D f1_p(NX, Vec2D(NY, Vec1D(lattice.Q, 0.)));
    Vec3D f2_p(NX, Vec2D(NY, Vec1D(lattice.Q, 0.)));
    Vec3D f1_m(NX, Vec2D(NY, Vec1D(lattice.Q, 0.)));
    Vec3D f2_m(NX, Vec2D(NY, Vec1D(lattice.Q, 0.)));

    // Electromagnetic fields
    Vec2D phi(NX, Vec1D(NY, 0.));
    Vec2D Ex(NX, Vec1D(NY, 0.));
    Vec2D Ey(NX, Vec1D(NY, 0.));
    Vec2D q_rho(NX, Vec1D(NY, 0.));
    Vec2D Bz(NX, Vec1D(NY, 0.));
    Vec2D Jx(NX, Vec1D(NY, 0.));
    Vec2D Jy(NX, Vec1D(NY, 0.));

    // No fixed background charge
    Vec2D q_rho_fixed(NX, Vec1D(NY, 0.));

    std::cout << "nu_phys: " << lattice.nu_phys
              << "  tau: "   << lattice.tau
              << "  q: "     << lattice.q_charge
              << "  N: "     << NX << std::endl;

    initialise(rho_p, rho_m, ux_p, uy_p, ux_m, uy_m,
               f1_p, f2_p, f1_m, f2_m, lattice, phi, Ex, Ey, q_rho);

    Poisson_E(phi, Ex, Ey, q_rho, lattice, &q_rho_fixed);
    Poisson_B(Bz, Jx, Jy, lattice);

    std::ofstream csv("probe_signal.csv");
    csv << "t,delta_rho,rho_p,rho_m\n";

    double max_delta_rho = 0.0;
    double max_rho_p = 0.0;
    double max_rho_m = 0.0;

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {

            double drho = std::abs(rho_p[i][j] - rho_m[i][j]);

            if (drho > max_delta_rho) {
                max_delta_rho = drho;
                max_rho_p = rho_p[i][j];
                max_rho_m = rho_m[i][j];
            }
        }
    }

    csv << 0 << ","
        << max_delta_rho << ","
        << max_rho_p << ","
        << max_rho_m << "\n";

    int n_steps  = 300;
    int out_freq = 1;

    for (int t = 1; t <= n_steps; t++) {

        int ret = algoLB(rho_p, rho_m, ux_p, uy_p, ux_m, uy_m,
                         f1_p, f2_p, f1_m, f2_m, lattice, Ex, Ey, q_rho, Bz);

        std::swap(f1_p, f2_p);
        std::swap(f1_m, f2_m);

        // Compute current density
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++) {
                Jx[i][j] = lattice.q_charge * (rho_p[i][j]*ux_p[i][j] - rho_m[i][j]*ux_m[i][j]);
                Jy[i][j] = lattice.q_charge * (rho_p[i][j]*uy_p[i][j] - rho_m[i][j]*uy_m[i][j]);
            }

        // Solve Poisson equations
        Poisson_E(phi, Ex, Ey, q_rho, lattice, &q_rho_fixed);
        Poisson_B(Bz, Jx, Jy, lattice);

        if (ret) {
            std::cout << "INSTABILITY at t=" << t << std::endl;
            break;
        }

        if (t % out_freq == 0) {

            max_delta_rho = 0.0;
            max_rho_p = 0.0;
            max_rho_m = 0.0;

            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {

                    double drho = std::abs(rho_p[i][j] - rho_m[i][j]);

                    if (drho > max_delta_rho) {
                        max_delta_rho = drho;
                        max_rho_p = rho_p[i][j];
                        max_rho_m = rho_m[i][j];
                    }
                }
            }

            csv << t << ","
                << max_delta_rho << ","
                << max_rho_p << ","
                << max_rho_m << "\n";

            std::cout << "t: " << t
                    << ", max_delta_rho: " << max_delta_rho
                    << std::endl;
        }
    }

    csv.close();

    return 0;
}
