#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include "BQBN-Reborn.h"

using namespace std;

static void InitializeLayerFlow()
{
    int width = 35;

    // Initialize rho_B to rho0_B for all cells
    std::fill(rho_R.begin(), rho_R.end(), 0.0);
    std::fill(rho_B.begin(), rho_B.end(), rho0_B);

    // Initialize rho_R, rho_B, rho, and distributions
    for (int x = 0; x < lx; x++)
    {
        for (int y = 0; y < ly; y++)
        {
            const size_t id = IDX0(x, y);
            if (obstacle[id] != 1 && y >= width && y < ly - width)
            {
                rho_R[id] = rho0_R;
                rho_B[id] = 0.0;
            }

            rho[id] = rho_R[id] + rho_B[id];

            const size_t id9_base = IDX9(x, y, 0);
            for (int i = 0; i < 9; i++)
            {
                const size_t id9 = id9_base + i;
                f_R[id9] = rho_R[id] * Ci_R[i];
                f_B[id9] = rho_B[id] * Ci_B[i];
                ff[id9] = f_R[id9] + f_B[id9];
            }
        }
    }
}

static void InitializeDrainage()
{
    // Set rho_B to rho0_B for all cells
    std::fill(rho_R.begin(), rho_R.end(), 0.0);
    std::fill(rho_B.begin(), rho_B.end(), rho0_B);

    // Initialize rho_R, rho_B, rho, and distributions
    for (int x = 0; x < lx; x++)
    {
        const bool left_region = (x < 5);
        const size_t base0_x = static_cast<size_t>(x) * static_cast<size_t>(ly);
        const size_t base9_x = base0_x * 9u;

        for (int y = 0; y < ly; y++)
        {
            const size_t ys = static_cast<size_t>(y);
            const size_t id = base0_x + ys;
            const size_t id9_base = base9_x + ys * 9u;

            if (left_region)
            {
                rho_R[id] = rho0_R;
                rho_B[id] = 0.0;
            }

            rho[id] = rho_R[id] + rho_B[id];

            for (int i = 0; i < 9; i++)
            {
                const size_t id9 = id9_base + static_cast<size_t>(i);
                f_R[id9] = rho_R[id] * Ci_R[i];
                f_B[id9] = rho_B[id] * Ci_B[i];
                ff[id9] = f_R[id9] + f_B[id9];
            }

            if (obstacle[id] == 1)
            {
                rho_R[id] = rho_wr;
                rho_B[id] = rho_wb;
            }
        }
    }
}

static void Stream()
{
    static std::vector<double> f_temp_R(lx * ly * 9, 0.0);
    static std::vector<double> f_temp_B(lx * ly * 9, 0.0);

    for (int x = 0; x < lx; x++)
    {
        const size_t base_x = static_cast<size_t>(x) * static_cast<size_t>(ly) * 9u;
        for (int y = 0; y < ly; y++)
        {
            const size_t base_y = base_x + static_cast<size_t>(y) * 9;

            for (int i = 0; i < 9; i++)
            {
                int x_temp = (x + ex[i] + lx) % lx;
                int y_temp = (y + ey[i] + ly) % ly;
                const size_t base_temp = static_cast<size_t>(x_temp) * ly * 9 + static_cast<size_t>(y_temp) * 9 + i;

                f_temp_R[base_temp] = f_R[base_y + i];
                f_temp_B[base_temp] = f_B[base_y + i];
            }
        }
    }

    f_R.swap(f_temp_R);
    f_B.swap(f_temp_B);
}

static void PressureBoundary()
{
    for (int y = 1; y < ly - 1; y++)
    {
        const size_t idw = static_cast<size_t>(y) * 9;
        const size_t ide = static_cast<size_t>(lx-1) * ly * 9 + static_cast<size_t>(y) * 9;
        
        for (int i = 0; i < 9; i++)
        {
            f_R[idw + i] = 0.0;
            f_B[idw + i] = 0.0;
            f_R[ide + i] = 0.0;
            f_B[ide + i] = 0.0;
        }
        
        f_R[idw + 1] = 2.0 / 3 * rho_in;
        f_R[idw + 5] = 1.0 / 6 * rho_in;
        f_R[idw + 8] = 1.0 / 6 * rho_in; 

        f_B[ide + 3] = -2.0 / 3 * rho_out;
        f_B[ide + 7] = -1.0 / 6 * rho_out;
        f_B[ide + 6] = -1.0 / 6 * rho_out;
    }
}

static void GetUnRho()
{
    for (int x = 0; x < lx; ++x)
    {
        const size_t base0_x = static_cast<size_t>(x) * ly;
        const size_t base9_x = base0_x * 9u;
        const size_t base2_x = base0_x * 2u;

        for (int y = 0; y < ly; ++y)
        {
            const size_t id0 = base0_x + static_cast<size_t>(y);
            const size_t id9 = base9_x + static_cast<size_t>(y) * 9u;
            const size_t id2 = base2_x + static_cast<size_t>(y) * 2u;

            if (obstacle[id0] == 0)
            {
                // Load f's once
                const double r0 = f_R[id9 + 0], r1 = f_R[id9 + 1], r2 = f_R[id9 + 2];
                const double r3 = f_R[id9 + 3], r4 = f_R[id9 + 4], r5 = f_R[id9 + 5];
                const double r6 = f_R[id9 + 6], r7 = f_R[id9 + 7], r8 = f_R[id9 + 8];

                const double b0 = f_B[id9 + 0], b1 = f_B[id9 + 1], b2 = f_B[id9 + 2];
                const double b3 = f_B[id9 + 3], b4 = f_B[id9 + 4], b5 = f_B[id9 + 5];
                const double b6 = f_B[id9 + 6], b7 = f_B[id9 + 7], b8 = f_B[id9 + 8];

                const double rhoR = r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8;
                const double rhoB = b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8;
                const double rhoT = rhoR + rhoB;

                rho_R[id0] = rhoR;
                rho_B[id0] = rhoB;
                rho[id0]   = rhoT;

                const double inv_rho = 1.0 / rhoT;

                // u_x
                u[id2 + 0] =
                    ((-r6 + r5 - r3 + r1 - r7 + r8)
                        + (-b6 + b5 - b3 + b1 - b7 + b8)) * inv_rho;

                // u_y
                u[id2 + 1] =
                    ((r6 + r2 + r5 - r7 - r4 - r8)
                        + (b6 + b2 + b5 - b7 - b4 - b8)) * inv_rho;
            }
            else
            {
                rho_R[id0] = rho_wr;
                rho_B[id0] = rho_wb;
            }
        }
    }
}

static void CollisionMRT()
{
    const double inv36 = 1.0 / 36.0;

    for (int x = 0; x < lx; ++x)
    {
        const size_t base0_x = static_cast<size_t>(x) * static_cast<size_t>(ly);
        const size_t base9_x = static_cast<size_t>(x) * static_cast<size_t>(ly) * 9u;
        const size_t base2_x = static_cast<size_t>(x) * static_cast<size_t>(ly) * 2u;

        for (int y = 0; y < ly; ++y)
        {
            const size_t ys = static_cast<size_t>(y);
            const size_t id0 = base0_x + ys;

            if (obstacle[id0] != 0) continue;

            const size_t id9 = base9_x + ys * 9u;
            const size_t id2 = base2_x + ys * 2u;

            const double Fx = df;
            const double Fy = 0.0;  // as in your code

            const double ux = u[id2 + 0];
            const double uy = u[id2 + 1];
            const double u_sqr = ux * ux + uy * uy;

            // --- load total f = f_R + f_B (only once) ---
            const double f0 = f_R[id9 + 0] + f_B[id9 + 0];
            const double f1 = f_R[id9 + 1] + f_B[id9 + 1];
            const double f2 = f_R[id9 + 2] + f_B[id9 + 2];
            const double f3 = f_R[id9 + 3] + f_B[id9 + 3];
            const double f4 = f_R[id9 + 4] + f_B[id9 + 4];
            const double f5 = f_R[id9 + 5] + f_B[id9 + 5];
            const double f6 = f_R[id9 + 6] + f_B[id9 + 6];
            const double f7 = f_R[id9 + 7] + f_B[id9 + 7];
            const double f8 = f_R[id9 + 8] + f_B[id9 + 8];

            const double rhoR = rho_R[id0];
            const double rhoB = rho_B[id0];
            const double rhoT = rhoR + rhoB;

            // --- build feq_total[i] directly (no need for feq_R + feq_B arrays) ---
            // feq_total = rhoR*Ci_R + rhoB*Ci_B + rhoT*W*(3 u·e + 4.5(u·e)^2 - 1.5 u^2)
            // Precompute u·e for D2Q9 explicitly:
            const double ue0 = 0.0;
            const double ue1 = ux;
            const double ue2 = uy;
            const double ue3 = -ux;
            const double ue4 = -uy;
            const double ue5 = ux + uy;
            const double ue6 = -ux + uy;
            const double ue7 = -ux - uy;
            const double ue8 = ux - uy;

            auto feq_total = [&](int i, double ue) -> double {
                const double common = (3.0 * ue + 4.5 * ue * ue - 1.5 * u_sqr);
                return rhoR * Ci_R[i] + rhoB * Ci_B[i] + rhoT * W[i] * common;
                };

            const double fe0 = feq_total(0, ue0);
            const double fe1 = feq_total(1, ue1);
            const double fe2 = feq_total(2, ue2);
            const double fe3 = feq_total(3, ue3);
            const double fe4 = feq_total(4, ue4);
            const double fe5 = feq_total(5, ue5);
            const double fe6 = feq_total(6, ue6);
            const double fe7 = feq_total(7, ue7);
            const double fe8 = feq_total(8, ue8);

            // --- moment transform m = M * f (UNROLLED using your M) ---
            const double m0 = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
            const double m1 = -4 * f0 - f1 - f2 - f3 - f4 + 2 * f5 + 2 * f6 + 2 * f7 + 2 * f8;
            const double m2 = 4 * f0 - 2 * f1 - 2 * f2 - 2 * f3 - 2 * f4 + f5 + f6 + f7 + f8;
            const double m3 = f1 - f3 + f5 - f6 - f7 + f8;
            const double m4 = -2 * f1 + 2 * f3 + f5 - f6 - f7 + f8;
            const double m5 = f2 - f4 + f5 + f6 - f7 - f8;
            const double m6 = -2 * f2 + 2 * f4 + f5 + f6 - f7 - f8;
            const double m7 = f1 - f2 + f3 - f4;
            const double m8 = f5 - f6 + f7 - f8;

            // --- moment transform meq = M * feq_total (UNROLLED) ---
            const double me0 = fe0 + fe1 + fe2 + fe3 + fe4 + fe5 + fe6 + fe7 + fe8;
            const double me1 = -4 * fe0 - fe1 - fe2 - fe3 - fe4 + 2 * fe5 + 2 * fe6 + 2 * fe7 + 2 * fe8;
            const double me2 = 4 * fe0 - 2 * fe1 - 2 * fe2 - 2 * fe3 - 2 * fe4 + fe5 + fe6 + fe7 + fe8;
            const double me3 = fe1 - fe3 + fe5 - fe6 - fe7 + fe8;
            const double me4 = -2 * fe1 + 2 * fe3 + fe5 - fe6 - fe7 + fe8;
            const double me5 = fe2 - fe4 + fe5 + fe6 - fe7 - fe8;
            const double me6 = -2 * fe2 + 2 * fe4 + fe5 + fe6 - fe7 - fe8;
            const double me7 = fe1 - fe2 + fe3 - fe4;
            const double me8 = fe5 - fe6 + fe7 - fe8;

            // --- force term in moment space (same as your code) ---
            double fm0 = 0.0;
            double fm1 = 6.0 * (Fx * ux + Fy * uy);
            double fm2 = -fm1;
            double fm3 = Fx;
            double fm4 = -Fx;
            double fm5 = Fy;
            double fm6 = -Fy;
            double fm7 = 2.0 * (Fx * ux - Fy * uy);
            double fm8 = Fx * ux + Fy * uy;

            // --- tau selection ---
            const double rhoN = (rhoR - rhoB) / rhoT;
            double tau;
            if (rhoN > delta) tau = tau_R;
            else if (rhoN > 0.0)  tau = st1 + s2 * rhoN + s3 * rhoN * rhoN;
            else if (rhoN > -delta) tau = st1 + t2 * rhoN + t3 * rhoN * rhoN;
            else tau = tau_B;

            const double inv_tau = 1.0 / tau;

            // S and Sf exactly as your code
            const double S0 = 0.0, S1 = inv_tau, S2 = 1.4, S3 = 0.0, S4 = 1.2,
                S5 = 0.0, S6 = 1.2, S7 = inv_tau, S8 = inv_tau;

            const double Sf0 = 1.0 - 0.5 * S0;
            const double Sf1 = 1.0 - 0.5 * S1;
            const double Sf2 = 1.0 - 0.5 * S2;
            const double Sf3 = 1.0 - 0.5 * S3;
            const double Sf4 = 1.0 - 0.5 * S4;
            const double Sf5 = 1.0 - 0.5 * S5;
            const double Sf6 = 1.0 - 0.5 * S6;
            const double Sf7 = 1.0 - 0.5 * S7;
            const double Sf8 = 1.0 - 0.5 * S8;

            // --- collide in moment space: m' = m - S*(m-meq) + Sf*fm ---
            const double mp0 = m0 - S0 * (m0 - me0) + Sf0 * fm0;
            const double mp1 = m1 - S1 * (m1 - me1) + Sf1 * fm1;
            const double mp2 = m2 - S2 * (m2 - me2) + Sf2 * fm2;
            const double mp3 = m3 - S3 * (m3 - me3) + Sf3 * fm3;
            const double mp4 = m4 - S4 * (m4 - me4) + Sf4 * fm4;
            const double mp5 = m5 - S5 * (m5 - me5) + Sf5 * fm5;
            const double mp6 = m6 - S6 * (m6 - me6) + Sf6 * fm6;
            const double mp7 = m7 - S7 * (m7 - me7) + Sf7 * fm7;
            const double mp8 = m8 - S8 * (m8 - me8) + Sf8 * fm8;

            // --- back transform: f' = M_inv * m' (UNROLLED using your M_inv_pre / 36) ---
            ff[id9 + 0] = inv36 * (4 * mp0 - 4 * mp1 + 4 * mp2);
            ff[id9 + 1] = inv36 * (4 * mp0 - 1 * mp1 - 2 * mp2 + 6 * mp3 - 6 * mp4 + 9 * mp7);
            ff[id9 + 2] = inv36 * (4 * mp0 - 1 * mp1 - 2 * mp2 + 6 * mp5 - 6 * mp6 - 9 * mp7);
            ff[id9 + 3] = inv36 * (4 * mp0 - 1 * mp1 - 2 * mp2 - 6 * mp3 + 6 * mp4 + 9 * mp7);
            ff[id9 + 4] = inv36 * (4 * mp0 - 1 * mp1 - 2 * mp2 - 6 * mp5 + 6 * mp6 - 9 * mp7);
            ff[id9 + 5] = inv36 * (4 * mp0 + 2 * mp1 + 1 * mp2 + 6 * mp3 + 3 * mp4 + 6 * mp5 + 3 * mp6 + 9 * mp8);
            ff[id9 + 6] = inv36 * (4 * mp0 + 2 * mp1 + 1 * mp2 - 6 * mp3 - 3 * mp4 + 6 * mp5 + 3 * mp6 - 9 * mp8);
            ff[id9 + 7] = inv36 * (4 * mp0 + 2 * mp1 + 1 * mp2 - 6 * mp3 - 3 * mp4 - 6 * mp5 - 3 * mp6 + 9 * mp8);
            ff[id9 + 8] = inv36 * (4 * mp0 + 2 * mp1 + 1 * mp2 + 6 * mp3 + 3 * mp4 - 6 * mp5 - 3 * mp6 - 9 * mp8);
        }
    }
}

static void Redistribution()
{
    for (int x = 1; x < lx - 1; x++)
    {
        const size_t base0_x = static_cast<size_t>(x) * static_cast<size_t>(ly);
        const size_t base9_x = base0_x * 9u;

        for (int y = 0; y < ly; y++)
        {
            const size_t ys = static_cast<size_t>(y);
            const size_t id0 = base0_x + ys;

            if (obstacle[id0] != 0) continue;

            const size_t id9 = base9_x + ys * 9u;

            const double rho_here = rho[id0];
            const double rhoR_here = rho_R[id0];
            const double rhoB_here = rho_B[id0];

            const double inv_rho = 1.0 / rho_here;
            const double ratioR = rhoR_here * inv_rho;
            const double ratioB = rhoB_here * inv_rho;

            // Pointers to the 9 populations of this cell
            double* const ff_cell = &ff[id9];
            double* const fR_cell = &f_R[id9];
            double* const fB_cell = &f_B[id9];

            // --- Compute Fx, Fy ---
            double Fx = 0.0, Fy = 0.0;
            for (int i = 0; i < 9; i++)
            {
                int xi = x + ex[i];

                int yi = y + ey[i];
                if (yi < 0) yi = ly - 1;
                else if (yi > ly-1) yi = 0;

                const size_t idn = static_cast<size_t>(xi) * static_cast<size_t>(ly) + static_cast<size_t>(yi);

                const double drho = rho_R[idn] - rho_B[idn];
                Fx += e[i][0] * drho;
                Fy += e[i][1] * drho;
            }

            const double Fm2 = Fx * Fx + Fy * Fy;

            if (Fm2 < 1e-16)
            {
                for (int i = 0; i < 9; ++i)
                {
                    const double base = ff_cell[i];
                    fR_cell[i] = ratioR * base;
                    fB_cell[i] = ratioB * base;
                }
                continue;
            }

            const double Fm = std::sqrt(Fm2);
            const double invFm = 1.0 / Fm;

            // Precompute dot-products once
            double ebyf[9]{};
            for (int i = 0; i < 9; i++)
                ebyf[i] = e[i][0] * Fx + e[i][1] * Fy;

            const double AFm = A * Fm;

            // ---- update ff ----
            // uses ebyf_over_Fm to avoid dividing twice
            for (int i = 0; i < 9; i++)
            {
                const double t = ebyf[i] * invFm; // ebyf/Fm
                ff_cell[i] += AFm * (W[i] * (t * t) - Bi[i]);
            }

            const double K = beta_ * ratioR * ratioB * rho_here;

            // i=0 is always just proportional split
            fR_cell[0] = ratioR * ff_cell[0];
            fB_cell[0] = ratioB * ff_cell[0];

            // i=1..8 includes perturbation term
            for (int i = 1; i < 9; i++)
            {
                const double cosphi = ebyf[i] * invFm * inv_sqrt_e_sqr[i];
                const double the_2nd_part = K * W[i] * cosphi;
                const double base = ff_cell[i]; //load ff only once

                fR_cell[i] = ratioR * base + the_2nd_part;
                fB_cell[i] = ratioB * base - the_2nd_part;
            }
        }
    }
}

static void BounceBack()
{
    for (int x = 0; x < lx; ++x)
    {
        const size_t base0_x = static_cast<size_t>(x) * static_cast<size_t>(ly);
        const size_t base9_x = static_cast<size_t>(x) * static_cast<size_t>(ly) * 9u;

        for (int y = 0; y < ly; ++y)
        {
            const size_t ys = static_cast<size_t>(y);
            const size_t id0 = base0_x + ys;
            if (obstacle[id0] != 1) continue;

            const size_t id9 = base9_x + ys * 9u;

            auto swapR = [&](size_t a, size_t b) {
                double t = f_R[id9 + a];
                f_R[id9 + a] = f_R[id9 + b];
                f_R[id9 + b] = t;
                };
            auto swapB = [&](size_t a, size_t b) {
                double t = f_B[id9 + a];
                f_B[id9 + a] = f_B[id9 + b];
                f_B[id9 + b] = t;
                };

            swapR(1u, 3u); swapB(1u, 3u);
            swapR(2u, 4u); swapB(2u, 4u);
            swapR(5u, 7u); swapB(5u, 7u);
            swapR(6u, 8u); swapB(6u, 8u);
            // 0 stays unchanged
        }
    }
}

void LBM()
{
    VoronoiObst();
    if (tunnel_r >0) CreateHorizontalTunnel(tunnel_r);
    if (mirror) MirrorTrans();
    SetBoundaryObstacles();

    InitializeDrainage();

    int step = 0;
    while (step < cycles)
    {
        Stream();

        PressureBoundary();

        GetUnRho();

        CollisionMRT();

        Redistribution();

        BounceBack();

        if ((step % checkpoint) == 0)
        {
            OutThreeKindom(step);
            OutStatus(step);
        }
        /*if ((step % frame_point) == 0)
        {
            FrameStatus(step, project);
        }*/
        step++;
    }
}