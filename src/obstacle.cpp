#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <random>

#include "BQBN-Reborn.h"

using namespace std;

static int Uniform(int MIN, int MAX);
static double p(int x);
static bool PointOfChoice(int x);
static double q(int x);

void SetBoundaryObstacles()
{
    for (int y = 0; y < ly; y++)
    {
        if (y == 0 || y == ly - 1)
        {
            for (int x = 0; x < lx; x++)
            {
                obstacle[IDX0(x, y)] = 1;
            }
        }
    }
}

void VoronoiObst()
{
    // 1. Generate seed cores
    std::vector<int> core_x(n), core_y(n);
    int obst = 0;

    while (obst < n)
    {
        int x = Uniform(pool_size, lx - pool_size - 1);
        if (PointOfChoice(x))
        {
            core_x[obst] = x;
            core_y[obst] = Uniform(1, ly - 2);
            obst++;
        }
    }

    // 2. Initialize obstacles to 1 (full obstacle)
    std::fill(obstacle.begin(), obstacle.end(), 1);

    // 3. Assign each cell to the two nearest seeds using squared distance
    std::vector<int> belong(lx * ly * 2, 0); // flattened 3D: x*ly*2 + y*2 + k
    for (int y = 0; y < ly; y++)
    {
        for (int x = 0; x < lx; x++)
        {
            int d2_min_1 = lx * lx + ly * ly;
            int d2_min_2 = lx * lx + ly * ly;

            for (int i = 0; i < n; i++)
            {
                int dx = x - core_x[i];
                int dy = y - core_y[i];
                int d2 = dx * dx + dy * dy; // squared distance

                if (d2 < d2_min_1)
                {
                    d2_min_2 = d2_min_1;
                    belong[IDX2(x, y, 1)] = belong[IDX2(x, y, 0)];

                    d2_min_1 = d2;
                    belong[IDX2(x, y, 0)] = i;
                }
                else if (d2 < d2_min_2)
                {
                    d2_min_2 = d2;
                    belong[IDX2(x, y, 1)] = i;
                }
            }
        }
    }

    // 4. Mark cells near borders between seeds as non-obstacle
    for (int x = 0; x < lx; x++)
    {
        for (int y = 0; y < ly; y++)
        {
            for (int i = 0; i < 9; i++)
            {
                int nx = (x + ex[i] + lx) % lx;
                int ny = (y + ey[i] + ly) % ly;
                if (belong[IDX2(nx, ny, 0)] != belong[IDX2(x, y, 0)])
                {
                    obstacle[IDX0(x, y)] = 0;
                    break; // exit the i-loop early
                }
            }
        }
    }

    // 5. Bolding: expand obstacles probabilistically
    std::vector<int> temp_obst = obstacle;
    std::vector<int> done(lx * ly, 0);
    int MAX = 12000 * (lx - pool_size * 2);
    int scale = 1;

    for (int x = 0; x < lx; x++)
    {
        int threshold = static_cast<int>(q(x) * MAX * scale);
        for (int y = 0; y < ly; y++)
        {
            bool if_not_obst = obstacle[IDX0(x, y)] == 0;
            bool if_not_done = done[IDX0(x, y)] == 0;

            if (if_not_obst && if_not_done)
            {
                if (Uniform(1, MAX) <= threshold)
                {
                    bool single = Uniform(1, 100) <= 30;

                    for (int yp = 0; yp < ly; yp++)
                    {
                        for (int xp = 0; xp < lx; xp++)
                        {
                            bool if_obst_p = obstacle[IDX0(xp, yp)] == 0;
                            bool if_same = belong[IDX2(xp, yp, 0)] == belong[IDX2(x, y, 0)]
                                && belong[IDX2(xp, yp, 1)] == belong[IDX2(x, y, 1)];

                            if (if_obst_p && if_same)
                            {
                                if (single)
                                {
                                    for (int i = 0; i < 9; i++)
                                    {
                                        int nx = (xp + ex[i] + lx) % lx;
                                        int ny = (yp + ey[i] + ly) % ly;
                                        temp_obst[IDX0(nx, ny)] = 0;
                                    }
                                    done[IDX0(xp, yp)] = 1;
                                }
                                else
                                {
                                    for (int i = 0; i < 9; i++)
                                    {
                                        int nx = (xp + ex[i] + lx) % lx;
                                        int ny = (yp + ey[i] + ly) % ly;
                                        for (int j = 0; j < 9; j++)
                                        {
                                            int nnx = (nx + ex[j] + lx) % lx;
                                            int nny = (ny + ey[j] + ly) % ly;
                                            temp_obst[IDX0(nnx, nny)] = 0;
                                        }
                                    }
                                    done[IDX0(xp, yp)] = 1;
                                }
                            }
                        }
                    }
                }
                else // mark done for remaining cells
                {
                    for (int yp = 0; yp < ly; yp++)
                    {
                        for (int xp = 0; xp < lx; xp++)
                        {
                            bool if_obst_p = obstacle[IDX0(xp, yp)] == 0;
                            bool if_same = belong[IDX2(xp, yp, 0)] == belong[IDX2(x, y, 0)]
                                && belong[IDX2(xp, yp, 1)] == belong[IDX2(x, y, 1)];
                            if (if_obst_p && if_same)
                                done[IDX0(xp, yp)] = 1;
                        }
                    }
                }
            }
            done[IDX0(x, y)] = 1;
        }
    }

    // 6. Copy temp_obst back to obstacle
    obstacle = temp_obst;

    // 7. Clear obstacles near left/right pools
    for (int y = 0; y < ly; y++)
    {
        for (int x = 0; x < pool_size; x++)
            obstacle[IDX0(x, y)] = 0;
        for (int x = lx - pool_size; x < lx; x++)
            obstacle[IDX0(x, y)] = 0;
    }
}

void CreateHorizontalTunnel(int radius)
{
    int mid = ly / 2;
    for (int x = 0; x < lx; x++)
    {
        for (int y = mid - radius; y <= mid + radius; y++)
        {
            obstacle[IDX0(x, y)] = 0;
        }
    }
}

void MirrorTrans()
{
    std::vector<int> temp(obstacle.size(), 0);

    for (int x = 0; x < lx; x++)
    {
        for (int y = 0; y < ly; y++)
        {
            temp[IDX0(lx - 1 - x, y)] = obstacle[IDX0(x, y)];
        }
    }

    obstacle = temp;
}

static int Uniform(int MIN, int MAX)
{
    return (rand() % (MAX - MIN + 1)) + MIN;
}

static double p(int x)
{
    double x_real = x - pool_size + 1;
    double arena = lx - pool_size * 2;
    double mid = arena - left_size - right_size;

    double dp = (LR_ratio - 1) / (LR_ratio + 1) / arena;

    double p0 = 1 / arena + dp;
    double pn = 1 / arena - dp;
    double k = (pn - p0) / (mid + 1);

    double p = 0;

    if (x_real <= left_size)
    {
        p = p0;	
    }
    else if (x_real <= arena - right_size)
    {
        p = p0 + k * (x_real - left_size);
    }
    else
    {
        p = pn;	
    }
    return p;
}

static bool PointOfChoice(int x)
{
    return (static_cast<double>(rand()) / RAND_MAX) < p(x);
}

static double q(int x)
{
    double x_real = x - pool_size + 1;
    double arena = lx - pool_size * 2;
    double mid = arena - left_size - right_size;

    double q0 = 0;
    double qn = 0;
    if (LR_ratio < 1)
    {
        q0 = 2 / arena;
        qn = 0;
    }
    else if (LR_ratio == 1)
    {
        q0 = 1 / arena;
        qn = 1 / arena;
    }
    else if (LR_ratio > 1)
    {
        q0 = 2 / arena;
        qn = 0;
    }

    double k = (qn - q0) / (mid + 1);

    double q = 0;

    if (x_real <= left_size)
    {
        q = q0;
    }
    else if (x_real <= arena - right_size)
    {
        q = q0 + k * (x_real - left_size);
    }
    else
    {
        q = qn;
    }
    return q;
}
