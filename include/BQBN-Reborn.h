#pragma once
#include <string>
#include <vector>
#include <filesystem>
// --- existing extern declarations ---
extern int lx;
extern int ly;
extern int n;

extern std::filesystem::path input_file;
extern std::filesystem::path output_folder;
extern std::string project;
extern int tunnel_r;
extern bool mirror;

extern int seed;
extern int cycles;
extern int checkpoint;

extern int pool_size;
extern int left_size;
extern int right_size;
extern double LR_ratio;

extern double A;
extern double df;

extern double wis_fold;
extern double wis_R;
extern double tau_R;
extern double tau_B;

extern double alpha_R;
extern double alpha_B;
extern double rho0_R;
extern double rho0_B;
extern double rho_in;
extern double rho_out;
extern double c_sqr_R;
extern double c_sqr_B;

extern double rho_wr;
extern double rho_wb;
extern double beta_;
extern double delta;
extern double sigma;

extern double st1;
extern double s2;
extern double s3;
extern double t2;
extern double t3;

extern const int ex[9];
extern const int ey[9];
extern const double e[9][2];

extern double Ci_R[9];
extern double Ci_B[9];

extern const double M[9][9];
extern double M_inv[9][9];

extern const double cc;
extern const double c_sqr;

extern const double W[9];
extern const double Bi[9];
extern const double e_sqr[9];
extern const double inv_sqrt_e_sqr[9];

extern std::vector<int> obstacle;

extern std::vector<double> f_R;
extern std::vector<double> f_B;
extern std::vector<double> ff;
extern std::vector<double> u;

extern std::vector<double> rho_R;
extern std::vector<double> rho_B;
extern std::vector<double> rho;

inline size_t IDX9(int x, int y, int i)
{
    return static_cast<size_t>(x) * static_cast<size_t>(ly) * 9u
        + static_cast<size_t>(y) * 9u
        + static_cast<size_t>(i);
}

inline size_t IDX2(int x, int y, int i)
{
    return static_cast<size_t>(x) * static_cast<size_t>(ly) * 2u
        + static_cast<size_t>(y) * 2u
        + static_cast<size_t>(i);
}

inline size_t IDX0(int x, int y)
{
    return static_cast<size_t>(x) * static_cast<size_t>(ly)
        + static_cast<size_t>(y);
}

inline size_t BELONG(int x, int y, int k)
{
    return static_cast<size_t>(x) * static_cast<size_t>(ly) * 2u
        + static_cast<size_t>(y) * 2u
        + static_cast<size_t>(k);
}

void LBM();

// --- new function to load parameters ---
void GlobalConst();
bool load_parameters(const std::string& filename);
void print_parameters();

void SetBoundaryObstacles();
void VoronoiObst();
void CreateHorizontalTunnel(int radius = 4);
void MirrorTrans();

void OutInfo(const std::optional<int>& seed);
void OutThreeKindom(int step);
void OutStatus(int step);
void OutSlice(int step, int x);