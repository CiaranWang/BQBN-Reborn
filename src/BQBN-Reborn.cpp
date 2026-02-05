#include <iostream>
#include <sstream>
#include <fstream>
#include <optional>
#include <filesystem>
#include <unordered_map>

#include "BQBN-Reborn.h"
using namespace std;

// ----------------------------
// Globals definition
// ----------------------------
int lx = 0;
int ly = 0;
int n = 0;

std::filesystem::path input_file;
std::filesystem::path output_folder;
std::string project;
int tunnel_r = 0;
bool mirror = false;

int cycles = 1000;
int checkpoint = 200;

int pool_size = 0;
int left_size = 0;
int right_size = 0;
double LR_ratio = 0.0;

double A = 0.0;
double df = 0.0;

double wis_fold = 0.0;
double wis_R = 0.0;
double tau_R = 0.0;
double tau_B = 0.0;

double alpha_R = 0.0;
double alpha_B = 0.0;
double rho0_R = 0.0;
double rho0_B = 0.0;
double rho_in = 0.0;
double rho_out = 0.0;
double c_sqr_R = 0.0;
double c_sqr_B = 0.0;

double rho_wr = 0.0;
double rho_wb = 0.0;
double beta_ = 0.0;
double delta = 0.0;
double sigma = 0.0;

double st1 = 0.0;
double s2 = 0.0;
double s3 = 0.0;
double t2 = 0.0;
double t3 = 0.0;

// Lattice constants
const int ex[9]{ 0,1,0,-1,0,1,-1,-1,1 };
const int ey[9]{ 0,0,1,0,-1,1,1,-1,-1 };
const double e[9][2] = {
    {  0.0,  0.0 },
    {  1.0,  0.0 },
    {  0.0,  1.0 },
    { -1.0,  0.0 },
    {  0.0, -1.0 },
    {  1.0,  1.0 },
    { -1.0,  1.0 },
    { -1.0, -1.0 },
    {  1.0, -1.0 }
};

double Ci_R[9]{};
double Ci_B[9]{};

const double M[9][9] = {
    {  1,  1,  1,  1,  1,  1,  1,  1,  1 },
    { -4, -1, -1, -1, -1,  2,  2,  2,  2 },
    {  4, -2, -2, -2, -2,  1,  1,  1,  1 },
    {  0,  1,  0, -1,  0,  1, -1, -1,  1 },
    {  0, -2,  0,  2,  0,  1, -1, -1,  1 },
    {  0,  0,  1,  0, -1,  1,  1, -1, -1 },
    {  0,  0, -2,  0,  2,  1,  1, -1, -1 },
    {  0,  1, -1,  1, -1,  0,  0,  0,  0 },
    {  0,  0,  0,  0,  0,  1, -1,  1, -1 }
};
double M_inv[9][9] = { 0 };

const double cc = 1.0;
const double c_sqr = cc * cc / 3.0;

const double W[9]{ 4.0 / 9,1.0 / 9,1.0 / 9,1.0 / 9,1.0 / 9,1.0 / 36,1.0 / 36,1.0 / 36,1.0 / 36 };
const double Bi[9]{ -4.0 / 27,2.0 / 27,2.0 / 27,2.0 / 27,2.0 / 27,5.0 / 108,5.0 / 108,5.0 / 108,5.0 / 108 };
const double e_sqr[9]{ 0.0,1.0,1.0,1.0,1.0,2.0,2.0,2.0,2.0 };
const double inv_sqrt_e_sqr[9]{
        0.0,                 // i=0 unused (cosphi not used)
        1.0, 1.0, 1.0, 1.0,  // |e| = 1
        0.7071067811865475,  // 1/sqrt(2)
        0.7071067811865475,
        0.7071067811865475,
        0.7071067811865475
};

std::vector<int> obstacle;
std::vector<double> f_R, f_B, ff, u;
std::vector<double> rho_R, rho_B, rho;

// ----------------------------
// Initialize M_inv
// ----------------------------
void GlobalConst()
{
    double M_inv_pre[9][9] = {
        {  4, -4,  4,  0,  0,  0,  0,  0,  0 },
        {  4, -1, -2,  6, -6,  0,  0,  9,  0 },
        {  4, -1, -2,  0,  0,  6, -6, -9,  0 },
        {  4, -1, -2, -6,  6,  0,  0,  9,  0 },
        {  4, -1, -2,  0,  0, -6,  6, -9,  0 },
        {  4,  2,  1,  6,  3,  6,  3,  0,  9 },
        {  4,  2,  1, -6, -3,  6,  3,  0, -9 },
        {  4,  2,  1, -6, -3, -6, -3,  0,  9 },
        {  4,  2,  1,  6,  3, -6, -3,  0, -9 }
    };

    double factor = 1.0 / 36.0;
    for (int i = 0; i < 9; i++)
        for (int j = 0; j < 9; j++)
            M_inv[i][j] = factor * M_inv_pre[i][j];
}

// ----------------------------
// Simplified .ini parser
// ----------------------------
bool load_parameters(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: cannot open " << filename << "\n";
        return false;
    }

    std::string line, section;
    while (std::getline(file, line)) {
        // remove comments
        if (auto pos = line.find('#'); pos != std::string::npos)
            line = line.substr(0, pos);

        // trim
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);
        if (line.empty()) continue;

        // section
        if (line.front() == '[' && line.back() == ']') {
            section = line.substr(1, line.size() - 2);
            continue;
        }

        // key = value
        if (auto pos = line.find('='); pos != std::string::npos) {
            std::string key = line.substr(0, pos);
            std::string val = line.substr(pos + 1);
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            val.erase(0, val.find_first_not_of(" \t"));
            val.erase(val.find_last_not_of(" \t") + 1);

            double dval = std::stod(val);
            int ival = static_cast<int>(dval);

            // map key to global
            if (section == "grid") {
                if (key == "lx") lx = ival;
                else if (key == "ly") ly = ival;
                else if (key == "n") n = ival;
                else if (key == "pool_size") pool_size = ival;
                else if (key == "LR_ratio") LR_ratio = dval;
            }
            else if (section == "constants") {
                if (key == "A") A = dval;
                else if (key == "df") df = dval;
                else if (key == "wis_fold") wis_fold = dval;
                else if (key == "wis_R") wis_R = dval;
                else if (key == "alpha_R") alpha_R = dval;
                else if (key == "alpha_B") alpha_B = dval;
                else if (key == "rho_wr") rho_wr = dval;
                else if (key == "rho_wb") rho_wb = dval;
                else if (key == "beta") beta_ = dval;
                else if (key == "delta") delta = dval;
                else if (key == "sigma") sigma = dval;
            }
        }
    }

    // derived quantities
    left_size = (lx - pool_size * 2) / 4;
    right_size = (lx - pool_size * 2) / 4;
    tau_R = 0.5 + wis_R;
    tau_B = 0.5 + wis_R * wis_fold;
    rho0_R = 1 - alpha_B;
    rho0_B = 1 - alpha_R;
    rho_in = rho0_R + 0.1;
    rho_out = rho0_B - 0.1;
    c_sqr_R = 0.6 * (1 - alpha_R);
    c_sqr_B = 0.6 * (1 - alpha_B);
    st1 = 2 * tau_R * tau_B / (tau_R + tau_B);
    s2 = 2 * (tau_R - st1) / delta;
    s3 = -s2 / (2 * delta);
    t2 = 2 * (st1 - tau_B) / delta;
    t3 = t2 / (2 * delta);

    Ci_R[0] = alpha_R;
    for (int i = 1; i <= 4; ++i) Ci_R[i] = (1 - alpha_R) / 5;
    for (int i = 5; i <= 8; ++i) Ci_R[i] = (1 - alpha_R) / 20;
    
    Ci_B[0] = alpha_B;
    for (int i = 1; i <= 4; ++i) Ci_B[i] = (1 - alpha_B) / 5;
    for (int i = 5; i <= 8; ++i) Ci_B[i] = (1 - alpha_B) / 20;

    obstacle.assign(lx * ly, 1);

    f_R.assign(lx * ly * 9, 0.0);
    f_B.assign(lx * ly * 9, 0.0);
    ff.assign(lx * ly * 9, 0.0);
    u.assign(lx * ly * 2, 0.0);

    rho_R.assign(lx * ly, 0.0);
    rho_B.assign(lx * ly, 0.0);
    rho.assign(lx * ly, 0.0);

    return true;
}

void print_parameters() {
    std::cout << "lx = " << lx << ", ly = " << ly << "\n";
    std::cout << "pool_size = " << pool_size << ", LR_ratio = " << LR_ratio << "\n";
    std::cout << "A = " << A << ", df = " << df << "\n";
    std::cout << "alpha_R = " << alpha_R << ", alpha_B = " << alpha_B << "\n";
    std::cout << "rho_wr = " << rho_wr << ", rho_wb = " << rho_wb << "\n";
}