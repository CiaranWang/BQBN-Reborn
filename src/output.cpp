#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <string>
//#include <opencv2/opencv.hpp>
#include "BQBN-Reborn.h"

using namespace std;
//using namespace cv;
extern std::filesystem::path output_folder; 

void OutInfo(const std::optional<int>& seed)
{
    // Build the file path properly
    std::filesystem::path file_path = output_folder / (project + ".info");

    // Open the file
    std::ofstream file(file_path);
    if (!file.is_open())
    {
        std::cerr << "Failed to open file: " << file_path << std::endl;
        return;
    }
    
    // Write seed safely
    file << "seed" << "\t";
    if (seed.has_value())
        file << seed.value() << "\n";
    else
        file << "none\n";

    file << "lx" << "\t" << lx << "\n";
    file << "ly" << "\t" << ly << "\n";
    file << "Obst_nr" << "\t" << n << "\n";
    file << "pool_size" << "\t" << pool_size << "\n";
    file << "left_size" << "\t" << left_size << "\n";
    file << "right_size" << "\t" << right_size << "\n";
    file << "L:R ratio\t1:" << 1 / LR_ratio << "\n";
    file << "A" << "\t" << A << "\n";
    file << "df" << "\t" << df << "\n";
    file << "wis_R" << "\t" << wis_R << "\n";
    file << "wis_fold" << "\t" << wis_fold << "\n";
    file << "tau_R" << "\t" << tau_R << "\n";
    file << "tau_B" << "\t" << tau_B << "\n";
    file << "alpha_R" << "\t" << alpha_R << "\n";
    file << "alpha_B" << "\t" << alpha_B << "\n";
    file << "rho_in" << "\t" << rho_in << "\n";
    file << "rho_out" << "\t" << rho_out << "\n";
    file << "rho_wr" << "\t" << rho_wr << "\n";
    file << "rho_wb" << "\t" << rho_wb << "\n";
    file << "beta" << "\t" << beta_ << "\n";
    file << "delta" << "\t" << delta << "\n";
    file << "sigma" << "\t" << sigma << "\n";

    file.close();
    std::cout << "Info saved to: " << file_path << std::endl;
}

void OutThreeKindom(int step)
{
    const std::filesystem::path file_name =
        output_folder / (project + "_" + std::to_string(step) + "_Saturation.plt");

    std::ofstream file(file_name);
    if (!file)
    {
        std::cerr << "Error opening file: " << file_name << "\n";
        return;
    }

    file << "variables = x, porosity, red_saturation, blue_saturation\n";

    const double inv_ly = 1.0 / static_cast<double>(ly);

    for (int x = 0; x < lx; ++x)
    {
        const size_t base = static_cast<size_t>(x) * static_cast<size_t>(ly);

        double obst = 0.0, red = 0.0, blue = 0.0;

        for (int y = 0; y < ly; ++y)
        {
            const size_t idx = base + static_cast<size_t>(y);

            if (obstacle[idx] == 1)
            {
                obst += 1.0;
            }
            else
            {
                // classify fluid cell
                red += (rho_R[idx] >= rho_B[idx]) ? 1.0 : 0.0;
                blue += (rho_R[idx] < rho_B[idx]) ? 1.0 : 0.0;
            }
        }

        const double denom = static_cast<double>(ly) - obst;

        const double porosity = 1.0 - obst * inv_ly;
        const double red_sat = (denom > 0.0) ? (red / denom) : 0.0;
        const double blue_sat = (denom > 0.0) ? (blue / denom) : 0.0;

        file << x << '\t' << porosity << '\t' << red_sat << '\t' << blue_sat << '\n';
    }

    std::cout << "Saved: " << file_name << "\n";
}

void OutStatus(int step)
{
    const std::filesystem::path file_name =
        output_folder / (project + "_" + std::to_string(step) + "_Status.plt");

    std::ofstream file(file_name);
    if (!file)
    {
        std::cerr << "Error opening file: " << file_name << '\n';
        return;
    }

    file << "variables = x, y, stat, u_x, u_y, rhoR, rhoB, rho, obst\n";
    file << "zone i=" << lx << ", j=" << ly << ", f=point\n";

    for (int x = 0; x < lx; x++)
    {
        const size_t base0_x = static_cast<size_t>(x) * static_cast<size_t>(ly);
        const size_t base2_x = base0_x * 2u;

        for (int y = 0; y < ly; y++)
        {
            const size_t ys = static_cast<size_t>(y);
            const size_t id0 = base0_x + ys;       // IDX0(x,y)
            const size_t id2 = base2_x + ys * 2u;  // IDX2(x,y,0) is id2+0, IDX2(x,y,1) is id2+1

            const int obst = obstacle[id0];

            int stat;
            if (obst == 1)          stat = 1;   // obstacle
            else if (rho_R[id0] >= rho_B[id0]) stat = 2; // red dominates
            else                    stat = 0;   // blue dominates

            file << x << '\t'
                << y << '\t'
                << stat << '\t'
                << u[id2 + 0u] << '\t'
                << u[id2 + 1u] << '\t'
                << rho_R[id0] << '\t'
                << rho_B[id0] << '\t'
                << rho[id0] << '\t'
                << obst << '\n';
        }
    }

    std::cout << file_name << '\n';
}

void OutSlice(int step, int x)
{
    // Safety check (optional but recommended)
    if (x < 0 || x >= lx)
    {
        std::cerr << "OutSlice: x out of range: " << x << std::endl;
        return;
    }

    // Construct output file path
    std::filesystem::path file_name =
        output_folder / (project + "_" + std::to_string(step) + "_" + std::to_string(x) + "_Slice.plt");

    std::ofstream file(file_name);
    if (!file)
    {
        std::cerr << "Error opening file: " << file_name << std::endl;
        return;
    }

    file << "variables = y, ux, uy, rho_R, rho_B, rho, obst\n";

    for (int y = 0; y < ly; y++)
    {
        size_t id = IDX0(x, y);

        file << y << "\t"
            << u[IDX2(x, y, 0)] << "\t"
            << u[IDX2(x, y, 1)] << "\t"
            << rho_R[id] << "\t"
            << rho_B[id] << "\t"
            << rho[id] << "\t"
            << obstacle[id] << "\n";
    }

    file.close();
    std::cout << file_name << std::endl;
}

/*void FrameStatus(int step, string project)
{
    string file_name = "";
    file_name = path + "\\" + project + "_" + to_string(step) + "_Shot" + ".png";
    Mat img(l_y, l_x, CV_8UC3, Scalar(0,0,0));

    for (int x = 0; x < l_x; x++)
    {
        for (int y = 0; y < l_y; y++)
        {
            if (obstacle[x][y] == 0 && rho_R[x][y] >= rho_B[x][y])
            {
                img.at<Vec3b>(y, x)[0] = 0;
                img.at<Vec3b>(y, x)[1] = 0;
                img.at<Vec3b>(y, x)[2] = 255;
            }
            else if (obstacle[x][y] == 0 && rho_R[x][y] < rho_B[x][y])
            {
                img.at<Vec3b>(y, x)[0] = 255;
                img.at<Vec3b>(y, x)[1] = 0;
                img.at<Vec3b>(y, x)[2] = 0;
            }
        }
    }
    resize(img, img, Size(l_x, l_y));
    imwrite(file_name, img);
    cout << file_name << endl;
}

void AV_Maker(string project,int frame_point) // Academic Video Maker, don't over-think!
{
    string file_name = "";
    file_name = path + "\\" + project + "_Video" + ".mp4";
    VideoWriter video(file_name, VideoWriter::fourcc('P','I','M','1'), 25.0, Size(l_x, l_y),1);

    for (int i = 0; i<1000000, i++;)
    {
        string img_path = path + "\\" + project + "_" + to_string(frame_point * i) + "_Shot" + ".png";
        Mat img = imread(img_path);
        if (img.empty())
        {
            break;
        }
        else
        {
            video << img;
            cout << "processing " << to_string(frame_point * i) << endl;
        }
        //resize(img, img, Size(1920, 1080));

    }
    cout << file_name << endl;
}*/