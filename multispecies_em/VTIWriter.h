#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#define _USE_MATH_DEFINES
#include <cstring>

using Vec1D = std::vector<double>;
using Vec2D = std::vector<Vec1D>;
using Vec3D = std::vector<Vec2D>;

class VTIWriter {
public:
    static void write2DVTI(const std::string &filename,
                           const Vec2D &Ex,
                           const Vec2D &Ey,
                           const Vec2D &Bz,
                           const Vec2D &phi,
                           const Vec2D &Jx,
                           const Vec2D &Jy,
                           const Vec2D &q_rho,
                           const Vec2D &rho_e,
                           const Vec2D &rho_i,
                           const Vec2D &ux_e,
                           const Vec2D &uy_e,
                           const Vec2D &ux_i,
                           const Vec2D &uy_i)
    {
        int NX = Ex.size();
        int NY = Ex[0].size();

        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error opening file " << filename << std::endl;
            return;
        }

        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        file << "  <ImageData WholeExtent=\"0 " << NX-1 << " 0 " << NY-1 << " 0 0\" "
             << "Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
        file << "    <Piece Extent=\"0 " << NX-1 << " 0 " << NY-1 << " 0 0\">\n";
        file << "      <PointData>\n";

        // --- Helpers ---
        auto writeScalar = [&](const Vec2D &field, const std::string &name){
            file << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
            for (int j = 0; j < NY; ++j) {
                for (int i = 0; i < NX; ++i) {
                    file << std::setprecision(12) << field[i][j] << " ";
                }
                file << "\n";
            }
            file << "        </DataArray>\n";
        };

        auto writeVector = [&](const Vec2D &vx, const Vec2D &vy, const std::string &name){
            file << "        <DataArray type=\"Float64\" Name=\"" << name 
                 << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
            for (int j = 0; j < NY; ++j) {
                for (int i = 0; i < NX; ++i) {
                    file << std::setprecision(12) << vx[i][j] << " " << vy[i][j] << " 0.0 ";
                }
                file << "\n";
            }
            file << "        </DataArray>\n";
        };

        // --- Scalars ---
        writeScalar(phi, "phi");
        writeScalar(q_rho, "q_rho");
        writeScalar(rho_e, "rho_m");
        writeScalar(rho_i, "rho_p");

        // --- Vectors ---
        writeVector(ux_e, uy_e, "U_m");  // electron velocity
        writeVector(ux_i, uy_i, "U_p");  // ion velocity
        writeVector(Ex, Ey, "E");        // electric field
        writeVector(Jx, Jy, "J");        // current density

        // Magnetic field B is only z-component
        file << "        <DataArray type=\"Float64\" Name=\"B\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int j = 0; j < NY; ++j) {
            for (int i = 0; i < NX; ++i) {
                file << "0.0 0.0 " << std::setprecision(12) << Bz[i][j] << " ";
            }
            file << "\n";
        }
        file << "        </DataArray>\n";

        file << "      </PointData>\n";
        file << "      <CellData>\n";
        file << "      </CellData>\n";
        file << "    </Piece>\n";
        file << "  </ImageData>\n";
        file << "</VTKFile>\n";

        file.close();
        std::cout << "Written VTI file: "
              << "\033[32m"   // switch to green
              << filename
              << "\033[0m"    // reset colour
              << std::endl;
    }
};
