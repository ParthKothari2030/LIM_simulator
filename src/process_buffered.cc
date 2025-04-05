#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <string>
#include <valarray>

#include "cic_buffered.hh"
#include "grid_halo.hh"
#include "lum_buffered.hh"
#include "map.hh"
#include "read_halo.hh"
#include "sfr_buffered.hh"

extern std::string Dir;
extern std::string Output;
extern float z;
extern float scaling_unit;
extern float h;
extern float Omega_m;
extern double Lambda;
extern float ext_alpha,ext_beta; // -> CO param
// extern float ext_Omega_HI; // -->  HI param

extern uint64_t N_cell_x_orig;
extern uint64_t N_cell_x;

extern uint32_t Buf_sz;

inline double Hubble (float zz) {

  return 100 * h * sqrt(Omega_m * powf(1.0 + zz, 3.0) + (1.0 - Omega_m));
}

inline uint64_t Ind (uint64_t xx, uint64_t yy, uint64_t zz) {

  return zz + N_cell_x * (yy + N_cell_x * xx);
}

void Intensity (std::valarray<float> &);
void Brightness_Temp (std::valarray<float> &);
void HI_Brightness_Temp_Navarro(std::valarray<float> &);
void
LIM_Buffered ()
{
  std::ifstream num_halo {Dir + "/num_halos.txt"};

  uint64_t N_halo;
  num_halo >> N_halo;
  num_halo.close();

  uint32_t loop = N_halo / Buf_sz;
  uint32_t rem = N_halo % Buf_sz;

  std::valarray<float> Pos_x(Buf_sz), Pos_y(Buf_sz), Pos_z(Buf_sz);
  std::valarray<float> Lum(Buf_sz);

  std::array<char *, 4> buff{
      reinterpret_cast<char *>(&Pos_x[0]), reinterpret_cast<char *>(&Pos_y[0]),
      reinterpret_cast<char *>(&Pos_z[0]), reinterpret_cast<char *>(&Lum[0])};

  std::valarray<float> Map(N_cell_x * N_cell_x * N_cell_x);

  Init_Map (Map);

  std::array<std::ifstream, 4> in;

  File_Handler (in, 0);

  for (uint32_t i = 0; i < loop; ++i) {
    Read_Halos_Buffered (in, buff, Buf_sz);
    Coarse_Grid_Buffered (Pos_x, Pos_y, Pos_z, Buf_sz);
    SFR_Buffered (Lum, Buf_sz, z);
    //Cii_Lum_Buffered(Lum, Buf_sz);
    CO_Lum_Buffered_1_new (Lum, Buf_sz,ext_alpha,ext_beta);
    // HI_temp_Buffered(Lum, Buf_sz, ext_Omega_HI,N_halo);
    Cloud_in_Cell_Buffered (Pos_x, Pos_y, Pos_z, Lum, Buf_sz, Map);

  }

  Read_Halos_Buffered (in, buff, rem);
  Coarse_Grid_Buffered (Pos_x, Pos_y, Pos_z, rem);
  SFR_Buffered (Lum, rem, z);
  //Cii_Lum_Buffered (Lum, rem);
  CO_Lum_Buffered_1 (Lum, rem,ext_alpha,ext_beta);
  // HI_temp_Buffered(Lum,rem,ext_Omega_HI,N_halo);  
  Cloud_in_Cell_Buffered (Pos_x, Pos_y, Pos_z, Lum, rem, Map);

  File_Handler (in, 1);

  //Intensity(Map);
  Brightness_Temp (Map) ; // --> for CO maps
  // HI_Brightness_Temp_Navarro(Map); // --> for HI maps

  Write_Map (Map);

} // End of LIM_Buffered()

void
Intensity (std::valarray<float> &Mapp)
{
  const float solar_lum = 3.828e26; // in units of W
  const double mpc3_to_m3 = pow (3.086e22, 3.0);
  const double lum_den_fac = solar_lum / mpc3_to_m3;
  const float jansky = 1e-26;
  const float Hubble_fac = 3.24e-20;

  float fac =
      lum_den_fac * Lambda / (4 * M_PI * Hubble_fac * Hubble(z)) / jansky;
  uint_fast64_t ii, jj, kk;

#pragma omp parallel for num_threads(4) collapse(2) private (ii, jj, kk)
  for (ii = 0; ii < N_cell_x; ++ii)
    for (jj = 0; jj < N_cell_x; ++jj)
      for (kk = 0; kk < N_cell_x; ++kk)
        Mapp[Ind(ii, jj, kk)] *= fac * scaling_unit;

} // End of Intensity()

void
Brightness_Temp (std::valarray<float> &Mapp)
{
  const float solar_lum = 3.828e26; // in units of W
  const double mpc3_to_m3 = pow (3.086e22, 3.0);
  const double lum_den_fac = solar_lum / mpc3_to_m3;
  const float k_B = 1.38e-23;
  const float Hubble_fac = 3.24e-20;

  const float fac =
      lum_den_fac * pow (Lambda, 3.0) * powf (1.0 + z, 2.0) /
      (8.0 * M_PI * k_B * Hubble_fac * Hubble(z)); // Convert to Kelvin

  uint_fast64_t ii, jj, kk;

#pragma omp parallel for num_threads(4) collapse(2) private (ii, jj, kk)
  for (ii = 0; ii < N_cell_x; ++ii)
    for (jj = 0; jj < N_cell_x; ++jj)
      for (kk = 0; kk < N_cell_x; ++kk)
        Mapp[Ind(ii, jj, kk)] *= fac * scaling_unit;

} // End of Brightness_Temp()

void HI_Brightness_Temp_Navarro(std::valarray<float> &Mapp)
// returns value in mK; Villaescusa-Navarro et al. 2018
// matches with Spinelli et al. 2020 power spectrum
{
  const float G = 6.67e-11; // gravitational constant
  const float Hubble_fac = 3.24e-20;  // conversion factor from km/s/Mpc to s^-1
  const float msun_to_kg = 1.99e30;
  const double mpc3_to_m3 = pow(3.086e22, 3.0);
  float rho_c = 3 * pow(Hubble(z)*Hubble_fac,2) / (8*M_PI*G); // critical density at redshift z (SI units)
  double fac = 18900 * pow(h,2) * pow(1.0+z,2) * msun_to_kg / mpc3_to_m3;
  fac /= Hubble(z) * rho_c;
  // std::cout<<"Factor: "<<fac<<std::endl;
  int_fast64_t ii, jj, kk;

  #pragma omp parallel for num_threads(4) collapse(2) private (ii, jj, kk)
    for (ii = 0; ii < N_cell_x; ++ii)
      for (jj = 0; jj < N_cell_x; ++jj)
        for (kk = 0; kk < N_cell_x; ++kk)
          Mapp[Ind(ii, jj, kk)] *= fac ; //* scaling_unit;
} // End of HI_brightness_temp_navarro