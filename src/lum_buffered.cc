#include <cmath>
#include <cstdint>
#include <math.h>
#include <omp.h>
#include <valarray>
#include <numeric>  // Include this header for std::accumulate
#include <iostream>

extern float h;
extern float Omega_m;
extern float z;

inline double
Hubble (float zz)
{
  return 100 * 3.24e-20 * h *
         sqrt(Omega_m * powf(1.0 + zz, 3.0) + (1.0 - Omega_m));
}

inline double
sum (float x, float y)
{
  return x + static_cast<double>(y);
}


void
Cii_Lum_Buffered (std::valarray<float> &Lum, uint32_t buff_sz)
{
  uint_fast64_t i;

  const float a = 0.8475, b = 7.2203; // Silva et al. 2015, model: m1

#pragma omp parallel for num_threads(4)
  for (i = 0; i < buff_sz; ++i)
    Lum[i] = pow (10.0, a * log10 (Lum[i]) + b);
}

void
CO_Lum_Buffered (std::valarray<float> &Lum, uint32_t buff_sz)
{
  uint_fast64_t i;

  const float del_mf = 1.0;
  const float del_mf_fac = 1.0 / (del_mf * 1e-10);

  const float alpha = 1.11;
  const float alpha_inv = 1.0 / alpha;
  const float beta = 0.6;

  const int8_t J = 2;
  const float fac = 4.95e-5 * powf (1.0 * J, 3.0);

#pragma omp parallel for num_threads(4)
  for (i = 0; i < buff_sz; ++i)
    {
      Lum[i] *= del_mf_fac; // SFR_to_L-IR
      Lum[i] = powf (10.0, (log10 (Lum[i]) - beta) * alpha_inv);
      // L-IR_to_L-CO_prime

      Lum[i] *= fac; // CO(2-1) line luminosity
    }
} // End of CO_Lum_Buffered()

void
CO_Lum_Buffered_1 (std::valarray<float> &Lum, uint32_t buff_sz)
{
  uint_fast64_t i;

  const float del_mf = 1.0;
  const float del_mf_fac = 1.0 / (del_mf * 1e-10);

  const float alpha = 1.37;
  const float beta = -1.74; // Carilli & Walter 2013, Li et al. 2016

  //const float alpha = 0.81;
  //const float beta = 0.54; // Sargent et al. 2014

  const float alpha_inv = 1.0 / alpha;

  const float line_ratio = 0.76; // L_CO(2-1) / L_CO(1-0) Fonseca et al. 2016

  const int8_t J = 2;
  const float fac = 4.95e-5 * powf (1.0 * J, 3.0);

#pragma omp parallel for num_threads(4)
  for (i = 0; i < buff_sz; ++i)
    {
      Lum[i] *= del_mf_fac; // SFR_to_L-IR
      Lum[i] = powf (10, (log10 (Lum[i]) - beta) * alpha_inv); // L-IR_to_L-CO_prime (1-0)
      Lum[i] *= line_ratio;
      Lum[i] *= fac;
    }
}

void
CO_Lum_Buffered_1_new (std::valarray<float> &Lum, uint32_t buff_sz,float ext_alpha, float ext_beta)
{
  uint_fast64_t i;

  const float del_mf = 1.0;
  const float del_mf_fac = 1.0 / (del_mf * 1e-10);

  const float alpha = ext_alpha;
  const float beta = ext_beta; // Carilli & Walter 2013, Li et al. 2016
  //const float alpha = 0.81;
  //const float beta = 0.54; // Sargent et al. 2014

  const float alpha_inv = 1.0 / alpha;

  const float line_ratio = 1.0; // L_CO(2-1) / L_CO(1-0) Fonseca et al. 2016

  const int8_t J = 1;
  const float fac = 4.95e-5 * powf (1.0 * J, 3.0);

#pragma omp parallel for num_threads(4)
  for (i = 0; i < buff_sz; ++i)
    {
      Lum[i] *= del_mf_fac; // SFR_to_L-IR
      Lum[i] = powf (10, (log10 (Lum[i]) - beta) * alpha_inv); // L-IR_to_L-CO_prime (1-0)
      Lum[i] *= line_ratio;
      Lum[i] *= fac;
    }
}


void 
HI_temp_Buffered(std::valarray<float> &Lum, u_int32_t buff_sz, float ext_omega_HI,uint64_t N_halo)
{
    uint_fast64_t i;
    //  Bagla et al. 2010 ; scheme three

    // Mass conversion
    Lum *= 1e10/h; 

    // printf("LUM_1 = %f\n", Lum[0]); // --> Conversion cross check

    // Velocity limits in km/s
    const float v_min = 30.0;
    const float v_max = 200.0;

    // Constants for mass calculations (Solar mass)
    float M_min = powf((v_min / 60.0), 3.0) * powf(( (1 + z) / 4.0 ), -1.5) * 1e10;
    float M_max = powf((v_max / 60.0), 3.0) * powf(( (1 + z) / 4.0 ), -1.5) * 1e10;
    
    // Simulation volume in Mpc^3
    const float Vol_Sim = powf(302.6,3);

    // Physical constants
    const float mpc_to_m = 3.086e22;       // Conversion factor from Mpc to meters
    const float H0 = 67740.0 / mpc_to_m;  // Hubble constant in s^-1
    const float G = 6.67430e-11;           // Gravitational constant in m^3/kg/s^2
    const float solar_mass = 1.989e30;     // Solar mass in kg
    

    // Critical density calculations
    const double rho_critical = (3 * pow(H0, 2)) / (8 * M_PI * G); // Critical density in kg/m^3
    const double rho_critical_new = (rho_critical / solar_mass) * pow(mpc_to_m, 3); // Critical density in Solar Mass/Mpc^3

    // Apply HI_recipe logic to the Lum array ; Bagla et al. 2010, Debanjan et al. 2016
    int Halo_num = 0;
    
    #pragma omp parallel for num_threads(4)  
      for (i = 0; i < buff_sz; i++) 
        {
            if (Lum[i] >= M_min) 
              {
                  Lum[i] = Lum[i] / (1 + (Lum[i] / M_max));
                  Halo_num ++;
              } 
            else 
              {
                  Lum[i] = 0.0; 
              }
        }
    
    
    // Sum the transformed Lum array  ; use accumulate function of C++
    float Lum_sum_Acu = std::accumulate(std::begin(Lum), std::end(Lum), 0.0f, sum);

    // std::cout << "Num halo in LUM : " << N_halo << std::endl;
    // printf("Lum_sum_acu = %f\n",Lum_sum_Acu);
    
    // Calculate scaling factor f3 (prop to Omega_HI)
    float f3 = (rho_critical_new * ext_omega_HI * Vol_Sim) / Lum_sum_Acu;  
    // std::cout<<"f3: "<<f3<<std::endl; // problem solved
    
    // Scale the Lum array by f3
    #pragma omp parallel for num_threads(4) // loop for f3 reassignment
      for (i = 0; i < buff_sz; ++i)
        {
          Lum[i] = Lum[i] * f3;
        }
}
