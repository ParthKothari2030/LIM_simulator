#ifndef _LUM_BUFFERED_HH
#define _LUM_BUFFERED_HH

#include <valarray>
#include  <cstdint>

void Cii_Lum_Buffered (std::valarray<float>&, uint32_t);
void CO_Lum_Buffered (std::valarray<float>&, uint32_t);
void CO_Lum_Buffered_1_new (std::valarray<float>&, uint32_t, float, float);
void HI_temp_Buffered(std::valarray<float>&, uint32_t,float, uint64_t);

#endif
