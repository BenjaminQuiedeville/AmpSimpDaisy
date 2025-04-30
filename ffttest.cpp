#include <stdint.h>
#include "shy_fft.h"

#include "pffft/pffft.c"

int main() {
    
    static const int fftsize = 32;
    
    PFFFT_Setup *setup = pffft_new_setup(fftsize, PFFFT_REAL);
    ShyFFT<float, fftsize> fft;
    
    alignas(64) float pffft_input[fftsize] = {0};
    pffft_input[0] = 1.0f; 
    pffft_input[1] = -1.0f; 
    
    alignas(64) float pffft_dft[fftsize+2] = {0};
    
    alignas(16) float input[fftsize] = {0};
    alignas(16) float dft[fftsize] = {0};
    input[0] = 1.0f;
    input[1] = -1.0f;
    // input[2] = 0.5f;
    // input[3] = -0.25f;
    

    
    fft.Init();
    
    fft.Direct(input, dft);
    pffft_transform(setup, pffft_input, pffft_dft, nullptr, PFFFT_FORWARD);

    fft.Inverse(dft, input);
    pffft_transform(setup, pffft_dft, pffft_input, nullptr, PFFFT_BACKWARD);
    
    for (int index = 0; index < fftsize; index++) {
        input[index] /= (float)fftsize;
        pffft_input[index] /= (float)fftsize;
    }

    return 0;
}
