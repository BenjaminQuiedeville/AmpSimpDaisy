#include <stdint.h>

#include "daisy_seed.h"

#include "daisysp.h"

#define PFFFT_SIMD_DISABLE
#include "pffft/pffft.c"

// #include "shy_fft.h"

#define local_const static const
#define global_const static const 
#define local_persist static
#define global static

typedef uint32_t u32;
typedef int32_t i32;
typedef uint16_t u16;
typedef int16_t i16;
typedef uint8_t u8;
typedef int8_t i8;


using namespace daisy;
using namespace daisysp;

global_const u16 BLOCK_SIZE = 16;
global_const u16 UPSAMPLE_FACTOR = 4;
global_const u16 UP_BLOCK_SIZE = BLOCK_SIZE * UPSAMPLE_FACTOR;

global_const float tube_gain = 100.0f;

// float ir[ir_size]
#include "data/base_ir.h"

global_const u16 npartitions = ir_size / BLOCK_SIZE;
global_const u16 fft_size = 32;


// global float signal_dft_buffer[fft_size] = {0};

global_const u16 dft_size = fft_size + 2;
alignas(16) global float signal_time_buffer[fft_size] = {0};
alignas(16) global float ols_buffer[fft_size] = {0};
alignas(16) global float *fdl_ptrs[npartitions] = {0}; 
alignas(16) global float fdl[npartitions][dft_size] = {0};
alignas(16) global float ir_parts_dft[npartitions][dft_size] = {0};
alignas(16) global float convolution_dft[dft_size] = {0};
// alignas(16) global float convolution_result[fft_size] = {0};


global u16 monitor_index = 0;
global_const u16 monitor_buffer_size = BLOCK_SIZE * 60;
global float monitor_buffer[monitor_buffer_size] = {0};


enum KnobsCodes {

    Knob1,  
    Knob2,  
    Knob3,  
    Knob4,  
    Knob5,  
    Knob6,  
    
    nKNOBS,
};


enum SwitchCodes {
    Toggle4,
    Toggle3,
    Toggle2,
    Toggle1,

    FS1,
    FS2,
    NSWITCHES,
};


DaisySeed  hardware;
Led led1;
Led led2;

Switch switches[NSWITCHES];


static inline float dbtoa(float db) { return powf(10.0f, db * 0.05f); }
static inline float atodb(float amp) { return 20.0f * log10f(amp); }

static inline float scale(float x, float min, float max, float newmin, float newmax, float curve) {
    return powf((x - min) / (max - min), curve) * (newmax - newmin) + newmin;
}

static inline float scale_linear(float x, float min, float max, float newmin, float newmax) {
    return (x - min) / (max - min) * (newmax - newmin) + newmin;
}

static inline float clip(float x, float min, float max) {
    return x > max ? max : x < min ? min : x;
}

static inline void apply_gain_linear(float gain, float *buffer, u32 nSamples) {
    for (u32 index = 0; index < nSamples; index++) {
        buffer[index] *= gain;
    }
}

// -------------------- Smooth parameters -------------------- 

struct SmoothParam {
    float current_value = 0.0f;
    float target = 0.0f;
    float b0 = 0.0f;
    float a1 = 0.0f;
    float y1 = 0.0f;
};

void smooth_param_init(SmoothParam *param, float init_value) {
    param->current_value = init_value;
    param->target = init_value;
    param->y1 = init_value;
    param->b0 = 1.0f;
    param->a1 = 0.0f;
}

void smooth_param_new_target(SmoothParam *param, float new_target, float tau_ms, float samplerate) {
    param->b0 = sinf(PI_F / (samplerate * tau_ms * 0.001f));
    param->a1 = param->b0 - 1.0f;
    param->target = new_target;
}

float smooth_param_next_value(SmoothParam *param) {
    param->current_value = param->target * param->b0 - param->y1 * param->a1;
    param->y1 = param->current_value;

    return param->current_value;
}


// -------------------- Biquad Filter -------------------- 

enum BiquadType : u8 {
    BIQUAD_LOWPASS = 0,
    BIQUAD_HIGHPASS,
    BIQUAD_PEAK,
    BIQUAD_LOWSHELF,
    BIQUAD_HIGHSHELF,
    BIQUAD_NFILTERTYPES,
};

struct Biquad {

    float b0 = 1.0;
    float b1 = 0.0;
    float b2 = 0.0;
    float a1 = 0.0;
    float a2 = 0.0;

    float w1 = 0.0f;
    float w2 = 0.0f;
};


void biquad_reset(Biquad *f) {
    f->w1 = 0.0f;
    f->w2 = 0.0f;
}

void biquad_set_coeffs(Biquad *f, float frequency, float Q, float gaindB, float samplerate, BiquadType type) {
        
    float w0 = 2.0f * PI_F / samplerate * frequency;
    float cosw0 = cosf(w0);
    float sinw0 = sinf(w0);

    float alpha = sinw0/(2.0f*Q);

    float A = 0.0f;
    float a0inv = 0.0f;
        
    switch (type) {
        case BIQUAD_LOWPASS: {
            a0inv = 1.0f/(1.0f + alpha);

            f->b0 = (1.0f - cosw0) * 0.5f * a0inv;
            f->b1 = 2.0f * f->b0;
            f->b2 = f->b0;
            f->a1 = -2.0f * cosw0 * a0inv;
            f->a2 = (1.0f - alpha) * a0inv;
            break;
        }

        case BIQUAD_HIGHPASS: {
            a0inv = 1.0f/(1.0f + alpha);

            f->b0 = (1.0f + cosw0) * 0.5f * a0inv;
            f->b1 = -2.0f * f->b0;
            f->b2 = f->b0;
            f->a1 = -2.0f * cosw0 * a0inv;
            f->a2 = (1.0f - alpha) * a0inv;
            break;
        }

        case BIQUAD_PEAK:  {
            A = powf(10.0f, gaindB/40.0f);
            
            a0inv = 1.0f/(1.0f + alpha/A);

            f->b0 = (1.0f + alpha * A) * a0inv;
            f->b1 = -2.0f * cosw0 * a0inv;
            f->b2 = (1.0f - alpha * A) * a0inv; 
            f->a1 = f->b1;
            f->a2 = (1.0f - alpha / A) * a0inv;
            break;
        }

        default: {
            f->b0 = 1.0f;
            f->b1 = 0.0f;
            f->b2 = 0.0f;
            f->a1 = 0.0f;
            f->a2 = 0.0f;
            break;
        }
    }
}

void biquad_process(Biquad *f, float *buffer, u16 n_samples) {

    for (u16 index = 0; index < n_samples; index++) {

        float w = buffer[index] - f->a1 * f->w1 - f->a2 * f->w2;
        buffer[index] = f->b0 * w + f->b1 * f->w1 + f->b2 * f->w2;
        
        f->w2 = f->w1;
        f->w1 = w;        
    }
}



// -------------------- Shelf Filter --------------------

enum ShelfType : u8 { lowshelf, highshelf };

struct ShelfFilter {

    float b0 = 1.0f;
    float b1 = 0.0f;
    float a1 = 0.0f;
    float x1 = 0.0f;
    float y1 = 0.0f;
    float out_gain = 1.0f;
};

void shelf_reset(ShelfFilter *f) {
    f->x1 = 0.0f;
    f->y1 = 0.0f;
}

void shelf_set_coeffs(ShelfFilter *f, float freq, float gain_db, float samplerate, ShelfType type) {
    float gain_linear = 1.0f;

    switch (type) {
        case lowshelf: {
            f->out_gain = (float)dbtoa(gain_db);
            gain_linear = (float)dbtoa(-gain_db);
            break;
        }
        
        case highshelf: {
            f->out_gain = 1.0f;
            gain_linear = (float)dbtoa(gain_db);
            break;
        }
        default: {
            f->b0 = 0.0f;
            f->b1 = 0.0f; 
            f->a1 = 0.0f;
            return;
        } 
    }

    float freqRadian = freq / samplerate * PI_F;

    local_const float PI_4 = HALFPI_F * 0.5f;

    float eta = (gain_linear + 1.0)/(gain_linear - 1.0);
    float rho = sinf(PI_F * freqRadian * 0.5 - PI_4) 
                / sinf(PI_F * freqRadian * 0.5 + PI_4);

    float etaSign = eta > 0.0f ? 1.0f : -1.0f;
    float alpha1 = gain_linear == 1.0f ? 0.0f : eta - etaSign*sqrtf(eta*eta - 1.0f);

    float beta0 = ((1.0f + gain_linear) + (1.0f - gain_linear) * alpha1) * 0.5f;
    float beta1 = ((1.0f - gain_linear) + (1.0f + gain_linear) * alpha1) * 0.5f;

    f->b0 = (beta0 + rho * beta1)/(1.0f + rho * alpha1);
    f->b1 = (beta1 + rho * beta0)/(1.0f + rho * alpha1);
    f->a1 = (rho + alpha1)/(1.0f + rho * alpha1);
}

void shelf_process(ShelfFilter *f, float *buffer, u16 n_samples) {
    
    for (u16 index = 0; index < n_samples; index++) {
        float out_sample = f->b0*buffer[index] + f->b1*f->x1 - f->a1*f->y1;
        f->x1 = buffer[index];
        f->y1 = out_sample;
        
        buffer[index] = out_sample;
    }
}


// ---------------- Main DSP data ----------------

struct DSPData {
    
    float samplerate = 0.0f;
    float upsamplerate = 0.0f;
    
    struct ParamStates {
        float gain1 = 0.0f;
        float gain2 = 0.0f; 
        float volume = 0.0f;
        
        float low = 0.5f;
        float mid = 0.5f;
        float trebble = 0.5f;
        
        bool bright = true;
        u8 channel = 0;
        bool do_gate = true;
        bool do_boost = false;

    } param_states;   

    
    struct NoiseGate {
        float threshold = 0.0;
        
        float *buffer = nullptr;
        u32 buffer_length = 0;
        u32 buffer_index = 0;
        
        float absolute_sum = 0.0;
        
        SmoothParam gain;
    
        bool is_open = false;
    
    } gate;
    
    OnePole tightFilter;
    Biquad boostFilter;
    
    struct Preamp {
        
        ShelfFilter cathode_bypass_filter0;
        OnePole input_filter;
        OnePole stage0LP;
        
        OnePole inputMudFilter;
        Biquad midBoost {BIQUAD_PEAK};
        ShelfFilter brightCapFilter;
        
        
        ShelfFilter cathode_bypass_filter1;
        OnePole couplingFilter1;
        OnePole stage1LP;


        ShelfFilter cathode_bypass_filter2;
        OnePole couplingFilter2;
        OnePole stage2LP;
        const float attenuation2 = 0.5f * tube_gain;

        ShelfFilter cathode_bypass_filter3;
        OnePole couplingFilter3;
        OnePole stage3LP;
        const float attenuation3 = 0.1f * tube_gain;

        ShelfFilter cathode_bypass_filter4;
        OnePole couplingFilter4;
        OnePole stage4LP;
    
        struct {
            Biquad upsample_filter1;
            Biquad upsample_filter2;
            
            Biquad downsample_filter1;
            Biquad downsample_filter2;
        } oversampler;
        
        float up_buffer[UP_BLOCK_SIZE] = {0};
    
        float stage0_bias[2] = {0};
        float stage1_bias[2] = {0};
        float stage2_bias[2] = {0};
        float stage3_bias[2] = {0};
        float stage4_bias[2] = {0};
    
        // float outputAttenuationdB = -34.0f;
    } preamp;
    
    ShelfFilter resonance;
    ShelfFilter presence;
    
    struct Tonestack {

        float b0 = 1.0;
        float b1 = 0.0;
        float b2 = 0.0;
        float b3 = 0.0;
        
        float a1 = 0.0;
        float a2 = 0.0;
        float a3 = 0.0;
        
        float x1 = 0.0f;
        float x2 = 0.0f;
        float x3 = 0.0f;
    
        float y1 = 0.0f;
        float y2 = 0.0f;
        float y3 = 0.0f;
        
    } tonestack;

    struct IRLoader {
        PFFFT_Setup* fft_setup = nullptr;
        // ShyFFT<float, fft_size> fft_setup;
        // arm_rfft_fast_instance_f32 fft_setup;
        
    } irloader;

    Biquad noise_filter1;
    Biquad noise_filter2;
} dsp;


global_const struct TonestackCtes {
    const double beta11 = 0.0001175;
    const double beta12 = 0.0005;
    const double beta13 = 0.02047;
    const double beta14 = 0.00051175;
    const double beta21 = 0.0000002209;
    const double beta22 = 0.000000255875;
    const double beta23 = 0.000000314625;
    const double beta24 = 0.0000032336;
    const double beta25 = 0.000010235;
    const double beta26 = 0.00000008084;
    const double beta31 = 0.0000000013959;
    const double beta32 = 0.0000000000348975;
    const double beta33 = 0.0000000000348975;
    const double beta34 = 0.000000000055225;
    const double beta35 = 0.000000000055225;
    const double beta36 = 0.000000002209;
    const double alpha11 = 0.00250925;
    const double alpha12 = 0.0005;
    const double alpha13 = 0.02047;
    const double alpha21 = -0.000000155375;
    const double alpha22 = 0.000010235;
    const double alpha23 = 0.000000255875;
    const double alpha24 = 0.0000220336;
    const double alpha25 = 0.00000077174;
    const double alpha31 = 0.0000000013959;
    const double alpha32 = 0.0000000000348975;
    const double alpha33 = -0.0000000000203275;
    const double alpha34 = 0.000000002209;
    const double alpha35 = 0.000000000055225;
} ctes;


void grid_conduction(float *buffer, u16 n_samples) {

    local_const float gridCondThresh = 1.0f;
    local_const float gridCondRatio  = 2.0f;
    // local_const float gridCondKnee   = 0.05f;
    local_const float gridCondKnee   = gridCondRatio / 4.0f;

    // @TODO: pour le bias, puis utiliser un Onepole pour smooth l'offset calculé sur ~10ms
    // refaire une fonction de calcul d'un seul echantillon pour les onepole

    for (u16 index = 0; index < n_samples; index++) {

        float sample = buffer[index];
    
        if (2.0f * (sample - gridCondThresh) > gridCondKnee) {
            sample = gridCondThresh + (sample - gridCondThresh)/gridCondRatio;
                    // sample = thresh + sample/ratio - thresh/ratio
                    // sample = thresh*ratio + sample - thresh
                    // sample = sample + thresh*ratio - thresh
        } else if (2.0f * abs(sample - gridCondThresh) <= gridCondKnee) {
            sample += ((1.0f/gridCondRatio - 1.0f) * powf(sample - gridCondThresh + gridCondKnee * 0.5f, 2))
                        /(2.0f * gridCondKnee);
        }
        
        buffer[index] = sample;
    }
}

void tube_sim(float *buffer, u16 n_samples, float pre_bias, float post_bias) {

    for (u16 index = 0; index < n_samples; index++) {

        float out_sample = buffer[index];
    
        out_sample *= -1.0f;
        out_sample += pre_bias;
    
        local_const float positiveLinRange = 0.2f;
        local_const float negClipPoint = 3.0f;
    
        if (out_sample > positiveLinRange) {
            // out_sample = tanh(out_sample - positiveLinRange) + positiveLinRange;
            out_sample = 1.0f - expf(positiveLinRange - out_sample) + positiveLinRange;
        }
        // else if (out_sample < -negClipPoint) {
        //     out_sample = -negClipPoint;
        // }
        else if (out_sample < 0.0f) {
            out_sample *= 2.0f/(3.0f*negClipPoint);
            out_sample = out_sample < -1.0f ? -2.0f/3.0f : out_sample - 1.0f/3.0f * out_sample*out_sample*out_sample;
            out_sample *= negClipPoint * 1.5f;
        }
        buffer[index] = out_sample - post_bias;
    }
}

float generate_output_bias(float input_bias) {
    local_const float positiveLinRange = 0.2f;
        
    if (input_bias > positiveLinRange) {
        return 1.0f - expf(positiveLinRange - input_bias) + positiveLinRange;
    } else {
        return input_bias;
    }
}


void update_tonestack() {
    float low = hardware.adc.GetFloat(Knob4);
    float mid = hardware.adc.GetFloat(Knob5);
    float trebble = hardware.adc.GetFloat(Knob6);

    if (low != dsp.param_states.low
        || mid != dsp.param_states.mid
        || trebble != dsp.param_states.trebble)    
    {
        // recompute tonestack
        
        dsp.param_states.low = low;
        dsp.param_states.mid = mid;
        dsp.param_states.trebble = trebble;
     
        low = scale_linear(low, 0.0f, 1.0f, 0.0f, 1.5f);
        mid = scale_linear(mid, 0.0f, 1.0f, 0.0f, 1.5f);
        trebble = scale_linear(trebble, 0.0f, 1.0f, 0.0f, 1.0f);
        
        double Low = exp((low-1.0)*3.4);
    
        double B1 = trebble*ctes.beta11 + mid*ctes.beta12 + Low*ctes.beta13 + ctes.beta14;
    
        double B2 = trebble*ctes.beta21 
                  - mid*mid*ctes.beta22
                  + mid*ctes.beta23
                  + Low*ctes.beta24
                  + Low*mid*ctes.beta25
                  + ctes.beta26;
    
        double B3 = Low*mid*ctes.beta31
                  - mid*mid*ctes.beta32
                  + mid*ctes.beta33
                  + trebble*ctes.beta34 - trebble*mid*ctes.beta35
                  + trebble*Low*ctes.beta36;
    
        double A0 = 1.0;
    
        double A1 = ctes.alpha11
                  + mid*ctes.alpha12 + Low*ctes.alpha13;
    
        double A2 = mid*ctes.alpha21
                  + Low*mid*ctes.alpha22
                  - mid*mid*ctes.alpha23
                  + Low*ctes.alpha24
                  + ctes.alpha25;
    
        double A3 = Low*mid*ctes.alpha31
                  - mid*mid*ctes.alpha32
                  + mid*ctes.alpha33
                  + Low*ctes.alpha34
                  + ctes.alpha35;
    
        double c = 2.0*dsp.samplerate;
    
    
        double a0 = -A0 - A1*c - A2 * pow(c, 2.0) - A3 * pow(c, 3.0);
    
        dsp.tonestack.b0 = (float)((-B1*c - B2 * pow(c, 2.0) - B3 * pow(c, 3.0))/a0);
        dsp.tonestack.b1 = (float)((-B1*c + B2 * pow(c, 2.0) + 3*B3 * pow(c, 3.0))/a0);
        dsp.tonestack.b2 = (float)((B1*c + B2 * pow(c, 2.0) - 3*B3 * pow(c, 3.0))/a0);
        dsp.tonestack.b3 = (float)((B1*c - B2 * pow(c, 2.0) + B3 * pow(c, 3.0))/a0);
        dsp.tonestack.a1 = (float)((-3*A0 - A1*c + A2 * pow(c, 2.0) + 3*A3 * pow(c, 3.0))/a0);
        dsp.tonestack.a2 = (float)((-3*A0 + A1*c + A2 * pow(c, 2.0) - 3*A3 * pow(c, 3.0))/a0);
        dsp.tonestack.a3 = (float)((-A0 + A1*c - A2 * pow(c, 2.0) + A3 * pow(c, 3.0))/a0);
    }
}




void AudioCallback(AudioHandle::InputBuffer  in,
                   AudioHandle::OutputBuffer out,
                   size_t nsamples)
{
    update_tonestack();

    float audio_buffer[BLOCK_SIZE] = {0};
    memcpy(audio_buffer, in[0], nsamples * sizeof(float));
    
    
    biquad_process(&dsp.noise_filter1, audio_buffer, nsamples);
    biquad_process(&dsp.noise_filter2, audio_buffer, nsamples);

    
    // process gate
    if (1) {
        DSPData::NoiseGate &gate = dsp.gate;
    
        local_const float attack_time_ms = 1.0f;
        local_const float release_time_ms = 15.0f;
        // local_const float hysteresis = 0.0f;
        // local_const float return_gain = 0.0f;
    
        for (u16 index = 0; index < nsamples; index++) {
            
            gate.absolute_sum -= abs(gate.buffer[gate.buffer_index]);
            
            gate.buffer[gate.buffer_index] = audio_buffer[index];
            gate.absolute_sum += abs(gate.buffer[gate.buffer_index]);
            
            gate.buffer_index++;
            if (gate.buffer_index == gate.buffer_length) { gate.buffer_index = 0; }

            // if close && > thresh -> open 
            // if close && < thresh -> close
            // if open && > thresh - hyst -> open 
            // if open && < thresh -hyst -> close
            
            bool should_open = false;
            float amplitude = gate.absolute_sum / gate.buffer_length;
            
            if (amplitude > gate.threshold) { 
                should_open = true;  
            }
            
            // if (gate.is_open && amplitude > return_gain) { 
            //     should_open = true;  
            // }
            
            // if (amplitude < gate.threshold) { 
            //     should_open = false; 
            // }
            
            if (should_open && !gate.is_open) {
                smooth_param_new_target(&gate.gain, 1.0f, attack_time_ms, dsp.samplerate);
                gate.is_open = true;
                
            } else if (!should_open && gate.is_open) {
                smooth_param_new_target(&gate.gain, 0.0f, release_time_ms, dsp.samplerate);
                gate.is_open = false;
            }
            
            float gain_value = smooth_param_next_value(&gate.gain); 
            
            audio_buffer[index] *= gain_value;            
        }
    }
    
    // tight boost
    
    if (dsp.param_states.do_boost) {
        dsp.tightFilter.ProcessBlock(audio_buffer, nsamples);
        biquad_process(&dsp.boostFilter, audio_buffer, nsamples);
    }
    
    // preamp
    {
        u16 up_nsamples = UPSAMPLE_FACTOR * nsamples;
        float *up_buffer = dsp.preamp.up_buffer;
        DSPData::Preamp &preamp = dsp.preamp;
    
        memset(up_buffer, 0, up_nsamples * sizeof(float));
        for (u16 index = 0; index < nsamples; index++) {
            up_buffer[UPSAMPLE_FACTOR * index] = audio_buffer[index]; 
        }
        
        biquad_process(&preamp.oversampler.upsample_filter1, up_buffer, up_nsamples);
        biquad_process(&preamp.oversampler.upsample_filter2, up_buffer, up_nsamples);
        apply_gain_linear(UPSAMPLE_FACTOR, up_buffer, up_nsamples);
    
        {
            // ------------ Stage 0 ------------
            grid_conduction(up_buffer, up_nsamples);
            shelf_process(&preamp.cathode_bypass_filter0, up_buffer, up_nsamples);
            tube_sim(up_buffer, up_nsamples, preamp.stage0_bias[0], preamp.stage0_bias[1]);
            
            preamp.input_filter.ProcessBlock(up_buffer, up_nsamples);
            preamp.stage0LP.ProcessBlock(up_buffer, up_nsamples);
            
            apply_gain_linear(dsp.param_states.gain1, up_buffer, up_nsamples);
            
            preamp.inputMudFilter.ProcessBlock(up_buffer, up_nsamples);
            biquad_process(&preamp.midBoost, up_buffer, up_nsamples);
            
            if (dsp.param_states.bright) {
                shelf_process(&preamp.brightCapFilter, up_buffer, up_nsamples);
            }
    
            // ------------ Stage 1 ------------
            grid_conduction(up_buffer, up_nsamples);
            shelf_process(&preamp.cathode_bypass_filter1, up_buffer, up_nsamples);
            tube_sim(up_buffer, up_nsamples, preamp.stage1_bias[0], preamp.stage1_bias[1]);
            
            preamp.couplingFilter1.ProcessBlock(up_buffer, up_nsamples);
            preamp.stage1LP.ProcessBlock(up_buffer, up_nsamples);
            
            apply_gain_linear(dsp.param_states.gain2, up_buffer, up_nsamples);


            // ------------ Stage 2 ------------
            grid_conduction(up_buffer, up_nsamples);
            shelf_process(&preamp.cathode_bypass_filter2, up_buffer, up_nsamples);
            tube_sim(up_buffer, up_nsamples, preamp.stage2_bias[0], preamp.stage2_bias[1]);
            
            preamp.couplingFilter2.ProcessBlock(up_buffer, up_nsamples);
            preamp.stage2LP.ProcessBlock(up_buffer, up_nsamples);
    
            if (dsp.param_states.channel == 0) {
                goto gain_stages_end_of_scope;
            }
    
            apply_gain_linear(preamp.attenuation2, up_buffer, up_nsamples);
            
            
            // ------------ Stage 3 ------------
            grid_conduction(up_buffer, up_nsamples);
            shelf_process(&preamp.cathode_bypass_filter3, up_buffer, up_nsamples);
            tube_sim(up_buffer, up_nsamples, preamp.stage3_bias[0], preamp.stage3_bias[1]);
            
            preamp.couplingFilter3.ProcessBlock(up_buffer, up_nsamples);
            preamp.stage3LP.ProcessBlock(up_buffer, up_nsamples);
    
    
            apply_gain_linear(preamp.attenuation3, up_buffer, up_nsamples);
    
    
            // ------------ Stage 4 ------------
            grid_conduction(up_buffer, up_nsamples);
            shelf_process(&preamp.cathode_bypass_filter4, up_buffer, up_nsamples);
            tube_sim(up_buffer, up_nsamples, preamp.stage4_bias[0], preamp.stage4_bias[1]);
            
            preamp.couplingFilter4.ProcessBlock(up_buffer, up_nsamples);
            preamp.stage4LP.ProcessBlock(up_buffer, up_nsamples);
        
            gain_stages_end_of_scope: {}
        }
        
        biquad_process(&preamp.oversampler.downsample_filter1, up_buffer, up_nsamples);
        biquad_process(&preamp.oversampler.downsample_filter2, up_buffer, up_nsamples);
        
        for (u16 index = 0; index < nsamples; index++) {
            audio_buffer[index] = up_buffer[index*UPSAMPLE_FACTOR];
        }
    }


    // tonestack
    {
        DSPData::Tonestack &ts = dsp.tonestack;
        
        for (u16 index = 0; index < nsamples; index++) {
            float out_sample = audio_buffer[index] * ts.b0
                                + ts.x1 * ts.b1
                                + ts.x2 * ts.b2
                                + ts.x3 * ts.b3
                                - ts.y1 * ts.a1
                                - ts.y2 * ts.a2
                                - ts.y3 * ts.a3;
            
            ts.x3 = ts.x2;
            ts.x2 = ts.x1;
            ts.x1 = audio_buffer[index];
            
            ts.y3 = ts.y2;
            ts.y2 = ts.y1;
            ts.y1 = out_sample;
            audio_buffer[index] = out_sample;
        }
    }
    
    // res / pres
    shelf_process(&dsp.resonance, audio_buffer, nsamples);
    shelf_process(&dsp.presence, audio_buffer, nsamples);
    
    
    local_const float inv_fft_size = 1.0f/fft_size;
    
    // irloader
    if (dsp.param_states.do_gate) {

        memcpy(ols_buffer, ols_buffer + nsamples, (fft_size - nsamples) * sizeof(float));
        memcpy(&ols_buffer[fft_size-nsamples], audio_buffer, nsamples * sizeof(float));
    
        float *temp = fdl_ptrs[npartitions - 1];
        for (u16 index = npartitions - 1; index > 0; index--) {
            fdl_ptrs[index] = fdl_ptrs[index-1];
        }
        fdl_ptrs[0] = temp;
    
        pffft_transform(dsp.irloader.fft_setup, ols_buffer, fdl_ptrs[0],  nullptr, PFFFT_FORWARD);
        
        memset(convolution_dft, 0, dft_size * sizeof(float));
    
        for (u16 part_index = 0; part_index < npartitions; part_index++) {
            pffft_zconvolve_accumulate(dsp.irloader.fft_setup, fdl_ptrs[part_index], ir_parts_dft[part_index], convolution_dft, 1.0f);
        }

        // local_const float band_gain = dbtoa(-36.0f);
        // convolution_dft[4] *= band_gain;
        // convolution_dft[4 + fft_size / 2] *= band_gain;

        pffft_transform(dsp.irloader.fft_setup, convolution_dft, signal_time_buffer, nullptr, PFFFT_BACKWARD);

        for (u32 index = 0; index < nsamples; index++) {
            audio_buffer[index] = signal_time_buffer[fft_size - nsamples + index] * inv_fft_size;
        }
        
    }

    for (u32 index = 0; index < nsamples; index++) {
        audio_buffer[index] *= dsp.param_states.volume;
    }

    if (monitor_index + BLOCK_SIZE > monitor_buffer_size) {
    
        u16 first_copy_size = monitor_buffer_size - monitor_index;
        u16 remainder = BLOCK_SIZE - first_copy_size;
    
        memcpy(&monitor_buffer[monitor_index], audio_buffer, first_copy_size * sizeof(float));
        monitor_index = 0;
        memcpy(&monitor_buffer[monitor_index], &audio_buffer[first_copy_size], remainder * sizeof(float));
    } else {
    
        memcpy(&monitor_buffer[monitor_index], audio_buffer, nsamples * sizeof(float));
        monitor_index += nsamples;
    }
    
    
    memcpy(out[0], audio_buffer, nsamples * sizeof(float));
    memset(out[1], 0, nsamples * sizeof(float));
}


int main(void) {
    hardware.Configure();
    hardware.Init();
    hardware.SetAudioBlockSize(BLOCK_SIZE);

    dsp.samplerate = hardware.AudioSampleRate();
    dsp.upsamplerate = UPSAMPLE_FACTOR * dsp.samplerate;

    //Create an ADC configuration
    AdcChannelConfig adc_configs[nKNOBS];
    
    adc_configs[Knob1].InitSingle(seed::D16);
    adc_configs[Knob2].InitSingle(seed::D17);
    adc_configs[Knob3].InitSingle(seed::D18);
    adc_configs[Knob4].InitSingle(seed::D19);
    adc_configs[Knob5].InitSingle(seed::D20);
    adc_configs[Knob6].InitSingle(seed::D21);

    led1.Init(seed::D22, false);
    led2.Init(seed::D23, false);
    
    led1.Set(0.0f);
    led1.Update();

    switches[Toggle4].Init(seed::D7);
    switches[Toggle3].Init(seed::D8);
    switches[Toggle2].Init(seed::D9);
    switches[Toggle1].Init(seed::D10);
    switches[FS1].Init(seed::D25);
    switches[FS2].Init(seed::D26);


    hardware.adc.Init(adc_configs, nKNOBS);
    //Start the adc
    hardware.adc.Start();    

    hardware.StartLog();

    {
        DSPData::NoiseGate &gate = dsp.gate;
    
        local_const float gate_buffer_length_sec = 0.01f;

        gate.buffer_length = (u32)(dsp.samplerate * gate_buffer_length_sec);
        gate.buffer_index = 0;
        gate.absolute_sum = 0.0;
        gate.threshold = dbtoa(-60.0);
        // gate.return_gain = gate.threshold;
        gate.is_open = false;

        gate.buffer = (float *)calloc(gate.buffer_length, sizeof(float));        
        smooth_param_init(&gate.gain, 0.0f);
    }


    biquad_reset(&dsp.noise_filter1);
    biquad_set_coeffs(&dsp.noise_filter1, 3000.0f, 
                        0.54119610f, 0.0f, dsp.samplerate, BIQUAD_LOWPASS);
    
    biquad_reset(&dsp.noise_filter2);
    biquad_set_coeffs(&dsp.noise_filter2, 3000.0f, 
                        1.3065630f, 0.0f, dsp.samplerate, BIQUAD_LOWPASS);

    dsp.tightFilter.Init();
    dsp.tightFilter.SetFrequency(400.0f/dsp.samplerate);
    dsp.tightFilter.SetFilterMode(OnePole::FilterMode::FILTER_MODE_HIGH_PASS);
    
    biquad_reset(&dsp.boostFilter);
    biquad_set_coeffs(&dsp.boostFilter, 1200.0f, 0.2, 8.0f, dsp.samplerate, BIQUAD_PEAK);
    
    {
        DSPData::Preamp &preamp = dsp.preamp;
        
        preamp.input_filter.Init();
        preamp.input_filter.SetFrequency(100.0f/dsp.upsamplerate); // 100.0
        preamp.input_filter.SetFilterMode(OnePole::FilterMode::FILTER_MODE_HIGH_PASS);
                
        preamp.inputMudFilter.Init();
        preamp.inputMudFilter.SetFrequency(100.0f/dsp.upsamplerate);
        preamp.inputMudFilter.SetFilterMode(OnePole::FilterMode::FILTER_MODE_HIGH_PASS);

        biquad_reset(&preamp.midBoost);
        biquad_set_coeffs(&preamp.midBoost, 1000.0f, 0.2f, 3.0f, dsp.upsamplerate, BIQUAD_PEAK);

        shelf_reset(&preamp.brightCapFilter);
        shelf_set_coeffs(&preamp.brightCapFilter, 550.0f, 0.0f, dsp.upsamplerate, lowshelf);


        preamp.couplingFilter1.Init();
        preamp.couplingFilter1.SetFrequency(15.0f / dsp.upsamplerate);
        preamp.couplingFilter1.SetFilterMode(OnePole::FilterMode::FILTER_MODE_HIGH_PASS);
        
        preamp.couplingFilter2.Init();
        preamp.couplingFilter2.SetFrequency(15.0f / dsp.upsamplerate);
        preamp.couplingFilter2.SetFilterMode(OnePole::FilterMode::FILTER_MODE_HIGH_PASS);
        
        preamp.couplingFilter3.Init();
        preamp.couplingFilter3.SetFrequency(15.0f / dsp.upsamplerate);
        preamp.couplingFilter3.SetFilterMode(OnePole::FilterMode::FILTER_MODE_HIGH_PASS);
        
        preamp.couplingFilter4.Init();
        preamp.couplingFilter4.SetFrequency(15.0f / dsp.upsamplerate);
        preamp.couplingFilter4.SetFilterMode(OnePole::FilterMode::FILTER_MODE_HIGH_PASS);
        
    
        preamp.stage0LP.Init();
        preamp.stage0LP.SetFrequency(12000.0f / dsp.upsamplerate);
        preamp.stage0LP.SetFilterMode(OnePole::FilterMode::FILTER_MODE_LOW_PASS);
        
        preamp.stage1LP.Init();
        preamp.stage1LP.SetFrequency(12000.0f / dsp.upsamplerate);
        preamp.stage1LP.SetFilterMode(OnePole::FilterMode::FILTER_MODE_LOW_PASS);
        
        preamp.stage2LP.Init();
        preamp.stage2LP.SetFrequency(12000.0f / dsp.upsamplerate);
        preamp.stage2LP.SetFilterMode(OnePole::FilterMode::FILTER_MODE_LOW_PASS);
        
        preamp.stage3LP.Init();
        preamp.stage3LP.SetFrequency(12000.0f / dsp.upsamplerate);
        preamp.stage3LP.SetFilterMode(OnePole::FilterMode::FILTER_MODE_LOW_PASS);
        
        preamp.stage4LP.Init();
        preamp.stage4LP.SetFrequency(12000.0f / dsp.upsamplerate);
        preamp.stage4LP.SetFilterMode(OnePole::FilterMode::FILTER_MODE_LOW_PASS);
        
    
        shelf_reset(&preamp.cathode_bypass_filter0);
        shelf_set_coeffs(&preamp.cathode_bypass_filter0, 280.0f, -4.0f, dsp.upsamplerate, lowshelf);
        
        shelf_reset(&preamp.cathode_bypass_filter1);
        shelf_set_coeffs(&preamp.cathode_bypass_filter1, 280.0f, -3.0f, dsp.upsamplerate, lowshelf);
        
        shelf_reset(&preamp.cathode_bypass_filter2);
        shelf_set_coeffs(&preamp.cathode_bypass_filter2, 280.0f, 0.0f, dsp.upsamplerate, lowshelf);
        
        shelf_reset(&preamp.cathode_bypass_filter3);
        shelf_set_coeffs(&preamp.cathode_bypass_filter3, 280.0f, -5.0f, dsp.upsamplerate, lowshelf);
        
        shelf_reset(&preamp.cathode_bypass_filter4);
        shelf_set_coeffs(&preamp.cathode_bypass_filter4, 280.0f, -3.0f, dsp.upsamplerate, lowshelf);
        
        
        biquad_reset(&preamp.oversampler.upsample_filter1);
        biquad_set_coeffs(&preamp.oversampler.upsample_filter1, dsp.samplerate/2.0f * 0.9f, 
                            0.54119610f, 0.0f, dsp.upsamplerate, BIQUAD_LOWPASS);
        
        biquad_reset(&preamp.oversampler.upsample_filter2);
        biquad_set_coeffs(&preamp.oversampler.upsample_filter2, dsp.samplerate/2.0f * 0.9f, 
                            1.3065630f, 0.0f, dsp.upsamplerate, BIQUAD_LOWPASS);
        
        biquad_reset(&preamp.oversampler.downsample_filter1);
        biquad_set_coeffs(&preamp.oversampler.downsample_filter1, dsp.samplerate/2.0f * 0.9f, 
                            0.54119610f, 0.0f, dsp.upsamplerate, BIQUAD_LOWPASS);
        
        biquad_reset(&preamp.oversampler.downsample_filter2);
        biquad_set_coeffs(&preamp.oversampler.downsample_filter2, dsp.samplerate/2.0f * 0.9f, 
                            1.3065630f, 0.0f, dsp.upsamplerate, BIQUAD_LOWPASS);                    
        
        
        preamp.stage0_bias[0] = 0.1f;
        preamp.stage0_bias[1] = generate_output_bias(preamp.stage0_bias[0]);

        preamp.stage1_bias[0] = 0.5f;
        preamp.stage1_bias[1] = generate_output_bias(preamp.stage1_bias[0]);

        preamp.stage2_bias[0] = 0.2f;
        preamp.stage2_bias[1] = generate_output_bias(preamp.stage2_bias[0]);

        preamp.stage3_bias[0] = 0.1f;
        preamp.stage3_bias[1] = generate_output_bias(preamp.stage3_bias[0]);

        preamp.stage4_bias[0] = 0.0f;
        preamp.stage4_bias[1] = generate_output_bias(preamp.stage4_bias[0]);
    }
    
    
    local_const float pres_amount = 0.6f;
    local_const float res_amount = 0.5f;
    
    shelf_reset(&dsp.resonance);
    shelf_set_coeffs(&dsp.resonance, 120.0f, scale_linear(res_amount, 0.0f, 1.0f, 0.0f, 24.0f), dsp.samplerate, lowshelf);
    
    shelf_reset(&dsp.presence);
    shelf_set_coeffs(&dsp.presence, 500.0f, scale_linear(pres_amount, 0.0f, 1.0f, 0.0f, 18.0f), dsp.samplerate, highshelf);

    {
        // dsp.irloader.fft_setup.Init();
        dsp.irloader.fft_setup = pffft_new_setup(fft_size, PFFFT_REAL);
        
        for (u16 part_index = 0; part_index < npartitions; part_index++) {
            fdl_ptrs[part_index] = fdl[part_index];
        }
        
        for (u16 part_index = 0; part_index < npartitions; part_index++) {
            memset(signal_time_buffer, 0, fft_size * sizeof(float));
            memcpy(signal_time_buffer, &ir[part_index * BLOCK_SIZE], BLOCK_SIZE * sizeof(float));
            
            // dsp.irloader.fft_setup.Direct(signal_time_buffer, ir_parts_dft[part_index]);            
            pffft_transform(dsp.irloader.fft_setup, signal_time_buffer, ir_parts_dft[part_index], nullptr, PFFFT_FORWARD);
        }
    }


    //Start calling the audio callback
    hardware.StartAudio(AudioCallback);

    // Loop forever
    while (1) {
        float gain1 = hardware.adc.GetFloat(Knob1);
        float gain2 = hardware.adc.GetFloat(Knob2);
        float volume = hardware.adc.GetFloat(Knob3);
        
        switches[Toggle4].Debounce();
        switches[Toggle3].Debounce();
        switches[Toggle2].Debounce();
        switches[Toggle1].Debounce();
        // switches[FS1].Debounce();
        // switches[FS2].Debounce();
            
        // 0 = 3 stages, 1 = 5 stages
        dsp.param_states.channel = (u8)switches[Toggle4].Pressed();
        dsp.param_states.do_gate = switches[Toggle2].Pressed();
        
        dsp.param_states.gain1 = tube_gain * scale(gain1, 0.0f, 1.0f, 0.0f, 1.0f, 2.5f);

        shelf_set_coeffs(&dsp.preamp.brightCapFilter, 550.0f, 
                        scale_linear(gain1, 0.0f, 1.0f, -15.0f, 0.0f),
                        dsp.upsamplerate, lowshelf);
        
        dsp.param_states.gain2 = tube_gain * scale(gain2, 0.0f, 1.0f, 0.0f, 1.0f, 2.5f);
        
        dsp.param_states.volume = scale(volume, 0.0f, 1.0f, 0.0f, 1.0f, 3.0f);
        
        dsp.param_states.bright = switches[Toggle3].Pressed();
        
        dsp.param_states.do_boost = switches[Toggle1].Pressed();
        // volatile bool switch1_state = switches[FS1].Pressed();
        // volatile bool switch2_state = switches[FS2].Pressed();
        
    
        dsp.tightFilter.SetFrequency(400.0f/dsp.samplerate);
        
        biquad_set_coeffs(&dsp.boostFilter, 1200.0f, 0.2, 8.0f, dsp.samplerate, BIQUAD_PEAK);
            
        shelf_set_coeffs(&dsp.resonance, 120.0f, scale_linear(res_amount, 0.0f, 1.0f, 0.0f, 24.0f), dsp.samplerate, lowshelf);        
        shelf_set_coeffs(&dsp.presence, 500.0f, scale_linear(pres_amount, 0.0f, 1.0f, 0.0f, 18.0f), dsp.samplerate, highshelf);
    
        {
            float sum = 0.0f;
            for (u16 index = 0; index < monitor_buffer_size; index++) {
                float sample = monitor_buffer[index];
                sum += sample * sample;
            }
            
            sum = sqrtf(sum / (float)monitor_buffer_size);
            
            float db = atodb(sum);
            db = clip(db, -60.0f, 0.0f);
            db = scale_linear(db, -60.0f, 0.0f, 0.0f, 1.0f);
            
            led1.Set(db);
            led1.Update();
        }
        
        // System::Delay(30);
    }
}
