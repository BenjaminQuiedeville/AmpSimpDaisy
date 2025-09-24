#include <stdint.h>

#include "daisy_seed.h"

// #include "daisysp.h"

#define PI_F 3.1415927410125732421875f
#define TWOPI_F (2.0f * PI_F)
#define HALFPI_F (PI_F * 0.5f)

#define PFFFT_SIMD_DISABLE
#include "pffft/pffft.c"

// #include "../../stmlib/fft/sht_fft.h"

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
// using namespace daisysp;

global_const u16 BLOCK_SIZE = 16;
global_const u16 UPSAMPLE_FACTOR = 4;
global_const u16 UP_BLOCK_SIZE = BLOCK_SIZE * UPSAMPLE_FACTOR;

global_const float tube_gain = 100.0f;
global_const float smooth_tau_ms = 100.0f/BLOCK_SIZE;

// float ir[ir_size]
#include "data/base_ir.h"

global_const u16 npartitions = ir_size / BLOCK_SIZE;
global_const u16 fft_size = 32;
global_const u16 dft_size = fft_size + 2;

global u16 monitor_index = 0;
global_const u16 monitor_buffer_size = BLOCK_SIZE * 60;
global float monitor_buffer[monitor_buffer_size] = {0};

global_const u32 gate_buffer_length = 96;
global float gate_buffer[gate_buffer_length] = {0};

global_const float gate_thresh[2][2] = {
    {-85.0f, -80.0f},
    {-70.0f, -60.0f}
};

global_const float channel_biases[4][5] = {
    {0.1f,  0.8f,  0.0f,  0.0f,  0.0f},
    {0.1f,  0.2f,  0.1f,  0.0f,  0.0f},
    {0.1f,  0.2f,  0.8f,  0.2f,  0.1f},
    {0.1f,  0.1f,  0.1f,  0.2f,  0.1f},
};

global_const float cathode_bp_amount[4][5] {
    {-2.0f,  0.0f, -1.0f,  0.0f,  0.0f},
    {-2.0f, -1.0f, -1.0f,  0.0f,  0.0f},
    {-3.0f, -2.0f,  0.0f, -3.0f, -1.0f},
    {-2.0f, -3.0f, -1.0f, -2.0f, -3.0f},
};

global_const float attenuations[4][4] = {
    {0.6f,  0.3f,  0.0f,  0.0f},
    {1.0f,  0.8f,  0.0f,  0.0f},
    {0.7f,  0.8f,  0.1f,  0.4f},
    {0.7f,  0.8f,  0.5f,  0.4f},
};

enum Channels {
    CLEAN, 
    OVERDRIVE,
};

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

static inline bool is_far(float y, float x, float tol) {
    return ((y + tol) <= x) || ((y - tol) >= x);
}

// -------------------- Smooth parameters --------------------

struct SmoothParam {
    float current_value = 0.0f;
    float target = 0.0f;
    float b0 = 0.0f;
    float a1 = 0.0f;
    float y1 = 0.0f;
    float value_buffer[BLOCK_SIZE] = {0};
};

static void smooth_param_init(SmoothParam *param, float init_value) {
    param->current_value = init_value;
    param->target = init_value;
    param->y1 = init_value;
    param->b0 = 1.0f;
    param->a1 = 0.0f;
}

static void smooth_param_new_target(SmoothParam *param, float new_target, float tau_ms, float samplerate) {
    param->b0 = sinf(PI_F / (samplerate * tau_ms * 0.001f));
    param->a1 = param->b0 - 1.0f;
    param->target = new_target;
}

static inline float smooth_param_next_value(SmoothParam *param) {
    param->current_value = param->target * param->b0 - param->y1 * param->a1;
    param->y1 = param->current_value;

    return param->current_value;
}

// static void smooth_param_compute_next_block(SmoothParam *param, u16 nvalues) {
//     for (u16 index = 0; index < nvalues; index++) {
//         param->value_buffer[index] = smooth_param_next_value(param);
//     }
// }

// -------------------- Onepole Filter --------------------

struct Onepole {
    float b0 = 0.0f;
    float a1 = 0.0f;
    float y1 = 0.0f;
};

// static inline void onepole_reset(Onepole *f) {
//     f->y1 = 0.0f;    
// }

static void onepole_set_coeffs(Onepole *f, float frequency, float samplerate) {
    f->b0 = sinf(PI_F / samplerate * frequency);
    f->a1 = 1.0f - f->b0;
}

static void onepole_process_lowpass(Onepole *f, float *buffer, u32 nsamples) {
    
    float b0 = f->b0;
    float a1 = f->a1;
    float y1 = f->y1;
    
    for (u32 index = 0; index < nsamples; index++) {
        float sample = buffer[index];
        float out_sample = sample*b0 + a1 *y1;
        y1 = out_sample;
        buffer[index] = out_sample;
    }
    
    f->y1 = y1;
}

static void onepole_process_highpass(Onepole *f, float *buffer, u32 nsamples) {
    
    float b0 = f->b0;
    float a1 = f->a1;
    float y1 = f->y1;
    
    for (u32 index = 0; index < nsamples; index++) {
        float sample = buffer[index];
        float lp_sample = sample*b0 + a1 *y1;
        y1 = lp_sample;
        buffer[index] = sample - lp_sample;
    }
    
    f->y1 = y1;
}

// -------------------- Biquad Filter --------------------

enum BiquadType : u8 {
    BIQUAD_LOWPASS = 0,
    BIQUAD_PEAK,
    BIQUAD_NFILTERTYPES,
};

struct Biquad {

    float b0 = 0.0;
    float b1 = 0.0;
    float b2 = 0.0;
    float a1 = 0.0;
    float a2 = 0.0;

    float w1 = 0.0f;
    float w2 = 0.0f;
};


// static void biquad_reset(Biquad *f) {
//     f->w1 = 0.0f;
//     f->w2 = 0.0f;
// }

static void biquad_set_coeffs(Biquad *f, float frequency, float Q, float gaindB, float samplerate, BiquadType type) {

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
            f->a1 = 2.0f * cosw0 * a0inv;
            f->a2 = -(1.0f - alpha) * a0inv;
            break;
        }

        case BIQUAD_PEAK:  {
            A = powf(10.0f, gaindB/40.0f);

            a0inv = 1.0f/(1.0f + alpha/A);

            f->b0 = (1.0f + alpha * A) * a0inv;
            f->b1 = -2.0f * cosw0 * a0inv;
            f->b2 = (1.0f - alpha * A) * a0inv;
            f->a1 = -f->b1;
            f->a2 = -(1.0f - alpha / A) * a0inv;
            break;
        }

        default: {
            f->b0 = 0.0f;
            f->b1 = 0.0f;
            f->b2 = 0.0f;
            f->a1 = 0.0f;
            f->a2 = 0.0f;
            break;
        }
    }
}

static void biquad_process(Biquad *f, float *buffer, u16 n_samples) {
    float b0 = f->b0;
    float b1 = f->b1;
    float b2 = f->b2;
    float a1 = f->a1;
    float a2 = f->a2;
    
    float w1 = f->w1;
    float w2 = f->w2;

    for (u16 index = 0; index < n_samples; index++) {
        
        float w = buffer[index] + a1 * w1 + a2 * w2;
        buffer[index] = b0 * w + b1 * w1 + b2 * w2;

        w2 = w1;
        w1 = w;
    }
    
    f->w1 = w1;
    f->w2 = w2;    
}



// -------------------- Shelf Filter --------------------

enum ShelfType : u8 { lowshelf, highshelf };

struct ShelfFilter {

    float b0 = 0.0f;
    float b1 = 0.0f;
    float a1 = 0.0f;
    float x1 = 0.0f;
    float y1 = 0.0f;
    float out_gain = 1.0f;
};

// static void shelf_reset(ShelfFilter *f) {
//     f->x1 = 0.0f;
//     f->y1 = 0.0f;
// }

static void shelf_set_coeffs(ShelfFilter *f, float freq, float gain_db, float samplerate, ShelfType type) {
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
    f->a1 = -(rho + alpha1)/(1.0f + rho * alpha1);
}

static void shelf_process(ShelfFilter *f, float *buffer, u16 n_samples) {

    float b0 = f->b0;
    float b1 = f->b1;
    float a1 = f->a1;
    float x1 = f->x1;
    float y1 = f->y1;

    for (u16 index = 0; index < n_samples; index++) {
        float input_sample = buffer[index];
        
        float out_sample = b0*input_sample + b1*x1 + a1*y1;
        x1 = input_sample;
        y1 = out_sample;

        buffer[index] = out_sample * f->out_gain; // * gain linear ? c'est pour ca que j'ai l'impression que ca fait rien ?
    }
    
    f->x1 = x1;
    f->y1 = y1;
}


// ---------------- Main DSP data ----------------

struct DSPData {

    float samplerate = 0.0f;
    float upsamplerate = 0.0f;
    float output_attenuation = dbtoa(-9.0f);

    struct {
        float gain1 = 0.0f;
        float gain2 = 0.0f;
        float volume = 0.0f;

        SmoothParam low;
        SmoothParam mid;
        SmoothParam trebble;

        u8 channel = CLEAN;
        u8 preamp_variant = 0;
    } params;


    struct NoiseGate {

        float threshold = 0.0;
        u32 buffer_index = 0;
        float absolute_sum = 0.0;
        SmoothParam gain;
        bool is_open = false;

    } gate;

    Onepole tight_filter;
    Biquad boost_filter;

    struct Preamp {

        ShelfFilter cathode_bypass_filter0;
        Onepole input_filter;
        Onepole stage0LP;

        Onepole input_mud_filter;
        Biquad mid_boost {BIQUAD_PEAK};
        ShelfFilter bright_cap_filter;

        ShelfFilter cathode_bypass_filter1;
        Onepole coupling_filter1;
        Onepole stage1LP;

        ShelfFilter cathode_bypass_filter2;
        Onepole coupling_filter2;
        Onepole stage2LP;

        ShelfFilter cathode_bypass_filter3;
        Onepole coupling_filter3;
        Onepole stage3LP;

        ShelfFilter cathode_bypass_filter4;
        Onepole coupling_filter4;
        Onepole stage4LP;

        struct {
            Biquad upsample_filter1;
            Biquad upsample_filter2;

            Biquad downsample_filter1;
            Biquad downsample_filter2;
        } oversampler;
        
        float attenuations[4] = {0};
        
        float stage0_bias[2] = {0};
        float stage1_bias[2] = {0};
        float stage2_bias[2] = {0};
        float stage3_bias[2] = {0};
        float stage4_bias[2] = {0};
    } preamp;

    ShelfFilter resonance;
    ShelfFilter presence;

    struct Tonestack {

        float b0 = 0.0f;
        float b1 = 0.0f;
        float b2 = 0.0f;
        float b3 = 0.0f;

        float a1 = 0.0f;
        float a2 = 0.0f;
        float a3 = 0.0f;

        float x1 = 0.0f;
        float x2 = 0.0f;
        float x3 = 0.0f;

        float y1 = 0.0f;
        float y2 = 0.0f;
        float y3 = 0.0f;

    } tonestack;

    struct IRLoader {
        PFFFT_Setup* fft_setup = nullptr;
        alignas(16) float signal_time_buffer[fft_size] = {0};
        alignas(16) float ols_buffer[fft_size] = {0};
        alignas(16) float *fdl_ptrs[npartitions] = {0};
        alignas(16) float fdl[npartitions][dft_size] = {0};
        alignas(16) float ir_parts_dft[npartitions][dft_size] = {0};
        alignas(16) float convolution_dft[dft_size] = {0};

    } irloader;

} dsp;


static void grid_conduction(float *buffer, u16 n_samples) {

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

static void tube_sim(float *buffer, u16 n_samples, float pre_bias, float post_bias) {

    for (u16 index = 0; index < n_samples; index++) {

        float sample = buffer[index];

        sample *= -1.0f;
        sample += pre_bias;

        constexpr float positiveLinRange = 0.2f;
        constexpr float negClipPoint = 3.0f;

        if (sample > positiveLinRange) {
            float temp = (sample-positiveLinRange) * 0.5f + 1.0f;
            sample = 1.0f - 1.0f/(temp*temp) + positiveLinRange;
        }
        // else if (sample < -negClipPoint) {
        //     sample = -negClipPoint;
        // }
        else if (sample < 0.0f) {
            sample *= 2.0f/(3.0f*negClipPoint);
            sample = sample < -1.0f ? -2.0f/3.0f : sample - 1.0f/3.0f * sample*sample*sample;
            sample *= negClipPoint * 1.5f;
        }
        buffer[index] = sample - post_bias;
    }
}

static float generate_output_bias(float input_bias) {
    local_const float positiveLinRange = 0.2f;

    if (input_bias > positiveLinRange) {
        return 1.0f - expf(positiveLinRange - input_bias) + positiveLinRange;
    } else {
        return input_bias;
    }
}

static void set_biases(DSPData::Preamp *preamp, u32 index) {
        
    preamp->stage0_bias[0] = channel_biases[index][0];
    preamp->stage0_bias[1] = generate_output_bias(preamp->stage0_bias[0]);

    preamp->stage1_bias[0] = channel_biases[index][1];
    preamp->stage1_bias[1] = generate_output_bias(preamp->stage1_bias[0]);

    preamp->stage2_bias[0] = channel_biases[index][2];
    preamp->stage2_bias[1] = generate_output_bias(preamp->stage2_bias[0]);

    preamp->stage3_bias[0] = channel_biases[index][3];
    preamp->stage3_bias[1] = generate_output_bias(preamp->stage3_bias[0]);

    preamp->stage4_bias[0] = channel_biases[index][4];
    preamp->stage4_bias[1] = generate_output_bias(preamp->stage4_bias[0]);
}


static void set_cathode_bypasses(DSPData::Preamp *preamp, u32 index) {
//     shelf_reset(&preamp->cathode_bypass_filter0);
    shelf_set_coeffs(&preamp->cathode_bypass_filter0, 280.0f, cathode_bp_amount[index][0], dsp.upsamplerate, lowshelf);

//     shelf_reset(&preamp->cathode_bypass_filter1);
    shelf_set_coeffs(&preamp->cathode_bypass_filter1, 280.0f, cathode_bp_amount[index][1], dsp.upsamplerate, lowshelf);

//     shelf_reset(&preamp->cathode_bypass_filter2);
    shelf_set_coeffs(&preamp->cathode_bypass_filter2, 280.0f, cathode_bp_amount[index][2], dsp.upsamplerate, lowshelf);

//     shelf_reset(&preamp->cathode_bypass_filter3);
    shelf_set_coeffs(&preamp->cathode_bypass_filter3, 280.0f, cathode_bp_amount[index][3], dsp.upsamplerate, lowshelf);

//     shelf_reset(&preamp->cathode_bypass_filter4);
    shelf_set_coeffs(&preamp->cathode_bypass_filter4, 280.0f, cathode_bp_amount[index][4], dsp.upsamplerate, lowshelf);
}

static void set_attenuations(DSPData::Preamp *preamp, u32 index) {
    for (u32 i = 0; i < 4; i++) {
        preamp->attenuations[i] = attenuations[index][i];
    }
}


void AudioCallback(AudioHandle::InputBuffer  in, AudioHandle::OutputBuffer out, size_t nsamples) {

    float gain1 = hardware.adc.GetFloat(Knob1);
    float gain2 = hardware.adc.GetFloat(Knob2);
    float volume = hardware.adc.GetFloat(Knob3);

    switches[Toggle4].Debounce();
    switches[Toggle3].Debounce();
    switches[Toggle2].Debounce();
    switches[Toggle1].Debounce();
    // switches[FS2].Debounce();

    bool do_boost = switches[Toggle1].Pressed();
    bool do_gate = switches[Toggle2].Pressed();
    
    u8 variant = (u8)switches[Toggle3].Pressed();
    
    bool update_biases_and_cathodes = false;
    
    switches[FS1].Debounce();
    if (switches[FS1].RisingEdge()) {
        dsp.params.channel = !dsp.params.channel;
        update_biases_and_cathodes = true;
        
        //@TODO update everything here
    }

    if (variant != dsp.params.preamp_variant) {
        dsp.params.preamp_variant = variant;
        update_biases_and_cathodes = true;
    }


    // 0 = 3 stages, 1 = 5 stages
    if (update_biases_and_cathodes) {
    
        // change biases        
        u16 index = dsp.params.channel*2 + variant;
        
        set_biases(&dsp.preamp, index);
        set_cathode_bypasses(&dsp.preamp, index);
        set_attenuations(&dsp.preamp, index);
    }

    dsp.params.gain1 = tube_gain * scale(gain1, 0.0f, 1.0f, 0.0f, 1.0f, 2.0f);

    shelf_set_coeffs(&dsp.preamp.bright_cap_filter, 550.0f,
                    scale_linear(gain1, 0.0f, 1.0f, -10.0f, 0.0f),
                    dsp.upsamplerate, lowshelf);

    dsp.params.gain2 = tube_gain * scale(gain2, 0.0f, 1.0f, 0.0f, 1.0f, 2.0f);

    dsp.params.volume = dsp.output_attenuation * scale(volume, 0.0f, 1.0f, 0.0f, 1.0f, 2.0f);

    // switches[Toggle3].Pressed();

    float knob4_value = hardware.adc.GetFloat(Knob4);
    float knob5_value = hardware.adc.GetFloat(Knob5);
    float knob6_value = hardware.adc.GetFloat(Knob6);

    // float low = hardware.adc.GetFloat(Knob4);
    // float mid = hardware.adc.GetFloat(Knob5);
    // float trebble = hardware.adc.GetFloat(Knob6);

    smooth_param_new_target(&dsp.params.low, knob4_value, smooth_tau_ms, dsp.samplerate);
    smooth_param_new_target(&dsp.params.mid, knob5_value, smooth_tau_ms, dsp.samplerate);
    smooth_param_new_target(&dsp.params.trebble, knob6_value, smooth_tau_ms, dsp.samplerate);

    {

        enum {
            beta11, beta12, beta13, beta14, 
            beta21, beta22, beta23, beta24, beta25, beta26, 
            beta31, beta32, beta33, beta34, beta35, beta36, 
            alpha11, alpha12, alpha13, 
            alpha21, alpha22, alpha23, alpha24, alpha25, 
            alpha31, alpha32, alpha33, alpha34, alpha35
        };

        // Soldano SLO
        local_const double ctes[29] = {
            0.0001175, 0.0005, 0.02047, 0.00051175,
            0.0000002209, 0.000000255875, 0.000000314625, 0.0000032336, 0.000010235, 0.00000008084,
            0.0000000013959, 0.0000000000348975, 0.0000000000348975, 0.000000000055225, 0.000000000055225, 0.000000002209,
            0.00250925, 0.0005, 0.02047,
            -0.000000155375, 0.000010235, 0.000000255875, 0.0000220336, 0.00000077174,
            0.0000000013959, 0.0000000000348975, -0.0000000000203275, 0.000000002209, 0.000000000055225,
        }; 

        // recompute tonestack
        float low     = smooth_param_next_value(&dsp.params.low);
        float mid     = smooth_param_next_value(&dsp.params.mid);
        float trebble = smooth_param_next_value(&dsp.params.trebble);

        double Low = scale_linear(low, 0.0f, 1.0f, 0.0f, 1.5f);
        double Mid = scale_linear(mid, 0.0f, 1.0f, 0.0f, 1.5f);
        double Trebble = scale_linear(trebble, 0.0f, 1.0f, 0.0f, 1.0f);

        Low = exp((Low-1.0)*3.4);

        double B1 = Trebble*ctes[beta11] + Mid*ctes[beta12] + Low*ctes[beta13] + ctes[beta14];

        double B2 = Trebble*ctes[beta21]
                  - Mid*Mid*ctes[beta22]
                  + Mid*ctes[beta23]
                  + Low*ctes[beta24]
                  + Low*Mid*ctes[beta25]
                  + ctes[beta26];

        double B3 = Low*Mid*ctes[beta31]
                  - Mid*Mid*ctes[beta32]
                  + Mid*ctes[beta33]
                  + Trebble*ctes[beta34] - Trebble*Mid*ctes[beta35]
                  + Trebble*Low*ctes[beta36];

        double A0 = 1.0;

        double A1 = ctes[alpha11]
                  + Mid*ctes[alpha12] + Low*ctes[alpha13];

        double A2 = Mid*ctes[alpha21]
                  + Low*Mid*ctes[alpha22]
                  - Mid*Mid*ctes[alpha23]
                  + Low*ctes[alpha24]
                  + ctes[alpha25];

        double A3 = Low*Mid*ctes[alpha31]
                  - Mid*Mid*ctes[alpha32]
                  + Mid*ctes[alpha33]
                  + Low*ctes[alpha34]
                  + ctes[alpha35];

        double c = 2.0*dsp.samplerate;


        double a0 = -A0 - A1*c - A2 * pow(c, 2.0) - A3 * pow(c, 3.0);

        dsp.tonestack.b0 = (float)((-B1*c - B2 * pow(c, 2.0) - B3 * pow(c, 3.0))/a0);
        dsp.tonestack.b1 = (float)((-B1*c + B2 * pow(c, 2.0) + 3*B3 * pow(c, 3.0))/a0);
        dsp.tonestack.b2 = (float)((B1*c + B2 * pow(c, 2.0) - 3*B3 * pow(c, 3.0))/a0);
        dsp.tonestack.b3 = (float)((B1*c - B2 * pow(c, 2.0) + B3 * pow(c, 3.0))/a0);
        dsp.tonestack.a1 = -(float)((-3*A0 - A1*c + A2 * pow(c, 2.0) + 3*A3 * pow(c, 3.0))/a0);
        dsp.tonestack.a2 = -(float)((-3*A0 + A1*c + A2 * pow(c, 2.0) - 3*A3 * pow(c, 3.0))/a0);
        dsp.tonestack.a3 = -(float)((-A0 + A1*c - A2 * pow(c, 2.0) + A3 * pow(c, 3.0))/a0);
    }

    // Boost section
    // @TODO do not recompute everytime 
    if (dsp.params.channel == CLEAN) {
        onepole_set_coeffs(&dsp.tight_filter, 300.0f, dsp.samplerate);
        biquad_set_coeffs(&dsp.boost_filter, 800.0f, 0.2f, 6.0f, dsp.samplerate, BIQUAD_PEAK);
    } else {
        onepole_set_coeffs(&dsp.tight_filter, 400.0f, dsp.samplerate);
        biquad_set_coeffs(&dsp.boost_filter, 1200.0f, 0.2f, 0.0f, dsp.samplerate, BIQUAD_PEAK);
    }

    // main process
    float audio_buffer[BLOCK_SIZE] = {0};
    memcpy(audio_buffer, in[0], nsamples * sizeof(float));


    // process gate
    if (1) {
        DSPData::NoiseGate &gate = dsp.gate;

        gate.threshold = dbtoa(gate_thresh[dsp.params.channel][do_gate]);

        local_const float attack_time_ms = 0.5f;
        local_const float release_time_ms = 15.0f;
        // local_const float hysteresis = 0.0f;
        // local_const float return_gain = 0.0f;
        
        float gate_absolute_sum = gate.absolute_sum;
        u32 gate_buffer_index = gate.buffer_index; 
        const float gate_threshold = gate.threshold;
                        
        for (u16 index = 0; index < nsamples; index++) {

            gate_absolute_sum -= abs(gate_buffer[gate_buffer_index]);

            gate_buffer[gate_buffer_index] = audio_buffer[index];
            gate_absolute_sum += abs(gate_buffer[gate_buffer_index]);

            gate_buffer_index++;
            if (gate_buffer_index == gate_buffer_length) { gate_buffer_index = 0; }

            // if close && > thresh -> open
            // if close && < thresh -> close
            // if open && > thresh - hyst -> open
            // if open && < thresh -hyst -> close

            bool should_open = false;
            float amplitude = gate_absolute_sum / (float)gate_buffer_length;

            if (amplitude > gate_threshold) {
                should_open = true;
            }

            // if (gate.is_open && amplitude > return_gain) {
            //     should_open = true;
            // }

            // if (amplitude < gate_threshold) {
            //     should_open = false;
            // }

            if (should_open && !gate.is_open) {
                smooth_param_new_target(&gate.gain, 1.0f, attack_time_ms, dsp.samplerate);
                gate.is_open = true;

            } else if (!should_open && gate.is_open) {
                smooth_param_new_target(&gate.gain, 0.0f, release_time_ms, dsp.samplerate);
                gate.is_open = false;
            }

            gate.gain.value_buffer[index] = smooth_param_next_value(&gate.gain);
        }
        
        for (u16 index = 0; index < nsamples; index++) {
            audio_buffer[index] *= gate.gain.value_buffer[index];
        }

        gate.absolute_sum = gate_absolute_sum;
        gate.buffer_index = gate_buffer_index; 
    }

    // tight boost

    if (do_boost) {
        onepole_process_highpass(&dsp.tight_filter, audio_buffer, nsamples);
        // dsp.tight_filter.ProcessBlock(audio_buffer, nsamples);
        biquad_process(&dsp.boost_filter, audio_buffer, nsamples);
    }

    // preamp
    {
        DSPData::Preamp &preamp = dsp.preamp;
        u16 up_nsamples = UPSAMPLE_FACTOR * nsamples;
        float up_buffer[UP_BLOCK_SIZE] = {0};

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

            onepole_process_highpass(&preamp.input_filter, up_buffer, up_nsamples);
            onepole_process_lowpass(&preamp.stage0LP, up_buffer, up_nsamples);

            apply_gain_linear(dsp.params.gain1 * dsp.preamp.attenuations[0], up_buffer, up_nsamples);

            onepole_process_highpass(&preamp.input_mud_filter, up_buffer, up_nsamples);
            biquad_process(&preamp.mid_boost, up_buffer, up_nsamples);


            shelf_process(&preamp.bright_cap_filter, up_buffer, up_nsamples);

            // ------------ Stage 1 ------------
            grid_conduction(up_buffer, up_nsamples);
            shelf_process(&preamp.cathode_bypass_filter1, up_buffer, up_nsamples);
            tube_sim(up_buffer, up_nsamples, preamp.stage1_bias[0], preamp.stage1_bias[1]);

            onepole_process_highpass(&preamp.coupling_filter1, up_buffer, up_nsamples);
            onepole_process_lowpass(&preamp.stage1LP, up_buffer, up_nsamples);

            apply_gain_linear(dsp.params.gain2 * dsp.preamp.attenuations[1], up_buffer, up_nsamples);


            // ------------ Stage 2 ------------
            grid_conduction(up_buffer, up_nsamples);
            shelf_process(&preamp.cathode_bypass_filter2, up_buffer, up_nsamples);
            tube_sim(up_buffer, up_nsamples, preamp.stage2_bias[0], preamp.stage2_bias[1]);

            onepole_process_highpass(&preamp.coupling_filter2, up_buffer, up_nsamples);
            onepole_process_lowpass(&preamp.stage2LP, up_buffer, up_nsamples);

            if (dsp.params.channel == CLEAN) {

                local_const float clean_compensation_gain = 0.4f;
                apply_gain_linear(clean_compensation_gain, up_buffer, up_nsamples);
                goto gain_stages_end_of_scope;
            }

            apply_gain_linear(dsp.preamp.attenuations[2] * tube_gain, up_buffer, up_nsamples);


            // ------------ Stage 3 ------------
            grid_conduction(up_buffer, up_nsamples);
            shelf_process(&preamp.cathode_bypass_filter3, up_buffer, up_nsamples);
            tube_sim(up_buffer, up_nsamples, preamp.stage3_bias[0], preamp.stage3_bias[1]);

            onepole_process_highpass(&preamp.coupling_filter3, up_buffer, up_nsamples);
            onepole_process_lowpass(&preamp.stage3LP, up_buffer, up_nsamples);

            apply_gain_linear(dsp.preamp.attenuations[3] * tube_gain, up_buffer, up_nsamples);


            // ------------ Stage 4 ------------
            grid_conduction(up_buffer, up_nsamples);
            shelf_process(&preamp.cathode_bypass_filter4, up_buffer, up_nsamples);
            tube_sim(up_buffer, up_nsamples, preamp.stage4_bias[0], preamp.stage4_bias[1]);

            onepole_process_highpass(&preamp.coupling_filter4, up_buffer, up_nsamples);
            onepole_process_lowpass(&preamp.stage4LP, up_buffer, up_nsamples);

            local_const float dist_compensation_gain = dbtoa(-27.0f);
            apply_gain_linear(dist_compensation_gain, up_buffer, up_nsamples);


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
        
        const float b0 = ts.b0;
        const float b1 = ts.b1;
        const float b2 = ts.b2;
        const float b3 = ts.b3;
        const float a1 = ts.a1;
        const float a2 = ts.a2;
        const float a3 = ts.a3;
        
        float x1 = ts.x1;
        float x2 = ts.x2;
        float x3 = ts.x3;
        float y1 = ts.y1;
        float y2 = ts.y2;
        float y3 = ts.y3;
        
        for (u16 index = 0; index < nsamples; index++) {
            float input_sample = audio_buffer[index];
            float out_sample = input_sample * b0
                                + x1 * b1
                                + x2 * b2
                                + x3 * b3
                                + y1 * a1
                                + y2 * a2
                                + y3 * a3;

            x3 = x2;
            x2 = x1;
            x1 = input_sample;

            y3 = y2;
            y2 = y1;
            y1 = out_sample;
            audio_buffer[index] = out_sample;
        }
        
        ts.x1 = x1;
        ts.x2 = x2;
        ts.x3 = x3;
        ts.y1 = y1;
        ts.y2 = y2;
        ts.y3 = y3;
    }

    // res / pres
    shelf_process(&dsp.resonance, audio_buffer, nsamples);
    shelf_process(&dsp.presence, audio_buffer, nsamples);


    local_const float inv_fft_size = 1.0f/fft_size;

    // irloader
    {
        DSPData::IRLoader &irloader = dsp.irloader;

        memcpy(irloader.ols_buffer, irloader.ols_buffer + nsamples, (fft_size - nsamples) * sizeof(float));
        memcpy(&irloader.ols_buffer[fft_size-nsamples], audio_buffer, nsamples * sizeof(float));

        float *temp = irloader.fdl_ptrs[npartitions - 1];
        for (u16 index = npartitions - 1; index > 0; index--) {
            irloader.fdl_ptrs[index] = irloader.fdl_ptrs[index-1];
        }
        irloader.fdl_ptrs[0] = temp;

        pffft_transform(irloader.fft_setup, irloader.ols_buffer, irloader.fdl_ptrs[0],  nullptr, PFFFT_FORWARD);

        memset(irloader.convolution_dft, 0, dft_size * sizeof(float));

        for (u16 part_index = 0; part_index < npartitions; part_index++) {
            pffft_zconvolve_accumulate(irloader.fft_setup, irloader.fdl_ptrs[part_index], irloader.ir_parts_dft[part_index], irloader.convolution_dft, 1.0f);
        }

        pffft_transform(irloader.fft_setup, irloader.convolution_dft, irloader.signal_time_buffer, nullptr, PFFFT_BACKWARD);

        for (u32 index = 0; index < nsamples; index++) {
            audio_buffer[index] = irloader.signal_time_buffer[fft_size - nsamples + index] * inv_fft_size;
        }

    }

    for (u32 index = 0; index < nsamples; index++) {
        audio_buffer[index] *= dsp.params.volume;
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
    switches[FS1].Init(seed::D25, 1000.0f);
    switches[FS2].Init(seed::D26);


    hardware.adc.Init(adc_configs, nKNOBS);
    hardware.adc.Start();

    smooth_param_init(&dsp.params.low, 0.5f);
    smooth_param_init(&dsp.params.mid, 0.5f);
    smooth_param_init(&dsp.params.trebble, 0.5f);

    {
        DSPData::NoiseGate &gate = dsp.gate;

        gate.buffer_index = 0;
        gate.absolute_sum = 0.0;
        gate.threshold = dbtoa(gate_thresh[CLEAN][0]);
        // gate.return_gain = gate.threshold;
        gate.is_open = false;

        smooth_param_init(&gate.gain, 0.0f);
    }

    onepole_set_coeffs(&dsp.tight_filter, 400.0f, dsp.samplerate);

//     biquad_reset(&dsp.boost_filter);
    biquad_set_coeffs(&dsp.boost_filter, 1200.0f, 0.2, 8.0f, dsp.samplerate, BIQUAD_PEAK);

    {
        DSPData::Preamp &preamp = dsp.preamp;

        onepole_set_coeffs(&preamp.input_filter, 100.0f, dsp.samplerate);
        onepole_set_coeffs(&preamp.input_mud_filter, 400.0f, dsp.upsamplerate);

//         biquad_reset(&preamp.mid_boost);
        biquad_set_coeffs(&preamp.mid_boost, 1000.0f, 0.2f, 3.0f, dsp.upsamplerate, BIQUAD_PEAK);

//         // shelf_reset(&preamp.bright_cap_filter);
        shelf_set_coeffs(&preamp.bright_cap_filter, 550.0f, 0.0f, dsp.upsamplerate, lowshelf);


        onepole_set_coeffs(&preamp.coupling_filter1, 15.0f, dsp.upsamplerate);
        onepole_set_coeffs(&preamp.coupling_filter2, 15.0f, dsp.upsamplerate);
        onepole_set_coeffs(&preamp.coupling_filter3, 15.0f, dsp.upsamplerate);
        onepole_set_coeffs(&preamp.coupling_filter4, 15.0f, dsp.upsamplerate);

        onepole_set_coeffs(&preamp.stage0LP, 10000.0f, dsp.upsamplerate);
        onepole_set_coeffs(&preamp.stage1LP, 12000.0f, dsp.upsamplerate);
        onepole_set_coeffs(&preamp.stage2LP, 12000.0f, dsp.upsamplerate);
        onepole_set_coeffs(&preamp.stage3LP, 12000.0f, dsp.upsamplerate);
        onepole_set_coeffs(&preamp.stage4LP, 12000.0f, dsp.upsamplerate);
        
        u32 config_index = dsp.params.channel * 2 + dsp.params.preamp_variant;
        set_biases(&dsp.preamp, config_index);
        set_cathode_bypasses(&dsp.preamp, config_index);
        set_attenuations(&dsp.preamp, config_index);

//         biquad_reset(&preamp.oversampler.upsample_filter1);
        biquad_set_coeffs(&preamp.oversampler.upsample_filter1, dsp.samplerate/2.0f * 0.9f,
                            0.54119610f, 0.0f, dsp.upsamplerate, BIQUAD_LOWPASS);

//         biquad_reset(&preamp.oversampler.upsample_filter2);
        biquad_set_coeffs(&preamp.oversampler.upsample_filter2, dsp.samplerate/2.0f * 0.9f,
                            1.3065630f, 0.0f, dsp.upsamplerate, BIQUAD_LOWPASS);

//         biquad_reset(&preamp.oversampler.downsample_filter1);
        biquad_set_coeffs(&preamp.oversampler.downsample_filter1, dsp.samplerate/2.0f * 0.9f,
                            0.54119610f, 0.0f, dsp.upsamplerate, BIQUAD_LOWPASS);

//         biquad_reset(&preamp.oversampler.downsample_filter2);
        biquad_set_coeffs(&preamp.oversampler.downsample_filter2, dsp.samplerate/2.0f * 0.9f,
                            1.3065630f, 0.0f, dsp.upsamplerate, BIQUAD_LOWPASS);
        
    }

//     shelf_reset(&dsp.resonance);
    shelf_set_coeffs(&dsp.resonance, 120.0f, 12.0f, dsp.samplerate, lowshelf);

//     shelf_reset(&dsp.presence);
    shelf_set_coeffs(&dsp.presence, 500.0f, 10.0f, dsp.samplerate, highshelf);

    {
        dsp.irloader.fft_setup = pffft_new_setup(fft_size, PFFFT_REAL);
    
        for (u16 part_index = 0; part_index < npartitions; part_index++) {
            dsp.irloader.fdl_ptrs[part_index] = dsp.irloader.fdl[part_index];
        }
    
        for (u16 part_index = 0; part_index < npartitions; part_index++) {
            memset(dsp.irloader.signal_time_buffer, 0, fft_size * sizeof(float));
            memcpy(dsp.irloader.signal_time_buffer, &ir[part_index * BLOCK_SIZE], BLOCK_SIZE * sizeof(float));
    
            pffft_transform(dsp.irloader.fft_setup, dsp.irloader.signal_time_buffer, dsp.irloader.ir_parts_dft[part_index], nullptr, PFFFT_FORWARD);
        }
    }


    //Start calling the audio callback
    hardware.StartAudio(AudioCallback);

    // Loop forever
    while (1) {

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
}
