#include "daisy_seed.h"

// #define USE_DAISYSP_LGPL
#include "daisysp.h"

using namespace daisy;
using namespace daisysp;

enum ButtonsCodes {

    Toggle4 = 7,
    Toggle3 = 8,
    Toggle2 = 9,
    Toggle1 = 10,

    Knob1 = 16, 
    Knob2 = 17, 
    Knob3 = 18, 
    Knob4 = 19, 
    Knob5 = 20, 
    Knob6 = 21, 

    LED1 = 22,
    LED2 = 23,
    
    FS1 = 25,
    FS2 = 26,
};

DaisySeed  hardware;

// OnePole dans DaisySP/Source/Filters
// Biquad dans


struct DSPData {

    // int placeholder
    
    struct {
    
    } gate;
    
    struct {
    
    } preamp;
    
    struct {
    
    } tonestack;

    struct {

    } irloader;

} dsp_data;

void AudioCallback(AudioHandle::InterleavingInputBuffer  in,
                   AudioHandle::InterleavingOutputBuffer out,
                   size_t                                size)
{
    //Fill the block with samples
    for(size_t i = 0; i < size; i += 2)
    {
 
    }
}


int main(void)
{
    // Configure and Initialize the Daisy Seed
    // These are separate to allow reconfiguration of any of the internal
    // components before initialization.
    hardware.Configure();
    hardware.Init();
    hardware.SetAudioBlockSize(4);

    //How many samples we'll output per second
    float samplerate = hardware.AudioSampleRate();

    // spawn un tableau de n AdcChannelConfig pour n entrées
    // les init indépendament et les envoyer dans hardware.adc.init(&adccfg, n)

    //Create an ADC configuration
    AdcChannelConfig adcConfig;
    //Add pin 21 as an analog input in this config. We'll use this to read the knob
    adcConfig.InitSingle(hardware.GetPin(16));

    //Initialize the button on pin 28
    // button1.Init(hardware.GetPin(25), samplerate / 48.f);
    
    //Set the ADC to use our configuration
    hardware.adc.Init(&adcConfig, 1);


    //Start the adc
    hardware.adc.Start();

    //Start calling the audio callback
    hardware.StartAudio(AudioCallback);

    // Loop forever
    for(;;) {}
}
