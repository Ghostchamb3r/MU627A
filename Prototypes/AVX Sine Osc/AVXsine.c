#include "csdl.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <immintrin.h>

/*
To Compile:
gcc -O2 -mavx -mavx2 -shared -o AVXhypersaw.dll AVXhypersaw2.c -DUSE_DOUBLE -I"C:\Program Files\Csound6_x64\include\csound"
icx -O2 -mavx -mavx2 -shared -o AVXsine.dll AVXsine.c -DUSE_DOUBLE -I"C:\Program Files\Csound6_x64\include\csound"

*/


// DATA
typedef struct _AVXsine{
OPDS h;
MYFLT *out;
MYFLT *in1, *in2, *in3;
MYFLT freq;
MYFLT amp;
MYFLT phase;
} AVXsine;

// INITIALIZE 

int AVXsine_init(CSOUND *csound, AVXsine *p){
     
    p->phase = 0.0f;
    // p->out = (MYFLT*)_aligned_malloc(sizeof(MYFLT) * CS_KSMPS, 32);  // Align for AVX


return OK;
}


// PROCESS

int AVXsine_process_audio(CSOUND *csound, AVXsine *p){
int numSamples = CS_KSMPS;
MYFLT freq = *p->in1;
MYFLT amp = *p->in2;
if (amp == 0) {
    amp = 1.0;  // Artificially boost it to confirm sound
}
// csound->Message(csound, "Amplitude Before Processing: %f\n", *p->in2);

MYFLT phase = *p->in3;

MYFLT sr = csound->GetSr(csound);
// float phaseIncrement = (2 * M_PI * freq) / sr;

// Compute Sine Wave
// Buffer to store phase increments per sample
float phaseIncBuffer[CS_KSMPS];
float sineBuffer[CS_KSMPS];  // Store sine values before AVX processing


// Generate sine wave using dynamic phase increments
__m256 sinewave;
// float* outputBuffer = (float*)p->out;
MYFLT *aout = p->out;

for (int i = 0; i < numSamples; i += 8) {

    phaseIncBuffer[i] = (2 * M_PI * freq) / sr;  // Ensure this updates per sample

    __m256 phaseOffsets = _mm256_set_ps(
        p->phase,
        p->phase + phaseIncBuffer[i],
        p->phase + phaseIncBuffer[i] * 2,
        p->phase + phaseIncBuffer[i] * 3,
        p->phase + phaseIncBuffer[i] * 4,
        p->phase + phaseIncBuffer[i] * 5,
        p->phase + phaseIncBuffer[i] * 6,
        p->phase + phaseIncBuffer[i] * 7
    );

    sinewave = _mm256_mul_ps(_mm256_set1_ps(amp), _mm256_sin_ps(phaseOffsets));
    _mm256_storeu_ps(&sineBuffer[i], sinewave);

// aout[i] = sineBuffer[i];
// csound->Message(csound, "Output Buffer[%d]: %f\n", i, p->out[i]);

float temp[8];  // Temporary array for AVX results
_mm256_storeu_ps(temp, sinewave);  // Store AVX values

for (int j = 0; j < 8; j++) {
    if (isnan(temp[j])) temp[j] = 0.0;  // Replace NaN with silence
aout[i + j] = (MYFLT)temp[j];

}
// csound->Message(csound, "Output Buffer after out[%d]: %f\n", i, p->out[i]);
// csound->Message(csound, "Phase Before Update: %f\n", p->phase);
    // Update phase correctly
    p->phase += phaseIncBuffer[i] * 8;
    if (p->phase > 2 * M_PI) {
        p->phase -= 2 * M_PI;
    }
        // csound->Message(csound, "Output Buffer after phase wrap[%d]: %f\n", i, p->out[i]);
// csound->Message(csound, "Phase After Wrap: %f\n", p->phase);
}




return OK;
}
//END PROCESSING


// REGISTER

static OENTRY localops[] = {
{ "AVXsine", sizeof(AVXsine), 0, 7, "a", "iii",(SUBR) AVXsine_init,
(SUBR) AVXsine_process_audio  }
};


LINKAGE

