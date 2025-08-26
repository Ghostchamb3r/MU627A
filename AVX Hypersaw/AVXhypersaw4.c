#include "csdl.h"
#include <math.h>
#define _USE_MATH_DEFINES
#include <emmintrin.h>  // For SSE2
#include <immintrin.h>  // For AVX

static inline __m256 approx_sin(__m256 x) {
    // So we have to use a polynomial approximation for sine (simplified) since there's no predefined sin/cos functions in the intrinsics
    __m256 x2 = _mm256_mul_ps(x, x);
    return _mm256_sub_ps(x, _mm256_mul_ps(x2, _mm256_set1_ps(0.16605)));  
}

static inline __m256 approx_cos(__m256 x) {
    __m256 x2 = _mm256_mul_ps(x, x);
    return _mm256_sub_ps(_mm256_set1_ps(1.0f), _mm256_mul_ps(x2, _mm256_set1_ps(0.5f))); 
}

/*
To Compile:

This is the correct command lol
gcc -O2 -dynamiclib -o hyperhermiteosc.dylib hyperhermiteosc.c -DUSE_DOUBLE -I/Library/Frameworks/CsoundLib64.framework/Versions/6.0/Headers -arch x86_64
gcc -O2 -mavx -mavx2 -shared -o AVXsaw.dll AVXhypersaw4.c -DUSE_DOUBLE -I"C:\Program Files\Csound6_x64\include\csound"

*/
/*
Might need to manually define pi
const float pi = 3.141592653589793;
*/
// implement tanh as a table lookup instead of calling tanh every sample
// Compute tanh

/*
Hermite interpolation reconfigured from:
https://people.sc.fsu.edu/~jburkardt/c_src/hermite_interpolant/hermite_interpolant.c
*/

#define TANH_TABLE_SIZE 1024  // Adjust for desired resolution
#define TANH_MIN -1.0  // Range mapping
#define TANH_MAX 1.0  

static MYFLT tanhTable[TANH_TABLE_SIZE];
static MYFLT tanhDerivTable[TANH_TABLE_SIZE]; // Store tanh' (x)

void generateTanhTable() { // we can possibly use AVX here to precompute the tables faster
    for (int i = 0; i < TANH_TABLE_SIZE; i++) {
        MYFLT x = TANH_MIN + (i * (TANH_MAX - TANH_MIN) / (TANH_TABLE_SIZE - 1));
        tanhTable[i] = tanh(x);
        tanhDerivTable[i] = 1 - tanhTable[i] * tanhTable[i]; // tanh'(x) = 1 - tanh^2(x)
    }
}


static MYFLT computeTanhHermite(MYFLT x) {
    if (x <= TANH_MIN) return tanhTable[0];
    if (x >= TANH_MAX) return tanhTable[TANH_TABLE_SIZE - 1];

    MYFLT normIndex = ((x - TANH_MIN) / (TANH_MAX - TANH_MIN)) * (TANH_TABLE_SIZE - 1); //2.98485588666341e-4
    int index = (int)normIndex; // 2
    MYFLT frac = normIndex - index; // .98485588666341e-4

    MYFLT f0 = tanhTable[index];
    MYFLT f1 = tanhTable[index + 1];
    MYFLT d0 = tanhDerivTable[index];
    MYFLT d1 = tanhDerivTable[index + 1];

    // Cubic Hermite interpolation formula
    MYFLT h00 = (1 + 2 * frac) * (1 - frac) * (1 - frac);  // 1.96971177332682e-4 * 0.9999015144113337 * 0.9999015144113337
    MYFLT h10 = frac * (1 - frac) * (1 - frac);
    MYFLT h01 = frac * frac * (3 - 2 * frac);
    MYFLT h11 = frac * frac * (frac - 1);

    return h00 * f0 + h10 * d0 + h01 * f1 + h11 * d1;
}




// DATA
typedef struct _AVXsaw{
OPDS h;
MYFLT *out;
MYFLT *in1, *in2, *in3, *in4;
MYFLT freq;
MYFLT amp;
MYFLT phase;
MYFLT m;
} AVXsaw;

// INITIALIZE 

int AVXsaw_init(CSOUND *csound, AVXsaw *p){
    p->phase = 0.0f;
    p->m = *p->in4;
    if (p->m < 0.0) p->m = 0.0;
    if (p->m > 1.0) p->m = 1.0;
    generateTanhTable();  // Precompute tanh values

return OK;
}

// PROCESS

int AVXsaw_process_audio(CSOUND *csound, AVXsaw *p){

int inNumSamples = CS_KSMPS;
MYFLT *aout = p->out;
MYFLT freq = *p->in1;
MYFLT amp = *p->in2;
MYFLT phase = *p->in3;
MYFLT m = *p->in4;
MYFLT sr = csound->GetSr(csound);
float phaseIncBuffer[CS_KSMPS];
float sinewavefin[CS_KSMPS]; // final sine buffer to pass temp to
float sineBuffer[CS_KSMPS];  // Store sine values before AVX processing
__m256 sinewave;

float phaseIncrement = (2 * M_PI * freq) / sr;

for(int i = 0; i < inNumSamples; i += 8){

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

    sinewave = _mm256_mul_ps(_mm256_set1_ps(amp), approx_sin(phaseOffsets));
    _mm256_storeu_ps(&sineBuffer[i], sinewave);

float temp[8];  // Temporary array for AVX results
_mm256_storeu_ps(temp, sinewave);  // Store AVX values
        // BEGIN CONDITION LOOP
            for (int j = 0; j < 8; j++) {
                if (isnan(temp[j])) temp[j] = 0.0;  // Replace NaN with silence
                sinewavefin[i + j] = (MYFLT)temp[j];


                __m256 k = _mm256_set1_ps(12000 / (freq * log10(freq)) );
                float shapedwave[8];  // Temporary storage
                _mm256_storeu_ps(shapedwave, _mm256_mul_ps(k, sinewave));
                    // BEGIN CONDITION LOOP
                    for (int j = 0; j < 8; j++) {
                    shapedwave[j] = computeTanhHermite(shapedwave[j]);
                    }
                    // END CONDITION LOOP
            float ringmodder[8];
            float sawtoothwave[8];
                    // BEGIN CONDITON LOOP
                    for (int j = 0; j < 8; j++) {
                    ringmodder[j] = shapedwave[j] * cos(p->phase);
                    sawtoothwave[j] = shapedwave[j] * (1 - p->m / 2) + ringmodder[j] * p->m;
                    }
                    // END CONDITION LOOP

            float sawtooth[8];
            _mm256_storeu_ps(sawtooth, _mm256_add_ps(
            _mm256_mul_ps(_mm256_set1_ps(1 - p->m / 2), _mm256_loadu_ps(shapedwave)),
            _mm256_mul_ps(_mm256_set1_ps(p->m), approx_cos(_mm256_set1_ps(p->phase)))
            ));
                    // BEGIN CONDITION LOOP
                    for (int j = 0; j < 8; j++) {
                    aout[i + j] = sawtooth[j];  
                    }
                    // END CONDITION LOOP
            } // END CONDITION LOOP
p->phase += phaseIncrement * 8;  
if (p->phase > 2 * M_PI) p->phase -= 2 * M_PI;

}// END INTIIAL CONDITION LOOP
return OK;

}//END PROCESSING

// REGISTER

static OENTRY localops[] = {
{ "AVXsaw", sizeof(AVXsaw), 0, 7, "a", "iiii",(SUBR) AVXsaw_init,
(SUBR) AVXsaw_process_audio }
};


LINKAGE

