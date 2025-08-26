#include "csdl.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <immintrin.h>


#define TANH_TABLE_SIZE 1024  // Adjust for desired resolution
#define TANH_MIN -1.0  // Range mapping
#define TANH_MAX 1.0  

static float tanhTable[TANH_TABLE_SIZE];
static float tanhDerivTable[TANH_TABLE_SIZE]; // Store tanh' (x)



void TanhTables() {
const __m256 Table_size = _mm256_set1_ps(TANH_TABLE_SIZE);
const __m256 Table_min = _mm256_set1_ps(TANH_MIN);
const __m256 Table_max = _mm256_set1_ps(TANH_MAX);
const __m256 one = _mm256_set1_ps(1.0f);
const __m256 two = _mm256_set1_ps(2.0f);
const __m256 three = _mm256_set1_ps(3.0f);

    for (int i = 0; i < TANH_TABLE_SIZE; i+=8){
        __m256 index = _mm256_set_ps(i+7, i+6, i+5, i+4, i+3, i+2, i+1, i);
        __m256 term1 = _mm256_mul_ps(index, _mm256_sub_ps(Table_max, Table_min) );
        __m256 term2 = _mm256_sub_ps(Table_size, one);
        __m256 term3 = _mm256_div_ps(term1, term2);
        __m256 x = _mm256_add_ps(Table_min, term3);
        __m256 tanhresult = _mm256_tanh_ps(x);
        _mm256_storeu_ps(&tanhTable[i], tanhresult);
        __m256 tanhderivs = _mm256_sub_ps(one,_mm256_mul_ps(tanhresult, tanhresult));
        _mm256_storeu_ps(&tanhDerivTable[i], tanhderivs);
    }
}

__m256 AVXtanhApprox(__m256 x){
const __m256 Table_size = _mm256_set1_ps(TANH_TABLE_SIZE);
const __m256 Table_min = _mm256_set1_ps(TANH_MIN);
const __m256 Table_max = _mm256_set1_ps(TANH_MAX);
const __m256 one = _mm256_set1_ps(1.0f);
const __m256 two = _mm256_set1_ps(2.0f);
const __m256 three = _mm256_set1_ps(3.0f);

__m256 term1 = _mm256_div_ps(_mm256_sub_ps(x, Table_min), _mm256_sub_ps(Table_max, Table_min));
__m256 term2 = _mm256_sub_ps(Table_size, one);
__m256 normindex = _mm256_mul_ps(term1, term2);

// truncate and cast to float
__m256i indextemp = _mm256_cvttps_epi32(normindex);
    //Clamp 
    __m256i maxIndex = _mm256_set1_epi32(TANH_TABLE_SIZE - 2);
    indextemp = _mm256_min_epi32(indextemp, maxIndex);
    //cast to float
// __m256 index = _mm256_castsi256_ps(indextemp);
   __m256 index = _mm256_cvtepi32_ps(indextemp);

__m256 frac = _mm256_sub_ps(normindex, index);


// index + 1
__m256i index1 = _mm256_add_epi32(indextemp, _mm256_set1_epi32(1));

// non-contigious memory access patterns; need to load data from scattered memory locations into a vector register.
__m256 f0 = _mm256_i32gather_ps(tanhTable, indextemp, 4);
__m256 f1 = _mm256_i32gather_ps(tanhTable, index1, 4);
__m256 d0 = _mm256_i32gather_ps(tanhDerivTable, indextemp, 4);
__m256 d1 = _mm256_i32gather_ps(tanhDerivTable, index1, 4);


__m256 h00 = _mm256_mul_ps(_mm256_add_ps(one, _mm256_mul_ps(two,frac)), _mm256_mul_ps(_mm256_sub_ps(one, frac), _mm256_sub_ps(one, frac)));
__m256 h10 = _mm256_mul_ps(frac, _mm256_mul_ps(_mm256_sub_ps(one, frac), _mm256_sub_ps(one, frac)));
__m256 h01 = _mm256_mul_ps(frac, _mm256_mul_ps(frac, _mm256_sub_ps(three, _mm256_mul_ps(two, frac))));
__m256 h11 = _mm256_mul_ps(frac, _mm256_mul_ps(frac, _mm256_sub_ps(frac, one)));


__m256 term3 = _mm256_mul_ps(h00, f0);
__m256 term4 = _mm256_mul_ps(h10, d0);
__m256 term5 = _mm256_mul_ps(h01, f1);
__m256 term6 = _mm256_mul_ps(h11, d1);

__m256 term7 = _mm256_add_ps(term3, term4);
__m256 term8 = _mm256_add_ps(term5, term6);

return _mm256_add_ps(term7,term8);
//__m256 result = _mm256_add_ps(term7,term8);
// return result;

}



/*
To Compile:

gcc -O2 -mavx -mavx2 -shared -o AVXhypersaw.dll AVXconfusion.c -DUSE_DOUBLE -I"C:\Program Files\Csound6_x64\include\csound"
icx -O2 -mavx -mavx2 -shared -o AVXhypersaw.dll AVXconfusion.c -DUSE_DOUBLE -I"C:\Program Files\Csound6_x64\include\csound"
icx -O2 -mavx -mavx2 -shared -o AVXhypersaw.dll AVXconfusion.c -I"C:\Program Files\Csound6_x64\include\csound"

*/
// DATA
typedef struct _AVXhypersaw{
OPDS h;
MYFLT *out;
MYFLT *in1, *in2, *in3, *in4;
MYFLT freq;
MYFLT amp;
MYFLT phase;
MYFLT m;
} AVXhypersaw;

// INITIALIZE 

int AVXhypersaw_init(CSOUND *csound, AVXhypersaw *p){
     
    p->phase = 0.0f;
    // p->out = (MYFLT*)_aligned_malloc(sizeof(MYFLT) * CS_KSMPS, 32);  // Align for AVX
    p->m = *p->in4;
    if (p->m < 0.0) p->m = 0.0;
    if (p->m > 1.0) p->m = 1.0;
    void TanhTables();

return OK;
}


// PROCESS

int AVXhypersaw_process_audio(CSOUND *csound, AVXhypersaw *p){

const __m256 Table_size = _mm256_set1_ps(TANH_TABLE_SIZE);
const __m256 Table_min = _mm256_set1_ps(TANH_MIN);
const __m256 Table_max = _mm256_set1_ps(TANH_MAX);
const __m256 one = _mm256_set1_ps(1.0f);
const __m256 two = _mm256_set1_ps(2.0f);
const __m256 three = _mm256_set1_ps(3.0f);

int numSamples = CS_KSMPS;
MYFLT freq = *p->in1;
MYFLT amp = *p->in2;
if (amp == 0) {
    amp = 1.0;  // Artificially boost it to confirm sound
}
// csound->Message(csound, "Amplitude Before Processing: %f\n", *p->in2);

MYFLT phase = *p->in3;
MYFLT m = *p->in4;

MYFLT sr = csound->GetSr(csound);
// float phaseIncrement = (2 * M_PI * freq) / sr;

// Compute Sine Wave
// Buffer to store phase increments per sample
float phaseIncBuffer[CS_KSMPS];


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

float phasef = (float)phase;
float mf = (float)m;

    sinewave = _mm256_mul_ps(_mm256_set1_ps(amp), _mm256_sin_ps(phaseOffsets));
    float k = 12000 / (freq * log10(freq));
    __m256 kVec = _mm256_set1_ps(k);
    __m256 weightedwave = _mm256_mul_ps(kVec, sinewave);
    __m256 shapedwave = AVXtanhApprox(weightedwave);
    // __m256 cosineVec = _mm256_cos_ps(p->phase);
    __m256 phasefVec = _mm256_set1_ps(phasef);
    __m256 cosineVec = _mm256_cos_ps(phasefVec);
    __m256 ringmod = _mm256_mul_ps(shapedwave, cosineVec);
    // __m256 emVec= _mm256_set1_ps(p->m);
    __m256 emVec= _mm256_set1_ps(mf);
    __m256 modterm1 = _mm256_add_ps(_mm256_sub_ps(one, _mm256_div_ps(emVec, two)), ringmod);
    __m256 modterm2 = _mm256_mul_ps(modterm1, emVec);
    __m256 sawtoothwave = _mm256_mul_ps(weightedwave, modterm2);


float temp[8];  // Temporary array for AVX results
_mm256_storeu_ps(temp, sawtoothwave);  // Store AVX values

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
{ "AVXhypersaw", sizeof(AVXhypersaw), 0, 7, "a", "iiii",(SUBR) AVXhypersaw_init,
(SUBR) AVXhypersaw_process_audio  }
};


LINKAGE

