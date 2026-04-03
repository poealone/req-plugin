/*
 * ReQ — Parametric EQ for PocketDAW
 * 4-band stereo biquad EQ: Low Shelf | Peak | Peak | High Shelf
 *
 * Parameters (0.0–1.0 normalized):
 *   0  Band 1 Freq    20–500 Hz     (low shelf)
 *   1  Band 1 Gain    -12 to +12 dB
 *   2  Band 2 Freq    80–2000 Hz    (peak)
 *   3  Band 2 Gain    -12 to +12 dB
 *   4  Band 2 Q       0.3–4.0
 *   5  Band 3 Freq    800–12000 Hz  (peak)
 *   6  Band 3 Gain    -12 to +12 dB
 *   7  Band 3 Q       0.3–4.0
 *   8  Band 4 Freq    2000–20000 Hz (high shelf)
 *   9  Band 4 Gain    -12 to +12 dB
 *  10  Master Gain    -6 to +6 dB
 *  11  (reserved)
 *
 * Build:
 *   gcc -shared -fPIC -O2 -o req.so req.c -lm
 *   aarch64-none-linux-gnu-gcc -shared -fPIC -O2 -o req.so req.c -lm
 *   x86_64-w64-mingw32-gcc -shared -o req.dll req.c -lm
 */

#include "pocketdaw.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ── Biquad filter state (Direct Form I) ── */
typedef struct {
    double b0, b1, b2, a1, a2; /* coefficients */
    double x1L, x2L, y1L, y2L; /* left  channel state */
    double x1R, x2R, y1R, y2R; /* right channel state */
} Biquad;

static void biquad_process(Biquad* f, const float* inL, const float* inR,
                            float* outL, float* outR, int n) {
    double b0=f->b0, b1=f->b1, b2=f->b2, a1=f->a1, a2=f->a2;
    double x1L=f->x1L, x2L=f->x2L, y1L=f->y1L, y2L=f->y2L;
    double x1R=f->x1R, x2R=f->x2R, y1R=f->y1R, y2R=f->y2R;
    for (int i = 0; i < n; i++) {
        double xL = (double)inL[i];
        double yL = b0*xL + b1*x1L + b2*x2L - a1*y1L - a2*y2L;
        x2L=x1L; x1L=xL; y2L=y1L; y1L=yL;
        outL[i] = (float)yL;

        double xR = (double)inR[i];
        double yR = b0*xR + b1*x1R + b2*x2R - a1*y1R - a2*y2R;
        x2R=x1R; x1R=xR; y2R=y1R; y1R=yR;
        outR[i] = (float)yR;
    }
    f->x1L=x1L; f->x2L=x2L; f->y1L=y1L; f->y2L=y2L;
    f->x1R=x1R; f->x2R=x2R; f->y1R=y1R; f->y2R=y2R;
}

/* ── Coefficient computation ── */

/* Low shelf — RBJ cookbook */
static void coeff_low_shelf(Biquad* f, double sr, double freq, double gainDb) {
    double A  = pow(10.0, gainDb / 40.0);
    double w0 = 2.0 * M_PI * freq / sr;
    double cw = cos(w0), sw = sin(w0);
    double S  = 1.0;
    double a  = sw / 2.0 * sqrt((A + 1.0/A) * (1.0/S - 1.0) + 2.0);
    double b0 =  A*((A+1.0) - (A-1.0)*cw + 2.0*sqrt(A)*a);
    double b1 =  2.0*A*((A-1.0) - (A+1.0)*cw);
    double b2 =  A*((A+1.0) - (A-1.0)*cw - 2.0*sqrt(A)*a);
    double a0 =  (A+1.0) + (A-1.0)*cw + 2.0*sqrt(A)*a;
    double a1 = -2.0*((A-1.0) + (A+1.0)*cw);
    double a2 =  (A+1.0) + (A-1.0)*cw - 2.0*sqrt(A)*a;
    f->b0=b0/a0; f->b1=b1/a0; f->b2=b2/a0;
    f->a1=a1/a0; f->a2=a2/a0;
}

/* Peaking EQ — RBJ cookbook */
static void coeff_peak(Biquad* f, double sr, double freq, double gainDb, double Q) {
    double A  = pow(10.0, gainDb / 40.0);
    double w0 = 2.0 * M_PI * freq / sr;
    double alpha = sin(w0) / (2.0 * Q);
    double b0 =  1.0 + alpha * A;
    double b1 = -2.0 * cos(w0);
    double b2 =  1.0 - alpha * A;
    double a0 =  1.0 + alpha / A;
    double a1 = -2.0 * cos(w0);
    double a2 =  1.0 - alpha / A;
    f->b0=b0/a0; f->b1=b1/a0; f->b2=b2/a0;
    f->a1=a1/a0; f->a2=a2/a0;
}

/* High shelf — RBJ cookbook */
static void coeff_high_shelf(Biquad* f, double sr, double freq, double gainDb) {
    double A  = pow(10.0, gainDb / 40.0);
    double w0 = 2.0 * M_PI * freq / sr;
    double cw = cos(w0), sw = sin(w0);
    double S  = 1.0;
    double a  = sw / 2.0 * sqrt((A + 1.0/A) * (1.0/S - 1.0) + 2.0);
    double b0 =  A*((A+1.0) + (A-1.0)*cw + 2.0*sqrt(A)*a);
    double b1 = -2.0*A*((A-1.0) + (A+1.0)*cw);
    double b2 =  A*((A+1.0) + (A-1.0)*cw - 2.0*sqrt(A)*a);
    double a0 =  (A+1.0) - (A-1.0)*cw + 2.0*sqrt(A)*a;
    double a1 =  2.0*((A-1.0) - (A+1.0)*cw);
    double a2 =  (A+1.0) - (A-1.0)*cw - 2.0*sqrt(A)*a;
    f->b0=b0/a0; f->b1=b1/a0; f->b2=b2/a0;
    f->a1=a1/a0; f->a2=a2/a0;
}

/* ── Plugin state ── */
#define NUM_BANDS 4
#define NUM_PARAMS 12

typedef struct {
    float sr;
    float p[NUM_PARAMS];    /* normalized params */
    Biquad bands[NUM_BANDS];
    int   dirty;            /* recompute coeffs next process */
    int   selectedNode;     /* for custom draw highlight */
    char  fmt[64];          /* format_param buffer */

    /* Computed real values (for draw) */
    float bandFreq[NUM_BANDS];
    float bandGain[NUM_BANDS];
    float bandQ[NUM_BANDS];
} ReQ;

/* Map normalized params to real values */
static void _update_values(ReQ* eq) {
    eq->bandFreq[0] = 20.0f   + eq->p[0] * 480.0f;   /* 20–500 Hz  */
    eq->bandGain[0] = -12.0f  + eq->p[1] * 24.0f;    /* ±12 dB     */
    eq->bandQ[0]    = 0.7f;                            /* fixed shelf slope */

    eq->bandFreq[1] = 80.0f   + eq->p[2] * 1920.0f;  /* 80–2000 Hz */
    eq->bandGain[1] = -12.0f  + eq->p[3] * 24.0f;
    eq->bandQ[1]    = 0.3f    + eq->p[4] * 3.7f;     /* 0.3–4.0    */

    eq->bandFreq[2] = 800.0f  + eq->p[5] * 11200.0f; /* 800–12k Hz */
    eq->bandGain[2] = -12.0f  + eq->p[6] * 24.0f;
    eq->bandQ[2]    = 0.3f    + eq->p[7] * 3.7f;

    eq->bandFreq[3] = 2000.0f + eq->p[8] * 18000.0f; /* 2k–20k Hz  */
    eq->bandGain[3] = -12.0f  + eq->p[9] * 24.0f;
    eq->bandQ[3]    = 0.7f;
}

static void _recompute(ReQ* eq) {
    _update_values(eq);
    double sr = (double)eq->sr;
    coeff_low_shelf (&eq->bands[0], sr, eq->bandFreq[0], eq->bandGain[0]);
    coeff_peak      (&eq->bands[1], sr, eq->bandFreq[1], eq->bandGain[1], eq->bandQ[1]);
    coeff_peak      (&eq->bands[2], sr, eq->bandFreq[2], eq->bandGain[2], eq->bandQ[2]);
    coeff_high_shelf(&eq->bands[3], sr, eq->bandFreq[3], eq->bandGain[3]);
    eq->dirty = 0;
}

/* ── Required exports ── */

int pdfx_api_version(void) { return 3; }
const char* pdfx_name(void) { return "ReQ"; }

PdFxInstance pdfx_create(float sr) {
    ReQ* eq = calloc(1, sizeof(ReQ));
    eq->sr = sr;
    /* Defaults: bands centered, all at 0 dB, Q=1 */
    eq->p[0]  = 0.2f;  /* B1 freq  ~116 Hz  */
    eq->p[1]  = 0.5f;  /* B1 gain   0 dB    */
    eq->p[2]  = 0.2f;  /* B2 freq  ~464 Hz  */
    eq->p[3]  = 0.5f;  /* B2 gain   0 dB    */
    eq->p[4]  = 0.2f;  /* B2 Q     ~1.0     */
    eq->p[5]  = 0.2f;  /* B3 freq  ~3.0 kHz */
    eq->p[6]  = 0.5f;  /* B3 gain   0 dB    */
    eq->p[7]  = 0.2f;  /* B3 Q     ~1.0     */
    eq->p[8]  = 0.5f;  /* B4 freq  ~11 kHz  */
    eq->p[9]  = 0.5f;  /* B4 gain   0 dB    */
    eq->p[10] = 0.5f;  /* master    0 dB    */
    eq->selectedNode = -1;
    eq->dirty = 1;
    _recompute(eq);
    return eq;
}

void pdfx_destroy(PdFxInstance inst) { free(inst); }

void pdfx_process(PdFxInstance inst, PdFxAudio* a) {
    ReQ* eq = (ReQ*)inst;
    if (eq->dirty) _recompute(eq);

    float master = powf(10.0f, (-6.0f + eq->p[10] * 12.0f) / 20.0f);

    /* Temp buffers on the stack — max 4096 frames */
    int n = a->bufferSize;
    float tL[4096], tR[4096];

    /* Chain bands */
    biquad_process(&eq->bands[0], a->inputL, a->inputR, tL, tR, n);
    biquad_process(&eq->bands[1], tL, tR, a->outputL, a->outputR, n);
    biquad_process(&eq->bands[2], a->outputL, a->outputR, tL, tR, n);
    biquad_process(&eq->bands[3], tL, tR, a->outputL, a->outputR, n);

    /* Apply master gain */
    for (int i = 0; i < n; i++) {
        a->outputL[i] *= master;
        a->outputR[i] *= master;
    }
}

/* ── Parameters ── */

int pdfx_param_count(void) { return NUM_PARAMS; }

const char* pdfx_param_name(int i) {
    const char* names[] = {
        "B1 Freq", "B1 Gain",
        "B2 Freq", "B2 Gain", "B2 Q",
        "B3 Freq", "B3 Gain", "B3 Q",
        "B4 Freq", "B4 Gain",
        "Master",  "---"
    };
    return (i >= 0 && i < NUM_PARAMS) ? names[i] : NULL;
}

const char* pdfx_param_group(int i) {
    if (i <= 1) return "Band 1 (Low Shelf)";
    if (i <= 4) return "Band 2 (Peak)";
    if (i <= 7) return "Band 3 (Peak)";
    if (i <= 9) return "Band 4 (High Shelf)";
    return "Output";
}

float pdfx_get_param(PdFxInstance inst, int i) {
    return (i >= 0 && i < NUM_PARAMS) ? ((ReQ*)inst)->p[i] : 0.0f;
}

void pdfx_set_param(PdFxInstance inst, int i, float v) {
    ReQ* eq = (ReQ*)inst;
    if (i >= 0 && i < NUM_PARAMS) {
        eq->p[i] = v;
        eq->dirty = 1;
    }
}

const char* pdfx_format_param(PdFxInstance inst, int i) {
    ReQ* eq = (ReQ*)inst;
    _update_values(eq);
    if (i == 0) { snprintf(eq->fmt, 64, "%.0f Hz", eq->bandFreq[0]); return eq->fmt; }
    if (i == 1) { snprintf(eq->fmt, 64, "%+.1f dB", eq->bandGain[0]); return eq->fmt; }
    if (i == 2) { snprintf(eq->fmt, 64, "%.0f Hz", eq->bandFreq[1]); return eq->fmt; }
    if (i == 3) { snprintf(eq->fmt, 64, "%+.1f dB", eq->bandGain[1]); return eq->fmt; }
    if (i == 4) { snprintf(eq->fmt, 64, "%.2f Q", eq->bandQ[1]); return eq->fmt; }
    if (i == 5) {
        float f = eq->bandFreq[2];
        if (f >= 1000.0f) snprintf(eq->fmt, 64, "%.1f kHz", f/1000.0f);
        else snprintf(eq->fmt, 64, "%.0f Hz", f);
        return eq->fmt;
    }
    if (i == 6) { snprintf(eq->fmt, 64, "%+.1f dB", eq->bandGain[2]); return eq->fmt; }
    if (i == 7) { snprintf(eq->fmt, 64, "%.2f Q", eq->bandQ[2]); return eq->fmt; }
    if (i == 8) {
        float f = eq->bandFreq[3];
        if (f >= 1000.0f) snprintf(eq->fmt, 64, "%.1f kHz", f/1000.0f);
        else snprintf(eq->fmt, 64, "%.0f Hz", f);
        return eq->fmt;
    }
    if (i == 9)  { snprintf(eq->fmt, 64, "%+.1f dB", eq->bandGain[3]); return eq->fmt; }
    if (i == 10) { snprintf(eq->fmt, 64, "%+.1f dB", -6.0f + eq->p[10]*12.0f); return eq->fmt; }
    return NULL;
}

void pdfx_get_accent_color(PdFxInstance inst, uint8_t* r, uint8_t* g, uint8_t* b) {
    *r = 0; *g = 255; *b = 65; /* PocketDAW green */
}

/* ── Presets ── */

int pdfx_preset_count(void) { return 5; }

const char* pdfx_preset_name(int i) {
    const char* n[] = {"Flat","Low Cut","Vocal Presence","Bass Boost","Air"};
    return (i >= 0 && i < 5) ? n[i] : NULL;
}

void pdfx_load_preset(PdFxInstance inst, int i) {
    ReQ* eq = (ReQ*)inst;
    float pre[5][NUM_PARAMS] = {
        /* Flat */
        {0.2f,0.5f, 0.2f,0.5f,0.2f, 0.2f,0.5f,0.2f, 0.5f,0.5f, 0.5f,0.5f},
        /* Low Cut — cut lows, boost presence */
        {0.05f,0.2f, 0.15f,0.5f,0.3f, 0.35f,0.65f,0.25f, 0.6f,0.5f, 0.5f,0.5f},
        /* Vocal Presence — dip 300Hz, boost 3kHz */
        {0.2f,0.45f, 0.15f,0.35f,0.3f, 0.22f,0.7f,0.4f, 0.6f,0.55f, 0.5f,0.5f},
        /* Bass Boost — low shelf up, high shelf neutral */
        {0.15f,0.7f, 0.1f,0.6f,0.2f, 0.2f,0.5f,0.2f, 0.5f,0.5f, 0.45f,0.5f},
        /* Air — high shelf boost */
        {0.2f,0.5f, 0.2f,0.5f,0.2f, 0.3f,0.5f,0.2f, 0.65f,0.7f, 0.5f,0.5f},
    };
    if (i >= 0 && i < 5) {
        for (int j = 0; j < NUM_PARAMS; j++) eq->p[j] = pre[i][j];
        eq->dirty = 1;
    }
}

/* ── Custom draw — EQ curve + draggable nodes ── */

/* Convert frequency to X pixel (log scale) */
static int freqToX(float freq, int x, int w) {
    float logMin = log10f(20.0f);
    float logMax = log10f(20000.0f);
    float t = (log10f(freq < 20.0f ? 20.0f : freq) - logMin) / (logMax - logMin);
    return x + (int)(t * w);
}

/* Convert dB gain to Y pixel */
static int gainToY(float dB, int y, int h) {
    float t = (dB + 12.0f) / 24.0f; /* 0=bottom(-12dB), 1=top(+12dB) */
    return y + h - (int)(t * h);
}

/* Compute magnitude response of all 4 bands at a given frequency (dB) */
static float responseDb(ReQ* eq, float freq) {
    double w = 2.0 * M_PI * freq / (double)eq->sr;
    double totalDb = 0.0;

    for (int b = 0; b < NUM_BANDS; b++) {
        Biquad* f = &eq->bands[b];
        /* H(z) = (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2) */
        /* At z = e^(jw): evaluate numerator and denominator */
        double re_n = f->b0 + f->b1*cos(w) + f->b2*cos(2*w);
        double im_n =        -f->b1*sin(w) - f->b2*sin(2*w);
        double re_d = 1.0  + f->a1*cos(w) + f->a2*cos(2*w);
        double im_d =        -f->a1*sin(w) - f->a2*sin(2*w);
        double mag2 = (re_n*re_n + im_n*im_n) / (re_d*re_d + im_d*im_d + 1e-30);
        totalDb += 10.0 * log10(mag2 < 1e-30 ? 1e-30 : mag2);
    }
    return (float)totalDb;
}

static const uint8_t NODE_COLORS[NUM_BANDS][3] = {
    {255, 160,  40},  /* Band 1 — amber      */
    {  0, 200, 255},  /* Band 2 — cyan       */
    {200,  80, 255},  /* Band 3 — purple     */
    {255,  60,  60},  /* Band 4 — red        */
};

int pdfx_draw(PdFxInstance inst, SDL_Renderer* r, PdDrawContext* ctx) {
    ReQ* eq = (ReQ*)inst;
    if (eq->dirty) _recompute(eq);

    int x = ctx->x, y = ctx->y, w = ctx->w, h = ctx->h;
    int cy = y + h / 2;

    /* Background */
    SDL_SetRenderDrawColor(r, 8, 8, 8, 255);
    SDL_Rect bg = {x, y, w, h};
    SDL_RenderFillRect(r, &bg);

    /* Grid lines — dB markers */
    for (int db = -9; db <= 9; db += 3) {
        int gy = gainToY((float)db, y, h);
        SDL_SetRenderDrawColor(r, db == 0 ? 50 : 22, db == 0 ? 50 : 22, db == 0 ? 50 : 22, 255);
        SDL_RenderDrawLine(r, x, gy, x+w, gy);
    }
    /* Grid lines — frequency markers */
    float freqMarks[] = {50,100,200,500,1000,2000,5000,10000};
    for (int f2 = 0; f2 < 8; f2++) {
        int gx = freqToX(freqMarks[f2], x, w);
        SDL_SetRenderDrawColor(r, 22, 22, 22, 255);
        SDL_RenderDrawLine(r, gx, y, gx, y+h);
    }

    /* EQ response curve */
    SDL_SetRenderDrawColor(r, 0, 255, 65, 200);
    int prevPy = -1;
    for (int px = 0; px < w; px++) {
        float t    = (float)px / (float)(w - 1);
        float freq = powf(10.0f, 1.301f + t * 2.699f); /* 20Hz–20kHz log */
        float db   = responseDb(eq, freq);
        db = db < -20.0f ? -20.0f : (db > 20.0f ? 20.0f : db);
        int py = gainToY(db, y, h);
        if (prevPy >= 0) SDL_RenderDrawLine(r, x+px-1, prevPy, x+px, py);
        prevPy = py;
    }

    /* Fill under curve (subtle) */
    for (int px = 0; px < w; px++) {
        float t    = (float)px / (float)(w - 1);
        float freq = powf(10.0f, 1.301f + t * 2.699f);
        float db   = responseDb(eq, freq);
        db = db < -20.0f ? -20.0f : (db > 20.0f ? 20.0f : db);
        int py = gainToY(db, y, h);
        SDL_SetRenderDrawColor(r, 0, 255, 65, 18);
        if (py < cy) SDL_RenderDrawLine(r, x+px, py, x+px, cy);
        else         SDL_RenderDrawLine(r, x+px, cy, x+px, py);
    }

    /* Band nodes */
    int nodeFreqParam[NUM_BANDS] = {0, 2, 5, 8};
    int nodeGainParam[NUM_BANDS] = {1, 3, 6, 9};

    for (int b = 0; b < NUM_BANDS; b++) {
        int nx = freqToX(eq->bandFreq[b], x, w);
        int ny = gainToY(eq->bandGain[b], y, h);

        uint8_t cr = NODE_COLORS[b][0];
        uint8_t cg = NODE_COLORS[b][1];
        uint8_t cb = NODE_COLORS[b][2];

        /* Outer ring (selection indicator) */
        int selected = (eq->selectedNode == b);
        SDL_SetRenderDrawColor(r, cr, cg, cb, selected ? 200 : 80);
        int rs = selected ? 12 : 10;
        for (int a = 0; a < 360; a += (selected ? 10 : 15)) {
            double rad = a * M_PI / 180.0;
            SDL_RenderDrawPoint(r, nx + (int)(cos(rad)*rs), ny + (int)(sin(rad)*rs));
        }

        /* Inner dot */
        SDL_SetRenderDrawColor(r, cr, cg, cb, 255);
        SDL_Rect dot = {nx-4, ny-4, 8, 8};
        SDL_RenderFillRect(r, &dot);

        /* Gain label */
        /* (Host renders param labels — we just draw the node) */
    }

    /* Mouse hover / drag (SDK v3.1 desktop) */
    if (ctx->mouseDown) {
        /* Find closest node to mouse on click */
        float mx = (float)ctx->mouseX;
        float my = (float)ctx->mouseY;
        float bestDist = 400.0f;
        int   bestNode = -1;
        for (int b = 0; b < NUM_BANDS; b++) {
            int nx = freqToX(eq->bandFreq[b], x, w);
            int ny = gainToY(eq->bandGain[b], y, h);
            float dx = mx - (float)nx;
            float dy = my - (float)ny;
            float dist = dx*dx + dy*dy;
            if (dist < bestDist) { bestDist = dist; bestNode = b; }
        }
        if (bestNode >= 0) {
            eq->selectedNode = bestNode;
            /* Map mouse X → frequency param */
            float t = (mx - x) / (float)w;
            t = t < 0.0f ? 0.0f : (t > 1.0f ? 1.0f : t);
            eq->p[nodeFreqParam[bestNode]] = t;
            /* Map mouse Y → gain param */
            float gt = 1.0f - (my - y) / (float)h;
            gt = gt < 0.0f ? 0.0f : (gt > 1.0f ? 1.0f : gt);
            eq->p[nodeGainParam[bestNode]] = gt;
            eq->dirty = 1;
        }
    }

    return 1;
}
