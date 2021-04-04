/*
* Bilateral filter - VapourSynth plugin
* Copyright (C) 2014  mawen1250
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cmath>
#include <cstring>
#include <algorithm>

#include "avisynth.h"

#ifdef _WIN32
#define VS_ALIGNED_MALLOC(pptr, size, alignment) do { *(pptr) = _aligned_malloc((size), (alignment)); } while (0)
#define VS_ALIGNED_FREE(ptr) do { _aligned_free((ptr)); } while (0)
#else
#define VS_ALIGNED_MALLOC(pptr, size, alignment) do { if(posix_memalign((void**)(pptr), (alignment), (size))) *((void**)pptr) = NULL; } while (0)
#define VS_ALIGNED_FREE(ptr) do { free((ptr)); } while (0)
#endif

template<typename T = void>
static inline T* vs_aligned_malloc(size_t size, size_t alignment)
{
#ifdef _WIN32
    return (T*)_aligned_malloc(size, alignment);
#else
    void* tmp = NULL;
    if (posix_memalign(&tmp, alignment, size))
        tmp = 0;
    return (T*)tmp;
#endif
}

static inline void vs_aligned_free(void* ptr)
{
    VS_ALIGNED_FREE(ptr);
}

static constexpr size_t Alignment = 32;

template < typename T >
static int stride_cal(int width)
{
    const size_t Alignment2 = Alignment / sizeof(T);
    return static_cast<int>(width % Alignment2 == 0 ? width : (width / Alignment2 + 1) * Alignment2);
}

template < typename T >
static void data2buff(T* dst, const T* src, int xoffset, int yoffset, int bufheight, int bufwidth, int bufstride, int height, int width, int stride)
{
    int x, y;
    T* dstp;
    const T* srcp;

    for (y = 0; y < height; ++y)
    {
        dstp = dst + (yoffset + static_cast<int64_t>(y)) * bufstride;
        srcp = src + static_cast<int64_t>(y) * stride;
        for (x = 0; x < xoffset; ++x)
            dstp[x] = srcp[0];
        memcpy(dstp + xoffset, srcp, sizeof(T) * width);
        for (x = xoffset + width; x < bufwidth; ++x)
            dstp[x] = srcp[width - 1];
    }

    srcp = dst + static_cast<int64_t>(yoffset) * bufstride;
    for (y = 0; y < yoffset; ++y)
    {
        dstp = dst + static_cast<int64_t>(y) * bufstride;
        memcpy(dstp, srcp, sizeof(T) * bufwidth);
    }

    srcp = dst + (static_cast<int64_t>(yoffset) + height - 1) * bufstride;
    for (y = yoffset + height; y < bufheight; ++y)
    {
        dstp = dst + static_cast<int64_t>(y) * bufstride;
        memcpy(dstp, srcp, sizeof(T) * bufwidth);
    }
}

template < typename T >
static T* newbuff(const T* src, int xoffset, int yoffset, int bufheight, int bufwidth, int bufstride, int height, int width, int stride)
{
    T* dst = vs_aligned_malloc<T>(sizeof(T) * bufheight * bufstride, Alignment);
    data2buff(dst, src, xoffset, yoffset, bufheight, bufwidth, bufstride, height, width, stride);
    return dst;
}

template < typename T >
static void freebuff(T* buff)
{
    vs_aligned_free(buff);
}

static constexpr double Pi = 3.1415926535897932384626433832795;
const double sqrt_2Pi = sqrt(2 * Pi);

static constexpr double sigmaSMul = 2.0;
static constexpr double sigmaRMul = 32.0;

static inline double Gaussian_Function_sqr_x(double sqr_x, double sigma)
{
    return exp(sqr_x / (sigma * sigma * -2));
}

static inline double Normalized_Gaussian_Function(double x, double sigma)
{
    x /= sigma;
    return exp(x * x / -2) / (sqrt_2Pi * sigma);
}

static inline double* Gaussian_Function_Spatial_LUT_Generation(const int xUpper, const int yUpper, const double sigmaS)
{
    int x, y;
    double* GS_LUT = new double[xUpper * static_cast<int64_t>(yUpper)];

    for (y = 0; y < yUpper; ++y)
    {
        for (x = 0; x < xUpper; ++x)
        {
            GS_LUT[y * xUpper + x] = static_cast<double>(Gaussian_Function_sqr_x(static_cast<double>(static_cast<int64_t>(x) * x + static_cast<int64_t>(y) * y), sigmaS));
        }
    }

    return GS_LUT;
}

static inline double Gaussian_Distribution2D_Spatial_LUT_Lookup(const double* GS_LUT, const int xUpper, const int x, const int y)
{
    return GS_LUT[y * xUpper + x];
}

static inline void Gaussian_Function_Spatial_LUT_Free(double* GS_LUT)
{
    delete[] GS_LUT;
}

static inline double* Gaussian_Function_Range_LUT_Generation(const int ValueRange, double sigmaR)
{
    int i;
    int Levels = ValueRange + 1;
    const int upper = std::min(ValueRange, static_cast<int>(sigmaR * sigmaRMul * ValueRange + 0.5));
    double* GR_LUT = new double[Levels];

    for (i = 0; i <= upper; ++i)
    {
        GR_LUT[i] = static_cast<double>(Normalized_Gaussian_Function(static_cast<double>(i) / ValueRange, sigmaR));
    }
    // For unknown reason, when more range weights are too small or equal 0, the runtime speed gets lower - mainly in function Recursive_Gaussian2D_Horizontal.
    // To avoid this issue, we set range weights whose range values are larger than sigmaR*sigmaRMul to the Gaussian function value at sigmaR*sigmaRMul.
    if (i < Levels)
    {
        const double upperLUTvalue = GR_LUT[upper];
        for (; i < Levels; ++i)
        {
            GR_LUT[i] = upperLUTvalue;
        }
    }

    return GR_LUT;
}

template < typename T >
static inline double Gaussian_Distribution2D_Range_LUT_Lookup(const double* GR_LUT, const T Value1, const T Value2)
{
    return GR_LUT[Value1 > Value2 ? Value1 - Value2 : Value2 - Value1];
}

static inline void Gaussian_Function_Range_LUT_Free(double* GR_LUT)
{
    delete[] GR_LUT;
}

static void Recursive_Gaussian_Parameters(const double sigma, double& B, double& B1, double& B2, double& B3)
{
    const double q = sigma < 2.5 ? 3.97156 - 4.14554 * sqrt(1 - 0.26891 * sigma) : 0.98711 * sigma - 0.96330;

    const double b0 = 1.57825 + 2.44413 * q + 1.4281 * q * q + 0.422205 * q * q * q;
    const double b1 = 2.44413 * q + 2.85619 * q * q + 1.26661 * q * q * q;
    const double b2 = -(1.4281 * q * q + 1.26661 * q * q * q);
    const double b3 = 0.422205 * q * q * q;

    B = static_cast<double>(1 - (b1 + b2 + b3) / b0);
    B1 = static_cast<double>(b1 / b0);
    B2 = static_cast<double>(b2 / b0);
    B3 = static_cast<double>(b3 / b0);
}

static void Recursive_Gaussian2D_Vertical(double* output, const double* input, int height, int width, int stride, const double B, const double B1, const double B2, const double B3)
{
    int i0, i1, i2, i3, j, lower, upper;
    double P0, P1, P2, P3;

    if (output != input)
    {
        memcpy(output, input, sizeof(double) * width);
    }

    for (j = 0; j < height; ++j)
    {
        lower = stride * j;
        upper = lower + width;

        i0 = lower;
        i1 = j < 1 ? i0 : i0 - stride;
        i2 = j < 2 ? i1 : i1 - stride;
        i3 = j < 3 ? i2 : i2 - stride;

        for (; i0 < upper; ++i0, ++i1, ++i2, ++i3)
        {
            P3 = output[i3];
            P2 = output[i2];
            P1 = output[i1];
            P0 = input[i0];
            output[i0] = B * P0 + B1 * P1 + B2 * P2 + B3 * P3;
        }
    }

    for (j = height - 1; j >= 0; --j)
    {
        lower = stride * j;
        upper = lower + width;

        i0 = lower;
        i1 = j >= height - 1 ? i0 : i0 + stride;
        i2 = j >= height - 2 ? i1 : i1 + stride;
        i3 = j >= height - 3 ? i2 : i2 + stride;

        for (; i0 < upper; ++i0, ++i1, ++i2, ++i3)
        {
            P3 = output[i3];
            P2 = output[i2];
            P1 = output[i1];
            P0 = output[i0];
            output[i0] = B * P0 + B1 * P1 + B2 * P2 + B3 * P3;
        }
    }
}

static void Recursive_Gaussian2D_Horizontal(double* output, const double* input, int height, int width, int stride, const double B, const double B1, const double B2, const double B3)
{
    int i, j, lower, upper;
    double P0, P1, P2, P3;

    for (j = 0; j < height; ++j)
    {
        lower = stride * j;
        upper = lower + width;

        i = lower;
        output[i] = P3 = P2 = P1 = input[i];

        for (++i; i < upper; ++i)
        {
            P0 = B * input[i] + B1 * P1 + B2 * P2 + B3 * P3;
            P3 = P2;
            P2 = P1;
            P1 = P0;
            output[i] = P0;
        }

        --i;
        P3 = P2 = P1 = output[i];

        for (--i; i >= lower; --i)
        {
            P0 = B * output[i] + B1 * P1 + B2 * P2 + B3 * P3;
            P3 = P2;
            P2 = P1;
            P1 = P0;
            output[i] = P0;
        }
    }
}

class Bilateral : public GenericVideoFilter
{
    PClip _ref;
    bool joint = false;
    double _sigmaS[3];
    double _sigmaR[3];
    int process[3];
    int _algorithm[3];
    bool has_at_least_v8;

    int _PBFICnum[3];

    int radius[3];
    int step[3];

    double* GS_LUT_[3];
    double* GR_LUT_[3];

    template<typename T>
    void Bilateral2D(PVideoFrame& dstf, PVideoFrame& srcf, PVideoFrame& reff, IScriptEnvironment* env);
    template<typename T>
    void Bilateral2D_1(T* dst, const T* src, const T* ref, int plane, int height, int width, int stride, int dst_stride);
    template<typename T>
    void Bilateral2D_2(T* dst, const T* src, const T* ref, int plane, int height, int width, int stride, int dst_stride);
    template<typename T>
    void Bilateral2D_2(T* dst, const T* src, int plane, int height, int width, int stride, int dst_stride);

public:
    Bilateral(PClip _child, PClip ref, double sigmaSY, double sigmaSU, double sigmaSV, double sigmaRY, double sigmaRU, double sigmaRV, int algorithmY, int algorithmU, int algorithmV, int PBFICnumY, int PBFICnumU, int PBFICnumV, int y, int u, int v, IScriptEnvironment* env);
    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
    int __stdcall SetCacheHints(int cachehints, int frame_range)
    {
        return cachehints == CACHE_GET_MTMODE ? MT_NICE_FILTER : 0;
    }
    ~Bilateral();
};

template<typename T>
void Bilateral::Bilateral2D(PVideoFrame& dstf, PVideoFrame& srcf, PVideoFrame& reff, IScriptEnvironment* env)
{
    const int planes_y[3] = { PLANAR_Y, PLANAR_U, PLANAR_V };
    const int planes_r[3] = { PLANAR_R, PLANAR_G, PLANAR_B };
    const int* planes = (vi.IsRGB()) ? planes_r : planes_y;
    const int planecount = std::min(vi.NumComponents(), 3);

    for (int pl = 0; pl < planecount; ++pl)
    {
        const int plane = planes[pl];
        const int height = srcf->GetHeight(plane);

        if (process[pl] == 3)
        {
            const int stride = srcf->GetPitch(plane) / sizeof(T);
            const int dst_stride = dstf->GetPitch(plane) / sizeof(T);
            const int width = srcf->GetRowSize(plane) / sizeof(T);
            const T* src = reinterpret_cast<const T*>(srcf->GetReadPtr(plane));
            const T* ref = reinterpret_cast<const T*>(reff->GetReadPtr(plane));
            T* dst = reinterpret_cast<T*>(dstf->GetWritePtr(plane));

            switch (_algorithm[pl])
            {
                case 1:
                    Bilateral2D_1(dst, src, ref, pl, height, width, stride, dst_stride);
                    break;
                case 2:
                    if (joint)
                        Bilateral2D_2(dst, src, ref, pl, height, width, stride, dst_stride);
                    else
                        Bilateral2D_2(dst, src, pl, height, width, stride, dst_stride);
                    break;
            }
        }
        else if (process[pl] == 2)
            env->BitBlt(dstf->GetWritePtr(plane), dstf->GetPitch(plane), srcf->GetReadPtr(plane), srcf->GetPitch(plane), srcf->GetRowSize(plane), height);
    }
}

// Implementation of O(1) cross/joint Bilateral filter algorithm from "Qingxiong Yang, Kar-Han Tan, Narendra Ahuja - Real-Time O(1) Bilateral Filtering"
template<typename T>
void Bilateral::Bilateral2D_1(T* dst, const T* src, const T* ref, int plane, int height, int width, int stride, int dst_stride)
{
    const int pcount = stride * height;
    const int PBFICnum = _PBFICnum[plane];

    const long long peak = (1 << vi.BitsPerComponent()) - 1;

    const double* GR_LUT = GR_LUT_[plane];

    // Generate quantized PBFICs' parameters
    T* PBFICk = new T[PBFICnum];

    for (int k = 0; k < PBFICnum; ++k)
    {
        PBFICk[k] = static_cast<T>(llrint(static_cast<double>(peak) * k / (PBFICnum - 1)));
    }

    // Generate recursive Gaussian parameters
    double B, B1, B2, B3;
    Recursive_Gaussian_Parameters(_sigmaS[plane], B, B1, B2, B3);

    // Generate quantized PBFICs
    double** PBFIC = new double * [PBFICnum];
    double* Wk = vs_aligned_malloc<double>(sizeof(double) * pcount, Alignment);
    double* Jk = vs_aligned_malloc<double>(sizeof(double) * pcount, Alignment);
    int k;

    for (k = 0; k < PBFICnum; ++k)
    {
        PBFIC[k] = vs_aligned_malloc<double>(sizeof(double) * pcount, Alignment);

        for (int j = 0; j < height; ++j)
        {
            int i = stride * j;

            for (int upper = i + width; i < upper; ++i)
            {
                Wk[i] = Gaussian_Distribution2D_Range_LUT_Lookup(GR_LUT, PBFICk[k], ref[i]);
                Jk[i] = Wk[i] * static_cast<double>(src[i]);
            }
        }

        Recursive_Gaussian2D_Horizontal(Wk, Wk, height, width, stride, B, B1, B2, B3);
        Recursive_Gaussian2D_Vertical(Wk, Wk, height, width, stride, B, B1, B2, B3);
        Recursive_Gaussian2D_Horizontal(Jk, Jk, height, width, stride, B, B1, B2, B3);
        Recursive_Gaussian2D_Vertical(Jk, Jk, height, width, stride, B, B1, B2, B3);

        for (int j = 0; j < height; ++j)
        {
            int i = stride * j;

            for (int upper = i + width; i < upper; ++i)
            {
                PBFIC[k][i] = Wk[i] == 0 ? 0 : Jk[i] / Wk[i];
            }
        }
    }

    // Generate filtered result from PBFICs using linear interpolation
    for (int j = 0; j < height; ++j)
    {
        int i = dst_stride * j;

        for (int upper = i + width; i < upper; ++i)
        {
            for (k = 0; k < PBFICnum - 2; ++k)
            {
                if (ref[i] < PBFICk[k + 1] && ref[i] >= PBFICk[k]) break;
            }

            dst[i] = static_cast<T>(std::clamp(llrint(((PBFICk[k + 1] - ref[i]) * PBFIC[k][i] + (ref[i] - PBFICk[k]) * PBFIC[k + 1][i]) / (PBFICk[k + 1] - PBFICk[k])), 0LL, peak));
        }
    }

    // Clear
    for (int k = 0; k < PBFICnum; ++k)
        vs_aligned_free(PBFIC[k]);

    vs_aligned_free(Jk);
    vs_aligned_free(Wk);
    delete[] PBFIC;
    delete[] PBFICk;
}

// Implementation of cross/joint Bilateral filter with truncated spatial window and sub-sampling
template<typename T>
void Bilateral::Bilateral2D_2(T* dst, const T* src, const T* ref, int plane, int height, int width, int stride, int dst_stride)
{
    const int radiusx = radius[plane];
    const int radiusy = radius[plane];
    const int buffnum = radiusy * 2 + 1;
    const int samplestep = step[plane];
    const int samplecenter = buffnum / 2;

    const int bufheight = height + radiusy * 2;
    const int bufwidth = width + radiusx * 2;
    const int bufstride = stride_cal<T>(bufwidth);

    const long long peak = (1 << vi.BitsPerComponent()) - 1;

    const T* srcp = src;
    const T* refp = ref;
    T* dstp = dst;

    const double* GS_LUT = GS_LUT_[plane];
    const double* GR_LUT = GR_LUT_[plane];

    // Allocate buffs
    T* srcbuff = newbuff(srcp, radiusx, radiusy, bufheight, bufwidth, bufstride, height, width, stride);
    T* refbuff = newbuff(refp, radiusx, radiusy, bufheight, bufwidth, bufstride, height, width, stride);
    T* srcbuffp1, * refbuffp1, * srcbuffp2, * refbuffp2;

    // Process
    double SWei, RWei1, RWei2, RWei3, RWei4, WeightSum, Sum;
    const int xUpper = radiusx + 1, yUpper = radiusy + 1;

    const int yoffset = samplecenter;

    for (int j = 0; j < height; ++j, srcp += stride, refp += stride, dstp += dst_stride)
    {
        srcbuffp1 = srcbuff + (yoffset + static_cast<int64_t>(j)) * bufstride;
        refbuffp1 = refbuff + (yoffset + static_cast<int64_t>(j)) * bufstride;

        for (int i = 0; i < width; ++i)
        {
            const int xoffset = samplecenter + i;
            srcbuffp2 = srcbuffp1 + xoffset;
            refbuffp2 = refbuffp1 + xoffset;

            WeightSum = GS_LUT[0] * GR_LUT[0];
            Sum = srcp[i] * WeightSum;

            for (int y = 1; y < yUpper; y += samplestep)
            {
                for (int x = 1; x < xUpper; x += samplestep)
                {
                    SWei = Gaussian_Distribution2D_Spatial_LUT_Lookup(GS_LUT, xUpper, x, y);
                    RWei1 = Gaussian_Distribution2D_Range_LUT_Lookup(GR_LUT, refp[i], refbuffp2[+y * bufstride + x]);
                    RWei2 = Gaussian_Distribution2D_Range_LUT_Lookup(GR_LUT, refp[i], refbuffp2[+y * bufstride - x]);
                    RWei3 = Gaussian_Distribution2D_Range_LUT_Lookup(GR_LUT, refp[i], refbuffp2[-y * bufstride - x]);
                    RWei4 = Gaussian_Distribution2D_Range_LUT_Lookup(GR_LUT, refp[i], refbuffp2[-y * bufstride + x]);

                    WeightSum += SWei * (RWei1 + RWei2 + RWei3 + RWei4);
                    Sum += SWei * (
                        srcbuffp2[+y * bufstride + x] * RWei1 +
                        srcbuffp2[+y * bufstride - x] * RWei2 +
                        srcbuffp2[-y * bufstride - x] * RWei3 +
                        srcbuffp2[-y * bufstride + x] * RWei4);
                }
            }

            dstp[i] = static_cast<T>(std::clamp(llrint(Sum / WeightSum), 0LL, peak));
        }
    }

    // Clear and output
    freebuff(srcbuff);
    freebuff(refbuff);
}

// Implementation of Bilateral filter with truncated spatial window and sub-sampling
template<typename T>
void Bilateral::Bilateral2D_2(T* dst, const T* src, int plane, int height, int width, int stride, int dst_stride)
{
    const int radiusx = radius[plane];
    const int radiusy = radius[plane];
    const int buffnum = radiusy * 2 + 1;
    const int samplestep = step[plane];
    const int samplecenter = buffnum / 2;

    const int bufheight = height + radiusy * 2;
    const int bufwidth = width + radiusx * 2;
    const int bufstride = stride_cal<T>(bufwidth);

    const long long peak = (1 << vi.BitsPerComponent()) - 1;

    const T* srcp = src;
    T* dstp = dst;

    const double* GS_LUT = GS_LUT_[plane];
    const double* GR_LUT = GR_LUT_[plane];

    // Allocate buffs
    T* srcbuff = newbuff(srcp, radiusx, radiusy, bufheight, bufwidth, bufstride, height, width, stride);
    T* srcbuffp1, * srcbuffp2;

    // Process
    double SWei, RWei1, RWei2, RWei3, RWei4, WeightSum, Sum;
    const int xUpper = radiusx + 1, yUpper = radiusy + 1;

    const int yoffset = samplecenter;

    for (int j = 0; j < height; ++j, srcp += stride, dstp += dst_stride)
    {
        srcbuffp1 = srcbuff + (yoffset + static_cast<int64_t>(j)) * bufstride;

        for (int i = 0; i < width; ++i)
        {
            const int xoffset = samplecenter + i;
            srcbuffp2 = srcbuffp1 + xoffset;

            WeightSum = GS_LUT[0] * GR_LUT[0];
            Sum = srcp[i] * WeightSum;

            for (int y = 1; y < yUpper; y += samplestep)
            {
                for (int x = 1; x < xUpper; x += samplestep)
                {
                    SWei = Gaussian_Distribution2D_Spatial_LUT_Lookup(GS_LUT, xUpper, x, y);
                    RWei1 = Gaussian_Distribution2D_Range_LUT_Lookup(GR_LUT, srcp[i], srcbuffp2[+y * bufstride + x]);
                    RWei2 = Gaussian_Distribution2D_Range_LUT_Lookup(GR_LUT, srcp[i], srcbuffp2[+y * bufstride - x]);
                    RWei3 = Gaussian_Distribution2D_Range_LUT_Lookup(GR_LUT, srcp[i], srcbuffp2[-y * bufstride - x]);
                    RWei4 = Gaussian_Distribution2D_Range_LUT_Lookup(GR_LUT, srcp[i], srcbuffp2[-y * bufstride + x]);

                    WeightSum += SWei * (RWei1 + RWei2 + RWei3 + RWei4);
                    Sum += SWei * (
                        srcbuffp2[+y * bufstride + x] * RWei1 +
                        srcbuffp2[+y * bufstride - x] * RWei2 +
                        srcbuffp2[-y * bufstride - x] * RWei3 +
                        srcbuffp2[-y * bufstride + x] * RWei4);
                }
            }

            dstp[i] = static_cast<T>(std::clamp(llrint(Sum / WeightSum), 0LL, peak));
        }
    }

    // Clear and output
    freebuff(srcbuff);
}

Bilateral::Bilateral(PClip _child, PClip ref, double sigmaSY, double sigmaSU, double sigmaSV, double sigmaRY, double sigmaRU, double sigmaRV, int algorithmY, int algorithmU, int algorithmV, int PBFICnumY, int PBFICnumU, int PBFICnumV, int y, int u, int v, IScriptEnvironment* env)
    : GenericVideoFilter(_child), _ref(ref)
{
    _sigmaS[0] = sigmaSY;
    _sigmaS[1] = sigmaSU;
    _sigmaS[2] = sigmaSV;
    _sigmaR[0] = sigmaRY;
    _sigmaR[1] = sigmaRU;
    _sigmaR[2] = sigmaRV;
    _algorithm[0] = algorithmY;
    _algorithm[1] = algorithmU;
    _algorithm[2] = algorithmV;
    _PBFICnum[0] = PBFICnumY;
    _PBFICnum[1] = PBFICnumU;
    _PBFICnum[2] = PBFICnumV;

    joint = (!_ref) ? false : true;

    const int process_planes[3] = { y, u, v };
    for (int i = 0; i < 3; ++i)
    {
        if (vi.IsRGB())
            process[i] = 3;
        else
        {
            switch (process_planes[i])
            {
                case 3: process[i] = 3; break;
                case 2: process[i] = 2; break;
                default: process[i] = 1; break;
            }
        }
    }

    int orad[3];
    int samples[3];

    for (int i = 0; i < 3; ++i)
    {
        GS_LUT_[i] = nullptr;
        GR_LUT_[i] = nullptr;
        step[i] = 0;
        radius[i] = 0;

        if (process[i] == 3)
        {
            if (_PBFICnum[i] == 0)
            {
                if (_sigmaR[i] >= 0.08)
                    _PBFICnum[i] = 4;
                else if (_sigmaR[i] >= 0.015)
                    _PBFICnum[i] = std::min(16, static_cast<int>(4 * 0.08 / _sigmaR[i] + 0.5));
                else
                    _PBFICnum[i] = std::min(32, static_cast<int>(16 * 0.015 / _sigmaR[i] + 0.5));

                if (i > 0 && !vi.IsRGB() && _PBFICnum[i] % 2 == 0 && _PBFICnum[i] < 256) // Set odd PBFIC number to chroma planes by default
                    _PBFICnum[i]++;
            }

            orad[i] = std::max(static_cast<int>(_sigmaS[i] * sigmaSMul + 0.5), 1);

            step[i] = orad[i] < 4 ? 1 : orad[i] < 8 ? 2 : 3;

            samples[i] = 1;
            radius[i] = 1 + (samples[i] - 1) * step[i];

            while (orad[i] * 2 > radius[i] * 3)
            {
                samples[i]++;
                radius[i] = 1 + (samples[i] - 1) * step[i];
                if (radius[i] >= orad[i] && samples[i] > 2)
                {
                    samples[i]--;
                    radius[i] = 1 + (samples[i] - 1) * step[i];
                    break;
                }
            }

            if (_algorithm[i] <= 0)
                _algorithm[i] = step[i] == 1 ? 2 : _sigmaR[i] < 0.08 && samples[i] < 5 ? 2
                : 4 * samples[i] * samples[i] <= 15 * _PBFICnum[i] ? 2 : 1;

            if (_algorithm[i] == 2)
                GS_LUT_[i] = Gaussian_Function_Spatial_LUT_Generation(radius[i] + 1, radius[i] + 1, _sigmaS[i]);

            GR_LUT_[i] = Gaussian_Function_Range_LUT_Generation((1 << vi.BitsPerComponent()) - 1, _sigmaR[i]);
        }
    }

    has_at_least_v8 = true;
    try { env->CheckVersion(8); }
    catch (const AvisynthError&) { has_at_least_v8 = false; };
}    

Bilateral::~Bilateral()
{
    for (int i = 0; i < 3; ++i)
    {
        if (process[i] == 3)
        {
            Gaussian_Function_Spatial_LUT_Free(GS_LUT_[i]);
            Gaussian_Function_Range_LUT_Free(GR_LUT_[i]);
        }
    }
}

PVideoFrame __stdcall Bilateral::GetFrame(int n, IScriptEnvironment* env)
{
    PVideoFrame src = child->GetFrame(n, env);
    PVideoFrame dst = (has_at_least_v8) ? env->NewVideoFrameP(vi, &src) : env->NewVideoFrame(vi);
    PVideoFrame ref = (joint) ? _ref->GetFrame(n, env) : src;

    if (vi.ComponentSize() == 1)
        Bilateral2D<uint8_t>(dst, src, ref, env);
    else
        Bilateral2D<uint16_t>(dst, src, ref, env);

    return dst;
}

AVSValue __cdecl Create_Bilateral(AVSValue args, void* user_data, IScriptEnvironment* env)
{
    PClip clip = args[0].AsClip();
    const VideoInfo& vi = clip->GetVideoInfo();

    if (!vi.IsPlanar())
        env->ThrowError("Bilateral: clip must be in planar format.");
    if (vi.ComponentSize() == 4)
        env->ThrowError("Bilateral: only 8..16-bit formats supported.");

    PClip ref = (args[1].Defined()) ? args[1].AsClip() : nullptr;
    if (ref)
    {
        PClip ref = args[1].AsClip();
        const VideoInfo& viref = ref->GetVideoInfo();

        if (!viref.IsPlanar())
            env->ThrowError("Bilateral: reference clip must be in planar format.");
        if (viref.ComponentSize() == 4)
            env->ThrowError("Bilateral: reference clip must be 8..16-bit format.");
        if (vi.width != viref.width || vi.height != viref.height)
            env->ThrowError("Bilateral: input clip and reference clip must be of the same size.");
        if (!vi.IsSameColorspace(viref))
            env->ThrowError("Bilateral: input clip and reference clip must be of the same color family.");
    }

    const float sigmaSY = args[2].AsFloatf(3.0f);
    const float sigmaSU = [&]()
    {
        if (args[3].Defined())
            return args[3].AsFloatf();
        else
        {
            if (vi.NumComponents() > 1 && !vi.IsRGB())
                return sigmaSY / std::sqrt(static_cast<float>(1 << vi.GetPlaneHeightSubsampling(PLANAR_U)) * (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)));
            else
                return sigmaSY;
        }
    }();
    const float sigmaSV = args[4].AsFloatf(sigmaSU);
    const float sigmaRY = args[5].AsFloatf(0.02f);
    const float sigmaRU = args[6].AsFloatf(sigmaRY);
    const float sigmaRV = args[7].AsFloatf(sigmaRU);

    if (sigmaSY < 0 || sigmaSU < 0 || sigmaSV < 0)
        env->ThrowError("Bilateral: sigma must be non-negative value.");
    if (sigmaRY < 0 || sigmaRU < 0 || sigmaRV < 0)
        env->ThrowError("Bilateral: sigmaV must be non-negative value.");
    
    const int algorithmY = args[8].AsInt(0);
    const int algorithmU = args[9].AsInt(algorithmY);
    const int algorithmV = args[10].AsInt(algorithmU);

    if (algorithmY < 0 || algorithmU < 0 || algorithmV < 0)
        env->ThrowError("Bilateral: algorithm must be non-negative value.");

    const int PBFICnumY = args[11].AsInt(0);
    const int PBFICnumU = args[12].AsInt(PBFICnumY);
    const int PBFICnumV = args[13].AsInt(PBFICnumU);

    if (PBFICnumY < 0 || PBFICnumU < 0 || PBFICnumV < 0
        || PBFICnumY == 1 || PBFICnumU == 1 || PBFICnumV == 1
        || PBFICnumY > 256 || PBFICnumU > 256 || PBFICnumV > 256)
        env->ThrowError("Bilateral: PBFICnum must be integer ranges in [0,256] except 1.");

    const int y = args[14].AsInt(3);
    const int u = args[15].AsInt(1);
    const int v = args[16].AsInt(1);

    if (y < 1 || y > 3)
        env->ThrowError("Bilateral: y must be between 1..3.");
    if (u < 1 || u > 3)
        env->ThrowError("Bilateral: u must be between 1..3.");
    if (v < 1 || v > 3)
        env->ThrowError("Bilateral: v must be between 1..3.");

    return new Bilateral(
        clip,
        ref,
        sigmaSY,
        sigmaSU,
        sigmaSV,
        sigmaRY,
        sigmaRU,
        sigmaRV,
        algorithmY,
        algorithmU,
        algorithmV,
        PBFICnumY,
        PBFICnumU,
        PBFICnumV,
        y,
        u,
        v,
        env);
}

const AVS_Linkage* AVS_linkage;

extern "C" __declspec(dllexport)
const char* __stdcall AvisynthPluginInit3(IScriptEnvironment * env, const AVS_Linkage* const vectors)
{
    AVS_linkage = vectors;

    env->AddFunction("Bilateral", "c[ref]c[sigmaSY]f[sigmaSU]f[sigmaSV]f[sigmaRY]f[sigmaRU]f[sigmaRV]f[algorithmY]i[algorithmU]i[algorithmV]i[PBFICnumY]i[PBFICnumU]i[PBFICnumV]i[y]i[u]i[v]i", Create_Bilateral, 0);

    return "Bilateral";
}
