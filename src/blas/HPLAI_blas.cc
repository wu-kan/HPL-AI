/*
 * MIT License
 * 
 * Copyright (c) 2021 WuK
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include "hplai.hh"

#if !defined(HPLAI_NO_HPL_BLASPP)

static enum HPL_ORDER blaspp2Hpl(blas::Layout layout)
{
    return layout == blas::Layout::RowMajor ? HplRowMajor : HplColumnMajor;
}

static enum HPL_SIDE blaspp2Hpl(blas::Side side)
{
    return side == blas::Side::Left ? HplLeft : HplRight;
}

static enum HPL_UPLO blaspp2Hpl(blas::Uplo uplo)
{
    return uplo == blas::Uplo::Upper ? HplUpper : HplLower;
}

static enum HPL_TRANS blaspp2Hpl(blas::Op trans)
{
    return trans == blas::Op::NoTrans ? HplNoTrans
           : trans == blas::Op::Trans ? HplTrans
                                      : HplConjTrans;
}

static enum HPL_DIAG blaspp2Hpl(blas::Diag diag)
{
    return diag == blas::Diag::Unit ? HplUnit : HplNonUnit;
}

template <>
int64_t blas::iamax<double>(
    int64_t n,
    double const *x,
    int64_t incx)
{
    return blas::iamax(n, x, incx);
}

template <>
void blas::axpy<double, double>(
    int64_t n,
    blas::scalar_type<double, double> alpha,
    double const *x,
    int64_t incx,
    double *y,
    int64_t incy)
{
    HPL_daxpy(n, alpha, x, incx, y, incy);
}

template <>
void blas::copy<double, double>(
    int64_t n,
    double const *x,
    int64_t incx,
    double *y,
    int64_t incy)
{
    HPL_dcopy(n, x, incx, y, incy);
}

template <>
void blas::gemm<double, double, double>(
    blas::Layout layout,
    blas::Op transA,
    blas::Op transB,
    int64_t m,
    int64_t n,
    int64_t k,
    blas::scalar_type<double, double, double> alpha,
    double const *A,
    int64_t lda,
    double const *B,
    int64_t ldb,
    blas::scalar_type<double, double, double> beta,
    double *C,
    int64_t ldc)
{
    HPL_dgemm(
        blaspp2Hpl(layout),
        blaspp2Hpl(transA),
        blaspp2Hpl(transB),
        m,
        n,
        k,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc);
}

template <>
void blas::gemv<double, double, double>(
    blas::Layout layout,
    blas::Op trans,
    int64_t m,
    int64_t n,
    blas::scalar_type<double, double, double> alpha,
    double const *A,
    int64_t lda,
    double const *x,
    int64_t incx,
    blas::scalar_type<double, double, double> beta,
    double *y,
    int64_t incy)
{
    HPL_dgemv(
        blaspp2Hpl(layout),
        blaspp2Hpl(trans),
        m,
        n,
        alpha,
        A,
        lda,
        x,
        incx,
        beta,
        y,
        incy);
}

template <>
void blas::ger<double, double, double>(
    blas::Layout layout,
    int64_t m,
    int64_t n,
    blas::scalar_type<double, double, double> alpha,
    double const *x,
    int64_t incx,
    double const *y,
    int64_t incy,
    double *A,
    int64_t lda)
{
    HPL_dger(
        blaspp2Hpl(layout),
        m,
        n,
        alpha,
        x,
        incx,
        const_cast<double *>(y),
        incy,
        A,
        lda);
}

template <>
void blas::scal<double>(
    int64_t n,
    double alpha,
    double *x,
    int64_t incx)
{
    HPL_dscal(n, alpha, x, incx);
}

template <>
void blas::trsm<double, double>(
    blas::Layout layout,
    blas::Side side,
    blas::Uplo uplo,
    blas::Op trans,
    blas::Diag diag,
    int64_t m,
    int64_t n,
    blas::scalar_type<double, double> alpha,
    double const *A,
    int64_t lda,
    double *B,
    int64_t ldb)
{
    HPL_dtrsm(
        blaspp2Hpl(layout),
        blaspp2Hpl(side),
        blaspp2Hpl(uplo),
        blaspp2Hpl(trans),
        blaspp2Hpl(diag),
        m,
        n,
        alpha,
        A,
        lda,
        B,
        ldb);
}

template <>
void blas::trsv<double, double>(
    blas::Layout layout,
    blas::Uplo uplo,
    blas::Op trans,
    blas::Diag diag,
    int64_t n,
    double const *A,
    int64_t lda,
    double *x,
    int64_t incx)
{
    HPL_dtrsv(
        blaspp2Hpl(layout),
        blaspp2Hpl(uplo),
        blaspp2Hpl(trans),
        blaspp2Hpl(diag),
        n,
        A,
        lda,
        x,
        incx);
}

#endif

#if defined(HPLAI_DEVICE_BLASPP_GEMM) || defined(HPLAI_DEVICE_BLASPP_TRSM)

static blas::Queue *HPLAI_DEVICE_BLASPP_QUEUE = NULL;
static int64_t HPLAI_DEVICE_BLASPP_BUFFER_SIZE = 0;
static HPLAI_T_AFLOAT *HPLAI_DEVICE_BLASPP_BUFFER = NULL;

static void HPLAI_DEVICE_BLASPP_BUFFER_RESIZE(int64_t NEW_SIZE)
{
    if (HPLAI_DEVICE_BLASPP_BUFFER != NULL)
    {
        blas::device_free(HPLAI_DEVICE_BLASPP_BUFFER);
        HPLAI_DEVICE_BLASPP_BUFFER = NULL;
    }
    HPLAI_DEVICE_BLASPP_BUFFER_SIZE = NEW_SIZE;
    if (HPLAI_DEVICE_BLASPP_BUFFER_SIZE > 0)
    {
        HPLAI_DEVICE_BLASPP_BUFFER = blas::device_malloc<HPLAI_T_AFLOAT>(
            HPLAI_DEVICE_BLASPP_BUFFER_SIZE);
        if (HPLAI_DEVICE_BLASPP_BUFFER == NULL)
            HPLAI_pabort(
                __LINE__,
                "HPLAI_DEVICE_BLASPP_BUFFER_RESIZE",
                "Memory allocation failed for HPLAI_DEVICE_BLASPP_BUFFER");
    }
}

#endif

template <>
int64_t blas::iamax<float>(
    int64_t n,
    float const *x,
    int64_t incx)
{
    return blas::iamax(n, x, incx);
}

template <>
void blas::axpy<float, float>(
    int64_t n,
    blas::scalar_type<float, float> alpha,
    float const *x,
    int64_t incx,
    float *y,
    int64_t incy)
{
    blas::axpy(n, alpha, x, incx, y, incy);
}

template <>
void blas::copy<float, float>(
    int64_t n,
    float const *x,
    int64_t incx,
    float *y,
    int64_t incy)
{
    blas::copy(n, x, incx, y, incy);
}

#if defined(HPLAI_DEVICE_BLASPP_GEMM)

#if defined(HPLAI_CUBLASGEMMEX_COMPUTETYPE)

static cublasOperation_t HPLAI_op2cublas(blas::Op trans)
{
    return trans == blas::Op::NoTrans ? CUBLAS_OP_N
           : trans == blas::Op::Trans ? CUBLAS_OP_T
                                      : CUBLAS_OP_C;
}

static cudaDataType_t HPLAI_GET_cudaDataType_t(float t)
{
    return CUDA_R_32F;
}

static cudaDataType_t HPLAI_GET_cudaDataType_t(double t)
{
    return CUDA_R_64F;
}

#endif

template <>
void blas::gemm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
    blas::Layout layout,
    blas::Op TRANSA,
    blas::Op TRANSB,
    int64_t M,
    int64_t N,
    int64_t K,
    blas::scalar_type<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT> ALPHA,
    HPLAI_T_AFLOAT const *A,
    int64_t LDA,
    HPLAI_T_AFLOAT const *B,
    int64_t LDB,
    blas::scalar_type<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT> BETA,
    HPLAI_T_AFLOAT *C,
    int64_t LDC)
{
    //HPLAI_pabort( __LINE__, "blas::gemm", "Use HPLAI_DEVICE_BLASPP_GEMM" );
    if (layout != blas::Layout::ColMajor)
    {
        blas::gemm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
            blas::Layout::ColMajor,
            TRANSB, TRANSA, N, M, K, ALPHA, B, LDB, A, LDA, BETA, C, LDC);
        return;
    }
    int64_t i, j;

    if ((M == 0) || (N == 0) ||
        (((ALPHA == HPLAI_rzero) || (K == 0)) &&
         (BETA == HPLAI_rone)))
        return;

    if (ALPHA == HPLAI_rzero && BETA == HPLAI_rzero)
    {
        for (j = 0; j < N; j++)
        {
            for (i = 0; i < M; i++)
                *(C + i + j * LDC) = HPLAI_rzero;
        }
        return;
    }

#if !defined(HPLAI_DEVICE_BLASPP_GEMM_USE_CPU)
#define HPLAI_DEVICE_BLASPP_GEMM_USE_CPU 64
#endif
    const int64_t P = HPLAI_DEVICE_BLASPP_GEMM_USE_CPU;
    if (M < P || N < P || K < P)
    {
        blas::gemm(
            layout,
            TRANSA,
            TRANSB,
            M,
            N,
            K,
            ALPHA,
            A,
            LDA,
            B,
            LDB,
            BETA,
            C,
            LDC);
        return;
    }
    const int64_t
        M2 = M % P,
        N2 = N % P,
        K2 = K % P,
        M1 = M - M2,
        N1 = N - N2,
        K1 = K - K2;

    blas::Op TRANSC = blas::Op::NoTrans;

    int64_t rC = TRANSC == blas::Op::NoTrans ? M1 : N1;
    int64_t cC = TRANSC == blas::Op::NoTrans ? N1 : M1;
    int64_t dLDC = rC;
    int64_t dsC = cC * dLDC;

    int64_t rB = TRANSB == blas::Op::NoTrans ? K1 : N1;
    int64_t cB = TRANSB == blas::Op::NoTrans ? N1 : K1;
    int64_t dLDB = rB;
    int64_t dsB = cB * dLDB;

    int64_t rA = TRANSA == blas::Op::NoTrans ? M1 : K1;
    int64_t cA = TRANSA == blas::Op::NoTrans ? K1 : M1;
    int64_t dLDA = rA;
    int64_t dsA = cA * dLDA;

    if (HPLAI_DEVICE_BLASPP_BUFFER_SIZE < dsC + dsB + dsA)
        HPLAI_DEVICE_BLASPP_BUFFER_RESIZE(dsC + dsB + dsA);

    HPLAI_T_AFLOAT *dC = HPLAI_DEVICE_BLASPP_BUFFER;
    HPLAI_T_AFLOAT *dB = dC + dsC;
    HPLAI_T_AFLOAT *dA = dB + dsB;
#if defined(HPLAI_DEVICE_BLASPP_GEMM_MULTISTREAM)
    HPLAI_DEVICE_BLASPP_QUEUE->fork();
#endif
    blas::device_setmatrix<HPLAI_T_AFLOAT>(rA, cA, A, LDA, dA, dLDA, *HPLAI_DEVICE_BLASPP_QUEUE);
#if defined(HPLAI_DEVICE_BLASPP_GEMM_MULTISTREAM)
    HPLAI_DEVICE_BLASPP_QUEUE->revolve();
#endif
    blas::device_setmatrix<HPLAI_T_AFLOAT>(rB, cB, B, LDB, dB, dLDB, *HPLAI_DEVICE_BLASPP_QUEUE);
#if defined(HPLAI_DEVICE_BLASPP_GEMM_MULTISTREAM)
    HPLAI_DEVICE_BLASPP_QUEUE->revolve();
#endif
    blas::gemm(
        layout,
        TRANSA,
        TRANSB,
        M1,
        N1,
        K2,
        ALPHA,
        TRANSA == blas::Op::NoTrans ? A + K1 * LDA : A + K1,
        LDA,
        TRANSB == blas::Op::NoTrans ? B + K1 : B + K1 * LDB,
        LDB,
        BETA,
        C,
        LDC);
    blas::device_setmatrix<HPLAI_T_AFLOAT>(rC, cC, C, LDC, dC, dLDC, *HPLAI_DEVICE_BLASPP_QUEUE);
#if defined(HPLAI_DEVICE_BLASPP_GEMM_MULTISTREAM)
    HPLAI_DEVICE_BLASPP_QUEUE->join();
#endif

#if defined(HPLAI_CUBLASGEMMEX_COMPUTETYPE)
    HPLAI_T_AFLOAT rone = HPLAI_rone;
    cublasGemmEx(
        HPLAI_DEVICE_BLASPP_QUEUE->handle(),
        HPLAI_op2cublas(TRANSA),
        HPLAI_op2cublas(TRANSB),
        M1,
        N1,
        K1,
        &ALPHA,
        dA,
        HPLAI_GET_cudaDataType_t(HPLAI_rzero),
        dLDA,
        dB,
        HPLAI_GET_cudaDataType_t(HPLAI_rzero),
        dLDB,
        &rone,
        dC,
        HPLAI_GET_cudaDataType_t(HPLAI_rzero),
        dLDC,
        HPLAI_CUBLASGEMMEX_COMPUTETYPE,
        CUBLAS_GEMM_DEFAULT);
#else
    blas::gemm(
        layout,
        TRANSA,
        TRANSB,
        M1,
        N1,
        K1,
        ALPHA,
        dA,
        dLDA,
        dB,
        dLDB,
        HPLAI_rone,
        dC,
        dLDC,
        *HPLAI_DEVICE_BLASPP_QUEUE);
#endif

    blas::device_getmatrix<HPLAI_T_AFLOAT>(rC, cC, dC, dLDC, C, LDC, *HPLAI_DEVICE_BLASPP_QUEUE);

    blas::gemm(
        layout,
        TRANSA,
        TRANSB,
        M,
        N2,
        K,
        ALPHA,
        A,
        LDA,
        TRANSB == blas::Op::NoTrans ? B + N1 * LDB : B + N1,
        LDB,
        BETA,
        C + N1 * LDC,
        LDC);

    blas::gemm(
        layout,
        TRANSA,
        TRANSB,
        M2,
        N1,
        K,
        ALPHA,
        TRANSA == blas::Op::NoTrans ? A + M1 : A + M1 * LDA,
        LDA,
        B,
        LDB,
        BETA,
        C + M1,
        LDC);

    HPLAI_DEVICE_BLASPP_QUEUE->sync();
}

#elif defined(HPLAI_ACL_BLASPP_GEMM)

#include <acl/acl.h>

#define ACLCHECK(cmd)              \
    do                             \
    {                              \
        int e = cmd;               \
        if (e != ACL_ERROR_NONE)   \
        {                          \
            char str[9];           \
            sprintf(str, "%d", e); \
            HPLAI_pabort(          \
                __LINE__,          \
                "aclrt",           \
                str);              \
        }                          \
    } while (0)

#ifndef HPLAI_ACL_BLASPP_STREAM_SIZE
#define HPLAI_ACL_BLASPP_STREAM_SIZE 1
#endif

#ifndef HPLAI_ACL_BLASPP_GEMM_MODEL_DIR
#define HPLAI_ACL_BLASPP_GEMM_MODEL_DIR "op_models"
#endif

static aclrtStream HPLAI_ACL_BLASPP_STREAM[HPLAI_ACL_BLASPP_STREAM_SIZE];
static aclrtContext HPLAI_ACL_BLASPP_CONTEXT;
static aclrtRunMode HPLAI_ACL_BLASPP_RUNMODE;
static int64_t HPLAI_ACL_BLASPP_HOST_BUFFER_SIZE = 0;
static void *HPLAI_ACL_BLASPP_HOST_BUFFER = NULL;
static int64_t HPLAI_ACL_BLASPP_DEVICE_BUFFER_SIZE = 0;
static void *HPLAI_ACL_BLASPP_DEVICE_BUFFER = NULL;

static void HPLAI_ACL_BLASPP_HOST_BUFFER_RESIZE(int64_t NEW_SIZE)
{
    if (HPLAI_ACL_BLASPP_HOST_BUFFER != NULL)
    {
        free(HPLAI_ACL_BLASPP_HOST_BUFFER);
        HPLAI_ACL_BLASPP_HOST_BUFFER = NULL;
    }
    HPLAI_ACL_BLASPP_HOST_BUFFER_SIZE = NEW_SIZE;
    if (HPLAI_ACL_BLASPP_HOST_BUFFER_SIZE > 0)
    {
        HPLAI_ACL_BLASPP_HOST_BUFFER = malloc(HPLAI_ACL_BLASPP_HOST_BUFFER_SIZE);
        if (HPLAI_ACL_BLASPP_HOST_BUFFER == NULL)
            HPLAI_pabort(
                __LINE__,
                "HPLAI_ACL_BLASPP_HOST_BUFFER_RESIZE",
                "Memory allocation failed for HPLAI_ACL_BLASPP_HOST_BUFFER");
    }
}

static void HPLAI_ACL_BLASPP_DEVICE_BUFFER_RESIZE(int64_t NEW_SIZE)
{
    if (HPLAI_ACL_BLASPP_DEVICE_BUFFER != NULL)
    {
        ACLCHECK(aclrtFree(HPLAI_ACL_BLASPP_DEVICE_BUFFER));
        HPLAI_ACL_BLASPP_DEVICE_BUFFER = NULL;
    }
    HPLAI_ACL_BLASPP_DEVICE_BUFFER_SIZE = NEW_SIZE;
    if (HPLAI_ACL_BLASPP_DEVICE_BUFFER_SIZE > 0)
    {
        ACLCHECK(aclrtMalloc(
            &HPLAI_ACL_BLASPP_DEVICE_BUFFER,
            HPLAI_ACL_BLASPP_DEVICE_BUFFER_SIZE,
            ACL_MEM_MALLOC_HUGE_FIRST));

        if (HPLAI_ACL_BLASPP_DEVICE_BUFFER == NULL)
            HPLAI_pabort(
                __LINE__,
                "HPLAI_ACL_BLASPP_DEVICE_BUFFER_RESIZE",
                "Memory allocation failed for HPLAI_ACL_BLASPP_DEVICE_BUFFER");
    }
}

static aclDataType HPLAI_GET_ACL_DataType(aclFloat16 a)
{
    return ACL_FLOAT16;
}

static aclDataType HPLAI_GET_ACL_DataType(float a)
{
    return ACL_FLOAT;
}

static aclDataType HPLAI_GET_ACL_DataType(double a)
{
    return ACL_DOUBLE;
}

#if defined(HPLAI_ACL_BLASPP_GEMM_DEBUG)

#include <sys/file.h>

static void HPLAI_ACL_Cast_JSON(
    int64_t trans,
    int64_t m,
    int64_t n,
    aclDataType dataTypeA,
    aclFormat formatA,
    aclDataType dataTypeC,
    aclFormat formatC)
{
    int64_t dim[] = {m, n};
    if (trans)
        std::swap(dim[0], dim[1]);
    char filename[999];
    sprintf(
        filename,
        "%s/0_Cast_%lld_2_%lld_%lld_%lld_2_%lld_%lld.json",
        HPLAI_ACL_BLASPP_GEMM_MODEL_DIR,
        dataTypeA,
        dim[0],
        dim[1],
        dataTypeC,
        dim[0],
        dim[1]);
    FILE *json = fopen(filename, "w");
    flock(json->_fileno, LOCK_EX);
    HPLAI_fprintf(
        json,
        "[{\"op\":\"Cast\",\"input_desc\":[{\"format\":\"ND\",\"shape\":[%lld,%lld],\"type\":\"%s\"}],\"output_desc\":[{\"format\":\"ND\",\"shape\":[%lld,%lld],\"type\":\"%s\"}],\"attr\":[{\"name\":\"dst_type\",\"type\":\"int\",\"value\":%lld}]}]",
        dim[0],
        dim[1],
        dataTypeA == ACL_FLOAT     ? "float"
        : dataTypeA == ACL_FLOAT16 ? "float16"
        : dataTypeA == ACL_DOUBLE  ? "double"
                                   : "undefined",
        dim[0],
        dim[1],
        dataTypeC == ACL_FLOAT     ? "float"
        : dataTypeC == ACL_FLOAT16 ? "float16"
        : dataTypeC == ACL_DOUBLE  ? "double"
                                   : "undefined",
        dataTypeC);
    fclose(json);
    flock(json->_fileno, LOCK_UN);
}

static void HPLAI_ACL_MatMul_JSON(
    int64_t transA,
    int64_t transB,
    int64_t m,
    int64_t n,
    int64_t k,
    aclDataType dataTypeA,
    aclFormat formatA,
    aclDataType dataTypeB,
    aclFormat formatB,
    aclDataType dataTypeC,
    aclFormat formatC)
{
    int64_t dima[] = {m, k}, dimb[] = {k, n}, dimc[] = {m, n};
    if (transA)
        std::swap(dima[0], dima[1]);
    if (transB)
        std::swap(dimb[0], dimb[1]);
    //atc --soc_version=Ascend910 --output=op_models --singleop=op_models/MatMul.json
    char filename[999];
    sprintf(
        filename,
        "%s/0_MatMulV2_1_2_%lld_%lld_1_2_%lld_%lld_1_2_%lld_%lld.json",
        HPLAI_ACL_BLASPP_GEMM_MODEL_DIR,
        dima[0],
        dima[1],
        dimb[0],
        dimb[1],
        dimc[0],
        dimc[1]);
    FILE *json = fopen(filename, "w");
    flock(json->_fileno, LOCK_EX);
    HPLAI_fprintf(
        json,
        "[{\"op\":\"MatMulV2\",\"input_desc\":[{\"format\":\"ND\",\"shape\":[%lld,%lld],\"type\":\"float16\"},{\"format\":\"ND\",\"shape\":[%lld,%lld],\"type\":\"float16\"}],\"output_desc\":[{\"format\":\"ND\",\"shape\":[%lld,%lld],\"type\":\"float16\"}],\"attr\":[{\"name\":\"transpose_x1\",\"type\":\"bool\",\"value\":%s},{\"name\":\"transpose_x2\",\"type\":\"bool\",\"value\":%s}]}]",
        dima[0],
        dima[1],
        dimb[0],
        dimb[1],
        dimc[0],
        dimc[1],
        transA ? "true" : "false",
        transB ? "true" : "false");
    fclose(json);
    flock(json->_fileno, LOCK_UN);
}

#endif

#if !defined(HPLAI_ACL_BLASPP_GEMM_NO_COMPILE)
#include <acl/acl_op_compiler.h>
#endif

static void HPLAI_ACL_Cast(
    int64_t trans,
    int64_t m,
    int64_t n,
    void *matrixA,
    aclDataType dataTypeA,
    aclFormat formatA,
    void *matrixC,
    aclDataType dataTypeC,
    aclFormat formatC,
    aclrtStream stream)
{
    aclopAttr *opAttr = aclopCreateAttr();
    ACLCHECK(aclopSetAttrInt(opAttr, "dst_type", dataTypeC));
    int64_t dim[] = {m, n};
    if (trans)
        std::swap(dim[0], dim[1]);
    aclTensorDesc *aDesc = aclCreateTensorDesc(
        dataTypeA,
        2,
        dim,
        formatA);
    aclTensorDesc *cDesc = aclCreateTensorDesc(
        dataTypeC,
        2,
        dim,
        formatC);
    aclTensorDesc *inputDesc[] = {aDesc};
    aclTensorDesc *outputDesc[] = {cDesc};
    aclDataBuffer
        *dataA = aclCreateDataBuffer(
            matrixA, dim[0] * dim[1] * sizeof(dataTypeA)),
        *dataC = aclCreateDataBuffer(
            matrixC, dim[0] * dim[1] * sizeof(dataTypeC));
    aclDataBuffer *
        inputs[] = {dataA};
    aclDataBuffer *
        outputs[] = {dataC};

#if !defined(HPLAI_ACL_BLASPP_GEMM_NO_COMPILE)
    ACLCHECK(aclopCompileAndExecute(
        "Cast",
        1,
        inputDesc,
        inputs,
        1,
        outputDesc,
        outputs,
        opAttr,
        ACL_ENGINE_SYS,
        ACL_COMPILE_SYS,
        NULL,
        stream));
#else
    ACLCHECK(aclopExecuteV2(
        "Cast",
        1,
        inputDesc,
        inputs,
        1,
        outputDesc,
        outputs,
        opAttr,
        stream));
#endif

    ACLCHECK(aclDestroyDataBuffer(dataA));
    ACLCHECK(aclDestroyDataBuffer(dataC));
    aclDestroyTensorDesc(aDesc);
    aclDestroyTensorDesc(cDesc);
    aclopDestroyAttr(opAttr);
}

static void HPLAI_ACL_MatMul(
    int64_t transA,
    int64_t transB,
    int64_t m,
    int64_t n,
    int64_t k,
    void *matrixA,
    aclDataType dataTypeA,
    aclFormat formatA,
    void *matrixB,
    aclDataType dataTypeB,
    aclFormat formatB,
    void *matrixC,
    aclDataType dataTypeC,
    aclFormat formatC,
    aclrtStream stream)
{
    aclopAttr *opAttr = aclopCreateAttr();
    ACLCHECK(aclopSetAttrBool(opAttr, "transpose_x1", transA));
    ACLCHECK(aclopSetAttrBool(opAttr, "transpose_x2", transB));
    int64_t dima[] = {m, k}, dimb[] = {k, n}, dimc[] = {m, n};
    if (transA)
        std::swap(dima[0], dima[1]);
    if (transB)
        std::swap(dimb[0], dimb[1]);
    aclTensorDesc *aDesc = aclCreateTensorDesc(
        dataTypeA,
        2,
        dima,
        formatA);
    aclTensorDesc *bDesc = aclCreateTensorDesc(
        dataTypeB,
        2,
        dimb,
        formatB);
    aclTensorDesc *cDesc = aclCreateTensorDesc(
        dataTypeC,
        2,
        dimc,
        formatC);
    aclTensorDesc *inputDesc[] = {aDesc, bDesc};
    aclTensorDesc *outputDesc[] = {cDesc};
    aclDataBuffer
        *dataA = aclCreateDataBuffer(
            matrixA, dima[0] * dima[1] * sizeof(dataTypeA)),
        *dataB = aclCreateDataBuffer(
            matrixB, dimb[0] * dimb[1] * sizeof(dataTypeB)),
        *dataC = aclCreateDataBuffer(
            matrixC, dimc[0] * dimc[1] * sizeof(dataTypeC));
    aclDataBuffer *
        inputs[] = {dataA, dataB};
    aclDataBuffer *
        outputs[] = {dataC};
#if !defined(HPLAI_ACL_BLASPP_GEMM_NO_COMPILE)
    do
    {
        int e = aclopCompileAndExecute(
            "MatMulV2",
            2,
            inputDesc,
            inputs,
            1,
            outputDesc,
            outputs,
            opAttr,
            ACL_ENGINE_SYS,
            ACL_COMPILE_SYS,
            NULL,
            stream);
        if (e != ACL_ERROR_NONE)
        {
            char str[99];
            sprintf(str, "%d %d %d %d", e, m, n, k);
            HPLAI_pabort(
                __LINE__,
                "aclrt",
                str);
        }
    } while (0);
#else
    ACLCHECK(aclopExecuteV2(
        "MatMulV2",
        2,
        inputDesc,
        inputs,
        1,
        outputDesc,
        outputs,
        opAttr,
        stream));
#endif
    ACLCHECK(aclDestroyDataBuffer(dataA));
    ACLCHECK(aclDestroyDataBuffer(dataB));
    ACLCHECK(aclDestroyDataBuffer(dataC));
    aclDestroyTensorDesc(aDesc);
    aclDestroyTensorDesc(bDesc);
    aclDestroyTensorDesc(cDesc);
    aclopDestroyAttr(opAttr);
}

template <>
void blas::gemm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
    blas::Layout layout,
    blas::Op TRANSA,
    blas::Op TRANSB,
    int64_t M,
    int64_t N,
    int64_t K,
    blas::scalar_type<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT> ALPHA,
    HPLAI_T_AFLOAT const *A,
    int64_t LDA,
    HPLAI_T_AFLOAT const *B,
    int64_t LDB,
    blas::scalar_type<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT> BETA,
    HPLAI_T_AFLOAT *C,
    int64_t LDC)
{
    if (layout != blas::Layout::RowMajor)
    {
        blas::gemm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
            blas::Layout::RowMajor,
            TRANSB, TRANSA, N, M, K, ALPHA, B, LDB, A, LDA, BETA, C, LDC);
        return;
    }
    if ((M == 0) || (N == 0) ||
        (((ALPHA == HPLAI_rzero) || (K == 0)) &&
         (BETA == HPLAI_rone)))
        return;

#if !defined(HPLAI_ACL_BLASPP_GEMM_USE_CPU)
#define HPLAI_ACL_BLASPP_GEMM_USE_CPU 2048
#endif
    const int64_t P = HPLAI_ACL_BLASPP_GEMM_USE_CPU;
    if (M < P || N < P || K < P)
    {
        blas::gemm(
            layout,
            TRANSA,
            TRANSB,
            M,
            N,
            K,
            ALPHA,
            A,
            LDA,
            B,
            LDB,
            BETA,
            C,
            LDC);
        return;
    }
    const int64_t
        M2 = M % P,
        N2 = N % P,
        K2 = K % P,
        M1 = M - M2,
        N1 = N - N2,
        K1 = K - K2;

#if defined(HPLAI_ACL_BLASPP_GEMM_DEBUG)
    HPLAI_ACL_Cast_JSON(
        false,
        M1,
        K1,
        HPLAI_GET_ACL_DataType(HPLAI_rzero),
        ACL_FORMAT_ND,
        ACL_FLOAT16,
        ACL_FORMAT_ND);
    HPLAI_ACL_Cast_JSON(
        false,
        K1,
        N1,
        HPLAI_GET_ACL_DataType(HPLAI_rzero),
        ACL_FORMAT_ND,
        ACL_FLOAT16,
        ACL_FORMAT_ND);
    HPLAI_ACL_MatMul_JSON(
        false,
        false,
        M1,
        N1,
        K1,
        ACL_FLOAT16,
        ACL_FORMAT_ND,
        ACL_FLOAT16,
        ACL_FORMAT_ND,
        ACL_FLOAT16,
        ACL_FORMAT_ND);
    HPLAI_ACL_Cast_JSON(
        false,
        M1,
        N1,
        ACL_FLOAT16,
        ACL_FORMAT_ND,
        HPLAI_GET_ACL_DataType(HPLAI_rzero),
        ACL_FORMAT_ND);
#endif

    const int64_t
        sAsize = (M1 * K1 * sizeof(HPLAI_T_AFLOAT) + 63) / 32 * 32,
        sBsize = (K1 * N1 * sizeof(HPLAI_T_AFLOAT) + 63) / 32 * 32,
        sCsize = (M1 * N1 * sizeof(HPLAI_T_AFLOAT) + 63) / 32 * 32,
        hAsize = (M1 * K1 * sizeof(aclFloat16) + 63) / 32 * 32,
        hBsize = (K1 * N1 * sizeof(aclFloat16) + 63) / 32 * 32,
        hCsize = (M1 * N1 * sizeof(aclFloat16) + 63) / 32 * 32;

    int64_t device_buffer_size = hCsize + hBsize + hAsize + sAsize + sBsize;
    if (device_buffer_size < hCsize + sCsize)
        device_buffer_size = hCsize + sCsize;
    if (HPLAI_ACL_BLASPP_DEVICE_BUFFER_SIZE < device_buffer_size)
        HPLAI_ACL_BLASPP_DEVICE_BUFFER_RESIZE(device_buffer_size);

    char *hCdevice = reinterpret_cast<char *>(HPLAI_ACL_BLASPP_DEVICE_BUFFER);
    char *hBdevice = hCdevice + hCsize;
    char *hAdevice = hBdevice + hBsize;
    char *sCdevice = hCdevice + hCsize;
    char *sAdevice = hAdevice + hAsize;
    char *sBdevice = sAdevice + sAsize;
    char *sChost = sCdevice;
    char *sAhost = sAdevice;
    char *sBhost = sBdevice;

    if (HPLAI_ACL_BLASPP_RUNMODE == ACL_HOST)
    {
        int64_t host_buffer_size = sAsize + sBsize;
        if (host_buffer_size < sCsize)
            host_buffer_size = sCsize;

        if (HPLAI_ACL_BLASPP_HOST_BUFFER_SIZE < device_buffer_size)
            HPLAI_ACL_BLASPP_HOST_BUFFER_RESIZE(device_buffer_size);
        sChost = reinterpret_cast<char *>(HPLAI_ACL_BLASPP_HOST_BUFFER),
        sAhost = sChost;
        sBhost = sAhost + sAsize;
    }

    if (TRANSA == blas::Op::NoTrans)
    {
        HPLAI_alacpy(
            K1,
            M1,
            A,
            LDA,
            reinterpret_cast<HPLAI_T_AFLOAT *>(sAhost),
            K1);
    }
    else
    {
        HPLAI_alatcpy(
            K1,
            M1,
            A,
            LDA,
            reinterpret_cast<HPLAI_T_AFLOAT *>(sAhost),
            K1);
    }

    if (HPLAI_ACL_BLASPP_RUNMODE == ACL_HOST)
        ACLCHECK(aclrtMemcpyAsync(
            reinterpret_cast<void *>(sAdevice),
            sAsize,
            reinterpret_cast<const void *>(sAhost),
            sAsize,
            ACL_MEMCPY_HOST_TO_DEVICE,
            HPLAI_ACL_BLASPP_STREAM[0]));

    HPLAI_ACL_Cast(
        false,
        M1,
        K1,
        sAdevice,
        HPLAI_GET_ACL_DataType(HPLAI_rzero),
        ACL_FORMAT_ND,
        hAdevice,
        ACL_FLOAT16,
        ACL_FORMAT_ND,
        HPLAI_ACL_BLASPP_STREAM[0]);

    if (TRANSB == blas::Op::NoTrans)
    {
        HPLAI_alacpy(
            N1,
            K1,
            B,
            LDB,
            reinterpret_cast<HPLAI_T_AFLOAT *>(sBhost),
            N1);
    }
    else
    {
        HPLAI_alatcpy(
            N1,
            K1,
            B,
            LDB,
            reinterpret_cast<HPLAI_T_AFLOAT *>(sBhost),
            N1);
    }

    if (HPLAI_ACL_BLASPP_RUNMODE == ACL_HOST)
        ACLCHECK(aclrtMemcpyAsync(
            reinterpret_cast<void *>(sBdevice),
            sBsize,
            reinterpret_cast<const void *>(sBhost),
            sBsize,
            ACL_MEMCPY_HOST_TO_DEVICE,
            HPLAI_ACL_BLASPP_STREAM[HPLAI_ACL_BLASPP_STREAM_SIZE - 1]));

    HPLAI_ACL_Cast(
        false,
        K1,
        N1,
        sBdevice,
        HPLAI_GET_ACL_DataType(HPLAI_rzero),
        ACL_FORMAT_ND,
        hBdevice,
        ACL_FLOAT16,
        ACL_FORMAT_ND,
        HPLAI_ACL_BLASPP_STREAM[HPLAI_ACL_BLASPP_STREAM_SIZE - 1]);

#if HPLAI_ACL_BLASPP_STREAM_SIZE > 1
    ACLCHECK(aclrtSynchronizeStream(HPLAI_ACL_BLASPP_STREAM[HPLAI_ACL_BLASPP_STREAM_SIZE - 1]));
#endif

    HPLAI_ACL_MatMul(
        false,
        false,
        M1,
        N1,
        K1,
        hAdevice,
        ACL_FLOAT16,
        ACL_FORMAT_ND,
        hBdevice,
        ACL_FLOAT16,
        ACL_FORMAT_ND,
        hCdevice,
        ACL_FLOAT16,
        ACL_FORMAT_ND,
        HPLAI_ACL_BLASPP_STREAM[0]);

    HPLAI_ACL_Cast(
        false,
        M1,
        N1,
        hCdevice,
        ACL_FLOAT16,
        ACL_FORMAT_ND,
        sCdevice,
        HPLAI_GET_ACL_DataType(HPLAI_rzero),
        ACL_FORMAT_ND,
        HPLAI_ACL_BLASPP_STREAM[0]);

    if (HPLAI_ACL_BLASPP_RUNMODE == ACL_HOST)
        ACLCHECK(aclrtMemcpyAsync(
            reinterpret_cast<void *>(sChost),
            sCsize,
            reinterpret_cast<const void *>(sCdevice),
            sCsize,
            ACL_MEMCPY_DEVICE_TO_HOST,
            HPLAI_ACL_BLASPP_STREAM[0]));

    blas::gemm(
        layout,
        TRANSA,
        TRANSB,
        M2,
        N,
        K,
        ALPHA,
        TRANSA == blas::Op::NoTrans ? A + LDA * M1 : A + M1,
        LDA,
        B,
        LDB,
        BETA,
        C + LDC * M1,
        LDC);

    blas::gemm(
        layout,
        TRANSA,
        TRANSB,
        M1,
        N2,
        K,
        ALPHA,
        A,
        LDA,
        TRANSB == blas::Op::NoTrans ? B + N1 : B + LDB * N1,
        LDB,
        BETA,
        C + N1,
        LDC);

    blas::gemm(
        layout,
        TRANSA,
        TRANSB,
        M1,
        N1,
        K2,
        ALPHA,
        TRANSA == blas::Op::NoTrans ? A + K1 : A + LDA * K1,
        LDA,
        TRANSB == blas::Op::NoTrans ? B + LDB * K1 : B + K1,
        LDB,
        BETA,
        C,
        LDC);

    ACLCHECK(aclrtSynchronizeStream(HPLAI_ACL_BLASPP_STREAM[0]));
    {
        HPLAI_T_AFLOAT *C0 = C, *mC0 = reinterpret_cast<HPLAI_T_AFLOAT *>(sChost);
        for (int64_t j = 0; j < M1; ++j, C0 += LDC, mC0 += N1)
            blas::axpy<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(N1, ALPHA, mC0, 1, C0, 1);
    }

    /*
    blas::gemm(
        layout,
        blas::Op::NoTrans,
        blas::Op::NoTrans,
        M1,
        N1,
        K1,
        ALPHA,
        reinterpret_cast<HPLAI_T_AFLOAT *>(sAhost),
        K1,
        reinterpret_cast<HPLAI_T_AFLOAT *>(sBhost),
        N1,
        HPLAI_rone,
        C,
        LDC);
*/
}

#else

template <>
void blas::gemm<float, float, float>(
    blas::Layout layout,
    blas::Op transA,
    blas::Op transB,
    int64_t m,
    int64_t n,
    int64_t k,
    blas::scalar_type<float, float, float> alpha,
    float const *A,
    int64_t lda,
    float const *B,
    int64_t ldb,
    blas::scalar_type<float, float, float> beta,
    float *C,
    int64_t ldc)
{
    blas::gemm(
        layout,
        transA,
        transB,
        m,
        n,
        k,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc);
}

#endif

template <>
void blas::gemv<float, float, float>(
    blas::Layout layout,
    blas::Op trans,
    int64_t m,
    int64_t n,
    blas::scalar_type<float, float, float> alpha,
    float const *A,
    int64_t lda,
    float const *x,
    int64_t incx,
    blas::scalar_type<float, float, float> beta,
    float *y,
    int64_t incy)
{
    blas::gemv(
        layout,
        trans,
        m,
        n,
        alpha,
        A,
        lda,
        x,
        incx,
        beta,
        y,
        incy);
}

template <>
void blas::ger<float, float, float>(
    blas::Layout layout,
    int64_t m,
    int64_t n,
    blas::scalar_type<float, float, float> alpha,
    float const *x,
    int64_t incx,
    float const *y,
    int64_t incy,
    float *A,
    int64_t lda)
{
    blas::ger(
        layout,
        m,
        n,
        alpha,
        x,
        incx,
        y,
        incy,
        A,
        lda);
}

template <>
void blas::scal<float>(
    int64_t n,
    float alpha,
    float *x,
    int64_t incx)
{
    blas::scal(n, alpha, x, incx);
}

#if defined(HPLAI_DEVICE_BLASPP_TRSM)

template <>
void blas::trsm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
    blas::Layout layout,
    blas::Side SIDE,
    blas::Uplo UPLO,
    blas::Op TRANS,
    blas::Diag DIAG,
    int64_t M,
    int64_t N,
    blas::scalar_type<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT> ALPHA,
    HPLAI_T_AFLOAT const *A,
    int64_t LDA,
    HPLAI_T_AFLOAT *B,
    int64_t LDB)
{
#if !defined(HPLAI_DEVICE_BLASPP_TRSM_USE_CPU)
#define HPLAI_DEVICE_BLASPP_TRSM_USE_CPU 256
#endif
    const int64_t P = HPLAI_DEVICE_BLASPP_TRSM_USE_CPU;
    if (M < P || N < P)
    {
        blas::trsm(
            layout,
            SIDE,
            UPLO,
            TRANS,
            DIAG,
            M,
            N,
            ALPHA,
            A,
            LDA,
            B,
            LDB);
        return;
    }

    //HPLAI_pabort( __LINE__, "blas::trsm", "Use HPLAI_DEVICE_BLASPP_TRSM" );
    if (layout != blas::Layout::ColMajor)
    {
        blas::trsm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
            blas::Layout::ColMajor,
            (SIDE == blas::Side::Right ? blas::Side::Left : blas::Side::Right),
            (UPLO == blas::Uplo::Lower ? blas::Uplo::Upper : blas::Uplo::Lower),
            TRANS, DIAG, N, M, ALPHA, A, LDA, B, LDB);
        return;
    }
    int64_t i, j;

    if ((M == 0) || (N == 0))
        return;

    if (ALPHA == HPLAI_rzero)
    {
        for (j = 0; j < N; j++)
        {
            for (i = 0; i < M; i++)
                *(B + i + j * LDB) = HPLAI_rzero;
        }
        return;
    }

    int64_t padding_size = 128 / sizeof(HPLAI_T_AFLOAT);
    if (padding_size < 1)
        padding_size = 1;

    int64_t rB = M;
    int64_t cB = N;
    int64_t dLDB = (rB + padding_size - 1) / padding_size * padding_size;
    int64_t dsB = cB * dLDB;

    int64_t rA = SIDE == blas::Side::Left ? M : N;
    int64_t cA = SIDE == blas::Side::Left ? M : N;
    int64_t dLDA = (rA + padding_size - 1) / padding_size * padding_size;
    int64_t dsA = cA * dLDA;

    if (HPLAI_DEVICE_BLASPP_BUFFER_SIZE < dsB + dsA)
        HPLAI_DEVICE_BLASPP_BUFFER_RESIZE(dsB + dsA);

    HPLAI_T_AFLOAT *dB = HPLAI_DEVICE_BLASPP_BUFFER;
    HPLAI_T_AFLOAT *dA = dB + dsB;

#if defined(HPLAI_DEVICE_BLASPP_TRSM_MULTISTREAM)
    HPLAI_DEVICE_BLASPP_QUEUE->fork();
#endif
    blas::device_setmatrix<HPLAI_T_AFLOAT>(rB, cB, B, LDB, dB, dLDB, *HPLAI_DEVICE_BLASPP_QUEUE);
#if defined(HPLAI_DEVICE_BLASPP_TRSM_MULTISTREAM)
    HPLAI_DEVICE_BLASPP_QUEUE->revolve();
#endif
    blas::device_setmatrix<HPLAI_T_AFLOAT>(rA, cA, A, LDA, dA, dLDA, *HPLAI_DEVICE_BLASPP_QUEUE);
#if defined(HPLAI_DEVICE_BLASPP_TRSM_MULTISTREAM)
    HPLAI_DEVICE_BLASPP_QUEUE->join();
#endif

    blas::trsm(
        layout,
        SIDE,
        UPLO,
        TRANS,
        DIAG,
        M,
        N,
        ALPHA,
        dA,
        dLDA,
        dB,
        dLDB,
        *HPLAI_DEVICE_BLASPP_QUEUE);

    blas::device_getmatrix<HPLAI_T_AFLOAT>(rB, cB, dB, dLDB, B, LDB, *HPLAI_DEVICE_BLASPP_QUEUE);

    HPLAI_DEVICE_BLASPP_QUEUE->sync();
}

#elif defined(HPLAI_GEN_BLASPP_TRSM)

template <typename T>
static void HPLAI_trsmLLNN(
    const int64_t M,
    const int64_t N,
    const T ALPHA,
    const T *A,
    const int64_t LDA,
    T *B,
    const int64_t LDB)
{
    int64_t i, iaik, ibij, ibkj, j, jak, jbj, k;

    for (j = 0, jbj = 0; j < N; j++, jbj += LDB)
    {
        for (i = 0, ibij = jbj; i < M; i++, ibij += 1)
        {
            B[ibij] *= ALPHA;
        }
        for (k = 0, jak = 0, ibkj = jbj; k < M; k++, jak += LDA, ibkj += 1)
        {
            B[ibkj] /= A[k + jak];
            for (i = k + 1, iaik = k + 1 + jak, ibij = k + 1 + jbj;
                 i < M; i++, iaik += 1, ibij += 1)
            {
                B[ibij] -= B[ibkj] * A[iaik];
            }
        }
    }
}

template <typename T>
static void HPLAI_trsmLLNU(
    const int64_t M,
    const int64_t N,
    const T ALPHA,
    const T *A,
    const int64_t LDA,
    T *B,
    const int64_t LDB)
{
    int64_t i, iaik, ibij, ibkj, j, jak, jbj, k;

    for (j = 0, jbj = 0; j < N; j++, jbj += LDB)
    {
        for (i = 0, ibij = jbj; i < M; i++, ibij += 1)
        {
            B[ibij] *= ALPHA;
        }
        for (k = 0, jak = 0, ibkj = jbj; k < M; k++, jak += LDA, ibkj += 1)
        {
            for (i = k + 1, iaik = k + 1 + jak, ibij = k + 1 + jbj;
                 i < M; i++, iaik += 1, ibij += 1)
            {
                B[ibij] -= B[ibkj] * A[iaik];
            }
        }
    }
}

template <typename T>
static void HPLAI_trsmLLTN(
    const int64_t M,
    const int64_t N,
    const T ALPHA,
    const T *A,
    const int64_t LDA,
    T *B,
    const int64_t LDB)
{
    register T t0;
    int64_t i, iaki, ibij, ibkj, j, jai, jbj, k;

    for (j = 0, jbj = 0; j < N; j++, jbj += LDB)
    {
        for (i = M - 1, jai = (M - 1) * LDA, ibij = M - 1 + jbj;
             i >= 0; i--, jai -= LDA, ibij -= 1)
        {
            t0 = ALPHA * B[ibij];
            for (k = i + 1, iaki = i + 1 + jai, ibkj = i + 1 + jbj;
                 k < M; k++, iaki += 1, ibkj += 1)
            {
                t0 -= A[iaki] * B[ibkj];
            }
            t0 /= A[i + jai];
            B[ibij] = t0;
        }
    }
}

template <typename T>
static void HPLAI_trsmLLTU(
    const int64_t M,
    const int64_t N,
    const T ALPHA,
    const T *A,
    const int64_t LDA,
    T *B,
    const int64_t LDB)
{
    register T t0;
    int64_t i, iaki, ibij, ibkj, j, jai, jbj, k;

    for (j = 0, jbj = 0; j < N; j++, jbj += LDB)
    {
        for (i = M - 1, jai = (M - 1) * LDA, ibij = M - 1 + jbj;
             i >= 0; i--, jai -= LDA, ibij -= 1)
        {
            t0 = ALPHA * B[ibij];
            for (k = i + 1, iaki = i + 1 + jai, ibkj = i + 1 + jbj;
                 k < M; k++, iaki += 1, ibkj += 1)
            {
                t0 -= A[iaki] * B[ibkj];
            }
            B[ibij] = t0;
        }
    }
}

template <typename T>
static void HPLAI_trsmLUNN(
    const int64_t M,
    const int64_t N,
    const T ALPHA,
    const T *A,
    const int64_t LDA,
    T *B,
    const int64_t LDB)
{
    int64_t i, iaik, ibij, ibkj, j, jak, jbj, k;

    for (j = 0, jbj = 0; j < N; j++, jbj += LDB)
    {
        for (i = 0, ibij = jbj; i < M; i++, ibij += 1)
        {
            B[ibij] *= ALPHA;
        }
        for (k = M - 1, jak = (M - 1) * LDA, ibkj = M - 1 + jbj;
             k >= 0; k--, jak -= LDA, ibkj -= 1)
        {
            B[ibkj] /= A[k + jak];
            for (i = 0, iaik = jak, ibij = jbj;
                 i < k; i++, iaik += 1, ibij += 1)
            {
                B[ibij] -= B[ibkj] * A[iaik];
            }
        }
    }
}

template <typename T>
static void HPLAI_trsmLUNU(
    const int64_t M,
    const int64_t N,
    const T ALPHA,
    const T *A,
    const int64_t LDA,
    T *B,
    const int64_t LDB)
{
    int64_t i, iaik, ibij, ibkj, j, jak, jbj, k;

    for (j = 0, jbj = 0; j < N; j++, jbj += LDB)
    {
        for (i = 0, ibij = jbj; i < M; i++, ibij += 1)
        {
            B[ibij] *= ALPHA;
        }
        for (k = M - 1, jak = (M - 1) * LDA, ibkj = M - 1 + jbj;
             k >= 0; k--, jak -= LDA, ibkj -= 1)
        {
            for (i = 0, iaik = jak, ibij = jbj;
                 i < k; i++, iaik += 1, ibij += 1)
            {
                B[ibij] -= B[ibkj] * A[iaik];
            }
        }
    }
}

template <typename T>
static void HPLAI_trsmLUTN(
    const int64_t M,
    const int64_t N,
    const T ALPHA,
    const T *A,
    const int64_t LDA,
    T *B,
    const int64_t LDB)
{
    int64_t i, iaki, ibij, ibkj, j, jai, jbj, k;
    register T t0;

    for (j = 0, jbj = 0; j < N; j++, jbj += LDB)
    {
        for (i = 0, jai = 0, ibij = jbj; i < M; i++, jai += LDA, ibij += 1)
        {
            t0 = ALPHA * B[ibij];
            for (k = 0, iaki = jai, ibkj = jbj; k < i; k++, iaki += 1, ibkj += 1)
            {
                t0 -= A[iaki] * B[ibkj];
            }
            t0 /= A[i + jai];
            B[ibij] = t0;
        }
    }
}

template <typename T>
static void HPLAI_trsmLUTU(
    const int64_t M,
    const int64_t N,
    const T ALPHA,
    const T *A,
    const int64_t LDA,
    T *B,
    const int64_t LDB)
{
    register T t0;
    int64_t i, iaki, ibij, ibkj, j, jai, jbj, k;

    for (j = 0, jbj = 0; j < N; j++, jbj += LDB)
    {
        for (i = 0, jai = 0, ibij = jbj; i < M; i++, jai += LDA, ibij += 1)
        {
            t0 = ALPHA * B[ibij];
            for (k = 0, iaki = jai, ibkj = jbj; k < i; k++, iaki += 1, ibkj += 1)
            {
                t0 -= A[iaki] * B[ibkj];
            }
            B[ibij] = t0;
        }
    }
}

template <typename T>
static void HPLAI_trsmRLNN(
    const int64_t M,
    const int64_t N,
    const T ALPHA,
    const T *A,
    const int64_t LDA,
    T *B,
    const int64_t LDB)
{
    int64_t i, iakj, ibij, ibik, j, jaj, jbj, jbk, k;

    for (j = N - 1, jaj = (N - 1) * LDA, jbj = (N - 1) * LDB;
         j >= 0; j--, jaj -= LDA, jbj -= LDB)
    {
        for (i = 0, ibij = jbj; i < M; i++, ibij += 1)
        {
            B[ibij] *= ALPHA;
        }
        for (k = j + 1, iakj = j + 1 + jaj, jbk = (j + 1) * LDB;
             k < N; k++, iakj += 1, jbk += LDB)
        {
            for (i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1)
            {
                B[ibij] -= A[iakj] * B[ibik];
            }
        }
        for (i = 0, ibij = jbj; i < M; i++, ibij += 1)
        {
            B[ibij] /= A[j + jaj];
        }
    }
}

template <typename T>
static void HPLAI_trsmRLNU(
    const int64_t M,
    const int64_t N,
    const T ALPHA,
    const T *A,
    const int64_t LDA,
    T *B,
    const int64_t LDB)
{
    int64_t i, iakj, ibij, ibik, j, jaj, jbj, jbk, k;

    for (j = N - 1, jaj = (N - 1) * LDA, jbj = (N - 1) * LDB;
         j >= 0; j--, jaj -= LDA, jbj -= LDB)
    {
        for (i = 0, ibij = jbj; i < M; i++, ibij += 1)
        {
            B[ibij] *= ALPHA;
        }
        for (k = j + 1, iakj = j + 1 + jaj, jbk = (j + 1) * LDB;
             k < N; k++, iakj += 1, jbk += LDB)
        {
            for (i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1)
            {
                B[ibij] -= A[iakj] * B[ibik];
            }
        }
    }
}

template <typename T>
static void HPLAI_trsmRLTN(
    const int64_t M,
    const int64_t N,
    const T ALPHA,
    const T *A,
    const int64_t LDA,
    T *B,
    const int64_t LDB)
{
    register T t0;
    int64_t i, iajk, ibij, ibik, j, jak, jbj, jbk, k;

    for (k = 0, jak = 0, jbk = 0; k < N; k++, jak += LDA, jbk += LDB)
    {
        for (i = 0, ibik = jbk; i < M; i++, ibik += 1)
        {
            B[ibik] /= A[k + jak];
        }
        for (j = k + 1, iajk = (k + 1) + jak, jbj = (k + 1) * LDB;
             j < N; j++, iajk += 1, jbj += LDB)
        {
            t0 = A[iajk];
            for (i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1)
            {
                B[ibij] -= t0 * B[ibik];
            }
        }
        for (i = 0, ibik = jbk; i < M; i++, ibik += 1)
        {
            B[ibik] *= ALPHA;
        }
    }
}

template <typename T>
static void HPLAI_trsmRLTU(
    const int64_t M,
    const int64_t N,
    const T ALPHA,
    const T *A,
    const int64_t LDA,
    T *B,
    const int64_t LDB)
{
    register T t0;
    int64_t i, iajk, ibij, ibik, j, jak, jbj, jbk, k;

    for (k = 0, jak = 0, jbk = 0; k < N; k++, jak += LDA, jbk += LDB)
    {
        for (j = k + 1, iajk = (k + 1) + jak, jbj = (k + 1) * LDB;
             j < N; j++, iajk += 1, jbj += LDB)
        {
            t0 = A[iajk];
            for (i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1)
            {
                B[ibij] -= t0 * B[ibik];
            }
        }
        for (i = 0, ibik = jbk; i < M; i++, ibik += 1)
        {
            B[ibik] *= ALPHA;
        }
    }
}

template <typename T>
static void HPLAI_trsmRUNN(
    const int64_t M,
    const int64_t N,
    const T ALPHA,
    const T *A,
    const int64_t LDA,
    T *B,
    const int64_t LDB)
{
    int64_t i, iakj, ibij, ibik, j, jaj, jbj, jbk, k;

    for (j = 0, jaj = 0, jbj = 0; j < N; j++, jaj += LDA, jbj += LDB)
    {
        for (i = 0, ibij = jbj; i < M; i++, ibij += 1)
        {
            B[ibij] *= ALPHA;
        }
        for (k = 0, iakj = jaj, jbk = 0; k < j; k++, iakj += 1, jbk += LDB)
        {
            for (i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1)
            {
                B[ibij] -= A[iakj] * B[ibik];
            }
        }
        for (i = 0, ibij = jbj; i < M; i++, ibij += 1)
        {
            B[ibij] /= A[j + jaj];
        }
    }
}

template <typename T>
static void HPLAI_trsmRUNU(
    const int64_t M,
    const int64_t N,
    const T ALPHA,
    const T *A,
    const int64_t LDA,
    T *B,
    const int64_t LDB)
{
    int64_t i, iakj, ibij, ibik, j, jaj, jbj, jbk, k;

    for (j = 0, jaj = 0, jbj = 0; j < N; j++, jaj += LDA, jbj += LDB)
    {
        for (i = 0, ibij = jbj; i < M; i++, ibij += 1)
        {
            B[ibij] *= ALPHA;
        }
        for (k = 0, iakj = jaj, jbk = 0; k < j; k++, iakj += 1, jbk += LDB)
        {
            for (i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1)
            {
                B[ibij] -= A[iakj] * B[ibik];
            }
        }
    }
}

template <typename T>
static void HPLAI_trsmRUTN(
    const int64_t M,
    const int64_t N,
    const T ALPHA,
    const T *A,
    const int64_t LDA,
    T *B,
    const int64_t LDB)
{
    register T t0;
    int64_t i, iajk, ibij, ibik, j, jak, jbj, jbk, k;

    for (k = N - 1, jak = (N - 1) * LDA, jbk = (N - 1) * LDB;
         k >= 0; k--, jak -= LDA, jbk -= LDB)
    {
        for (i = 0, ibik = jbk; i < M; i++, ibik += 1)
        {
            B[ibik] /= A[k + jak];
        }
        for (j = 0, iajk = jak, jbj = 0; j < k; j++, iajk += 1, jbj += LDB)
        {
            t0 = A[iajk];
            for (i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1)
            {
                B[ibij] -= t0 * B[ibik];
            }
        }
        for (i = 0, ibik = jbk; i < M; i++, ibik += 1)
        {
            B[ibik] *= ALPHA;
        }
    }
}

template <typename T>
static void HPLAI_trsmRUTU(
    const int64_t M,
    const int64_t N,
    const T ALPHA,
    const T *A,
    const int64_t LDA,
    T *B,
    const int64_t LDB)
{
    register T t0;
    int64_t i, iajk, ibij, ibik, j, jak, jbj, jbk, k;

    for (k = N - 1, jak = (N - 1) * LDA, jbk = (N - 1) * LDB;
         k >= 0; k--, jak -= LDA, jbk -= LDB)
    {
        for (j = 0, iajk = jak, jbj = 0; j < k; j++, iajk += 1, jbj += LDB)
        {
            t0 = A[iajk];
            for (i = 0, ibij = jbj, ibik = jbk; i < M; i++, ibij += 1, ibik += 1)
            {
                B[ibij] -= t0 * B[ibik];
            }
        }
        for (i = 0, ibik = jbk; i < M; i++, ibik += 1)
        {
            B[ibik] *= ALPHA;
        }
    }
}

template <>
void blas::trsm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
    blas::Layout layout,
    blas::Side SIDE,
    blas::Uplo UPLO,
    blas::Op TRANS,
    blas::Diag DIAG,
    int64_t M,
    int64_t N,
    blas::scalar_type<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT> ALPHA,
    HPLAI_T_AFLOAT const *A,
    int64_t LDA,
    HPLAI_T_AFLOAT *B,
    int64_t LDB)
{
    //HPLAI_pabort( __LINE__, "blas::trsm", "Use HPLAI_GEN_BLASPP_TRSM" );
    if (layout != blas::Layout::ColMajor)
    {
        blas::trsm<HPLAI_T_AFLOAT, HPLAI_T_AFLOAT>(
            blas::Layout::ColMajor,
            (SIDE == blas::Side::Right ? blas::Side::Left : blas::Side::Right),
            (UPLO == blas::Uplo::Lower ? blas::Uplo::Upper : blas::Uplo::Lower),
            TRANS, DIAG, N, M, ALPHA, A, LDA, B, LDB);
        return;
    }
    int64_t i, j;

    if ((M == 0) || (N == 0))
        return;

    if (ALPHA == HPLAI_rzero)
    {
        for (j = 0; j < N; j++)
        {
            for (i = 0; i < M; i++)
                *(B + i + j * LDB) = HPLAI_rzero;
        }
        return;
    }

    if (SIDE == blas::Side::Left)
    {
        if (UPLO == blas::Uplo::Upper)
        {
            if (TRANS == blas::Op::NoTrans)
            {
                if (DIAG == blas::Diag::NonUnit)
                {
                    HPLAI_trsmLUNN(M, N, ALPHA, A, LDA, B, LDB);
                }
                else
                {
                    HPLAI_trsmLUNU(M, N, ALPHA, A, LDA, B, LDB);
                }
            }
            else
            {
                if (DIAG == blas::Diag::NonUnit)
                {
                    HPLAI_trsmLUTN(M, N, ALPHA, A, LDA, B, LDB);
                }
                else
                {
                    HPLAI_trsmLUTU(M, N, ALPHA, A, LDA, B, LDB);
                }
            }
        }
        else
        {
            if (TRANS == blas::Op::NoTrans)
            {
                if (DIAG == blas::Diag::NonUnit)
                {
                    HPLAI_trsmLLNN(M, N, ALPHA, A, LDA, B, LDB);
                }
                else
                {
                    HPLAI_trsmLLNU(M, N, ALPHA, A, LDA, B, LDB);
                }
            }
            else
            {
                if (DIAG == blas::Diag::NonUnit)
                {
                    HPLAI_trsmLLTN(M, N, ALPHA, A, LDA, B, LDB);
                }
                else
                {
                    HPLAI_trsmLLTU(M, N, ALPHA, A, LDA, B, LDB);
                }
            }
        }
    }
    else
    {
        if (UPLO == blas::Uplo::Upper)
        {
            if (TRANS == blas::Op::NoTrans)
            {
                if (DIAG == blas::Diag::NonUnit)
                {
                    HPLAI_trsmRUNN(M, N, ALPHA, A, LDA, B, LDB);
                }
                else
                {
                    HPLAI_trsmRUNU(M, N, ALPHA, A, LDA, B, LDB);
                }
            }
            else
            {
                if (DIAG == blas::Diag::NonUnit)
                {
                    HPLAI_trsmRUTN(M, N, ALPHA, A, LDA, B, LDB);
                }
                else
                {
                    HPLAI_trsmRUTU(M, N, ALPHA, A, LDA, B, LDB);
                }
            }
        }
        else
        {
            if (TRANS == blas::Op::NoTrans)
            {
                if (DIAG == blas::Diag::NonUnit)
                {
                    HPLAI_trsmRLNN(M, N, ALPHA, A, LDA, B, LDB);
                }
                else
                {
                    HPLAI_trsmRLNU(M, N, ALPHA, A, LDA, B, LDB);
                }
            }
            else
            {
                if (DIAG == blas::Diag::NonUnit)
                {
                    HPLAI_trsmRLTN(M, N, ALPHA, A, LDA, B, LDB);
                }
                else
                {
                    HPLAI_trsmRLTU(M, N, ALPHA, A, LDA, B, LDB);
                }
            }
        }
    }
}

#else

template <>
void blas::trsm<float, float>(
    blas::Layout layout,
    blas::Side side,
    blas::Uplo uplo,
    blas::Op trans,
    blas::Diag diag,
    int64_t m,
    int64_t n,
    blas::scalar_type<float, float> alpha,
    float const *A,
    int64_t lda,
    float *B,
    int64_t ldb)
{
    blas::trsm(
        layout,
        side,
        uplo,
        trans,
        diag,
        m,
        n,
        alpha,
        A,
        lda,
        B,
        ldb);
}

#endif

template <>
void blas::trsv<float, float>(
    blas::Layout layout,
    blas::Uplo uplo,
    blas::Op trans,
    blas::Diag diag,
    int64_t n,
    float const *A,
    int64_t lda,
    float *x,
    int64_t incx)
{
    blas::trsv(
        layout,
        uplo,
        trans,
        diag,
        n,
        A,
        lda,
        x,
        incx);
}

#ifdef __cplusplus
extern "C"
{
#endif

    MPI_Datatype HPLAI_MPI_AFLOAT;

#ifdef STDC_HEADERS
    void HPLAI_blas_init(
        const int RANK,
        const int SIZE)
#else
void HPLAI_blas_init(RANK, SIZE)
    const int RANK,
    SIZE;
#endif
    {
        MPI_Type_contiguous(sizeof(HPLAI_T_AFLOAT), MPI_BYTE, &HPLAI_MPI_AFLOAT);
        MPI_Type_commit(&HPLAI_MPI_AFLOAT);

#if defined(HPLAI_DEVICE_BLASPP_GEMM) || defined(HPLAI_DEVICE_BLASPP_TRSM)
        // https://github.com/NVIDIA/multi-gpu-programming-models/blob/master/mpi/jacobi.cpp
        MPI_Comm local_comm;
        MPI_Comm_split_type(
            MPI_COMM_WORLD,
            MPI_COMM_TYPE_SHARED,
            RANK,
            MPI_INFO_NULL,
            &local_comm);
        int local_rank = -1;
        MPI_Comm_rank(local_comm, &local_rank);
        MPI_Comm_free(&local_comm);
#if BLASPP_VERSION >= 20210400
        local_rank %= blas::get_device_count();
#endif
        HPLAI_DEVICE_BLASPP_QUEUE = new blas::Queue(local_rank, 0); // batch_limit_ = batch_size = 0 // 无需 blas::batch
#elif defined(HPLAI_ACL_BLASPP_GEMM)
    MPI_Comm local_comm;
    MPI_Comm_split_type(
        MPI_COMM_WORLD,
        MPI_COMM_TYPE_SHARED,
        RANK,
        MPI_INFO_NULL,
        &local_comm);
    int local_rank = -1;
    MPI_Comm_rank(local_comm, &local_rank);
    MPI_Comm_free(&local_comm);
    ACLCHECK(aclInit(NULL));
    uint32_t count = 0;
    aclrtGetDeviceCount(&count);
    ACLCHECK(aclrtCreateContext(&HPLAI_ACL_BLASPP_CONTEXT, local_rank % count));
    ACLCHECK(aclrtSetCurrentContext(HPLAI_ACL_BLASPP_CONTEXT));
    ACLCHECK(aclrtGetRunMode(&HPLAI_ACL_BLASPP_RUNMODE));
    ACLCHECK(aclopSetModelDir(HPLAI_ACL_BLASPP_GEMM_MODEL_DIR));
    for (int i = 0; i < HPLAI_ACL_BLASPP_STREAM_SIZE; ++i)
        ACLCHECK(aclrtCreateStream(HPLAI_ACL_BLASPP_STREAM + i));
#endif

#ifdef HPL_CALL_VSIPL
        vsip_init((void *)0);
#endif
    }

    void HPLAI_blas_finalize()
    {
#ifdef HPL_CALL_VSIPL
        vsip_finalize((void *)0);
#endif

#if defined(HPLAI_DEVICE_BLASPP_GEMM) || defined(HPLAI_DEVICE_BLASPP_TRSM)
        HPLAI_DEVICE_BLASPP_BUFFER_RESIZE(0);
        delete HPLAI_DEVICE_BLASPP_QUEUE;
#elif defined(HPLAI_ACL_BLASPP_GEMM)
    HPLAI_ACL_BLASPP_HOST_BUFFER_RESIZE(0);
    HPLAI_ACL_BLASPP_DEVICE_BUFFER_RESIZE(0);
    for (int i = 0; i < HPLAI_ACL_BLASPP_STREAM_SIZE; ++i)
        ACLCHECK(aclrtDestroyStream(HPLAI_ACL_BLASPP_STREAM[i]));
    ACLCHECK(aclrtDestroyContext(HPLAI_ACL_BLASPP_CONTEXT));
    ACLCHECK(aclFinalize());
#endif
        MPI_Type_free(&HPLAI_MPI_AFLOAT);
    }

#ifdef __cplusplus
}
#endif
