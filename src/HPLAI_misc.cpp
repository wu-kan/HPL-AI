#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

MPI_Datatype HPLAI_MPI_AFLOAT;

#ifdef STDC_HEADERS
    void HPLAI_init(
        const int RANK,
        const int SIZE)
#else
void HPLAI_init(RANK, SIZE)
    const int RANK,
    SIZE;
#endif
    {
        MPI_Type_contiguous(sizeof(HPLAI_T_AFLOAT), MPI_BYTE, &HPLAI_MPI_AFLOAT);
        MPI_Type_commit(&HPLAI_MPI_AFLOAT);
#ifdef HPL_CALL_VSIPL
        vsip_init((void *)0);
#endif
    }

    void HPLAI_finalize()
    {
#ifdef HPL_CALL_VSIPL
        vsip_finalize((void *)0);
#endif
        MPI_Type_free(&HPLAI_MPI_AFLOAT);
    }

#ifdef __cplusplus
}
#endif