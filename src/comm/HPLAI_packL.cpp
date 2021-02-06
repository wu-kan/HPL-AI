/*
 * Include files
 */
#include "hplai.h"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
   int HPLAI_packL(
       HPLAI_T_panel *PANEL,
       const int INDEX,
       const int LEN,
       const int IBUF)
#else
int HPLAI_packL(PANEL, INDEX, LEN, IBUF)
    HPLAI_T_panel *PANEL;
const int INDEX;
const int LEN;
const int IBUF;
#endif
   {
/* 
 * Purpose
 * =======
 *
 * HPLAI_packL forms  the MPI data type for the panel to be broadcast.
 * Successful  completion  is  indicated  by  the  returned  error  code
 * MPI_SUCCESS.
 *
 * Arguments
 * =========
 *
 * PANEL   (input/output)                HPLAI_T_panel *
 *         On entry,  PANEL  points to the  current panel data structure
 *         being broadcast.
 *
 * INDEX   (input)                       const int
 *         On entry,  INDEX  points  to  the  first entry of the  packed
 *         buffer being broadcast.
 *
 * LEN     (input)                       const int
 *         On entry, LEN is the length of the packed buffer.
 *
 * IBUF    (input)                       const int
 *         On entry, IBUF  specifies the panel buffer/count/type entries
 *         that should be initialized.
 *
 * ---------------------------------------------------------------------
 */
#ifdef HPL_USE_MPI_DATATYPE
/*
 * .. Local Variables ..
 */
#ifndef HPL_COPY_L
      MPI_Datatype *type = NULL;
      void ***bufs = NULL;
      HPLAI_T_AFLOAT *A;
      int *blen = NULL;
      MPI_Aint *disp = NULL;
      int curr, i, i1, ibuf, ierr = MPI_SUCCESS, j1,
                             jb, jbm, jbp1, lda, len, m, m1, nbufs;
#else
      int ierr;
#endif
/* ..
 * .. Executable Statements ..
 */
#ifdef HPL_COPY_L
      /*
 * Panel + L1 + DPIV  have been copied into a contiguous buffer - Create
 * and commit a contiguous data type
 */
      PANEL->buffers[IBUF] = (void ***)(PANEL->L2 + INDEX);
      PANEL->counts[IBUF] = 1;

      ierr = MPI_Type_contiguous(LEN * sizeof(HPLAI_T_AFLOAT), MPI_BYTE, &PANEL->dtypes[IBUF]);
      if (ierr == MPI_SUCCESS)
         ierr = MPI_Type_commit(&PANEL->dtypes[IBUF]);

      return (ierr);
#else
      /*
 * Panel is not contiguous (because of LDA and also L1 + DPIV) -  Create
 * and commit a struct data type
 */
      jbp1 = (jb = PANEL->jb) + 1;
      /*
 * Temporaries to create the type struct.
 */
      bufs = (void ***)malloc(jbp1 * sizeof(void **));
      blen = (int *)malloc(jbp1 * sizeof(int));
      disp = (MPI_Aint *)malloc(jbp1 * sizeof(MPI_Aint));
      type = (MPI_Datatype *)malloc(jbp1 * sizeof(MPI_Datatype));

      if ((bufs != NULL) && (blen != NULL) &&
          (disp != NULL) && (type != NULL))
      {
         m = PANEL->mp;
         curr = (int)(PANEL->grid->myrow == PANEL->prow);
         if (curr != 0)
            m -= jb;

         len = LEN;
         ibuf = INDEX;
         nbufs = 0;
         jbm = jb * m;

         if ((m > 0) && (ibuf < jbm))
         {
            /*
 * Retrieve proper pointers depending on process row and column
 */
            if (PANEL->grid->mycol == PANEL->pcol)
            {
               lda = PANEL->lda;
               if (curr != 0)
               {
                  A = Mptr(PANEL->A, jb, -jb, lda);
               }
               else
               {
                  A = Mptr(PANEL->A, 0, -jb, lda);
               }
            }
            else
            {
               lda = PANEL->ldl2;
               A = PANEL->L2;
            }
            /*
 * Pack the first (partial) column of L
 */
            m1 = m - (i1 = ibuf - (j1 = ibuf / m) * m);
            m1 = Mmin(len, m1);

            bufs[nbufs] = (void **)(Mptr(A, i1, j1, lda));
            type[nbufs] = MPI_BYTE;
            blen[nbufs] = m1 * sizeof(HPLAI_T_AFLOAT);
            if (ierr == MPI_SUCCESS)
               ierr = MPI_Get_address(bufs[nbufs], &disp[nbufs]);

            nbufs++;
            len -= m1;
            j1++;
            ibuf += m1;
            /*
 * Pack the remaining columns of L
 */
            while ((len > 0) && (j1 < jb))
            {
               m1 = Mmin(len, m);

               bufs[nbufs] = (void **)(Mptr(A, 0, j1, lda));
               type[nbufs] = MPI_BYTE;
               blen[nbufs] = m1 * sizeof(HPLAI_T_AFLOAT);
               if (ierr == MPI_SUCCESS)
                  ierr = MPI_Get_address(bufs[nbufs], &disp[nbufs]);

               nbufs++;
               len -= m1;
               j1++;
               ibuf += m1;
            }
         }
         /*
 * Pack L1, DPIV, DINFO
 */
         if (len > 0)
         { /* L1, DPIV, DINFO */
            bufs[nbufs] = (void **)(PANEL->L1 + ibuf - jbm);
            type[nbufs] = MPI_BYTE;
            blen[nbufs] = len * sizeof(HPLAI_T_AFLOAT);
            if (ierr == MPI_SUCCESS)
               ierr = MPI_Get_address(bufs[nbufs], &disp[nbufs]);
            nbufs++;
         }

         for (i = 1; i < nbufs; i++)
            disp[i] -= disp[0];
         disp[0] = 0;

         PANEL->buffers[IBUF] = (void ***)bufs[0];
         PANEL->counts[IBUF] = 1;
         /*
 * construct the struct type 
 */
         if (ierr == MPI_SUCCESS)
            ierr = MPI_Type_create_struct(nbufs, blen, disp, type,
                                          &PANEL->dtypes[IBUF]);
         /*
 * release temporaries
 */
         if (bufs)
            free(bufs);
         if (blen)
            free(blen);
         if (disp)
            free(disp);
         if (type)
            free(type);
         /*
 * commit the type 
 */
         if (ierr == MPI_SUCCESS)
            ierr = MPI_Type_commit(&PANEL->dtypes[IBUF]);

         return (ierr);
      }
      else
      {
         /*
 * Memory allocation failed -> abort
 */
         HPL_pabort(__LINE__, "HPLAI_packL", "Memory allocation failed");
         return (MPI_SUCCESS); /* never executed (hopefully ...) */
      }
#endif
#else
   /* HPL_USE_MPI_DATATYPE not defined - Oops, there is a bug
             somewhere, so, just in case  and until I find it ... */
   return (MPI_SUCCESS);
#endif
      /*
 * End of HPLAI_packL
 */
   }

#ifdef __cplusplus
}
#endif
