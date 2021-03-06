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
/*
 * Include files
 */
#include "hplai.hh"

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef STDC_HEADERS
    int main(
        int ARGC,
        char **ARGV)
#else
int main(ARGC, ARGV)
    /*
 * .. Scalar Arguments ..
 */
    int ARGC;
/*
 * .. Array Arguments ..
 */
char **ARGV;
#endif
    {
        /* 
 * Purpose
 * =======
 *
 * main is the main driver program for testing the HPL routines.
 * This  program is  driven  by  a short data file named  "HPL.dat".
 *
 * ---------------------------------------------------------------------
 */
        /*
 * .. Local Variables ..
 */
        int nval[HPLAI_MAX_PARAM],
            nbval[HPLAI_MAX_PARAM],
            pval[HPLAI_MAX_PARAM],
            qval[HPLAI_MAX_PARAM],
            nbmval[HPLAI_MAX_PARAM],
            ndvval[HPLAI_MAX_PARAM],
            ndhval[HPLAI_MAX_PARAM];

        HPLAI_T_FACT pfaval[HPLAI_MAX_PARAM],
            rfaval[HPLAI_MAX_PARAM];

        HPLAI_T_TOP topval[HPLAI_MAX_PARAM];

        HPLAI_T_grid grid;
        HPLAI_T_palg algo;
        HPLAI_T_test test;
        int L1notran, Unotran, align, equil, in, inb,
            inbm, indh, indv, ipfa, ipq, irfa, itop,
            mycol, myrow, ns, nbs, nbms, ndhs, ndvs,
            npcol, npfs, npqs, nprow, nrfs, ntps,
            rank, size, tswap;
        HPLAI_T_ORDER pmapping;
        HPLAI_T_FACT rpfa;
        HPLAI_T_SWAP fswap;
        /* ..
 * .. Executable Statements ..
 */
        MPI_Init(&ARGC, &ARGV);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        HPLAI_blas_init(rank, size);
        /*
 * Read and check validity of test parameters from input file
 *
 * HPL Version 1.0, Linpack benchmark input file
 * Your message here
 * HPL.out      output file name (if any)
 * 6            device out (6=stdout,7=stderr,file)
 * 4            # of problems sizes (N)
 * 29 30 34 35  Ns
 * 4            # of NBs
 * 1 2 3 4      NBs
 * 0            PMAP process mapping (0=Row-,1=Column-major)
 * 3            # of process grids (P x Q)
 * 2 1 4        Ps
 * 2 4 1        Qs
 * 16.0         threshold
 * 3            # of panel fact
 * 0 1 2        PFACTs (0=left, 1=Crout, 2=Right)
 * 2            # of recursive stopping criterium
 * 2 4          NBMINs (>= 1)
 * 1            # of panels in recursion
 * 2            NDIVs
 * 3            # of recursive panel fact.
 * 0 1 2        RFACTs (0=left, 1=Crout, 2=Right)
 * 1            # of broadcast
 * 0            BCASTs (0=1rg,1=1rM,2=2rg,3=2rM,4=Lng,5=LnM)
 * 1            # of lookahead depth
 * 0            DEPTHs (>=0)
 * 2            SWAP (0=bin-exch,1=long,2=mix)
 * 4            swapping threshold
 * 0            L1 in (0=transposed,1=no-transposed) form
 * 0            U  in (0=transposed,1=no-transposed) form
 * 1            Equilibration (0=no,1=yes)
 * 8            memory alignment in double (> 0)
 */
        HPLAI_pdinfo(&test, &ns, nval, &nbs, nbval, &pmapping, &npqs, pval, qval,
                     &npfs, pfaval, &nbms, nbmval, &ndvs, ndvval, &nrfs, rfaval,
                     &ntps, topval, &ndhs, ndhval, &fswap, &tswap, &L1notran,
                     &Unotran, &equil, &align);
        /*
 * Loop over different process grids - Define process grid. Go to bottom
 * of process grid loop if this case does not use my process.
 */
        for (ipq = 0; ipq < npqs; ipq++)
        {
            (void)HPLAI_grid_init(MPI_COMM_WORLD, pmapping, pval[ipq], qval[ipq],
                                  &grid);
            (void)HPLAI_grid_info(&grid, &nprow, &npcol, &myrow, &mycol);

            if ((myrow < 0) || (myrow >= nprow) ||
                (mycol < 0) || (mycol >= npcol))
                goto label_end_of_npqs;

            for (in = 0; in < ns; in++)
            { /* Loop over various problem sizes */
                for (inb = 0; inb < nbs; inb++)
                { /* Loop over various blocking factors */
                    for (indh = 0; indh < ndhs; indh++)
                    { /* Loop over various lookahead depths */
                        for (itop = 0; itop < ntps; itop++)
                        { /* Loop over various broadcast topologies */
                            for (irfa = 0; irfa < nrfs; irfa++)
                            { /* Loop over various recursive factorizations */
                                for (ipfa = 0; ipfa < npfs; ipfa++)
                                { /* Loop over various panel factorizations */
                                    for (inbm = 0; inbm < nbms; inbm++)
                                    { /* Loop over various recursive stopping criteria */
                                        for (indv = 0; indv < ndvs; indv++)
                                        { /* Loop over various # of panels in recursion */
                                            /*
 * Set up the algorithm parameters
 */
                                            algo.btopo = topval[itop];
                                            algo.depth = ndhval[indh];
                                            algo.nbmin = nbmval[inbm];
                                            algo.nbdiv = ndvval[indv];

                                            algo.pfact = rpfa = pfaval[ipfa];

                                            if (L1notran != 0)
                                            {
                                                if (rpfa == HPLAI_LEFT_LOOKING)
                                                    algo.pffun = HPLAI_papanllN;
                                                else if (rpfa == HPLAI_CROUT)
                                                    algo.pffun = HPLAI_papancrN;
                                                else
                                                    algo.pffun = HPLAI_papanrlN;

                                                algo.rfact = rpfa = rfaval[irfa];
                                                if (rpfa == HPLAI_LEFT_LOOKING)
                                                    algo.rffun = HPLAI_parpanllN;
                                                else if (rpfa == HPLAI_CROUT)
                                                    algo.rffun = HPLAI_parpancrN;
                                                else
                                                    algo.rffun = HPLAI_parpanrlN;

                                                if (Unotran != 0)
                                                    algo.upfun = HPLAI_paupdateNN;
                                                else
                                                    algo.upfun = HPLAI_paupdateNT;
                                            }
                                            else
                                            {
                                                if (rpfa == HPLAI_LEFT_LOOKING)
                                                    algo.pffun = HPLAI_papanllT;
                                                else if (rpfa == HPLAI_CROUT)
                                                    algo.pffun = HPLAI_papancrT;
                                                else
                                                    algo.pffun = HPLAI_papanrlT;

                                                algo.rfact = rpfa = rfaval[irfa];
                                                if (rpfa == HPLAI_LEFT_LOOKING)
                                                    algo.rffun = HPLAI_parpanllT;
                                                else if (rpfa == HPLAI_CROUT)
                                                    algo.rffun = HPLAI_parpancrT;
                                                else
                                                    algo.rffun = HPLAI_parpanrlT;

                                                if (Unotran != 0)
                                                    algo.upfun = HPLAI_paupdateTN;
                                                else
                                                    algo.upfun = HPLAI_paupdateTT;
                                            }

                                            algo.fswap = fswap;
                                            algo.fsthr = tswap;
                                            algo.equil = equil;
                                            algo.align = align;

                                            HPLAI_pdtest(&test, &grid, &algo, nval[in], nbval[inb]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            (void)HPLAI_grid_exit(&grid);
        label_end_of_npqs:;
        }
        /*
 * Print ending messages, close output file, exit.
 */
        if (rank == 0)
        {
            test.ktest = test.kpass + test.kfail + test.kskip;
#ifndef HPL_DETAILED_TIMING
            HPLAI_fprintf(test.outfp, "%s%s\n",
                          "========================================",
                          "========================================");
#else
        if (test.thrsh > HPLAI_rzero)
            HPLAI_fprintf(test.outfp, "%s%s\n",
                          "========================================",
                          "========================================");
#endif

            HPLAI_fprintf(test.outfp, "\n%s %6d %s\n", "Finished", test.ktest,
                          "tests with the following results:");
            if (test.thrsh > HPLAI_rzero)
            {
                HPLAI_fprintf(test.outfp, "         %6d %s\n", test.kpass,
                              "tests completed and passed residual checks,");
                HPLAI_fprintf(test.outfp, "         %6d %s\n", test.kfail,
                              "tests completed and failed residual checks,");
                HPLAI_fprintf(test.outfp, "         %6d %s\n", test.kskip,
                              "tests skipped because of illegal input values.");
            }
            else
            {
                HPLAI_fprintf(test.outfp, "         %6d %s\n", test.kpass,
                              "tests completed without checking,");
                HPLAI_fprintf(test.outfp, "         %6d %s\n", test.kskip,
                              "tests skipped because of illegal input values.");
            }

            HPLAI_fprintf(test.outfp, "%s%s\n",
                          "----------------------------------------",
                          "----------------------------------------");
            HPLAI_fprintf(test.outfp, "\nEnd of Tests.\n");
            HPLAI_fprintf(test.outfp, "%s%s\n",
                          "========================================",
                          "========================================");

            if ((test.outfp != stdout) && (test.outfp != stderr))
                (void)fclose(test.outfp);
        }
        HPLAI_blas_finalize();
        MPI_Finalize();
        exit(0);

        return (0);
        /*
 * End of main
 */
    }

#ifdef __cplusplus
}
#endif
