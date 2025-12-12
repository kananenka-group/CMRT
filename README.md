# CMRT

**Simulation of Excitonic Dynamics Using Coherent Modified Redfield
Theory**

Author: Aarti Sindhu, University of Delaware

This repository contains Fortran code for simulating dynamics using
**Coherent Modified Redfield Theory (CMRT)**.

## Repository Structure

-   **module:** `mod_spectra_cmrt.f90`
-   **input file:** `spectra_cmrt.inp`
    -   Contains all parameters except bath parameters (found in the
        `setup_parameters` subroutine).
-   **scripts:**
    -   `makefile` --- build configuration
    -   `run` --- execution script
    -   `submit.qs` --- job submission script

## How to Read the Code

To understand the flow of the program:

1.  Search for the **`setup`** routines.
2.  Then locate and follow the **`main`** subroutine.

## How to Run

### Step 1: Initial Calculation

1.  Set `iflow = 2` in `spectra_cmrt.inp`.

2.  Run:

    ``` bash
    ./run
    ```

3.  After completion, the following files appear in each generated
    folder\
    (number of folders = value of `ipara` in the input file):

    -   `response.out`
    -   `se.out`
    -   `gsb.out`
    -   `esa.out`

### Step 2: Combine Results

4.  Set `iflow = 3` in `spectra_cmrt.inp`.

5.  Run:

    ``` bash
    ./run
    ```

    This merges all `response.out`, `se.out`, `gsb.out`, and `esa.out`
    into the working directory.

### Step 3: Fourier Transform

6.  Copy `response.out` to the `fourier_trans` folder.

7.  Enter the folder:

    ``` bash
    cd fourier_trans
    ```

8.  Compile and run:

    ``` bash
    make
    ./aout
    ```

This creates `spectra.out`, which contains the 2D spectrum:

  Column   Meaning
  -------- ----------------------------
  1        w1
  2        w2
  3        Spectrum (arbitrary units)

## Additional Options

-   Generate **rephasing** or **non-rephasing** spectra by editing the
    `rep/nonrep` flag in the `fft_2d` subroutine.
-   Generate **SE**, **GSB**, or **ESA** spectra by changing the input
    file name in the `main` subroutine.
