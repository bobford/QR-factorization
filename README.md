# QR-factorization
This is a QR factorization of a complex matrix using assembly language on the Armv8.2-A architecture.

        This is a QR factorization of a complex matrix using the Householder technique.  It is
        a generalization to complex numbers of a previous QR factorization done in the early
        stages of the COVID-19 pandemic.

        This only implements the assembly language version on the Armv8.2-A architecture.  Note that
        v8.2 does not include implementation of complex arithmetic in hardware.

        If the matrix is small enough, e.g. 192x120, the efficiency is quite high, about 50-70%,
        depending on what you believe the speed of the cpu is, whether the advertised rate of
        2.5 GHz or the measured rate of 1.8 GHz.

        NB  This version runs the R calculation on the user interface because of the interest in
        timing for benchmark performance but runs the calculation of Q (and repeats R), along with
        the verification, as a background thread.
