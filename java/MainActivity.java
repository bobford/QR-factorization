package com.bob.complexqr;

/*
        This is a QR factorization of a complex matrix using the Householder technique.  It is
        a generalization to complex numbers of a previous QR factorization done in the early
        stages of the COVID-19 pandemic.

        This only implements the assembly language version on the Armv8.2-A architecture.  Note that
        v8.2 does not include implementation of complex arithmetic in hardware.

        If the matrix is small enough, e.g. 192x120, the efficiency is quite high, about 50-70%,
        depending on what you believe the speed of the cpu is, whether the advertised rate of
        2.5 GHz or the measured rate of 1.8 GHz.
        The limitations on the matrix size given below were chosen because it would take too long.

        NB  This version runs the R calculation on the user interface because of the interest in
        timing for benchmark performance but runs the calculation of Q (and repeats R), along with
        the verification, as a background thread.
*/

import static android.os.Process.getThreadPriority;
import static android.os.Process.myTid;

import androidx.appcompat.app.AppCompatActivity;

import android.app.Activity;
import android.os.Bundle;
import android.view.View;
import android.widget.ProgressBar;
import android.widget.TextView;

import android.graphics.Typeface;
import android.os.Build;
import android.text.method.ScrollingMovementMethod;
import android.util.Log;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Random;
import java.util.Scanner;


public class MainActivity extends AppCompatActivity {

    public static final String TAG = "bob";

    private static final String     PARAMETER_FILE_NAME     = "params.txt";     // to change number of rows or columns

    private static final String     ARM_ARCHITECTURE        = "arm64-v8a";      // does not run on 32 bit architecture

    private static final int        rowsMax     = 400;
    private static final int        colsMax     = 240;

    TextView tv;
    TextView tv0, tv1, tv2, tv3, tv4, tv5;                  // to print out columns of data, right justified
    ProgressBar spinner;

    ComplexMatrix displ = null;
    double[][] AAr;                                         // elements of the A matrix, random numbers [-1.0, +1.0]
    double[][] AAi;                                         // imaginary part

//      This section for permissions        * * * * * * * * * * * * * * * * * * * * * * * *
    private boolean ready = false;                                                      // *

    private static final String[]   PERMISSIONS = {
//            Manifest.permission.WRITE_EXTERNAL_STORAGE,
//            Manifest.permission.READ_EXTERNAL_STORAGE
    };
//      This section for permissions        * * * * * * * * * * * * * * * * * * * * * * * *

    // Used to load the 'complexqr' library on application startup.
    static {
        System.loadLibrary("complexqr");
    }

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        ready = true;                       // no permissions needed

    }


    public void onResume() {
        super.onResume();
        Log.d(TAG, "in onResume");

        if (ready) {
            initialize();
            process();
        } else {
            Log.d(TAG, "sleeping");
        }
    }

    private void initialize() {
        tv = findViewById(R.id.tv);
        tv.setMovementMethod(new ScrollingMovementMethod());
        tv.setTextIsSelectable(true);
        tv.setText(getResources().getString(R.string.instructions));
        tv0 = findViewById(R.id.tv0);
        tv1 = findViewById(R.id.tv1);
        tv2 = findViewById(R.id.tv2);
        tv3 = findViewById(R.id.tv3);
        tv4 = findViewById(R.id.tv4);
        tv5 = findViewById(R.id.tv5);
        spinner = findViewById(R.id.spinner);
        spinner.setVisibility(View.INVISIBLE);
        Log.d(TAG, "priority 0 = " + Thread.currentThread().getPriority());
        int tid = myTid();
        Log.d(TAG, "priority 0 = " + getThreadPriority(tid));
    }


    public void process() {

        if (!getArchitecture().equals(ARM_ARCHITECTURE)) {
            tv.append(getResources().getString(R.string.architecture));
            return;
        }

        List<String> list = getWatermarkParams();

        //      default values
        int nRows = 192;
        int nCols = 120;

        if (list != null) {             // if PARAMETER_FILE_NAME     = "params.txt" exists

            for (String item : list) {
                String[] params = item.split("=");
                if (params.length > 1) {
                    if ("rows".equals(params[0].trim())) {
                        int n = Integer.parseInt(params[1].trim());
                        if (n > 1 && n <= rowsMax) {
                            nRows = n;
                        } else {
                            tv.append(getResources().getString(R.string.badParam2));
                            return;
                        }
                    }
                    if ("cols".equals(params[0].trim())) {
                        int n = Integer.parseInt(params[1].trim());
                        if (n > 1 && n <= colsMax) {
                            nCols = n;
                        } else {
                            tv.append(getResources().getString(R.string.badParam2));
                            return;
                        }
                    }
                }
            }
        }

        if (nRows < nCols) {
            tv.append(getResources().getString(R.string.badParam));
            return;
        }

        generateAmatrix(nRows, nCols);

        double[] Amatrix = convertToComplexSingle(AAr, AAi, nRows, nCols);

        writeBytesToFile("A_matrix.dat", convertDoubleArrayToByteArray(Amatrix));

        long assemblyTime = executeQR(Amatrix, nRows, nCols);

        doSomeTaskAsync(nRows, nCols, assemblyTime);                          // this will compute R, Q and verify in background
    }

    private long executeQR(double[] A, int nRows, int nCols) {

        Log.d(TAG, "testComplexHouseholder");
        int computeQ = 0;

        double[][] QQ = identity(nRows);

        double[] Q = convertToComplexSingle(QQ, nRows, nRows);
        long start, end, assemblyTime;
        start = System.currentTimeMillis();

        try {
            assemblyTime = complexHouseholder(A, Q, nRows, nCols, computeQ);
        }  catch (Exception e) {
            tv.append(getResources().getString(R.string.badParam3));
            return 0L;
        }

        end = System.currentTimeMillis();
        Log.d(TAG, "execution time java " + (end-start) + " msec and, assembly only, " + assemblyTime + " microseconds");
        Log.d(TAG, "Mflops per second = " + flopCount(nRows, nCols) / assemblyTime);

        return assemblyTime;
    }

    private void generateAmatrix(int nRows, int nCols) {

        tv.append(getResources().getString(R.string.intro1) + nRows + "x" + nCols + getResources().getString(R.string.intro2));

        AAr = generate(nRows, nCols);
        AAi = generate(nRows, nCols);
    }


    private ComplexMatrix evaluateQR(int nRows, int nCols, long assemblyTime) {

/*
        240x240 18.9 msec, R only; 44.5 msec both R and Q
        192x120  3.7 msec, R only; 24.7 msec both R and Q
        200x96   3.3 msec, R only; 16.7 msec both R and Q
        300x60   1.9 msec, R only; 22.7 msec both R and Q
        384x240 45.6 msec, R only; GcSupervisor: GC congestion
 */

        int computeQ = 1;               // this computes Q and R in the background

        double[] B = convertToComplexSingle(AAr, AAi, nRows, nCols);

        double[][] QQ = identity(nRows);

        double[] Q = convertToComplexSingle(QQ, nRows, nRows);

        long start, end, assemblyTimeQ;

        start = System.currentTimeMillis();

        try {
            assemblyTimeQ = complexHouseholder(B, Q, nRows, nCols, computeQ);
        }  catch (Exception e) {
            tv.append(getResources().getString(R.string.badParam3));
            return null;
        }

        end = System.currentTimeMillis();

        Log.d(TAG, "execution time R & Q, assembly, " + assemblyTimeQ + " microseconds");

        ComplexMatrix displ = verify(AAr, AAi, B, Q);
        tv.append(getResources().getString(R.string.arm) + assemblyTime + getResources().getString(R.string.micro));
        tv.append(getResources().getString(R.string.forQ) + (end-start) + getResources().getString(R.string.sec));

        return displ;
    }

    private ComplexMatrix verify(double[][] ar, double[][] ai, double[] A, double[] qq) {

        writeBytesToFile("R_matrix.dat", convertDoubleArrayToByteArray(A));
        writeBytesToFile("Q_matrix.dat", convertDoubleArrayToByteArray(qq));

        ComplexMatrix matrix = new ComplexMatrix(ar, ai);

        int m = matrix.getNrow();
        int n = matrix.getNcol();

        double[][] Qreal = new double[m][m];
        double[][] Qimag = new double[m][m];
        double[][] Rreal = new double[m][n];
        double[][] Rimag = new double[m][n];
        for (int i=0; i<m; i++) {                               // row
            for (int j = 0; j < m; j++) {                       // column
                Qreal[i][j] = qq[2*m * i + 2 * j];
                Qimag[i][j] = qq[2*m * i + 2 * j + 1];
            }
            for (int j = 0; j < n; j++) {
                Rreal[i][j] = A[2*n * i + 2 * j];
                Rimag[i][j] = A[2*n * i + 2 * j + 1];
            }
        }

        ComplexMatrix Qmat = new ComplexMatrix(Qreal, Qimag);
        ComplexMatrix Rmat = new ComplexMatrix(Rreal, Rimag);

        ComplexMatrix Qcon = Qmat.conjugate();
        ComplexMatrix Qher = Qcon.transpose();
        ComplexMatrix Aest = Qher.times(Rmat);
        ComplexMatrix Qeye = Qher.times(Qmat);
        ComplexMatrix Err  = Aest.minus(matrix);

        Log.d(TAG, "norm of error is " + Err.norm());

        Complex det = Qeye.determinant();
        Log.d(TAG, "determinent of Qeye is " + det.toString());

        if ((Err.norm() < 1.E-8) && (Math.abs(det.getReal() - 1.0) < 1.E-10) && (Math.abs(det.getImag()) < 1.E-10)) {
            tv.append(getResources().getString(R.string.verify));
        }
        return Rmat;
    }

    public void display(ComplexMatrix Rmat) {

        if (Rmat.getNcol() < 3 || Rmat.getNrow() < 4) {
            return;
        }
        tv.append("\n");
        tv0.setText("Real\n");
        tv1.setText("Imag\n");
        tv2.setText("Real\n");
        tv3.setText("Imag\n");
        tv4.setText("Real\n");
        tv5.setText("Imag\n");
        Typeface font = Typeface.createFromAsset(getAssets(), "fonts/consolas.ttf");
        tv.setTypeface(font);
        for (int i=0; i<4; i++) {
            double[] x = new double[6];
            x[0] = Rmat.getElementCopy(i,0).getReal();
            x[1] = Rmat.getElementCopy(i,0).getImag();
            x[2] = Rmat.getElementCopy(i,1).getReal();
            x[3] = Rmat.getElementCopy(i,1).getImag();
            x[4] = Rmat.getElementCopy(i,2).getReal();
            x[5] = Rmat.getElementCopy(i,2).getImag();

            for (int j=0; j<6; j++) {
                switch (j) {
                    case 0:
                        tv0.append(String.format(Locale.US, "%+7.2f", x[0]) + "\n");
                        break;
                    case 1:
                        tv1.append(String.format(Locale.US, "%+7.2fi", x[1]) + "\n");
                        break;
                    case 2:
                        tv2.append(String.format(Locale.US, "%+7.2f", x[2]) + "\n");
                        break;
                    case 3:
                        tv3.append(String.format(Locale.US, "%+7.2fi", x[3]) + "\n");
                        break;
                    case 4:
                        tv4.append(String.format(Locale.US, "%+7.2f", x[4]) + "\n");
                        break;
                    case 5:
                        tv5.append(String.format(Locale.US, "%+7.2fi", x[5]) + "\n");
                        break;
                    default:
                        break;
                }
            }
        }
        tv.setTypeface(Typeface.DEFAULT);

    }

    public double[][] generate(int nRows, int nCols) {
//          generates a random 2-d array
        Random r = new Random();
        double min = -1.0;
        double max = 1.;

        double[][] A = new double[nRows][nCols];
        for (int i=0; i<nRows; i++) {
            for (int j=0; j<nCols; j++) {
                A[i][j] = min + r.nextDouble() * (max - min);
            }
        }
        return A;
    }


    public double[] convertToComplexSingle(double[][] A, int m, int n) {
        double[] a = new double[2*m*n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                a[2*(n * i + j)] = A[i][j];
                a[2*(n * i + j) + 1] = 0.0;
            }
        }
        return a;
    }

    public double[] convertToComplexSingle(double[][] A, double[][] B, int m, int n) {
        double[] a = new double[2*m*n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                a[2*(n * i + j)] = A[i][j];
                a[2*(n * i + j) + 1] = B[i][j];
            }
        }
        return a;
    }

    private double[][] identity(int n) {
        double[][] x = new double[n][n];
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                if (i == j)
                    x[i][j] = 1.;
            }
        }
        return x;
    }


    public List<String> getWatermarkParams() {
        /*
            Format of the parameter file is:
            String  =  Integer
            where "String" is the name of the parameter and Integer is the value, separated by an = character
            e.g.:
            rows = 384
         */
        List<String> list = new ArrayList<>();

        File parameters;
        try {
            File dir = getExternalFilesDir(null);
//            Log.d(TAG, dir.getAbsolutePath());             //  /storage/emulated/0/Android/data/com.bob.complexqr/files/
            parameters = new File(dir, PARAMETER_FILE_NAME);
        } catch (NullPointerException e) {
            return null;
        }
        if (parameters.exists()) {                     // got it
            // read it
            Log.d(TAG, "parameters");

            try {
                BufferedReader br = new BufferedReader(new FileReader(parameters));
                String line;
                while ((line = br.readLine()) != null) {
                    Log.d(TAG, line);

                    if (line.indexOf("//") != 0 && line.length() != 0) {
                        list.add(line);
                    }
                }
            } catch (IOException e) {
                Log.d(TAG, "impossible error");
            }
        } else {
            Log.d(TAG, "parameter file missing, using defaults");
            return null;
        }

        return list;
    }


    private String getArchitecture() {

        Log.d(TAG, "OS.ARCH : " + System.getProperty("os.arch"));

        Log.d(TAG, "SUPPORTED_ABIS : " + Arrays.toString(Build.SUPPORTED_ABIS));
        Log.d(TAG, "SUPPORTED_32_BIT_ABIS : " + Arrays.toString(Build.SUPPORTED_32_BIT_ABIS));
        Log.d(TAG, "SUPPORTED_64_BIT_ABIS : " + Arrays.toString(Build.SUPPORTED_64_BIT_ABIS));

        return Build.SUPPORTED_64_BIT_ABIS[0];
    }


    public String getAndroidVersion() {
//        String release = Build.VERSION.RELEASE;
        int sdkVersion = Build.VERSION.SDK_INT;
        return "Android SDK: " + sdkVersion;
    }

//      https://stackoverflow.com/questions/1995439/get-android-phone-model-programmatically

    public String getDeviceName() {
        String manufacturer = Build.MANUFACTURER;
        String model = Build.MODEL;
        if (model.toLowerCase().startsWith(manufacturer.toLowerCase())) {
            return capitalize(model);
        } else {
            return capitalize(manufacturer) + " " + model;
        }
    }


    private String capitalize(String s) {
        if (s == null || s.length() == 0) {
            return "";
        }
        char first = s.charAt(0);
        if (Character.isUpperCase(first)) {
            return s;
        } else {
            return Character.toUpperCase(first) + s.substring(1);
        }
    }

    public int flopCount(int m, int n) {
        int k, count;
        k = (m * n * n - n * n * n / 3) / 2;      // for the multiplies and the same for the adds
        count = 8 * k;                      // complex dot product have 8 floating point operations each
        count += 6 * k;                     // plus the complex outer product
        count += 2 * k;                     // plus the complex addition
        return count;
    }

    public void writeBytesToFile(String name, byte[] data) {
        String dir = getExternalFilesDir(null).getPath() + File.separator + name;
        Log.d(TAG, dir);
        File file = new File(dir);
        try {
            OutputStream os = new FileOutputStream(file);

            os.write(data);

            // Close the file
            os.close();
        } catch (Exception e) {
            Log.d(TAG,"Exception: " + e);
        }
    }
    private byte[] convertDoubleArrayToByteArray(double[] data) {
        if (data == null) return null;
        // ----------
        byte[] byts = new byte[data.length * Double.BYTES];
        for (int i = 0; i < data.length; i++)
            System.arraycopy(convertDoubleToByteArray(data[i]), 0, byts, i * Double.BYTES, Double.BYTES);
        return byts;
    }
    private byte [] convertDoubleToByteArray(double number) {
        ByteBuffer byteBuffer = ByteBuffer.allocate(Double.BYTES);
        byteBuffer.putDouble(number);
        return byteBuffer.array();
    }
    public double toDouble(byte[] bytes) {
        return ByteBuffer.wrap(bytes).getDouble();
    }
    public byte[] toByteArray(double value) {
        byte[] bytes = new byte[8];
        ByteBuffer.wrap(bytes).putDouble(value);
        return bytes;
    }

    /**
     * https://android.stackexchange.com/questions/19810/how-can-i-determine-max-cpu-speed-at-runtime
     *
     * Get max cpu rate.
     *
     * This works by examining the list of CPU frequencies in the pseudo file
     * "/sys/devices/system/cpu/cpu0/cpufreq/stats/time_in_state" and how much time has been spent
     * in each. It finds the highest non-zero time and assumes that is the maximum frequency (note
     * that sometimes frequencies higher than that which was designed can be reported.) So it is not
     * impossible that this method will return an incorrect CPU frequency.
     *
     * Also note that (obviously) this will not reflect different CPU cores with different
     * maximum speeds.
     *
     * @return cpu frequency in MHz
     */
    public static int getMaxCPUFreqMHz() {

        int maxFreq = -1;
        try {

            RandomAccessFile reader = new RandomAccessFile( "/sys/devices/system/cpu/cpu0/cpufreq/stats/time_in_state", "r" );

//            RandomAccessFile reader = new RandomAccessFile( "/sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_max_freq", "r" ); max freq direct

            while (true) {
                String line = reader.readLine();
                if ( null == line ) {
                    break;
                }
                String[] splits = line.split( "\\s+" );
                assert ( splits.length == 2 );
                int timeInState = Integer.parseInt( splits[1] );
                if ( timeInState > 0 ) {
                    int freq = Integer.parseInt( splits[0] ) / 1000;
                    if ( freq > maxFreq ) {
                        maxFreq = freq;
                    }
                }
            }

        } catch ( IOException ex ) {
            Log.d(TAG, ex.toString());
        }

        return maxFreq;
    }

    /*
    http://www.java2s.com/Code/Android/Hardware/GetCPUFrequencyCurrent.htm
    gives same result as other
     */
    public static int getCPUFrequencyCurrent() throws Exception {
        return readSystemFileAsInt();
    }

    private static int readSystemFileAsInt() throws Exception {
        InputStream in;
        try {
            final Process process = new ProcessBuilder(new String[] { "/system/bin/cat", "/sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq"}).start();

            in = process.getInputStream();
            final String content = readFully(in);
            return Integer.parseInt(content);
        } catch (final Exception e) {
            throw new Exception(e);
        }
    }
    public static String readFully(final InputStream pInputStream) {
        final StringBuilder sb = new StringBuilder();
        final Scanner sc = new Scanner(pInputStream);
        while(sc.hasNextLine()) {
            sb.append(sc.nextLine());
        }
        return sb.toString();
    }

    private void postProcess(int nRows, int nCols, long assemblyTime) {
        tv.append(getResources().getString(R.string.fini));

        String answer = getDeviceName();
        tv.append("\n" + answer + "\n");
        String version = getAndroidVersion();
        tv.append(version + "\n\n");
        tv.append(getResources().getString(R.string.epilogue));
        Log.d(TAG, "all done");

        if (displ != null) {
            display(displ);
        } else {
            Log.d(TAG, "displ is null");
        }
        int freq = getMaxCPUFreqMHz();
        Log.d(TAG, "getMaxCPUFreqMHz = " + freq);
        try {
            freq = getCPUFrequencyCurrent() / 1000;
        } catch (Exception e) {
            Log.d(TAG, "getCPUFrequencyCurrent " + e.toString());
        }
        Log.d(TAG, "getCPUFrequencyCurrent = " + freq);
        int eff = (int) (100 * flopCount(nRows, nCols) / assemblyTime / (4 * freq));
        Log.d(TAG, "efficiency = " + eff + "% based on clock frequency of " + freq);
    }

    private void doSomeTaskAsync(int nRows, int nCols, long assemblyTime) {
        Activity currActivity = MainActivity.this;
        spinner.setVisibility(View.VISIBLE);
        new BackgroundTask(currActivity) {

            String allDone;
            @Override
            public void doInBackground() {

                displ = evaluateQR(nRows, nCols, assemblyTime);

                allDone = "Thank you for your attention.\n";

            }

            @Override
            public void onPostExecute() {
                spinner.setVisibility(View.INVISIBLE);
                postProcess(nRows, nCols, assemblyTime);
            }
        }.execute();
    }

    /*
    private void doSomeTaskAsync(int nRows, int nCols, long assemblyTime) {
        Activity currActivity = MainActivity.this;
        spinner.setVisibility(View.VISIBLE);
        new BackgroundTask(currActivity) {

            String allDone;
            @Override
            public void doInBackground() {
                
                displ = evaluateQR(nRows, nCols, assemblyTime);

                allDone = "Thank you for your attention.\n";

            }

            @Override
            public void onPostExecute() {
                spinner.setVisibility(View.INVISIBLE);
                // here is result part same
                // same like post execute
                // UI Thread(update your UI widget)
                tv.append(getResources().getString(R.string.fini));

                String answer = getDeviceName();
                tv.append("\n" + answer + "\n");
                String version = getAndroidVersion();
                tv.append(version + "\n\n");
                tv.append(getResources().getString(R.string.epilogue));
                Log.d(TAG, "all done");

                if (displ != null) {
                    display(displ);
                } else {
                    Log.d(TAG, "displ is null");
                }
                int freq = getMaxCPUFreqMHz();
                Log.d(TAG, "getMaxCPUFreqMHz = " + freq);
                try {
                    freq = getCPUFrequencyCurrent() / 1000;
                } catch (Exception e) {
                    Log.d(TAG, "getCPUFrequencyCurrent " + e.toString());
                }
                Log.d(TAG, "getCPUFrequencyCurrent = " + freq);
                int eff = (int) (100 * flopCount(nRows, nCols) / assemblyTime / (4 * freq));
                Log.d(TAG, "efficiency = " + eff + "% based on clock frequency of " + freq);
            }
        }.execute();
    }
*/



    /**
     * A native method that is implemented by the 'complexqr' native library,
     * which is packaged with this application.
     */
    public native long complexHouseholder(double[] a, double[] qq, int m, int n, int q);
}