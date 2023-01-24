package com.bob.complexqr;

public class Matrix {

    private int numberOfRows = 0;                   // number of rows
    private int numberOfColumns = 0;                // number of columns
    private double matrix[][] = null; 	            // 2-D  Matrix
    private double hessenberg[][] = null; 	        // 2-D  Hessenberg equivalent
    private boolean hessenbergDone = false;         // = true when Hessenberg matrix calculated
    private int permutationIndex[] = null;          // row permutation index
    private double rowSwapIndex = 1.0D;             // row swap index
    private double[] eigenValues = null;            // eigen values of the matrix
    private double[][] eigenVector = null;          // eigen vectors of the matrix
    private double[] sortedEigenValues = null;      // eigen values of the matrix sorted into descending order
    private double[][] sortedEigenVector = null;    // eigen vectors of the matrix sorted to matching descending eigen value order
    private int numberOfRotations = 0;              // number of rotations in Jacobi transformation
    private int[] eigenIndices = null;              // indices of the eigen values before sorting into descending order
    private int maximumJacobiIterations = 100;      // maximum number of Jacobi iterations
    private boolean eigenDone = false;              // = true when eigen values and vectors calculated
    private boolean matrixCheck = true;             // check on matrix status
    // true - no problems encountered in LU decomposition
    // false - attempted a LU decomposition on a singular matrix

    private boolean supressErrorMessage = false;    // true - LU decompostion failure message supressed

    private double tiny = 1.0e-100;                 // small number replacing zero in LU decomposition

    // CONSTRUCTORS
    // Construct a numberOfRows x numberOfColumns matrix of variables all equal to zero
    public Matrix(int numberOfRows, int numberOfColumns){
        this.numberOfRows = numberOfRows;
        this.numberOfColumns = numberOfColumns;
        this.matrix = new double[numberOfRows][numberOfColumns];
        this.permutationIndex = new int[numberOfRows];
        for(int i=0;i<numberOfRows;i++)this.permutationIndex[i]=i;
    }

    // Construct a numberOfRows x numberOfColumns matrix of variables all equal to the number const
    public Matrix(int numberOfRows, int numberOfColumns, double constant){
        this.numberOfRows = numberOfRows;
        this.numberOfColumns = numberOfColumns;
        this.matrix = new double[numberOfRows][numberOfColumns];
        for(int i=0;i<numberOfRows;i++){
            for(int j=0;j<numberOfColumns;j++)this.matrix[i][j]=constant;
        }
        this.permutationIndex = new int[numberOfRows];
        for(int i=0;i<numberOfRows;i++)this.permutationIndex[i]=i;
    }

    // Construct matrix with a copy of an existing numberOfRows x numberOfColumns 2-D array of variables
    public Matrix(double[][] twoD){
        this.numberOfRows = twoD.length;
        this.numberOfColumns = twoD[0].length;
        this.matrix = new double[this.numberOfRows][this.numberOfColumns];
        for(int i=0; i<numberOfRows; i++){
            if(twoD[i].length!=numberOfColumns)throw new IllegalArgumentException("All rows must have the same length");
            for(int j=0; j<numberOfColumns; j++){
                this.matrix[i][j]=twoD[i][j];
            }
        }
        this.permutationIndex = new int[numberOfRows];
        for(int i=0;i<numberOfRows;i++)this.permutationIndex[i]=i;
    }

    // Construct matrix with a copy of an existing numberOfRows x numberOfColumns 2-D array of floats
    public Matrix(float[][] twoD){
        this.numberOfRows = twoD.length;
        this.numberOfColumns = twoD[0].length;
        for(int i=1; i<numberOfRows; i++){
            if(twoD[i].length!=numberOfColumns)throw new IllegalArgumentException("All rows must have the same length");
        }
        this.matrix = new double[this.numberOfRows][this.numberOfColumns];
        for(int i=0; i<numberOfRows; i++){
            for(int j=0; j<numberOfColumns; j++){
                this.matrix[i][j] = (double)twoD[i][j];
            }
        }
        this.permutationIndex = new int[numberOfRows];
        for(int i=0;i<numberOfRows;i++)this.permutationIndex[i]=i;
    }

    // Construct matrix with a copy of an existing numberOfRows x numberOfColumns 2-D array of longs
    public Matrix(long[][] twoD){
        this.numberOfRows = twoD.length;
        this.numberOfColumns = twoD[0].length;
        for(int i=1; i<numberOfRows; i++){
            if(twoD[i].length!=numberOfColumns)throw new IllegalArgumentException("All rows must have the same length");
        }
        this.matrix = new double[this.numberOfRows][this.numberOfColumns];
        for(int i=0; i<numberOfRows; i++){
            for(int j=0; j<numberOfColumns; j++){
                this.matrix[i][j] = (double)twoD[i][j];
            }
        }
        this.permutationIndex = new int[numberOfRows];
        for(int i=0;i<numberOfRows;i++)this.permutationIndex[i]=i;
    }


    // Construct matrix with a copy of an existing numberOfRows x numberOfColumns 2-D array of ints
    public Matrix(int[][] twoD){
        this.numberOfRows = twoD.length;
        this.numberOfColumns = twoD[0].length;
        for(int i=1; i<numberOfRows; i++){
            if(twoD[i].length!=numberOfColumns)throw new IllegalArgumentException("All rows must have the same length");
        }
        this.matrix = new double[this.numberOfRows][this.numberOfColumns];
        for(int i=0; i<numberOfRows; i++){
            for(int j=0; j<numberOfColumns; j++){
                this.matrix[i][j] = (double)twoD[i][j];
            }
        }
        this.permutationIndex = new int[numberOfRows];
        for(int i=0;i<numberOfRows;i++)this.permutationIndex[i]=i;
    }


    // SPECIAL MATRICES
    // Construct an identity matrix
    public static Matrix identityMatrix(int numberOfRows){
        Matrix special = new Matrix(numberOfRows, numberOfRows);
        for(int i=0; i<numberOfRows; i++){
            special.matrix[i][i]=1.0;
        }
        return special;
    }



    // GET VALUES
    // Return the number of rows
    public int getNumberOfRows(){
        return this.numberOfRows;
    }

    // Return the number of rows
    public int getNrow(){
        return this.numberOfRows;
    }

    // Return the number of columns
    public int getNumberOfColumns(){
        return this.numberOfColumns;
    }

    // Return the number of columns
    public int getNcol(){
        return this.numberOfColumns;
    }


    // Return a reference to the internal 2-D array
    public double[][] getArrayReference(){
        return this.matrix;
    }

    // Return a reference to the internal 2-D array
    // included for backward compatibility with incorrect earlier documentation
    public double[][] getArrayPointer(){
        return this.matrix;
    }

    // Return a copy of the internal 2-D array
    public double[][] getArrayCopy(){
        double[][] c = new double[this.numberOfRows][this.numberOfColumns];
        for(int i=0; i<numberOfRows; i++){
            for(int j=0; j<numberOfColumns; j++){
                c[i][j]=this.matrix[i][j];
            }
        }
        return c;
    }

    // Return a copy of a column
    public double[] getColumnCopy(int ii){
        if(ii>=this.numberOfColumns)throw new IllegalArgumentException("Column index, " + ii + ", must be less than the number of columns, " + this.numberOfColumns);
        if(ii<0)throw new IllegalArgumentException("column index, " + ii + ", must be zero or positive");
        double[] col = new double[this.numberOfRows];
        for(int i=0; i<numberOfRows; i++){
            col[i]=this.matrix[i][ii];
        }
        return col;
    }


    // Return  a single element of the internal 2-D array
    public double getElement(int i, int j){
        return this.matrix[i][j];
    }

    // Return a single element of the internal 2-D array
    // included for backward compatibility with incorrect earlier documentation
    public double getElementCopy(int i, int j){
        return this.matrix[i][j];
    }

    // Return a single element of the internal 2-D array
    // included for backward compatibility with incorrect earlier documentation
    public double getElementPointer(int i, int j){
        return this.matrix[i][j];
    }


    // Return a copy of the permutation index array
    public int[]  getIndexCopy(){
        int[] indcopy = new int[this.numberOfRows];
        for(int i=0; i<this.numberOfRows; i++){
            indcopy[i]=this.permutationIndex[i];
        }
        return indcopy;
    }

    // Return the row swap index
    public double getSwap(){
        return this.rowSwapIndex;
    }


}
