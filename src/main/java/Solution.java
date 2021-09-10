import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class Solution {
    /*
    * Use the markov chain theroem to solve this problem
    * 1. Identify terminate states and non-terminate states
    * 2. Calculate the I, R and Q matrices based on the non-terminate states
    * 3. Calculate the F matrix (i.e., F = (I - Q)^-1)
    * 4. Calculate the FR matrix
    * 5. The first row of the FR matrix are the probabilities
     */

    public static int[] solution(int[][] m) {
        if(m[0][0]==0 && m.length==1) {
            return new int[] {1, 1};
        }

        //categorize terminate and non-terminate states and collect the denominators for non-terminate states
        List<Integer> terminateStates = new ArrayList<>();
        List<Integer> nonTerminateStates = new ArrayList<>();
        List<Integer> denominators = new ArrayList<>();
        categorizeStates(m, terminateStates, nonTerminateStates, denominators);

        //if s0 is a terminate state, return [1, 0, ... , 0, 1] since it will not get to other states
        if (terminateStates.get(0) == 0) {
            int[] firstRowTerminateResult = new int[m.length + 1];
            firstRowTerminateResult[0] = 1;
            firstRowTerminateResult[firstRowTerminateResult.length - 1] = 1;
            return firstRowTerminateResult;
        }

        //calculate the I, Q and R matrices
        Matrix iMatrix = buildIMatrix(nonTerminateStates.size());
        Matrix qMatrix = buildQMatrix(m, nonTerminateStates, denominators);
        Matrix rMatrix = buildRMatrix(m, nonTerminateStates, denominators);

        //calculate the F matrix (i.e., F = (I - Q)^-1) -> FR matrix
        Matrix iMinusQ = iMatrix.minus(qMatrix);
        Matrix fMatrix = iMinusQ.getInverseMatrix();
        Matrix frMatrix = fMatrix.multiply(rMatrix);

        //the first row of the FR matrix are the probabilities
        List<Fraction> probabilities = frMatrix.getRow(0);
        List<Fraction> numerators = new ArrayList<>();
        List<Integer> denominatorList = new ArrayList<>();

        // find the numerators and denominators
        for (int i = 0; i < probabilities.size(); i++) {
            denominatorList.add(probabilities.get(i).getDenominator());
            numerators.add(probabilities.get(i));
        }
        int lcm = findLCM(denominatorList);
        int[] result = new int[probabilities.size() + 1];
        for (int i = 0; i < result.length - 1; i++) {
            numerators.set(i, numerators.get(i).multiply(new Fraction(lcm)));
            result[i] = numerators.get(i).getNumerator();
        }
        result[result.length - 1] = lcm;
        return result;
    }

    private static void categorizeStates(int[][] m, List<Integer> terminateStates, List<Integer> nonTerminateStates, List<Integer> denominators) {
        for (int i = 0; i < m.length; i++) {
            boolean isTerminateState = true;
            int denominator = 0;
            for (int j = 0; j < m[i].length; j++) {
                denominator +=m[i][j];
                if (m[i][j] != 0) {
                    isTerminateState = false;
                }
            }
            if (isTerminateState) {
                terminateStates.add(i);
            } else {
                nonTerminateStates.add(i);
                denominators.add(denominator);
            }
        }
    }

    private static Matrix buildIMatrix(int size) {
        List<List<Fraction>> iList = new ArrayList<>();
        for (int i = 0; i < size; i++) {
            List<Fraction> iRow = new ArrayList<>();
            for (int j = 0; j < size; j++) {
                if (i==j) {
                    iRow.add(new Fraction(1));
                } else {
                    iRow.add(new Fraction(0));
                }
            }
            iList.add(iRow);
        }
        return new Matrix(iList, size, size);
    }

    private static Matrix buildQMatrix(int[][] m, List<Integer> nonTerminateStates, List<Integer> denominators) {
        int size = nonTerminateStates.size();
        List<List<Fraction>> qList = new ArrayList<>();
        for (int i = 0; i < size; i++) {
            List<Fraction> qRow = new ArrayList<>();
            for (int j = 0; j < size; j++) {
                qRow.add(new Fraction(m[nonTerminateStates.get(i)][nonTerminateStates.get(j)], denominators.get(i)));
            }
            qList.add(qRow);
        }
        return new Matrix(qList, size, size);
    }

    private static Matrix buildRMatrix(int[][] m, List<Integer> nonTerminateStates, List<Integer> denominators) {
        int numOfRows = nonTerminateStates.size();
        int numOfCols = m.length - numOfRows;
        List<List<Fraction>> rList = new ArrayList<>();
        for (int i = 0; i < numOfRows; i++)
        {
            ArrayList<Fraction> RRow = new ArrayList<>();
            for (int j = 0; j < numOfCols; j++)
            {
                RRow.add(new Fraction(m[nonTerminateStates.get(i)][j + numOfRows], denominators.get(i)));
            }
            rList.add(RRow);
        }

        return new Matrix(rList, numOfRows, numOfCols);
    }

    private static int findLCM(List<Integer> list) {
        Collections.sort(list);
        int lcm = 1;
        for (int n : list) {
            lcm = lcm(lcm, n);
        }
        return lcm;
    }

    private static int lcm(int number1, int number2) {
        if (number1 == 0 || number2 == 0) {
            return 0;
        }
        int absNumber1 = Math.abs(number1);
        int absNumber2 = Math.abs(number2);
        int absHigherNumber = Math.max(absNumber1, absNumber2);
        int absLowerNumber = Math.min(absNumber1, absNumber2);
        int lcm = absHigherNumber;
        while (lcm % absLowerNumber != 0) {
            lcm += absHigherNumber;
        }
        return lcm;
    }

    private static class Matrix {

        private final int M;
        private final int N;
        private final Fraction det;
        private final List<List<Fraction>> matrix;
        private final List<List<Fraction>> inverseMatrix;

        public Matrix(List<List<Fraction>> matrix, int m, int n) {
            this.matrix = matrix;
            this.M = m;
            this.N = n;
            this.det = this.determinant(matrix, n);
            this.inverseMatrix = this.inverse();
        }

        private void getCofactor(List<List<Fraction>> matrix, List<List<Fraction>> tempMat, int p, int q, int n) {
            int i = 0;
            int j = 0;
            for (int row = 0; row < n; row++) {
                for (int col = 0; col < n; col++) {
                    if (row != p && col != q) {
                        tempMat.get(i).set(j++, matrix.get(row).get(col));
                        if (j == n - 1) {
                            j = 0;
                            i++;
                        }
                    }
                }
            }
        }

        private Fraction determinant(List<List<Fraction>> matrix, int n) {
            Fraction ans = new Fraction(0, 1);
            if (this.M != this.N) {
                return ans;
            }
            if (n == 1) {
                return matrix.get(0).get(0);
            }
            List<List<Fraction>> tempMat = new ArrayList<>();
            // Init 2d fraction arraylist
            for (int i = 0; i < this.M; i++) {
                List<Fraction> tempMatRow = new ArrayList<>();
                for (int j = 0; j < this.N; j++) {
                    tempMatRow.add(new Fraction(0, 1));
                }
                tempMat.add(tempMatRow);
            }

            int sign = 1;
            Fraction signFraction = new Fraction(sign, 1);
            for (int k = 0; k < n; k++) {
                this.getCofactor(matrix, tempMat, 0, k, n);
                ans = ans.plus(signFraction.multiply(matrix.get(0).get(k).multiply(determinant(tempMat, n - 1))));
                sign = -sign;
                signFraction = new Fraction(sign, 1);
            }
            return ans;
        }

        private void adjoint(List<List<Fraction>> matrix, List<List<Fraction>> adj) {
            if (this.N == 1) {
                adj.get(0).set(0, new Fraction(1, 1));
                return;
            }
            int sign = 1;

            List<List<Fraction>> tempMat = new ArrayList<>();
            for (int i = 0; i < this.N; i++) {
                ArrayList<Fraction> tempMatRow = new ArrayList<>();
                for (int j = 0; j < this.N; j++) {
                    tempMatRow.add(new Fraction(0, 1));
                }
                tempMat.add(tempMatRow);
            }

            for (int p = 0; p < this.N; p++) {
                for (int q = 0; q < this.N; q++) {
                    this.getCofactor(matrix, tempMat, p, q, this.N);
                    sign = ((p + q) % 2 == 0) ? 1 : -1;
                    Fraction signFraction = new Fraction(sign, 1);
                    adj.get(q).set(p, signFraction.multiply((this.determinant(tempMat, this.N - 1))));
                }
            }
        }

        private List<List<Fraction>> inverse() {
            List<List<Fraction>> inv = new ArrayList<>();
            for (int i = 0; i < this.M; i++) {
                ArrayList<Fraction> invRow = new ArrayList<>();
                for (int j = 0; j < this.N; j++) {
                    invRow.add(new Fraction(0, 1));
                }
                inv.add(invRow);
            }

            if (this.det.equals(new Fraction(0))) {
                return inv;
            }

            List<List<Fraction>> adj = new ArrayList<>();
            for (int i = 0; i < this.M; i++) {
                List<Fraction> adjRow = new ArrayList<>();
                for (int j = 0; j < this.N; j++) {
                    adjRow.add(new Fraction(0, 1));
                }
                adj.add(adjRow);
            }

            adjoint(this.matrix, adj);
            for (int p = 0; p < this.N; p++) {
                for (int q = 0; q < this.N; q++) {
                    Fraction temp = adj.get(p).get(q).dividedBy(this.det);
                    inv.get(p).set(q, temp);
                }
            }
            return inv;
        }

        public Matrix getInverseMatrix() {
            if (this.M != this.N) {
                return null;
            }
            return new Matrix(this.inverseMatrix, this.M, this.N);
        }

        public Fraction getElement(int m, int n) {
            return this.matrix.get(m).get(n);
        }

        public List<Fraction> getRow(int m) {
            if (m <= this.M) {
                return this.matrix.get(m);
            }
            return new ArrayList<>();
        }

        public Matrix minus(Matrix mat) {
            int M_m = mat.getDimension()[0];
            int N_m = mat.getDimension()[1];
            if (this.M != M_m || this.N != N_m) {
                return mat;
            } else {
                List<List<Fraction>> difference = new ArrayList<>();
                for (int i = 0; i < this.M; i++) {
                    List<Fraction> differenceRow = new ArrayList<>();
                    for (int j = 0; j < this.N; j++) {
                        differenceRow.add(new Fraction(0, 1));
                    }
                    difference.add(differenceRow);
                }
                for (int i = 0; i < this.M; i++) {
                    for (int j = 0; j < this.N; j++) {
                        difference.get(i).set(j, this.matrix.get(i).get(j).minus(mat.getElement(i, j)));
                    }
                }
                return new Matrix(difference, this.M, this.N);
            }
        }

        public Matrix multiply(Matrix matrix) {
            int M_m = matrix.getDimension()[0];
            int p_m = matrix.getDimension()[1];
            if (this.N != M_m) {
                return matrix;
            } else {
                List<List<Fraction>> product = new ArrayList<>();
                for (int i = 0; i < this.M; i++) {
                    List<Fraction> productRow = new ArrayList<>();
                    for (int j = 0; j < p_m; j++) {
                        productRow.add(new Fraction(0, 1));
                    }
                    product.add(productRow);
                }
                for (int i = 0; i < this.M; i++) {
                    for (int j = 0; j < p_m; j++) {
                        for (int k = 0; k < this.N; k++) {
                            Fraction temp = product.get(i).get(j);
                            product.get(i).set(j, temp.plus(this.matrix.get(i).get(k).multiply(matrix.getElement(k, j))));
                        }
                    }
                }
                return new Matrix(product, this.M, p_m);
            }

        }

        public int[] getDimension() {
            return new int[] { this.M, this.N };
        }

    }

    private static class Fraction {

        private int numerator;
        private int denominator = 1;
        private boolean sign = false; // true = negative, false = positive

        public Fraction(int num, int denom) {
            this.numerator = num;
            if (denom == 0) {
                System.out.println("Error: denominator cannot be 0");
            } else {
                this.denominator = denom;
            }
            this.simplify();
        }

        public Fraction(int num) {
            this.numerator = num;
            this.simplify();
        }

        private int getGcm(int num1, int num2) {
            return num2 == 0 ? num1 : this.getGcm(num2, num1 % num2);
        }

        // Simplify fraction to simplest form, runs in constructor
        public void simplify() {
            this.sign = !(this.numerator <= 0 && this.denominator <= 0) && !(this.numerator >= 0 && this.denominator >= 0);

            this.numerator = Math.abs(this.numerator);
            this.denominator = Math.abs(this.denominator);

            int gcm = this.getGcm(this.numerator, this.denominator);
            this.numerator = this.numerator / gcm;
            this.denominator = this.denominator / gcm;
            if (this.numerator == 0 && this.denominator != 0) {
                this.denominator = 1;
                this.sign = false;
            }
        }

        public Fraction plus(Fraction f1) {
            int num = 0;
            if (this.sign) { // fraction is negative
                if (f1.getSign()) { // f1 is negative
                    num = (-1) * this.numerator * f1.denominator + this.denominator * (-1) * f1.numerator;
                } else { // f1 is positive
                    num = (-1) * this.numerator * f1.denominator + this.denominator * f1.numerator;
                }
            } else { // fraction is positive
                if (f1.getSign()) { // f1 is negative
                    num = this.numerator * f1.denominator + this.denominator * (-1) * f1.numerator;
                } else { // f1 is positive
                    num = this.numerator * f1.denominator + this.denominator * f1.numerator;
                }
            }
            int denom = this.denominator * f1.getDenominator();
            return new Fraction(num, denom);
        }

        public Fraction minus(Fraction f1) {
            int num = 0;
            if (this.sign) { // fraction is negative
                if (f1.getSign()) { // f1 is negative
                    num = (-1) * this.numerator * f1.denominator + this.denominator * f1.numerator;
                } else { // f1 is positive
                    num = (-1) * this.numerator * f1.denominator - this.denominator * f1.numerator;
                }
            } else { // fraction is positive
                if (f1.getSign()) { // f1 is negative
                    num = this.numerator * f1.denominator + this.denominator * f1.numerator;
                } else { // f1 is positive
                    num = this.numerator * f1.denominator - this.denominator * f1.numerator;
                }
            }
            int denom = this.denominator * f1.getDenominator();
            return new Fraction(num, denom);
        }

        public Fraction multiply(Fraction f1) {
            int signInt = 1;
            // Either one fraction is negative will make the product fraction negative, but not for both fractions are negative.
            if (this.sign && !f1.getSign() || !this.sign && f1.getSign()) {
                signInt = -1;
            }
            return new Fraction(signInt * this.numerator * f1.getNumerator(), this.denominator * f1.getDenominator());
        }

        public Fraction dividedBy(Fraction f1) {
            int signInt = 1;
            // Either one fraction is negative will make the product fraction negative, but not for both fractions are negative.
            if (this.sign && !f1.getSign() || !this.sign && f1.getSign()) {
                signInt = -1;
            }
            return new Fraction(signInt *this.numerator * f1.getDenominator(), this.denominator * f1.getNumerator());
        }

        public boolean equals(Fraction f1) {
            return this.numerator == f1.getNumerator() && this.denominator == f1.getDenominator() && this.sign == f1.getSign();
        }

        public int getNumerator() {
            return this.numerator;
        }

        public int getDenominator() {
            return this.denominator;
        }

        public boolean getSign() {
            return this.sign;
        }
    }
}