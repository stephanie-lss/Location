import java.math.BigDecimal;

import static java.lang.Math.pow;
import static java.math.BigDecimal.ROUND_HALF_DOWN;
import static java.math.BigDecimal.ROUND_HALF_UP;

public class HeightFitting {
    /**
     * �������ܣ�GPS�߳��쳣��Ϸ�������������Ϸ�������6����ĸ�˹�������Ϣ��ͨ�������������ض���������ߣ�
     *
     * @param x1 ��֪���ƽ������x1
     * @param x2 ��֪���ƽ������x2
     * @param x3 ��֪���ƽ������x3
     * @param x4 ��֪���ƽ������x4
     * @param x5 ��֪���ƽ������x5
     * @param x6 ��֪���ƽ������x6
     * @param h1 ��֪��ĸ߳��쳣h1
     * @param h2 ��֪��ĸ߳��쳣h2
     * @param h3 ��֪��ĸ߳��쳣h3
     * @param h4 ��֪��ĸ߳��쳣h4
     * @param h5 ��֪��ĸ߳��쳣h5
     * @param h6 ��֪��ĸ߳��쳣h6
     * @return a0��a1,a2,a3,a4,a5
     **/

    public static double[][] deltax(double x1, double x2, double x3, double x4, double x5, double x6, double h1, double h2, double h3, double h4, double h5, double h6) {
        double[][] NBB1 = NBB(x1, x2, x3, x4, x5, x6);
        double[][] W1 = W(x1, x2, x3, x4, x5, x6, h1, h2, h3, h4, h5, h6);
        double[][] deltax = multiplyAB(NBB1, W1);
        return deltax;
    }

    public static double[][] NBB(double x1, double x2, double x3, double x4, double x5, double x6) {
        double[][] BB = new double[6][6];
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                if (j == 0) {
                    BB[i][j] = 1;
                }
                if (i == 0 && j != 0) {
                    BB[i][j] = pow(x1, j);
                }
                if (i == 1 && j != 0) {
                    BB[i][j] = pow(x2, j);
                }
                if (i == 2 && j != 0) {
                    BB[i][j] = pow(x3, j);
                }
                if (i == 3 && j != 0) {
                    BB[i][j] = pow(x4, j);
                }
                if (i == 4 && j != 0) {
                    BB[i][j] = pow(x5, j);
                }
                if (i == 5 && j != 0) {
                    BB[i][j] = pow(x6, j);
                }
            }
        }

        double[][] NBB = multiplyAB(getA_T(BB), BB);
        return NBB;

    }

    public static double[][] W(double x1, double x2, double x3, double x4, double x5, double x6, double h1, double h2, double h3, double h4, double h5, double h6) {
        double[][] BB = new double[6][6];
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                if (j == 0) {
                    BB[i][j] = 1;
                }
                if (i == 0 && j != 0) {
                    BB[i][j] = pow(x1, j);
                }
                if (i == 1 && j != 0) {
                    BB[i][j] = pow(x2, j);
                }
                if (i == 2 && j != 0) {
                    BB[i][j] = pow(x3, j);
                }
                if (i == 3 && j != 0) {
                    BB[i][j] = pow(x4, j);
                }
                if (i == 4 && j != 0) {
                    BB[i][j] = pow(x5, j);
                }
                if (i == 5 && j != 0) {
                    BB[i][j] = pow(x6, j);
                }
            }
        }
        double[][] LL = {{h1}, {h2}, {h3}, {h4}, {h5}, {h6}};

        double[][] W = multiplyAB(getA_T(BB), LL);
        return W;
    }


    /**
     * �������ܣ��̶������������2����ĸ�˹�������Ϣ�����δ֪�������ߣ�
     *
     * @param H1 ��֪��1��ظ�
     * @param h1 ��֪��1������
     * @param H2 ��֪��2��ظ�
     * @param h2 ��֪��2������
     * @param H  δ֪���ظ�
     * @return δ֪��������
     */
    public static double FixDiffMod(double H1, double h1, double H2, double h2, double H) {
        double Hh1 = H1 - h1;
        double Hh2 = H2 - h2;
        double average = getAverage(Hh1, Hh2);
        double h = H - average;
        return h;
    }

    public static double getAverage(double a, double b) {
        double sum = 0;
        sum = (a + b) / 2;
        return sum;
    }

    /**
     * �������ܣ�GPS�߳��쳣��Ϸ�������ƽ����Ϸ���������ڵ���3����ĸ�˹�������Ϣ��ͨ��ƽ���������ض����������
     *
     * @param Know_Point_Information   ��֪����Ϣ��ÿ��Ϊһ����֪����Ϣ���ֱ�Ϊ��˹����x��y����ظ�H��������h�����������ڵ���3������Ϊ4��
     * @param Unknow_Point_Information δ֪����Ϣ��ÿ��Ϊһδ֪����Ϣ���ֱ�Ϊ��˹����x��y����ظ�H�����������⣬����Ϊ3��
     * @return δ֪��������h��
     */
    public static double[] Plane_Height_Fitting(double[][] Know_Point_Information, double[][] Unknow_Point_Information) {
        int Know_Point_Number = Know_Point_Information.length;
        double[][] Know_Point_XY = new double[Know_Point_Number][3];
        double[] Know_Point_Height_GPS = new double[Know_Point_Number];
        double[] Know_Point_Height_Normal = new double[Know_Point_Number];
        for (int i = 0; i < Know_Point_Number; i++) {
            Know_Point_XY[i][0] = 1;
            Know_Point_XY[i][1] = Know_Point_Information[i][0];
            Know_Point_XY[i][2] = Know_Point_Information[i][1];
            Know_Point_Height_GPS[i] = Know_Point_Information[i][2];
            Know_Point_Height_Normal[i] = Know_Point_Information[i][3];
        }
        double[][] pass0 = getA_T(Know_Point_XY);
        double[][] pass1 = multiplyAB(pass0, Know_Point_XY);
        double[][] pass2 = multiplyAB(getN(pass1), getA_T(Know_Point_XY));
        double[] pass3 = minusArry(Know_Point_Height_GPS, Know_Point_Height_Normal);
        double[] Optimal_Parameters = multiplyAB(pass2, pass3);
        int Unknow_Point_Number = Unknow_Point_Information.length;
        double[][] Unknow_Point_XY = new double[Unknow_Point_Number][3];
        double[] Unknow_Point_Height_GPS = new double[Unknow_Point_Number];
        double[] Unknow_Point_Height_Normal = new double[Unknow_Point_Number];
        for (int i = 0; i < Unknow_Point_Number; i++) {
            Unknow_Point_XY[i][0] = 1;
            Unknow_Point_XY[i][1] = Unknow_Point_Information[i][0];
            Unknow_Point_XY[i][2] = Unknow_Point_Information[i][1];
            Unknow_Point_Height_GPS[i] = Unknow_Point_Information[i][2];
        }
        double Unknow_Height_Unnormal[] = multiplyAB(Unknow_Point_XY, Optimal_Parameters);
        for (int i = 0; i < Unknow_Point_Number; i++) {
            Unknow_Point_Height_Normal[i] = Unknow_Point_Height_GPS[i] - Unknow_Height_Unnormal[i];
        }
        return Unknow_Point_Height_Normal;
    }

    public static double[] multiplyAB(double[][] matrixA, double[] arryB) {
        double[] matrixAB = new double[matrixA.length];

        if (matrixA[0].length != arryB.length) {
            System.out.println("MatrixOperate");
            System.exit(0);
        } else {
            for (int i = 0; i < matrixA.length; i++) {
                double sum = 0;
                for (int j = 0; j < matrixA[0].length; j++) {
                    sum = sum + matrixA[i][j] * arryB[j];
                }
                matrixAB[i] = sum;
            }
        }
        return matrixAB;
    }

    public static double[][] getN(double[][] data) {
        double[][] newData = new double[data.length][data.length];
        if (data.length == 2 && data[0].length == 2) {
            newData = getN2m2(data);
        } else {
            if (data.length != data[0].length) {
            }
            double A = getHL(data);
            for (int i = 0; i < data.length; i++) {
                for (int j = 0; j < data.length; j++) {
                    double num;
                    if ((i + j) % 2 == 0) {
                        num = getHL(getDY(data, i + 1, j + 1));
                    } else {
                        num = -getHL(getDY(data, i + 1, j + 1));
                    }
                    newData[i][j] = num / A;
                }
            }
            newData = getA_T(newData);
        }
        return newData;
    }

    public static double[][] getA_T(double[][] A) {
        int h = A.length;
        int v = A[0].length;
        double[][] A_T = new double[v][h];
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < v; j++) {
                A_T[j][i] = A[i][j];
            }
        }
        return A_T;
    }

    public static double[][] getDY(double[][] data, int h, int v) {
        int H = data.length;
        int V = data[0].length;
        double[][] newData = new double[H - 1][V - 1];
        if (H != V) {
        } else {
            for (int i = 0; i < newData.length; i++) {
                if (i < h - 1) {
                    for (int j = 0; j < newData[i].length; j++) {
                        if (j < v - 1) {
                            newData[i][j] = data[i][j];
                        } else {
                            newData[i][j] = data[i][j + 1];
                        }
                    }
                } else {
                    for (int j = 0; j < newData[i].length; j++) {
                        if (j < v - 1) {
                            newData[i][j] = data[i + 1][j];
                        } else {
                            newData[i][j] = data[i + 1][j + 1];
                        }
                    }
                }
            }
        }
        return newData;
    }

    public static double getHL(double[][] data) {
        double total = 0;
        int num = data.length;
        double[] nums = new double[num];
        if (data.length != data[0].length) {
        } else {
            if (data.length == 2) {
                return data[0][0] * data[1][1] - data[0][1] * data[1][0];
            }
            for (int i = 0; i < num; i++) {
                if (i % 2 == 0) {
                    nums[i] = data[0][i] * getHL(getDY(data, 1, i + 1));
                } else {
                    nums[i] = -data[0][i] * getHL(getDY(data, 1, i + 1));
                }
            }
            for (int i = 0; i < num; i++) {
                total += nums[i];
            }
        }
        return total;
    }

    public static double[][] getN2m2(double[][] data) {
        double[][] newmatrix = new double[2][2];
        newmatrix[0][0] = data[1][1];
        newmatrix[0][1] = -1 * data[1][0];
        newmatrix[1][0] = -1 * data[0][1];
        newmatrix[1][1] = 1 * data[0][0];
        return multiplyAb(1 / (data[0][0] * data[1][1] - data[1][0] * data[0][1]), newmatrix);
    }

    public static double[][] multiplyAb(double a, double[][] A) {
        double[][] aA = new double[A.length][A[0].length];

        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[0].length; j++) {
                aA[i][j] = a * A[i][j];
            }
        }
        return aA;
    }

    public static double[] minusArry(double[] arryA, double[] arryB) {
        double[] arryMinus = new double[arryA.length];
        if (arryA.length != arryB.length) {
            System.out.println("MatrixOperate");
            System.exit(0);
        } else {
            for (int i = 0; i < arryA.length; i++) {
                arryMinus[i] = arryA[i] - arryB[i];
            }
        }
        return arryMinus;
    }

    public static double[][] multiplyAB(double[][] matrixA, double[][] matrixB) {
        double[][] matrixAB = new double[matrixA.length][matrixB[0].length];
        if (matrixA[0].length != matrixB.length) {
            System.out.println("MatrixOperate");
            System.exit(0);
        } else {
            for (int i = 0; i < matrixA.length; i++) {
                for (int j = 0; j < matrixB[0].length; j++) {
                    double sum = 0;
                    for (int k = 0; k < matrixB.length; k++) {
                        sum = sum + matrixA[i][k] * matrixB[k][j];
                    }
                    matrixAB[i][j] = sum;
                }
            }
        }
        return matrixAB;
    }

    public static double[][] NBB(double x, double y) {
        double[][] BB = new double[1][6];
        BB[0][0] = 1;
        BB[0][1] = x;
        BB[0][2] = y;
        BB[0][3] = x * x;
        BB[0][4] = y * y;
        BB[0][5] = x * y;
        double[][] NBB = multiplyAB(getA_T(BB), BB);
        return NBB;
    }

    public static double[][] W(double x, double y, double H, double h) {
        double[][] BB = new double[1][6];
        BB[0][0] = 1;
        BB[0][1] = x;
        BB[0][2] = y;
        BB[0][3] = x * x;
        BB[0][4] = y * y;
        BB[0][5] = x * y;
        double[][] LL = new double[1][1];
        LL[0][0] = H - h;
        double[][] W = multiplyAB(getA_T(BB), LL);
        return W;
    }

    /**
     * �������ܣ�GPS�߳��쳣��Ϸ�������������Ϸ�������6����ĸ�˹�������Ϣ��ͨ�����������������������ߵ�����������
     *
     * @param x1 ��֪���ƽ������x
     * @param x2 ��֪���ƽ������x
     * @param x3 ��֪���ƽ������x
     * @param x4 ��֪���ƽ������x
     * @param x5 ��֪���ƽ������x
     * @param x6 ��֪���ƽ������x
     * @param y1 ��֪���ƽ������y
     * @param y2 ��֪���ƽ������y
     * @param y3 ��֪���ƽ������y
     * @param y4 ��֪���ƽ������y
     * @param y5 ��֪���ƽ������y
     * @param y6 ��֪���ƽ������y
     * @param h1 ��֪���������h
     * @param h2 ��֪���������h
     * @param h3 ��֪���������h
     * @param h4 ��֪���������h
     * @param h5 ��֪���������h
     * @param h6 ��֪���������h
     * @param H1 ��֪���������H
     * @param H2 ��֪���������H
     * @param H3 ��֪���������H
     * @param H4 ��֪���������H
     * @param H5 ��֪���������H
     * @param H6 ��֪���������H
     * @return a0��a1,a2,a3,a4,a5
     */
    public static double[][] deltax(double x1, double x2, double x3, double x4, double x5, double x6, double y1, double y2, double y3, double y4, double y5, double y6,
                                    double h1, double h2, double h3, double h4, double h5, double h6, double H1, double H2, double H3, double H4, double H5, double H6) {
        double[][] NBB = addAB(addAB(addAB(addAB(addAB(NBB(x1, y1), NBB(x2, y2)), NBB(x3, y3)), NBB(x4, y4)), NBB(x5, y5)), NBB(x6, y6));
        double[][] W = addAB(addAB(addAB(addAB(addAB(W(x1, y1, H1, h1), W(x2, y2, H2, h2)), W(x3, y3, H3, h3)), W(x4, y4, H4, h4)), W(x5, y5, H5, h5)), W(x6, y6, H6, h6));
        double[][] NBB1 = getN(NBB);
        for (int i = 0; i < NBB.length; i++) {
            for (int j = 0; j < NBB[0].length; j++) {
                System.out.print(NBB1[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println();
        System.out.println();
        System.out.println();
        for (int i = 0; i < NBB.length; i++) {
            for (int j = 0; j < NBB[0].length; j++) {
                System.out.print(NBB[i][j] + " ");
            }
            System.out.println();
        }
        double[][] deltax = multiplyAB(inverse(NBB), W);
        return deltax;
    }

    public static double div(double v1, double v2) {
        return div(v1, v2, 15);
    }

    public static double div(double v1, double v2, int scale) {
        if (scale < 0) {

            throw new IllegalArgumentException(
                    "The scale must be a positive integer or zero");
        }
        BigDecimal b1 = new BigDecimal(Double.valueOf(v1));
        BigDecimal b2 = new BigDecimal(Double.valueOf(v2));
        return b1.divide(b2, scale, ROUND_HALF_UP)
                .doubleValue();
    }

    public static double[][] addAB(double[][] A, double[][] B) {
        double[][] ABadd = new double[A.length][A[0].length];
        if (A.length != B.length || A[0].length != B[0].length) {
            System.out.println("MatrixOperate");
            System.exit(0);
        } else {
            for (int i = 0; i < A.length; i++) {
                for (int j = 0; j < A[0].length; j++) {
                    ABadd[i][j] = A[i][j] + B[i][j];
                }
            }
        }
        return ABadd;
    }

    public static double[][] inverse(double[][] A) {
        int i = 0, j = 0, n = A[0].length;
        double[][] inv = new double[n][n];
        double[][] B = new double[n][n];
        for (i = 0; i < n; i++) {
            B[i][i] = 1;
        }
        for (i = 0; i < n; i++) {
            double[][] L = new double[n][n];
            double[][] U = new double[n][n];
            int[] P = new int[n];
            double[][] copyOfA = copy(A);
            LUP_Decomposition(copyOfA, L, U, P);
            inv[i] = LUP_Solve(L, U, P, B[i]);
        }
        transposition(inv);
        return inv;
    }

    public static double[][] copy(double[][] A) {
        int len = A[0].length;
        double[][] copy = new double[len][len];
        for (int i = 0; i < len; i++) {
            for (int j = 0; j < len; j++) {
                copy[i][j] = A[i][j];
            }
        }
        return copy;
    }

    public static double[][] transposition(double[][] A) {
        int i = 0, j = 0, n = A[0].length;
        double tmp = 0f;
        for (i = 0; i < n; i++) {
            for (j = 0; j < i; j++) {
                tmp = A[i][j];
                A[i][j] = A[j][i];
                A[j][i] = tmp;
            }
        }
        return A;
    }

    public static double[] LUP_Solve(double[][] L, double[][] U, int[] P, double[] b) {
        int n = L[0].length, i = 0, j = 0;
        double[] x = new double[n];
        double[] y = new double[n];
        for (i = 0; i < n; i++) {
            y[i] = b[P[i]];
            for (j = 0; j < i; j++) {
                BigDecimal bigDecimal = new BigDecimal(Double.toString(y[i]));
                BigDecimal bigDecima2 = new BigDecimal(Double.toString(y[j]));
                BigDecimal bigDecima3 = new BigDecimal(Double.toString(L[i][j]));
                bigDecima2.setScale(33, ROUND_HALF_UP);
                bigDecimal.setScale(33, ROUND_HALF_UP);
                bigDecima3.setScale(33, ROUND_HALF_UP);
                y[i] = (bigDecimal.subtract(bigDecima3.multiply(bigDecima2))).doubleValue();
            }
        }
        for (i = n - 1; i >= 0; i--) {
            x[i] = y[i];
            for (j = n - 1; j > i; j--) {
                BigDecimal bigDecimal = new BigDecimal(Double.toString(x[i]));
                BigDecimal bigDecima2 = new BigDecimal(Double.toString(x[j]));
                BigDecimal bigDecima3 = new BigDecimal(Double.toString(U[i][j]));
                bigDecimal.setScale(33, ROUND_HALF_UP);
                bigDecima2.setScale(33, ROUND_HALF_UP);
                bigDecima3.setScale(33, ROUND_HALF_UP);
                x[i] = (bigDecimal.subtract(bigDecima3.multiply(bigDecima2))).doubleValue();
            }
            x[i] = div(x[i], U[i][i], 36);
            System.out.println("U[i][i]" + U[i][i]);
        }
        return x;
    }

    public static void LUP_Decomposition(double[][] A, double[][] L, double[][] U, int[] P) {
        int n = A[0].length;
        int i = 0, j = 0, k = 0, row = 0;
        for (i = 0; i < n; i++) P[i] = i;
        for (i = 0; i < n - 1; i++) {
            double p = 0f;
            for (j = i; j < n; j++) {
                if (Math.abs(A[j][i]) > p) {
                    p = Math.abs(A[j][i]);
                    row = j;
                }
            }
            if (p == 0) {
                System.err.println("singular matrix");
                return;
            }
            int tmp = P[i];
            P[i] = P[row];
            P[row] = tmp;

            double tmp2 = 0f;
            for (j = 0; j < n; j++) {
                tmp2 = A[i][j];
                A[i][j] = A[row][j];
                A[row][j] = tmp2;
            }
            double u = A[i][i], l = 0f;
            for (j = i + 1; j < n; j++) {
                BigDecimal bigDecimal11 = new BigDecimal(Double.toString(A[j][i]));
                BigDecimal bigDecima2 = new BigDecimal(Double.toString(A[i][i]));
                l = bigDecimal11.divide(bigDecima2, 36, ROUND_HALF_DOWN).doubleValue();
                System.out.println("u:" + u);
                A[j][i] = l;
                for (k = i + 1; k < n; k++) {
                    BigDecimal bigDecimal = new BigDecimal(Double.toString(A[j][k]));
                    BigDecimal bigDecimal1 = new BigDecimal(Double.toString(A[i][k]));
                    BigDecimal bigDecimal2 = new BigDecimal(Double.toString(l));
                    bigDecimal.setScale(33, ROUND_HALF_UP);
                    bigDecimal1.setScale(33, ROUND_HALF_UP);
                    bigDecimal2.setScale(33, ROUND_HALF_UP);
                    A[j][k] = (bigDecimal.subtract(bigDecimal1.multiply(bigDecimal2))).doubleValue();
                }
            }
        }
        for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++) {
                if (i != j)
                    L[i][j] = A[i][j];
                else L[i][j] = 1;
            }
            for (k = i; k < n; k++) {
                U[i][k] = A[i][k];
            }
        }
    }
}
