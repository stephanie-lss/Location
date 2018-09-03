import static java.lang.Math.*;
import static java.lang.Math.cos;
import static java.lang.Math.pow;


import static java.lang.Math.*;
public class ParameterModel {
    static double PI = 3.1415926;
    static double yjcs=500000;
    static double a80=6378140.0;
    static double e80_pow=0.006694384999588;
    static double e80_pow_1=0.006739501819473;
    static double a84=6378137.0;
    static double e84_pow=0.00669437999013;
    static double L0=120.0000;
    static double rou=206264.806247096355;
    static double k = Math.pow(10, 6);

    /**
     * Ī���˹���߲������
     * @param Bo1  ԭʼ����ϵ�µ�һ����֪���γ��
     * @param Lo1  ԭʼ����ϵ�µ�һ����֪��ľ���
     * @param Ho1  ԭʼ����ϵ�µ�һ����֪��ĸ߶�
     * @param Bo2  ԭʼ����ϵ�µڶ�����֪���γ��
     * @param Lo2  ԭʼ����ϵ�µڶ�����֪��ľ���
     * @param Ho2  ԭʼ����ϵ�µڶ�����֪��ĸ߶�
     * @param Bo3  ԭʼ����ϵ�µ�������֪���γ��
     * @param Lo3  ԭʼ����ϵ�µ�������֪��ľ���
     * @param Ho3  ԭʼ����ϵ�µ�������֪��ĸ߶�
     * @param xd1  Ŀ������ϵ�µ�һ����֪���xd1
     * @param xd2  Ŀ������ϵ�µڶ�����֪���xd2
     * @param xd3  Ŀ������ϵ�µ�������֪���xd3
     * @param yd1  Ŀ������ϵ�µ�һ����֪���yd1
     * @param yd2  Ŀ������ϵ�µڶ�����֪���yd2
     * @param yd3  Ŀ������ϵ�µ�������֪���yd3
     * @param Hd1  Ŀ������ϵ�µ�һ����֪���Hd1
     * @param Hd2  Ŀ������ϵ�µڶ�����֪���Hd2
     * @param Hd3  Ŀ������ϵ�µ�������֪���Hd3
     * @return     7*1�Ķ�ά����deltax1���ܹ��߸���Ī���˹���߲�������Ϊ��X0(m) ��Y0(m) ��Z0(m) ��x(s) ��y(s) ��z(s)  m(ppm)
     */
    public static double[][] deltax1(double Bo1,double Lo1,double Ho1,double Bo2,double Lo2,double Ho2,double Bo3,double Lo3,double Ho3,
                                     double xd1,double xd2,double xd3,double yd1,double yd2,double yd3,double Hd1,double Hd2,double Hd3){
        double[] XYZo1 = BLHXYZ(Bo1,Lo1,Ho1, a84, e84_pow);
        double[] XYZo2 = BLHXYZ(Bo2,Lo2,Ho2, a84, e84_pow);
        double[] XYZo3 = BLHXYZ(Bo3,Lo3,Ho3, a84, e84_pow);
        double[] BL1 = gaussfs(xd1,yd1 - yjcs, L0, a80, e80_pow, e80_pow_1);
        double[] BL2 = gaussfs(xd2,yd2 - yjcs, L0, a80, e80_pow, e80_pow_1);
        double[] BL3 = gaussfs(xd3,yd3 - yjcs, L0, a80, e80_pow, e80_pow_1);
        double[] XYZd1 = BLHXYZ(BL1[0], BL1[1], Hd1, a80, e80_pow);
        double[] XYZd2 = BLHXYZ(BL2[0], BL2[1], Hd2, a80, e80_pow);
        double[] XYZd3 = BLHXYZ(BL3[0], BL3[1], Hd3, a80, e80_pow);
        double X0 = getAverage(XYZo1[0], XYZo2[0], XYZo3[0]);//Xo1,Xo2,Xo3
        double Y0 = getAverage(XYZo1[1], XYZo2[1], XYZo3[1]);
        double Z0 = getAverage(XYZo1[2], XYZo2[2], XYZo3[2]);
        double[][] NBB1 = molodenskyNBB(XYZo1[0], XYZo1[1], XYZo1[2], XYZd1[0], XYZd1[1], XYZd1[2], X0, Y0, Z0);
        double[][] NBB2 = molodenskyNBB(XYZo2[0], XYZo2[1], XYZo2[2], XYZd2[0], XYZd2[1], XYZd2[2], X0, Y0, Z0);
        double[][] NBB3 = molodenskyNBB(XYZo3[0], XYZo3[1], XYZo3[2], XYZd3[0], XYZd3[1], XYZd3[2], X0, Y0, Z0);
        double[][] w1 = molodenskyw(XYZo1[0], XYZo1[1], XYZo1[2], XYZd1[0], XYZd1[1], XYZd1[2], X0, Y0, Z0);
        double[][] w2 = molodenskyw(XYZo2[0], XYZo2[1], XYZo2[2], XYZd2[0], XYZd2[1], XYZd2[2], X0, Y0, Z0);
        double[][] w3 = molodenskyw(XYZo3[0], XYZo3[1], XYZo3[2], XYZd3[0], XYZd3[1], XYZd3[2], X0, Y0, Z0);
        double[][] sum1 = addAB(NBB1, NBB2);
        double[][] NBBsum = addAB(sum1, NBB3);
        double[][] sum2 = addAB(w1, w2);
        double[][] wsum = addAB(sum2, w3);
        double[][] deltax1 = multiplyAB(getN(NBBsum), wsum);
        return deltax1;
    }
    public static double[] gaussfs(double x, double y, double L0, double a, double e_pow, double e_pow_1) {
        double[] result = new double[2];
        double rL0 = dhhd(L0);
        double A0 = 1 + 3 * e_pow / 4 + 45 * (pow(e_pow, 2)) / 64 + 350 * (pow(e_pow, 3)) / 512 + 11025 * (pow(e_pow, 4)) / 16384;
        double B0 = x / (a * (1 - e_pow) * A0);
        double K0 = (3 * e_pow / 4 + 45 * (pow(e_pow, 2)) / 64 + 350 * (pow(e_pow, 3)) / 512 + 11025 * (pow(e_pow, 4)) / 16384) / 2;
        double K2 = -(63 * (pow(e_pow, 2)) / 64 + 1108 * (pow(e_pow, 3)) / 512 + 58239 * (pow(e_pow, 4)) / 16384) / 3;
        double K4 = (604 * (pow(e_pow, 3)) / 512 + 68484 * (pow(e_pow, 4)) / 16384) / 3;
        double K6 = -(26328 * (pow(e_pow, 4)) / 16384) / 3;
        double sB0 = pow(sin(B0), 2);
        double Bf = B0 + sin(2 * B0) * (K0 + sB0 * (K2 + sB0 * (K4 + K6 * sB0)));
        double n_powf = e_pow_1 * (pow(cos(Bf), 2));
        double tf = tan(Bf);
        double tf2 = pow(tf, 2);
        double tf4 = pow(tf, 4);
        double Wf = sqrt(1 - e_pow * (pow(sin(Bf), 2)));
        double Nf = a / Wf;
        double yNf = y / Nf;
        double Vf2 = 1 + n_powf;
        //B=Bf-0.5*Vf2*tf*(yNf^2-(5+3*tf2+n_powf-9*n_powf*tf2)*(yNf^4)/12+(61+90*tf2+45*tf4)*(yNf^6)/360);
        double B = Bf - 0.5 * Vf2 * tf * (pow(yNf, 2) - (5 + 3 * tf2 + n_powf - 9 * n_powf * tf2) * (pow(yNf, 4)) / 12 + (61 + 90 * tf2 + 45 * tf4) * (pow(yNf, 6)) / 360);
        //l=(yNf-(1+2*tf2+n_powf)*(yNf^3)/6+(5+28*tf2+24*tf4+6*n_powf+8*n_powf*tf2)*(yNf^5)/120)/cos(Bf);   %�Ի���Ϊ��λ�����
        double l = (yNf - (1 + 2 * tf2 + n_powf) * (pow(yNf, 3)) / 6 + (5 + 28 * tf2 + 24 * tf4 + 6 * n_powf + 8 * n_powf * tf2) * (pow(yNf, 5)) / 120) / cos(Bf);
        double L = rL0 + l;
        B = hdhd(B);
        L = hdhd(L);
        result[0] = B;
        result[1] = L;
        return result;
    }
    /**
     * ����ת��Ϊ����
     *
     * @param alfa ����alfa������ ����
     * @return ��Ӧ�����Ļ���ֵ
     */
    public static double dhhd(double alfa) {
        alfa = alfa + Math.pow(10, -10);
        double alfa1 = Math.floor(alfa) + Math.floor((alfa - Math.floor(alfa)) * 100) / 60;
        double alfa2 = (alfa * 100 - Math.floor(alfa * 100)) / 36;
        alfa = (alfa1 + alfa2) * PI / 180;
        return alfa;
    }

    /**
     * ������ת��Ϊ��
     *
     * @param alfa ����alfa������  ����
     * @return ��Ӧ����ֵ�Ķ���
     */
    public static double hdhd(double alfa) {
        alfa = alfa + Math.pow(10, -10);
        alfa = alfa * 180 / PI;
        double alfa1 = Math.floor(alfa) + Math.floor((alfa - Math.floor(alfa)) * 60) / 100;
        double alfa2 = (alfa * 60 - Math.floor(alfa * 60)) * 0.006;
        alfa = alfa1 + alfa2;
        return alfa;
    }
    /**
     * ���������ת��Ϊֱ������
     *
     * @param B     ����ֵ
     * @param L     γ��ֵ
     * @param H     �߶�ֵ
     * @param a     WGS84�������
     * @param e_pow WGS84�������
     * @return �������ת�����Ӧ�Ŀռ�ֱ������X Y Z
     */
    public static double[] BLHXYZ(double B, double L, double H, double a, double e_pow) {
        B = dhhd(B);
        L = dhhd(L);
        double[] result = new double[3];
        double w = Math.sqrt(1 - e_pow * (Math.sin(B) * Math.sin(B)));
        double N = a / w;
        double X = (N + H) * Math.cos(B) * Math.cos(L);
        double Y = (N + H) * Math.cos(B) * Math.sin(L);
        double Z = (N * (1 - e_pow) + H) * Math.sin(B);
        result[0] = X;
        result[1] = Y;
        result[2] = Z;
        return result;
    }
    public static double getAverage(double a, double b, double c) {
        double sum = 0;
        sum = (a + b + c) / 3;
        return sum;
    }
    /**
     * �������NBB��������С���˷���Ҫ�Ĳ���
     *
     * @param Xo Դ�����ڿռ�ֱ������ϵ�µ�x������
     * @param Yo Դ�����ڿռ�ֱ������ϵ�µ�y������
     * @param Zo Դ�����ڿռ�ֱ������ϵ�µ�z�������
     * @param Xd Ŀ�������ڿռ�ֱ������ϵ�µ�x������
     * @param Yd Ŀ�������ڿռ�ֱ������ϵ�µ�x������
     * @param Zd Ŀ�������ڿռ�ֱ������ϵ�µ�x������
     * @param X0 Դ����Xo��ƽ��ֵ
     * @param Y0 Դ����Yo��ƽ��ֵ
     * @param Z0 Դ����Zo��ƽ��ֵ
     * @return NBB����
     */
    public static double[][] molodenskyNBB(double Xo, double Yo, double Zo, double Xd, double Yd, double Zd, double X0, double Y0, double Z0) {
        double[][] BB = new double[3][7];
        BB[0][0] = 1;
        BB[0][1] = 0;
        BB[0][2] = 0;
        BB[0][3] = 0;
        BB[0][4] = -(Zo - Z0) / rou;
        BB[0][5] = (Yo - Y0) / rou;
        BB[0][6] = (Xo - X0) / k;
        BB[1][0] = 0;
        BB[1][1] = 1;
        BB[1][2] = 0;
        BB[1][3] = (Zo - Z0) / rou;
        BB[1][4] = 0;
        BB[1][5] = -(Xo - X0) / rou;
        BB[1][6] = (Yo - Y0) / k;
        BB[2][0] = 0;
        BB[2][1] = 0;
        BB[2][2] = 1;
        BB[2][3] = -(Yo - Y0) / rou;
        BB[2][4] = (Xo - X0) / rou;
        BB[2][5] = 0;
        BB[2][6] = (Zo - Z0) / k;
        double[][] NBB = multiplyAB(getA_T(BB), BB);
        return NBB;
    }
    public static double[][] molodenskyw(double Xo, double Yo, double Zo, double Xd, double Yd, double Zd, double X0, double Y0, double Z0) {
        double[][] BB = new double[3][7];
        double[][] LL = new double[3][1];
        BB[0][0] = 1;
        BB[0][1] = 0;
        BB[0][2] = 0;
        BB[0][3] = 0;
        BB[0][4] = -(Zo - Z0) / rou;
        BB[0][5] = (Yo - Y0) / rou;
        BB[0][6] = (Xo - X0) / k;
        BB[1][0] = 0;
        BB[1][1] = 1;
        BB[1][2] = 0;
        BB[1][3] = (Zo - Z0) / rou;
        BB[1][4] = 0;
        BB[1][5] = -(Xo - X0) / rou;
        BB[1][6] = (Yo - Y0) / k;
        BB[2][0] = 0;
        BB[2][1] = 0;
        BB[2][2] = 1;
        BB[2][3] = -(Yo - Y0) / rou;
        BB[2][4] = (Xo - X0) / rou;
        BB[2][5] = 0;
        BB[2][6] = (Zo - Z0) / k;
        LL[0][0] = Xd - Xo;
        LL[1][0] = Yd - Yo;
        LL[2][0] = Zd - Zo;
        double[][] w = multiplyAB(getA_T(BB), LL);
        return w;
    }
    public static double[][] multiplyAB(double[][] matrixA, double[][] matrixB) {
        double[][] matrixAB = new double[matrixA.length][matrixB[0].length];
        if (matrixA[0].length != matrixB.length) {
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
    public static double[][] addAB(double[][] A, double[][] B) {
        double[][] ABadd = new double[A.length][A[0].length];
        if (A.length != B.length || A[0].length != B[0].length) {
        } else {
            for (int i = 0; i < A.length; i++) {
                for (int j = 0; j < A[0].length; j++) {
                    ABadd[i][j] = A[i][j] + B[i][j];
                }
            }
        }
        return ABadd;
    }
    public static double[][] getN(double[][] data) {
        double[][] newData = new double[data.length][data.length];
        if (data.length == 2 && data[0].length == 2) {
            newData = getN2m2(data);

        }
        else {
            if (data.length != data[0].length) {
            }
            double A = getHL(data); //
            System.out.println("A Gethl");
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

    /**
     * ����ɯ�߲������
     * @param B1  ԭʼ����ϵ�µ�һ����֪���γ��
     * @param L1   ԭʼ����ϵ�µ�һ����֪��ľ���
     * @param H1   ԭʼ����ϵ�µ�һ����֪��ĸ߶�
     * @param B2   ԭʼ����ϵ�µڶ�����֪���γ��
     * @param L2   ԭʼ����ϵ�µڶ�����֪��ľ���
     * @param H2   ԭʼ����ϵ�µڶ�����֪��ĸ߶�
     * @param B3   ԭʼ����ϵ�µ�������֪���γ��
     * @param L3   ԭʼ����ϵ�µ�������֪��ľ���
     * @param H3   ԭʼ����ϵ�µ�������֪��ĸ߶�
     * @param X1   Ŀ������ϵ�µ�һ����֪���X
     * @param X2   Ŀ������ϵ�µ�һ����֪���Y
     * @param X3   Ŀ������ϵ�µ�һ����֪���Z
     * @param Y1   Ŀ������ϵ�µڶ�����֪���X
     * @param Y2   Ŀ������ϵ�µڶ�����֪���Y
     * @param Y3   Ŀ������ϵ�µڶ�����֪���Z
     * @param Z1   Ŀ������ϵ�µ�������֪���X
     * @param Z2   Ŀ������ϵ�µ�������֪���Y
     * @param Z3   Ŀ������ϵ�µ�������֪���Z
     * @return     7*1�Ķ�ά����deltax���ܹ��߸�������ɯ�߲�������Ϊ��X0(m) ��Y0(m) ��Z0(m) ��x(s) ��y(s) ��z(s)  m(ppm)
     */
    public static  double[][] deltax3(double B1,double L1,double H1,double B2,double L2,double H2,double B3,double L3,double H3,
                                     double X1,double X2,double X3,double Y1,double Y2,double Y3,double Z1,double Z2,double Z3)
    {
        double[] R0=BLHXYZ(B1,L1,H1,a84,e84_pow);
        double Xo1=R0[0];
        double Yo1=R0[1];
        double Zo1=R0[2];
        double[] R1=BLHXYZ(B2,L2,H2,a84,e84_pow);
        double Xo2=R1[0];
        double Yo2=R1[1];
        double Zo2=R1[2];
        double[] R2=BLHXYZ(B3,L3,H3,a84,e84_pow);
        double Xo3=R2[0];
        double Yo3=R2[1];
        double Zo3=R2[2];
        double[] E1=gaussfs(X1, Y1-yjcs, L0, a80,e80_pow,e80_pow_1);
        double Bd1=E1[0];
        double Ld1=E1[1];
        double[] E2=gaussfs(X2, Y2-yjcs, L0, a80,e80_pow, e80_pow_1);
        double Bd2=E2[0];
        double Ld2=E2[1];
        double[] E3=gaussfs(X3, Y3-yjcs, L0, a80,e80_pow, e80_pow_1);
        double Bd3=E3[0];
        double Ld3=E3[1];
        double[] T0=BLHXYZ(Bd1,Ld1,Z1,a80,e80_pow);
        double Xd1=T0[0];
        double Yd1=T0[1];
        double Zd1=T0[2];
        double[] T1=BLHXYZ(Bd2,Ld2,Z2,a80,e80_pow);
        double Xd2=T1[0];
        double Yd2=T1[1];
        double Zd2=T1[2];
        double[] T2=BLHXYZ(Bd3,Ld3,Z3,a80,e80_pow);
        double Xd3=T2[0];
        double Yd3=T2[1];
        double Zd3=T2[2];
        double[][] NBB123=addAB(addAB(scsqj1(Xo1,Yo1,Zo1),scsqj1(Xo2,Yo2,Zo2)),scsqj1(Xo3,Yo3,Zo3));
        double[][] W123=addAB(addAB(scsqj2(Xo1,Yo1,Zo1,Xd1,Yd1,Zd1),scsqj2(Xo2,Yo2,Zo2,Xd2,Yd2,Zd2)),scsqj2(Xo3,Yo3,Zo3,Xd3,Yd3,Zd3));
        double[][] deltax=multiplyAB(getN(NBB123),W123);
        return deltax;
    }
    /**
     * ����ɯ��������е��м�ֵNBB
     * @param Xo  ԭʼ����ϵ�µ�X
     * @param Yo  ԭʼ����ϵ�µ�Y
     * @param Zo  ԭʼ����ϵ�µ�Z
     * @return  NBB
     */
    public  static double[][] scsqj1(double Xo,double Yo,double Zo){

        double[][] BB=new double[3][7];
        BB[0][0]=1;
        BB[0][1]=0;
        BB[0][2]=0;
        BB[0][3]=0;
        BB[0][4]=-Zo/rou;
        BB[0][5]=Yo/rou;
        BB[0][6]=Xo/k;
        BB[1][0]=0;
        BB[1][1]=1;
        BB[1][2]=0;
        BB[1][3]=Zo/rou;
        BB[1][4]=0;
        BB[1][5]=-Xo/rou;
        BB[1][6]=Yo/k;
        BB[2][0]=0;
        BB[2][1]=0;
        BB[2][2]=1;
        BB[2][3]=-Yo/rou;
        BB[2][4]=Xo/rou;
        BB[2][5]=0;
        BB[2][6]=Zo/k;
        double[][] NBB=multiplyAB(getA_T(BB), BB);
        return NBB;
    }
    public static double[][] scsqj2(double Xo,double Yo,double Zo,double Xd,double Yd,double Zd){
        double[][] BB=new double[3][7];
        BB[0][0]=1;
        BB[0][1]=0;
        BB[0][2]=0;
        BB[0][3]=0;
        BB[0][4]=-Zo/rou;
        BB[0][5]=Yo/rou;
        BB[0][6]=Xo/k;
        BB[1][0]=0;
        BB[1][1]=1;
        BB[1][2]=0;
        BB[1][3]=Zo/rou;
        BB[1][4]=0;
        BB[1][5]=-Xo/rou;
        BB[1][6]=Yo/k;
        BB[2][0]=0;
        BB[2][1]=0;
        BB[2][2]=1;
        BB[2][3]=-Yo/rou;
        BB[2][4]=Xo/rou;
        BB[2][5]=0;
        BB[2][6]=Zo/k;
        double[][] LL=new double[3][1];
        LL[0][0]=Xd-Xo;
        LL[1][0]=Yd-Yo;
        LL[2][0]=Zd-Zo;
        double[][] W=multiplyAB(getA_T(BB), LL);
        return W;
    }

    /**
     * ����ʽ�������
     * @param Xo1  ��֪�غϵ�X ������
     * @param Yo1  ��֪�غϵ�Y ������
     * @param Zo1  ��֪�غϵ�Z ������
     * @param Xo2  ��֪�غϵ�X ������
     * @param Yo2  ��֪�غϵ�Y ������
     * @param Zo2  ��֪�غϵ�Z ������
     * @param Xo3  ��֪�غϵ�X ������
     * @param Yo3  ��֪�غϵ�Y ������
     * @param Zo3  ��֪�غϵ�Z ������
     * @param Xo4  ��֪�غϵ�X ������
     * @param Yo4  ��֪�غϵ�Y ������
     * @param Zo4  ��֪�غϵ�Z ������
     * @param Xo5  ��֪�غϵ�X ������
     * @param Yo5  ��֪�غϵ�Y ������
     * @param Zo5  ��֪�غϵ�Z ������
     * @param Xo6  ��֪�غϵ�X ������
     * @param Yo6  ��֪�غϵ�Y ������
     * @param Zo6  ��֪�غϵ�Z ������
     * @param Xo7  ��֪�غϵ�X ������
     * @param Yo7  ��֪�غϵ�Y ������
     * @param Zo7  ��֪�غϵ�Z ������
     * @param Xo8  ��֪�غϵ�X ������
     * @param Yo8  ��֪�غϵ�Y ������
     * @param Zo8  ��֪�غϵ�Z ������
     * @param Xd1  ��Xo1���غϵ�X ������
     * @param Yd1  ��Yo1���غϵ�Y ������
     * @param Zd1  ��Zo1���غϵ�Z ������
     * @param Xd2  ��Xo2���غϵ�X ������
     * @param Yd2  ��Yo2���غϵ�Y ������
     * @param Zd2  ��Zo2���غϵ�Z ������
     * @param Xd3  ��Xo3���غϵ�X ������
     * @param Yd3  ��Yo3���غϵ�Y ������
     * @param Zd3  ��Zo3���غϵ�Z ������
     * @param Xd4  ��Xo4���غϵ�X ������
     * @param Yd4  ��Yo4���غϵ�Y ������
     * @param Zd4  ��Zo4���غϵ�Z ������
     * @param Xd5  ��Xo1���غϵ�X ������
     * @param Yd5  ��Yo5���غϵ�Y ������
     * @param Zd5  ��Zo5���غϵ�Z ������
     * @param Xd6  ��Xo6���غϵ�X ������
     * @param Yd6  ��Yo6���غϵ�Y ������
     * @param Zd6  ��Zo6���غϵ�Z ������
     * @param Xd7  ��Xo7���غϵ�X ������
     * @param Yd7  ��Yo7���غϵ�Y ������
     * @param Zd7  ��Zo7���غϵ�Z ������
     * @param Xd8  ��Xo8���غϵ�X ������
     * @param Yd8  ��Yo8���غϵ�Y ������
     * @param Zd8  ��Zo8���غϵ�Z ������
     * @return     ת������
     */
    public double[][] deltax2(double Xo1,double Yo1,double Zo1,double Xo2,double Yo2,double Zo2,double Xo3,double Yo3,double Zo3,double Xo4,double Yo4,double Zo4,double Xo5,double Yo5,double Zo5,double Xo6,double Yo6,double Zo6,double Xo7,double Yo7,double Zo7,double Xo8,double Yo8,double Zo8,
                             double Xd1,double Yd1,double Zd1,double Xd2,double Yd2,double Zd2,double Xd3,double Yd3,double Zd3,double Xd4,double Yd4,double Zd4,double Xd5,double Yd5,double Zd5,double Xd6,double Yd6,double Zd6,double Xd7,double Yd7,double Zd7,double Xd8,double Yd8,double Zd8) {

        double[][] NBB1 = polyNBB(Xo1, Yo1, Zo1);
        double[][] NBB2 = polyNBB(Xo2, Yo2, Zo2);
        double[][] NBB3 = polyNBB(Xo3, Yo3, Zo3);
        double[][] NBB4 = polyNBB(Xo4, Yo4, Zo4);
        double[][] NBB5 = polyNBB(Xo5, Yo5, Zo5);
        double[][] NBB6 = polyNBB(Xo6, Yo6, Zo6);
        double[][] NBB7 = polyNBB(Xo7, Yo7, Zo7);
        double[][] NBB8 = polyNBB(Xo8, Yo8, Zo8);
        double[][] w1 = polyw(Xo1, Yo1, Zo1, Xd1, Yd1, Zd1);
        double[][] w2 = polyw(Xo2, Yo2, Zo2, Xd2, Yd2, Zd2);
        double[][] w3 = polyw(Xo3, Yo3, Zo3, Xd3, Yd3, Zd3);
        double[][] w4 = polyw(Xo4, Yo4, Zo4, Xd4, Yd4, Zd4);
        double[][] w5 = polyw(Xo5, Yo5, Zo5, Xd5, Yd5, Zd5);
        double[][] w6 = polyw(Xo6, Yo6, Zo6, Xd6, Yd6, Zd6);
        double[][] w7 = polyw(Xo7, Yo7, Zo7, Xd7, Yd7, Zd7);
        double[][] w8 = polyw(Xo8, Yo8, Zo8, Xd8, Yd8, Zd8);
        double[][] NBBsum = addAB(addAB(addAB(addAB(addAB(addAB(addAB(NBB1, NBB2), NBB3), NBB4), NBB5), NBB6), NBB7),NBB8);
        double[][] wsum = addAB(addAB(addAB(addAB(addAB(addAB(addAB(w1, w2), w3), w4), w5), w6), w7),w8);
        double[][] result = multiplyAB(getN(NBBsum), wsum);
        return result;
    }
    public double[][] polyNBB(double Xo,double Yo,double Zo){
        double[][] BB = new double[3][9];
        BB[0][0] = Xo;
        BB[0][1] = 1;
        BB[0][2] = Xo;
        BB[0][3] = Yo;
        BB[0][4] = Zo;
        BB[0][5] = Xo*Yo;
        BB[0][6] = Zo*Xo;
        BB[0][7] = Zo*Yo;
        BB[0][8] = Zo*Yo*Xo;
        BB[1][0] = Yo;
        BB[1][1] = 1;
        BB[1][2] = Xo;
        BB[1][3] = Yo;
        BB[1][4] = Zo;
        BB[1][5] = Xo*Yo;
        BB[1][6] = Zo*Xo;
        BB[1][7] = Zo*Yo;
        BB[1][8] = Zo*Yo*Xo;
        BB[2][0] = Zo;
        BB[2][1] = 1;
        BB[2][2] = Xo;
        BB[2][3] = Yo;
        BB[2][4] = Zo;
        BB[2][5] = Xo*Yo;
        BB[2][6] = Zo*Xo;
        BB[2][7] = Zo*Yo;
        BB[2][8] = Zo*Yo*Xo;
        double[][] NBB = multiplyAB(getA_T(BB),BB);
        return NBB;
    }

    public double[][] polyw(double Xo,double Yo,double Zo,double Xd,double Yd,double Zd){
        double[][] BB = new double[3][9];
        double[][] LL = new double[3][3];
        BB[0][0] = Xo;
        BB[0][1] = 1;
        BB[0][2] = Xo;
        BB[0][3] = Yo;
        BB[0][4] = Zo;
        BB[0][5] = Xo*Yo;
        BB[0][6] = Zo*Xo;
        BB[0][7] = Zo*Yo;
        BB[0][8] = Zo*Yo*Xo;
        BB[1][0] = Yo;
        BB[1][1] = 1;
        BB[1][2] = Xo;
        BB[1][3] = Yo;
        BB[1][4] = Zo;
        BB[1][5] = Xo*Yo;
        BB[1][6] = Zo*Xo;
        BB[1][7] = Zo*Yo;
        BB[1][8] = Zo*Yo*Xo;
        BB[2][0] = Zo;
        BB[2][1] = 1;
        BB[2][2] = Xo;
        BB[2][3] = Yo;
        BB[2][4] = Zo;
        BB[2][5] = Xo*Yo;
        BB[2][6] = Zo*Xo;
        BB[2][7] = Zo*Yo;
        BB[2][8] = Zo*Yo*Xo;
        LL[0][0] = Xd;
        LL[1][1] = Yd;
        LL[2][2] = Zd;
        double[][] w = multiplyAB(getA_T(BB),LL);
        return w;
    }
}


