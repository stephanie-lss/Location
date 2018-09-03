import static java.lang.Math.*;

public class CoordinateTransform {
    static double PI = 3.1415926;

    /**
     * ������������
     *
     * @param flatting_backward ���ʵ���
     * @return flatting:0.003324449296662885; eccentricity:0.0814729809826527; second_eccentricity:0.08174473721679898
     */
    public static double[] element(double flatting_backward) {
        double parameter[] = new double[3];
        double flatting = 1 / flatting_backward;
        double eccentricity = sqrt(2 * flatting - flatting * flatting);
        double second_eccentricity = sqrt((eccentricity * eccentricity) / (1 - eccentricity * eccentricity));
        parameter[0] = flatting;
        parameter[1] = eccentricity;
        parameter[2] = second_eccentricity;
        return parameter;
    }

    /**
     * ��˹-����-3�ȴ�
     *
     * @param theat             �����γ��
     * @param lamd              ����㾭��
     * @param semi_major        ������
     * @param flatting_backward ���ʵ���
     * @param theat_orignal     ��ʼԭ��γ��
     * @param lamd_orignal      ��ʼԭ�㾭��
     * @param E_offset          ��αƫ��
     * @param N_offset          ��αƫ��
     * @param scale_factor      �߶�����
     * @return E��N ����
     */
    public static double[] Gauss_Kruger_3(double theat, double lamd, double semi_major, double flatting_backward, double theat_orignal, double lamd_orignal, double E_offset, double N_offset, double scale_factor) {
        int bandwith = 3;
        double[] parameter = element(flatting_backward);
        double T = tan(theat) * tan(theat);
        double C = parameter[1] * parameter[1] * cos(theat) * cos(theat) / (1 - parameter[1] * parameter[1]);
        double A = (lamd - lamd_orignal) * cos(theat);
        double V = semi_major / (sqrt(1 - parameter[1] * parameter[1] * sin(theat) * sin(theat)));
        double M = semi_major * ((1 - parameter[1] * parameter[1] / 4 - 3 * pow(parameter[1], 4) / 64 - 5 * pow(parameter[1], 6) / 256) * theat - (3 * parameter[1] * parameter[1] / 8 + 3 * pow(parameter[1], 4) / 32 + 45 * pow(parameter[1], 6) / 1024) * sin(2 * theat) + (15 * pow(parameter[1], 4) / 256 + 45 * pow(parameter[1], 6) / 1024) * sin(4 * theat) - (35 * pow(parameter[1], 6) / 3072) * sin(6 * theat));
        double Mo = semi_major * ((1 - parameter[1] * parameter[1] / 4 - 3 * pow(parameter[1], 4) / 64 - 5 * pow(parameter[1], 6) / 256) * theat_orignal - (3 * parameter[1] * parameter[1] / 8 + 3 * pow(parameter[1], 4) / 32 + 45 * pow(parameter[1], 6) / 1024) * sin(2 * theat_orignal) + (15 * pow(parameter[1], 4) / 256 + 45 * pow(parameter[1], 6) / 1024) * sin(4 * theat_orignal) - (35 * pow(parameter[1], 6) / 3072) * sin(6 * theat_orignal));
        double E = E_offset + scale_factor * V * (A + (1 - T + C) * pow(A, 3) / 6 + (5 - 18 * T + T * T + 72 * C - 58 * pow(parameter[2], 2)) * pow(A, 5) / 120);
        double N = N_offset + scale_factor * (M - Mo + V * tan(theat) * (A * A / 2 + (5 - T + 9 * C + 4 * pow(C, 2)) * pow(A, 4) / 24 + (61 - 58 * T + pow(T, 2) + 600 * C - 330 * pow(parameter[2], 2)) * pow(A, 6) / 720));
        double Zone = Math.floor(((lamd_orignal + bandwith) / 3));
        double[] parameter2 = new double[3];
        parameter2[0] = E;
        parameter2[1] = N;
        parameter2[2] = Zone;
        return parameter2;
    }

    /**
     * ��˹-����-6�ȴ�
     *
     * @param theat             �����γ��
     * @param lamd              ����㾭��
     * @param semi_major        ������
     * @param flatting_backward ���ʵ���
     * @param theat_orignal     ��ʼԭ��γ��
     * @param lamd_orignal      ��ʼԭ�㾭��
     * @param E_offset          ��αƫ��
     * @param N_offset          ��αƫ��
     * @param scale_factor      �߶�����
     * @return E, N ����
     */
    public static double[] Gauss_Kruger_6(double theat, double lamd, double semi_major, double flatting_backward, double theat_orignal, double lamd_orignal, double E_offset, double N_offset, double scale_factor) {
        int bandwith = 6;
        double[] parameter = element(flatting_backward);
        double T = tan(theat) * tan(theat);
        double C = parameter[1] * parameter[1] * cos(theat) * cos(theat) / (1 - parameter[1] * parameter[1]);
        double A = (lamd - lamd_orignal) * cos(theat);
        double V = semi_major / (sqrt(1 - parameter[1] * parameter[1] * sin(theat) * sin(theat)));
        double M = semi_major * ((1 - parameter[1] * parameter[1] / 4 - 3 * pow(parameter[1], 4) / 64 - 5 * pow(parameter[1], 6) / 256) * theat - (3 * parameter[1] * parameter[1] / 8 + 3 * pow(parameter[1], 4) / 32 + 45 * pow(parameter[1], 6) / 1024) * sin(2 * theat) + (15 * pow(parameter[1], 4) / 256 + 45 * pow(parameter[1], 6) / 1024) * sin(4 * theat) - (35 * pow(parameter[1], 6) / 3072) * sin(6 * theat));
        double Mo = semi_major * ((1 - parameter[1] * parameter[1] / 4 - 3 * pow(parameter[1], 4) / 64 - 5 * pow(parameter[1], 6) / 256) * theat_orignal - (3 * parameter[1] * parameter[1] / 8 + 3 * pow(parameter[1], 4) / 32 + 45 * pow(parameter[1], 6) / 1024) * sin(2 * theat_orignal) + (15 * pow(parameter[1], 4) / 256 + 45 * pow(parameter[1], 6) / 1024) * sin(4 * theat_orignal) - (35 * pow(parameter[1], 6) / 3072) * sin(6 * theat_orignal));
        double E = E_offset - scale_factor * V * (A + (1 - T + C) * pow(A, 3) / 6 + (5 - 18 * T + T * T + 72 * C - 58 * pow(parameter[2], 2)) * pow(A, 5) / 120);
        double N = N_offset - scale_factor * (M - Mo + V * tan(theat) * (A * A / 2 + (5 - T + 9 * C + 4 * pow(C, 2)) * pow(A, 4) / 24 + (61 - 58 * T + pow(T, 2) + 600 * C - 330 * pow(parameter[2], 2)) * pow(A, 6) / 720));
        double Zone = Math.floor(((lamd_orignal + bandwith) / 3));
        double[] parameter3 = new double[3];
        parameter3[0] = E;
        parameter3[1] = N;
        parameter3[2] = Zone;
        return parameter3;
    }

    /**
     * ������б��ī����ͶӰ
     *
     * @param a                 ������
     * @param theat             �����γ��
     * @param lamd              ����㾭��
     * @param flatting_backward ���ʵ���
     * @param theat_orignal     ͶӰ����γ��
     * @param lamd_orignal      ͶӰ���ľ���
     * @param Angle_grid        ����б��λ��
     * @param Azimuth           ��ʼ��λ��
     * @param E_offset          ��αƫ��
     * @param N_offset          ��αƫ��
     * @param scale_factor      �߶�����
     * @return ͶӰ������� E��N
     */
    public static double[] Hotine_Oblique_Mercator(double a,double theat, double lamd, double flatting_backward, double theat_orignal, double lamd_orignal, double Azimuth, double Angle_grid, double E_offset, double N_offset, double scale_factor) {
        double F;
        double[] parameter = element(flatting_backward);
        double B = pow((1 + (parameter[1] * parameter[1] * pow(cos(theat_orignal), 4) / (1 - pow(parameter[1], 2)))), 0.5);
        double A = a * B * scale_factor * pow((1 - parameter[1] * parameter[1]), 0.5) / (1 - pow(parameter[1], 2) * pow(sin(theat_orignal), 2));
        double to = tan(PI / 4 - theat_orignal / 2) / pow(((1 - parameter[1] * sin(theat_orignal)) / (1 + parameter[1] * sin(theat_orignal))), (parameter[1] / 2));
        double D = B * pow((1 - parameter[1] * parameter[1]), 0.5) / (cos(theat_orignal) * pow((1 - pow(parameter[1], 2) * pow(sin(theat_orignal), 2)), 0.5));
        if (D < 1) {
            double D2 = 1;
            F = D2 + pow((D2 - 1), 0.5) * signum(theat_orignal);
        } else {
            F = D + pow((pow(D, 2) - 1), 0.5) * signum(theat_orignal);
        }
        double H = F * pow(to, B);
        double G = (F - 1 / F) / 2;
        double gamo = asin(sin(Azimuth) / D);
        double lamdo = lamd_orignal - (asin(G * tan(gamo))) / B;
        double t = tan(PI / 4 - theat / 2) / pow(((1 - parameter[1] * sin(theat)) / (1 + parameter[1] * sin(theat))), (parameter[1] / 2));
        double Q = H / pow(t, B);
        double S = (Q - 1 / Q) / 2;
        double T = (Q + 1 / Q) / 2;
        double V = sin(B * (lamd - lamdo));
        double U = (-1 * V * cos(gamo) + S * sin(gamo)) / T;
        double v = A * log((1 - U) / (1 + U)) / (2 * B);
        double u = (A * atan((S * cos(gamo) + V * sin(gamo)) / cos(B * (lamd - lamdo))) / B);
        double E = E_offset + v * cos(Angle_grid) + u * sin(Angle_grid);
        double N = N_offset + u * cos(Angle_grid) - v * sin(Angle_grid);
        double[] parameter4 = {E, N};
        return parameter4;
    }

    /**
     * ī����ͶӰ
     *
     * @param theat             �����γ��
     * @param lamd              ����㾭��
     * @param semi_major        ������
     * @param flatting_backward ���ʵ���
     * @param lamd_orignal      ��һ��׼γ��γ��
     * @param theat1            ��ʼԭ�㾭��
     * @param E_offset          ��αƫ��
     * @param N_offset          ��αƫ��
     * @return ͶӰ������� E��N
     */
    public static double[] Mercator(double theat, double lamd, double semi_major, double flatting_backward, double lamd_orignal, double theat1, double E_offset, double N_offset) {
        double[] parameter = element(flatting_backward);
        double ko = cos(theat1) / sqrt(1 - parameter[1] * parameter[1] * sin(theat1) * sin(theat1));
        double E = E_offset + semi_major * ko * (lamd - lamd_orignal);
        double N = N_offset + semi_major * ko * log(tan(PI / 4 + theat / 2) * pow(((1 - parameter[1] * sin(theat)) / (1 + parameter[1] * sin(theat))), (parameter[1] / 2)));
        double[] parameter5 = {E, N};
        return parameter5;
    }

    /**
     * б��ī����ͶӰ
     *
     * @param a                 ������
     * @param theat             ���������γ��
     * @param lamd              ��������꾭��
     * @param flatting_backward ���ʵ���
     * @param theat_orignal     ͶӰ����γ��
     * @param lamd_orignal      ͶӰ���ľ���
     * @param Azimuth           ��ʼ�߷�λ
     * @param Angle_grid        ����б��λ��
     * @param E_offset          ͶӰ���Ķ�αƫ��
     * @param N_offset          ͶӰ���ı�αƫ��
     * @param scale_factor      ��ʼ�߱�������
     * @return ͶӰ������� E��N
     */
    public static double[] Oblique_Mercator(double a,double theat, double lamd, double flatting_backward, double theat_orignal, double lamd_orignal, double Azimuth, double Angle_grid, double E_offset, double N_offset, double scale_factor) {
        double F;
        double[] parameter = element(flatting_backward);
        double B = pow((1 + (parameter[1] * parameter[1] * pow(cos(theat_orignal), 4) / (1 - pow(parameter[1], 2)))), 0.5);
        double A = a * B * scale_factor * pow((1 - parameter[1] * parameter[1]), 0.5) / (1 - pow(parameter[1], 2) * pow(sin(theat_orignal), 2));
        double to = tan(PI / 4 - theat_orignal / 2) / pow(((1 - parameter[1] * sin(theat_orignal)) / (1 + parameter[1] * sin(theat_orignal))), (parameter[1] / 2));
        double D = B * pow((1 - parameter[1] * parameter[1]), 0.5) / (cos(theat_orignal) * pow((1 - pow(parameter[1], 2) * pow(sin(theat_orignal), 2)), 0.5));
        if (D < 1) {
            double D2 = 1;
            F = D2 + pow((D2 - 1), 0.5) * signum(theat_orignal);
        } else {
            F = D + pow((pow(D, 2) - 1), 0.5) * signum(theat_orignal);
        }
        double H = F * pow(to, B);
        double G = (F - 1 / F) / 2;
        double gamo = asin(sin(Azimuth) / D);
        double lamdo = lamd_orignal - (asin(G * tan(gamo))) / B;
        double Uc;
        if (Azimuth == PI / 2) {
            Uc = A * (lamd_orignal - lamdo);
        } else {
            Uc = (A / B) * atan(pow((D * D - 1), 0.5) / cos(Azimuth)) * signum(theat_orignal);
        }
        double t = tan(PI / 4 - theat / 2) / pow(((1 - parameter[1] * sin(theat)) / (1 + parameter[1] * sin(theat))), (parameter[1] / 2));
        double Q = H / pow(t, B);
        double S = (Q - 1 / Q) / 2;
        double T = (Q + 1 / Q) / 2;
        double V = sin(B * (lamd - lamdo));
        double U = (-1 * V * cos(gamo) + S * sin(gamo)) / T;
        double v = A * Math.log((1 - U) / (1 + U)) / (2 * B);
        double u = (A * atan((S * cos(gamo) + V * sin(gamo)) / cos(B * (lamd - lamdo))) / B) - (Math.abs(Uc) * signum(theat_orignal));
        double E = E_offset + v * cos(Angle_grid) + u * sin(Angle_grid);
        double N = N_offset + u * cos(Angle_grid) - v * sin(Angle_grid);
        double[] parameter6 = {E, N};
        return parameter6;
    }

    /**
     * ������ͶӰ�ͳ������Ͷ
     *
     * @param theat             �����γ��
     * @param lamd              ����㾭��
     * @param semi_major        ������
     * @param flatting_backward ���ʵ���
     * @param theat_orignal     ��ʼԭ��γ��
     * @param lamd_orignal      ��ʼԭ�㾭��
     * @param E_offset          ��αƫ��
     * @param N_offset          ��αƫ��
     * @return ͶӰ������� E��N
     */
    public static double[] Oblique_Stereographic(double theat, double lamd, double semi_major, double flatting_backward, double theat_orignal, double lamd_orignal, double E_offset, double N_offset) {
        double[] parameter = element(flatting_backward);
        double R = semi_major * sqrt(1 - parameter[1] * parameter[1]) / (1 - parameter[1] * parameter[1] * sin(theat_orignal) * sin(theat_orignal));
        double n = sqrt(1 + (pow((parameter[1] * cos(theat_orignal)), 4) / (1 - parameter[1] * parameter[1])));
        double S1 = (1 + sin(theat_orignal)) / (1 - sin(theat_orignal));
        double S2 = (1 - parameter[1] * sin(theat_orignal)) / (1 + parameter[1] * sin(theat_orignal));
        double w1 = pow((S1 * pow((S2), parameter[1])), n);
        double sin_axo = (w1 - 1) / (w1 + 1);
        double c = ((n + sin(theat_orignal)) * (1 - sin_axo)) / ((n - sin(theat_orignal)) * (1 + sin_axo));
        double w2 = c * w1;
        double axo = asin((w2 - 1) / (w2 + 1));
        double lamd_o = lamd_orignal;
        double lamd_A = n * (lamd - lamd_o) + lamd_o;
        double Sa = (1 + sin(theat)) / (1 - sin(theat));
        double Sb = (1 - parameter[1] * sin(theat)) / (1 + parameter[1] * sin(theat));
        double w = c * pow((Sa * pow((Sb), parameter[1])), n);
        double ax = asin((w - 1) / (w + 1));
        double B = (1 + sin(ax) * sin(axo) + cos(ax) * cos(axo) * cos(lamd_A - lamd_o));
        double ko = 0.994;
        double E = E_offset + 2 * R * ko * cos(ax) * sin(lamd_A - lamd_o) / B;
        double N = N_offset + 2 * R * ko * (sin(ax) * cos(axo) - cos(ax) * sin(axo) * cos(lamd_A - lamd_o)) / B;
        double[] parameter7 = {E, N};
        return parameter7;
    }

    /**
     * ����ī����ͶӰ
     *
     * @param theat             �����γ��
     * @param lamd              ����㾭��
     * @param semi_major        ������
     * @param flatting_backward ���ʵ���
     * @param theat_orignal     ��ʼԭ��γ��
     * @param lamd_orignal      ��ʼԭ�㾭��
     * @param E_offset          ��αƫ��
     * @param N_offset          ��αƫ��
     * @param scale_factor      ��ʼԭ���������
     * @return ͶӰ������� E��N
     */
    		public static double[] TransverseMercator(double theat,double lamd,double semi_major,double flatting_backward,double theat_orignal,double lamd_orignal,double E_offset,double N_offset,double scale_factor){
		double[] ele=element(flatting_backward);
		double T=tan(theat)*tan(theat);
		double C=ele[1]*ele[1]*cos(theat)*cos(theat)/(1-ele[1]*ele[1]);
		double A=(lamd-lamd_orignal)*cos(theat);
		double V=semi_major/(sqrt(1-ele[1]*ele[1]*sin(theat)*sin(theat)));
		double M=semi_major*( (1-ele[1]*ele[1]/4-3*pow(ele[1], 4)/64-5*pow(ele[1], 6)/256)*theat-(3*ele[1]*ele[1]/8+3*pow(ele[1], 4)/32+45*pow(ele[1], 6)/1024)*sin(2*theat)+(15*pow(ele[1], 4)/256+45*pow(ele[1], 6)/1024)*sin(4*theat)-(35*pow(ele[1], 6)/3072)*sin(6*theat) );
		double Mo=semi_major*( (1-ele[1]*ele[1]/4-3*pow(ele[1], 4)/64-5*pow(ele[1], 6)/256)*theat_orignal-(3*ele[1]*ele[1]/8+3*pow(ele[1], 4)/32+45*pow(ele[1], 6)/1024)*sin(2*theat_orignal)+(15*pow(ele[1], 4)/256+45*pow(ele[1], 6)/1024)*sin(4*theat_orignal)-(35*pow(ele[1], 6)/3072)*sin(6*theat_orignal) );
		double E=E_offset+scale_factor*V*(A+(1-T+C)*pow(A, 3)/6+(5-18*T+T*T+72*C-58*pow(ele[2], 2))*pow(A, 5)/120);	
		double N=N_offset+scale_factor*( M-Mo+V*tan(theat)*(A*A/2+(5-T+9*C+4*pow(C, 2))*pow(A, 4)/24+(61-58*T+pow(T, 2)+600*C-330*pow(ele[2], 2))*pow(A, 6)/720));
		double[] TransverseMercator=new double[2];
		TransverseMercator[0]=E;
		TransverseMercator[1]=N;
		return TransverseMercator;
			    
	}

    /**
     * ��б��ƽͶӰ
     *
     * @param theat             ���������γ�ȣ��Ի���Ϊ��λ��
     * @param lamd              ��������꾫�ȣ��Ի���Ϊ��λ��
     * @param semi_major        ���򳤰���
     * @param flatting_backward ���ʵ���
     * @param theat_orignal     ��ʼԭ��γ��
     * @param lamd_orignal      ��ʼԭ�㾭��
     * @param E_offset          ��αƫ��
     * @param N_offset          ��αƫ��
     * @param scale_factor      ��ʼԭ���������
     * @return ͶӰ������� E��N
     */
    public static double[] WebMercatot(double theat, double lamd, double semi_major, double flatting_backward, double theat_orignal, double lamd_orignal, double E_offset, double N_offset) {
        //double[] ele=element(flatting_backward);
        double R = semi_major;
        double E = E_offset + R * (lamd - lamd_orignal);
        double N = N_offset + R * log(tan(PI / 4 + theat / 2));
        double[] WebMercatot = new double[2];
        WebMercatot[0] = E;
        WebMercatot[1] = N;
        return WebMercatot;
    }

    /**
     * ����ī����ͶӰ���ϲ�ר�ð汾)
     *
     * @param theat             ���������γ�ȣ��Ի���Ϊ��λ��
     * @param lamd              ��������꾭�ȣ��Ի���Ϊ��λ��
     * @param semi_major        ���򳤰���
     * @param flatting_backward ���ʵ���
     * @param theat_orignal     ��ʼԭ��γ�ȣ��Ի���Ϊ��λ��
     * @param lamd_orignal      ��ʼԭ�㾭�ȣ��Ի���Ϊ��λ��
     * @param E_offset          ��αƫ��
     * @param N_offset          ��αƫ��
     * @param scale_factor      ��ʼԭ���������
     * @return ͶӰ������� E��N
     */
    public static double[] SouthOrientatedTransverseMercator(double theat, double lamd, double semi_major, double flatting_backward, double theat_orignal, double lamd_orignal, double E_offset, double N_offset, double scale_factor) {
        double[] ele = element(flatting_backward);
        double T = tan(theat) * tan(theat);
        double C = ele[1] * ele[1] * cos(theat) * cos(theat) / (1 - ele[1] * ele[1]);
        double A = (lamd - lamd_orignal) * cos(theat);
        double V = semi_major / (sqrt(1 - ele[1] * ele[1] * sin(theat) * sin(theat)));
        double M = semi_major * ((1 - ele[1] * ele[1] / 4 - 3 * pow(ele[1], 4) / 64 - 5 * pow(ele[1], 6) / 256) * theat - (3 * ele[1] * ele[1] / 8 + 3 * pow(ele[1], 4) / 32 + 45 * pow(ele[1], 6) / 1024) * sin(2 * theat) + (15 * pow(ele[1], 4) / 256 + 45 * pow(ele[1], 6) / 1024) * sin(4 * theat) - (35 * pow(ele[1], 6) / 3072) * sin(6 * theat));
        double Mo = semi_major * ((1 - ele[1] * ele[1] / 4 - 3 * pow(ele[1], 4) / 64 - 5 * pow(ele[1], 6) / 256) * theat_orignal - (3 * ele[1] * ele[1] / 8 + 3 * pow(ele[1], 4) / 32 + 45 * pow(ele[1], 6) / 1024) * sin(2 * theat_orignal) + (15 * pow(ele[1], 4) / 256 + 45 * pow(ele[1], 6) / 1024) * sin(4 * theat_orignal) - (35 * pow(ele[1], 6) / 3072) * sin(6 * theat_orignal));
        double E = E_offset + scale_factor * V * (A + (1 - T + C) * pow(A, 3) / 6 + (5 - 18 * T + T * T + 72 * C - 58 * pow(ele[2], 2)) * pow(A, 5) / 120);
        double N = N_offset + scale_factor * (M - Mo + V * tan(theat) * (A * A / 2 + (5 - T + 9 * C + 4 * pow(C, 2)) * pow(A, 4) / 24 + (61 - 58 * T + pow(T, 2) + 600 * C - 330 * pow(ele[2], 2)) * pow(A, 6) / 720));
        double[] SouthOrientatedTransverseMercator = new double[2];
        SouthOrientatedTransverseMercator[0] = E;
        SouthOrientatedTransverseMercator[1] = N;
        return SouthOrientatedTransverseMercator;
    }

    /**
     * ˫��׼γ�������صȽ�Բ׶ͶӰ
     *
     * @param theat             ���������γ�ȣ��Ի���Ϊ��λ��
     * @param lamd              ��������꾫�ȣ��Ի���Ϊ��λ��
     * @param semi_major        ����������򳤰���
     * @param flatting_backward ����������ʵ���
     * @param theat1            ͶӰ���� ��һ��׼γ��γ��
     * @param theat2            ͶӰ���� �ڶ���׼γ��γ��
     * @param theat_f           ͶӰ����αԭ��γ��
     * @param lamd_f            ͶӰ����αԭ�㾭��
     * @param E_offset          ͶӰ����αԭ�㶫ƫ
     * @param N_offset          ͶӰ����αԭ�㱱ƫ
     * @return ͶӰ������� E��N
     */
    public static double[] LambertConicConformal_2SP(double theat, double lamd, double semi_major, double flatting_backward, double theat1, double theat2, double theat_f, double lamd_f, double E_offset, double N_offset) {
        double[] ele = element(flatting_backward);
        double m1 = cos(theat1) / sqrt(1 - ele[1] * ele[1] * sin(theat1) * sin(theat1));
        double m2 = cos(theat2) / sqrt(1 - ele[1] * ele[1] * sin(theat2) * sin(theat2));
        double t1 = tan(PI / 4 - theat1 / 2) / pow(((1 - ele[1] * sin(theat1)) / (1 + ele[1] * sin(theat1))), ele[1] / 2);
        double t2 = tan(PI / 4 - theat2 / 2) / pow(((1 - ele[1] * sin(theat2)) / (1 + ele[1] * sin(theat2))), ele[1] / 2);
        double tf = tan(PI / 4 - theat_f / 2) / pow(((1 - ele[1] * sin(theat_f)) / (1 + ele[1] * sin(theat_f))), ele[1] / 2);
        double t = tan(PI / 4 - theat / 2) / pow(((1 - ele[1] * sin(theat)) / (1 + ele[1] * sin(theat))), ele[1] / 2);
        double n = (log(m1) - log(m2)) / (log(t1) - log(t2));
        double F = m1 / (n * pow(t1, n));
        double r = semi_major * F * pow(t, n);
        double rf = semi_major * F * pow(tf, n);
        double phi = n * (lamd - lamd_f);
        double E = E_offset + r * sin(phi);
        double N = N_offset + rf - r * cos(phi);
        double[] LambertConicConformal_2SP = new double[2];
        LambertConicConformal_2SP[0] = E;
        LambertConicConformal_2SP[1] = N;
        return LambertConicConformal_2SP;
    }

    /**
     * ����׼γ�������صȽ�Բ׶ͶӰ
     *
     * @param theat             ���������γ�ȣ��Ի���Ϊ��λ��
     * @param lamd              ��������꾫�ȣ��Ի���Ϊ��λ��
     * @param semi_major        ���򳤰���
     * @param flatting_backward ���ʵ���
     * @param theat_original    ��ʼԭ��γ��
     * @param lamd_original     ��ʼԭ�㾭��
     * @param E_offset          ��αƫ��
     * @param N_offset          ��αƫ��
     * @param scale_factor      ��ʼԭ���������
     * @return ͶӰ������� E��N
     */
    public static double[] LambertConicConformal_1SP(double theat, double lamd, double semi_major, double flatting_backward, double theat_original, double lamd_original, double E_offset, double N_offset, double scale_factor) {
        double[] ele = element(flatting_backward);
        double mo = cos(theat_original) / sqrt(1 - ele[1] * ele[1] * sin(theat_original) * sin(theat_original));
        double to = tan(PI / 4 - theat_original / 2) / pow(((1 - ele[1] * sin(theat_original)) / (1 + ele[1] * sin(theat_original))), ele[1] / 2);
        double t = tan(PI / 4 - theat / 2) / pow(((1 - ele[1] * sin(theat)) / (1 + ele[1] * sin(theat))), ele[1] / 2);
        double n = sin(theat_original);
        double F = mo / (n * pow(to, n));
        double r = semi_major * F * pow(t, n) * scale_factor;
        double ro = semi_major * F * pow(to, n) * scale_factor;
        double phi = n * (lamd - lamd_original);
        double E = E_offset + r * sin(phi);
        double N = N_offset + ro - r * cos(phi);
        double[] LambertConicConformal_1SP = new double[2];
        LambertConicConformal_1SP[0] = E;
        LambertConicConformal_1SP[1] = N;
        return LambertConicConformal_1SP;
    }

    /**
     * ������-��������ͶӰ
     *
     * @param theat             ���������γ�ȣ��Ի���Ϊ��λ��
     * @param lamd              ��������꾫�ȣ��Ի���Ϊ��λ��
     * @param semi_major        ���򳤰���
     * @param flatting_backward ���ʵ���
     * @param theat_orignal     ��ʼԭ��γ��
     * @param lamd_orignal      ��ʼԭ�㾭��
     * @param E_offset          ��αƫ��
     * @param N_offset          ��αƫ��
     * @return E: 66644.73292887252;N: 2.8278444756823532E7
     */
    public static double[] CassiniSoldner(double theat, double lamd, double semi_major, double flatting_backward, double theat_orignal, double lamd_orignal, double E_offset, double N_offset) {
        double[] ele = element(flatting_backward);
        double M = semi_major * ((1 - ele[1] * ele[1] / 4 - 3 * pow(ele[1], 4) / 64 - 5 * pow(ele[1], 6) / 256) * theat - (3 * pow(ele[1], 2) / 8 + 3 * pow(ele[1], 4) / 32 + 45 * pow(ele[1], 6) / 1024) * sin(2 * theat) + (15 * pow(ele[1], 4) / 256 + 45 * pow(ele[1], 6) / 1024) * sin(4 * theat) - (35 * pow(ele[1], 6) / 3072) * sin(6 * theat));
        double Mo = semi_major * ((1 - ele[1] * ele[1] / 4 - 3 * pow(ele[1], 4) / 64 - 5 * pow(ele[1], 6) / 256) * theat_orignal - (3 * pow(ele[1], 2) / 8 + 3 * pow(ele[1], 4) / 32 + 45 * pow(ele[1], 6) / 1024) * sin(2 * theat_orignal) + (15 * pow(ele[1], 4) / 256 + 45 * pow(ele[1], 6) / 1024) * sin(4 * theat_orignal) - (35 * pow(ele[1], 6) / 3072) * sin(6 * theat_orignal));
        double A = (lamd - lamd_orignal) * cos(theat);
        double T = tan(theat) * tan(theat);
        double C = ele[1] * ele[1] * cos(theat) * cos(theat) * (1 - ele[1] * ele[1]);
        double v = semi_major / pow(1 - ele[1] * ele[1] * sin(theat) * sin(theat), 0.5);
        double E = E_offset + v * (A - T * pow(A, 3) / 6 - (8 - T + 8 * C) * T * pow(A, 5) / 120);
        double X = M - Mo + v * tan(theat) * (A * A / 2 + (5 - T + 6 * C) * pow(A, 4) / 24);
        double N = N_offset + X;
        double[] CassiniSoldner = new double[2];
        CassiniSoldner[0] = E;
        CassiniSoldner[1] = N;
        return CassiniSoldner;
    }
}


