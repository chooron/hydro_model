package tranditionModel;

import optimization.OptimizerByGA;
import tranditionModel.parent.parentModel;

import java.util.ArrayList;

/**
 * 新安江水文模型
 *
 * @author jing xin(增加应用于实际的功能)
 * @cite Wenxuan in github
 * @Email 1303541772@qq.com
 * @phone 17629150947
 */
public class XajModel implements parentModel {
    /*土壤蓄水量参数*/
    private double WUM;     // 上层土壤蓄水容量
    private double WLM;     // 下层土壤蓄水容量
    private double WDM;     // 底层土壤蓄水容量
    private double WM;      // 土壤总蓄水容量

    /*蒸散发参数*/
    private double K;       // 蒸发皿折算系数
    private double C;       // 深层蒸散发系数

    /*产流参数*/
    private double B;       // 蓄水容量曲线的指数
    private double Imp;     // 不透水面积与全流域面积的比值

    /*水源划分参数*/
    private double SM;      // 流域平均自由水蓄水容量
    private double EX;      // 自由水蓄水容量曲线指数
    private double KSS;     // 自由水蓄水库对壤中流的出流系数
    private double KG;      // 自由水蓄水库对地下径流出流系数

    /*汇流参数*/
    private double KKSS;    // 地下水库消退系数
    private double KKG;     // 壤中流水库消退系数
    private double Area;    // 水文单元面积

    private RunoffGenerationResult runoffGenerationResult;          // 产流计算结果
    private SourcePartitionResult sourcePartitionResult;            // 水源划分计算结果
    private RunoffConcentrationResult runoffConcentrationResult;    // 汇流计算结果

    ArrayList<state> states;
    double[] others;

    public ArrayList<state> getStates() {
        return states;
    }

    public void setStates(ArrayList<state> states) {
        this.states = states;
    }

    public double[] getOthers() {
        return others;
    }

    public void setOthers(double[] others) {
        this.others = others;
    }

    public void setArea(double area) {
        Area = area;
    }

    /**
     * 获取模型的计算结果
     *
     * @return 模型的最终计算结果（即汇流计算结果）
     */
    public double[] GetResult() {
        return runoffConcentrationResult.Q;
    }

    public static double[] getBestPara(parentModel model, int populationSize, int evolveTime, double[] target, double[] initPara, double threshold){
        OptimizerByGA ga = new OptimizerByGA(model, populationSize, evolveTime,target, initPara, threshold);
        return ga.GAmain();
    }

    public double[] paraOptimize(double[] paras){
        XajModel model = this;
        model = model.SetSoilWaterStorageParam(paras[0], paras[1], paras[2]);
        model = model.SetEvapotranspirationParam(paras[3], paras[4]);
        model = model.SetRunoffGenerationParam(paras[5], paras[6]);
        model = model.SetSourcePartitionParam(paras[7], paras[8], paras[9], paras[10]);
        model = model.SetRunoffConcentrationParam(paras[11], paras[12], Area);
        double[] P = new double[model.states.size()];
        double[] EI = new double[model.states.size()];
        for (int i = 0; i < model.states.size(); i++) {
            P[i] = model.states.get(i).P;
            EI[i] = model.states.get(i).E;
        }
        model = model.ComputeRunoffGeneration(P, EI, model.others[0], model.others[1], model.others[2]);
        model = model.ComputeSourcePartition(model.others[3], model.others[4]);
        model = model.ComputeRunoffConcentration(model.others[5], model.others[6], model.others[4]);
        return model.GetResult();
    }

    public static class state{
        double P;
        double E;
        double R;
        double[] other;

        public void setP(double p) {
            P = p;
        }

        public void setE(double e) {
            E = e;
        }

        public void setR(double r) {
            R = r;
        }

        public void setOther(double[] other) {
            this.other = other;
        }

        public double getP() {
            return P;
        }

        public double getE() {
            return E;
        }

        public double getR() {
            return R;
        }

        public double[] getOther() {
            return other;
        }
    }

    /**
     * 产流计算结果
     */
    public static class RunoffGenerationResult {
        /**
         * 序列长度
         */
        public final int Length;

        /**
         * 序列长度
         */
        public double[] PE;

        /**
         * 产流量序列
         */
        public double[] R;

        /**
         * 土壤含水量序列
         */
        public double[] W;

        /**
         * @param n 序列长度
         */
        public RunoffGenerationResult(int n) {
            PE = new double[n];
            R = new double[n];
            W = new double[n];
            Length = n;
        }

        /**
         * 设置序列指定位置处的产流计算结果
         *
         * @param idx 序列位置索引
         * @param pe  净雨量
         * @param r   产流量
         * @param w   土壤含水量
         */
        public void Set(int idx, double pe, double r, double w) {
            PE[idx] = pe;
            R[idx] = r;
            W[idx] = w;
        }
    }

    /**
     * 水源划分计算结果
     */
    static class SourcePartitionResult {
        /**
         * 序列长度
         */
        public final int Length;

        /**
         * 不透水面积产流序列
         */
        public double[] RIMP;

        /**
         * 地表径流产流序列
         */
        public double[] RS;

        /**
         * 壤中流产流序列
         */
        public double[] RSS;

        /**
         * 地下径流产流序列
         */
        public double[] RG;

        /**
         * 流域平均自由含水量序列
         */
        public double[] S;

        /**
         * @param n 序列长度
         */
        public SourcePartitionResult(int n) {
            RIMP = new double[n];
            RS = new double[n];
            RSS = new double[n];
            RG = new double[n];
            S = new double[n];
            Length = n;
        }

        /**
         * 设置序列指定位置处的水源划分计算结果
         *
         * @param idx  序列位置索引
         * @param rimp 不透水面积产流
         * @param rs   地表径流产流
         * @param rss  壤中流产流
         * @param rg   地下径流产流
         * @param s    流域平均自由含水量
         */
        public void Set(int idx, double rimp, double rs, double rss, double rg, double s) {
            RIMP[idx] = rimp;
            RS[idx] = rs;
            RSS[idx] = rss;
            RG[idx] = rg;
            S[idx] = s;
        }
    }

    /**
     * 汇流计算结果
     */
    public static class RunoffConcentrationResult {
        /**
         * 序列长度
         */
        public final int Length;

        /**
         * 地表径流序列
         */
        public double[] QRS;

        /**
         * 壤中流序列
         */
        public double[] QRSS;

        /**
         * 地下径流序列
         */
        public double[] QRG;

        /**
         * 总径流序列
         */
        public double[] Q;

        /**
         * @param QRS  地表径流序列
         * @param QRSS 壤中流序列
         * @param QRG  地下径流序列
         * @param Q    总径流序列
         */
        public RunoffConcentrationResult(double[] QRS, double[] QRSS,
                                         double[] QRG, double[] Q) {
            this.Length = QRS.length;

            this.QRS = QRS;
            this.QRSS = QRSS;
            this.QRG = QRG;
            this.Q = Q;
        }
    }

    /**
     * 执行产流计算
     *
     * @param P   降雨量序列
     * @param EI  蒸发皿蒸发量序列
     * @param wu0 初始时刻上层土壤含水量
     * @param wl0 初始时刻下层土壤含水量
     * @param wd0 初始时刻深层土壤含水量
     * @return 完成产流计算的模型实例
     */
    public XajModel ComputeRunoffGeneration(double[] P, double[] EI, double wu0, double wl0, double wd0) {
        int n = P.length;
        runoffGenerationResult = new RunoffGenerationResult(n);
        LayeredSoilParam w0 = new LayeredSoilParam(wu0, wl0, wd0);  // 时段初土壤分层含水量
        double w0Sum = w0.getSum();                   // 时段初土壤总含水量
        double wmmax = WM * (1 + B) / (1 - Imp);      // 流域点最大蓄水容量
        for (int i = 0; i < n; ++i) {
            double p = P[i];                          // 时段内降雨量
            double ep = EI[i] * K;                    // 计算时段内实际蒸发量
            double e = ComputeE(ep, p, w0.U, w0.L);   // 计算时段内土壤蒸发量

            double pe = p - e;                        // 计算时段内净雨量
            double r = ComputeR(pe, wmmax, w0Sum);    // 计算时段内产流量
            w0 = ComputeLayeredW(w0, pe, r);          // 计算时段末土壤分层含水量
            w0Sum = w0.getSum();                      // 计算时段末土壤总含水量

            runoffGenerationResult.Set(i, pe, r, w0Sum);
        }
        return this;
    }

    /**
     * 执行水源划分计算（需要先完成产流计算）
     *
     * @param s0 初始时刻流域平均自由含水量
     * @param dt 时段长度
     * @return 完成水源划分计算的模型实例
     */
    public XajModel ComputeSourcePartition(double s0, double dt) {

        int n = runoffGenerationResult.Length;
        sourcePartitionResult = new SourcePartitionResult(n);

        double KSSD = (1 - Math.pow(1 - (KG + KSS), dt / 24)) / (1 + KG / KSS);     // 转化后的壤中流出流系数
        double KGD = KSSD * KG / KSS;     // 转化后的地下水出流系数
        double Smax = (1 + EX) * SM;      // 流域点自由蓄水量最大值
        for (int i = 0; i < n; ++i) {
            double pe = runoffGenerationResult.PE[i];
            double w = runoffGenerationResult.W[i];
            double r = runoffGenerationResult.R[i];

            double au = Smax * (1 - Math.pow(1 - s0 / SM, 1 / (1 + EX)));
            double FR;              // 产流面积
            double coef;
            if (pe > 0) {
                FR = r / pe - Imp;
                coef = (pe + au < Smax) ? SM - SM * Math.pow(1 - (pe + au) / Smax, 1 + EX) : SM;
            } else {
                FR = 1 - Imp - Math.pow(1 - w / WM, B / (1 + B)) * (1 - Imp);
                coef = s0;
            }

            double rimp = pe * Imp;
            double rs = Math.max((pe + s0 - coef) * FR, 0.0);
            double rss = coef * KSSD * FR;
            double rg = coef * KGD * FR;
            s0 = coef * (1 - KSSD - KGD);

            sourcePartitionResult.Set(i, rimp, rs, rss, rg, s0);
        }

        return this;
    }

    /**
     * 执行汇流计算（需要先完成水源划分计算）
     *
     * @param qrss0 初始时刻壤中流流量
     * @param qrg0  初始时刻地下径流流量
     * @param dt    时段长度
     * @return 完成汇流计算的模型实例
     */
    public XajModel ComputeRunoffConcentration(double qrss0, double qrg0, double dt) {

        double[] QRS = UnitHydrograph(sourcePartitionResult.RS, sourcePartitionResult.RIMP, dt);    // 计算地表径流
        double[] QRSS = LinearReservoir(sourcePartitionResult.RSS, KKSS, qrss0, dt);    // 计算壤中流
        double[] QRG = LinearReservoir(sourcePartitionResult.RG, KKG, qrg0, dt);        // 计算地下径流
        double[] Q = ComputeQ(QRS, QRSS, QRG);                                          // 计算总径流
        runoffConcentrationResult = new RunoffConcentrationResult(QRS, QRSS, QRG, Q);

        return this;
    }

    /**
     * 设置流域土壤含水量参数
     *
     * @param wum 上层土壤含水量
     * @param wlm 下层土壤含水量
     * @param wdm 底层土壤含水量
     * @return 完成流域土壤含水量参数设置的模型实例
     */
    public XajModel SetSoilWaterStorageParam(double wum, double wlm, double wdm) {
        WUM = wum;
        WLM = wlm;
        WDM = wdm;
        WM = WUM + WLM + WDM;
        return this;
    }

    /**
     * 设置流域蒸散发参数
     *
     * @param k 蒸发皿折算系数
     * @param c 深层蒸散发系数
     * @return 完成流域蒸散发参数设置的模型实例
     */
    public XajModel SetEvapotranspirationParam(double k, double c) {
        K = k;
        C = c;
        return this;
    }

    /**
     * 设置流域产流计算参数
     *
     * @param b   蓄水容量曲线的指数
     * @param imp 不透水面积与全流域面积的比值
     * @return 完成流域产流计算参数设置的模型实例
     */
    public XajModel SetRunoffGenerationParam(double b, double imp) {
        B = b;
        Imp = imp;
        return this;
    }

    /**
     * 设置流域水源划分参数
     *
     * @param sm  流域平均自由水蓄水容量
     * @param ex  自由水蓄水容量曲线指数
     * @param kss 自由水蓄水库对壤中流的出流系数
     * @param kg  自由水蓄水库对地下径流的出流系数
     * @return 完成流域水源划分参数设置的模型实例
     */
    public XajModel SetSourcePartitionParam(double sm, double ex, double kss, double kg) {
        SM = sm;
        EX = ex;
        KSS = kss;
        KG = kg;
        return this;
    }

    /**
     * 设置流域汇流计算参数
     *
     * @param kkss 地下水库消退系数
     * @param kkg  壤中流水库消退系数
     * @param area 水文单元面积
     * @return 完成流域汇流计算参数设置的模型实例
     */
    public XajModel SetRunoffConcentrationParam(double kkss, double kkg, double area) {
        KKSS = kkss;
        KKG = kkg;
        Area = area;
        return this;
    }

    /**
     * 计算时段内的土壤蒸发量
     *
     * @param ep  时段内实际蒸发量
     * @param p   时段内降雨量
     * @param wu0 时段初上层土壤含水量
     * @param wl0 时段初下层土壤含水量
     * @return 时段内的土壤蒸发量
     */
    private double ComputeE(double ep, double p, double wu0, double wl0) {
        double eu, el, ed;
        if (p + wu0 >= ep) {
            eu = ep;
            el = ed = 0;
        } else {
            eu = p + wu0;
            if (wl0 >= C * WLM) {
                el = (ep - eu) * wl0 / WLM;
                ed = 0;
            } else if (wl0 >= C * (ep - eu)) {
                el = C * (ep - eu);
                ed = 0;
            } else {
                el = C * wl0;
                ed = C * (ep - eu) - el;
            }
        }
        return eu + el + ed;
    }

    /**
     * 计算时段内的产流量
     *
     * @param pe    时段内净雨量
     * @param wmmax 流域点最大蓄水容量
     * @param w0    时段初土壤总含水量
     * @return 时段内产流量
     */
    private double ComputeR(double pe, double wmmax, double w0) {
        double a = wmmax * (1 - Math.pow(1 - w0 / WM, 1 / (1 + B)));  // 时段初流域蓄水量w相应的纵坐标

        double r;  // 产流量
        if (pe <= 0)
            r = 0;
        else if (pe + a < wmmax)
            r = pe - WM + w0 + WM * Math.pow(1 - (pe + a) / wmmax, 1 + B);
        else
            r = pe - (WM - w0);

        return r;
    }

    /**
     * 计算时段末的土壤含水量
     *
     * @param w0 时段初土壤含水量
     * @param pe 时段内净雨量
     * @param r  时段内产流量
     * @return 时段末的土壤含水量
     */
    private LayeredSoilParam ComputeLayeredW(LayeredSoilParam w0, double pe, double r) {
        double wu0 = w0.U;
        double wl0 = w0.L;
        double wd0 = w0.D;
        double wu, wl, wd;
        double dw = pe - r;
        if (dw > 0) {
            if (wu0 + dw < WUM) {
                wu = wu0 + dw;
                wl = wl0;
                wd = wd0;
            } else {
                wu = WUM;
                wl = wl0 + dw - (WUM - wu0);
                if (wl < WLM)
                    wd = wd0;
                else {
                    wl = WLM;
                    wd = wd0 + dw - (WUM - wu0) - (WLM - wl0);
                }
            }
        } else {
            if (wu0 + dw > 0) {
                wu = wu0 + dw;
                wl = wl0;
                wd = wd0;
            } else {
                wu = 0;
                wl = wu0 + dw + wl0;
                if (wl > 0)
                    wd = wd0;
                else {
                    wl = 0;
                    wd = wu0 + wl0 + dw + wd0;
                }
            }
        }
        return new LayeredSoilParam(wu, wl, wd);
    }

    /**
     * 汇流阶段用以计算壤中流或地下径流的线性水库法
     *
     * @param R   产流序列
     * @param KK  消退系数
     * @param qr0 初始时刻产流量
     * @param dt  时段长度
     * @return 径流序列
     */
    private double[] LinearReservoir(double[] R, double KK, double qr0, double dt) {
        double[] QR = new double[R.length];
        double KKD = Math.pow(KK, dt / 24);
        for (int i = 0; i < R.length; ++i) {
            qr0 = qr0 * KKD + R[i] * Area / (3.6 * dt) * (1 - KKD);
            QR[i] = qr0;
        }
        return QR;
    }

    /**
     * 汇流阶段用以计算地表径流的单位线法
     *
     * @param RS   地表径流产流序列
     * @param RIMP 不透水面积产流序列
     * @param dt   时段长度
     * @return 径流序列
     */
    private double[] UnitHydrograph(double[] RS, double[] RIMP, double dt) {
        int n = RS.length;
        double[] result = new double[n];
        double[] UH = new double[]{0.3, 0.6, 0.1};
        for (int i = 0; i < n; ++i) {
            double rs = RS[i];
            double rimp = RIMP[i];
            for (int j = 0; j < UH.length && i + j < n; ++j)
                result[i + j] += (rs + rimp) * UH[j] * Area / (3.6 * dt);
        }

        return result;
    }

    /**
     * 计算径流总量序列
     *
     * @param QRS  地表径流序列
     * @param QRSS 壤中流序列
     * @param QRG  地下径流序列
     * @return 径流总量序列
     */
    private static double[] ComputeQ(double[] QRS, double[] QRSS, double[] QRG) {
        int n = QRS.length;
        double[] Q = new double[n];
        for (int i = 0; i < n; ++i)
            Q[i] = QRS[i] + QRSS[i] + QRG[i];

        return Q;
    }
}

/**
 * 土壤分层参数
 */
class LayeredSoilParam {
    /**
     * @param u 上层参数值
     * @param l 下层参数值
     * @param d 深层参数值
     */
    public LayeredSoilParam(double u, double l, double d) {
        Set(u, l, d);
    }

    /**
     * @param u 上层参数值
     * @param l 下层参数值
     * @param d 深层参数值
     */
    public void Set(double u, double l, double d) {
        U = u;
        L = l;
        D = d;
    }

    /**
     * @return 分层参数值总和
     */
    public double getSum() {
        return U + L + D;
    }

    /**
     * 上层参数值
     */
    public double U;

    /**
     * 下层参数值
     */
    public double L;

    /**
     * 深层参数值
     */
    public double D;
}



