package tranditionModel;

import optimization.OptimizerByGA;
import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import tranditionModel.parent.parentModel;
import tranditionModel.privateUtil.BeanRefUtil;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * @author jing xin
 * 该版为Sacramento模型的解耦版本，已基本通过数据测试，但相比matlab程序有误差有较小的
 */
public class SACmodel implements parentModel {

    public static class SACstate {
        double P;//降雨
        double EP;//蒸发

        // 含水量
        double AUZTW;
        double ALZTW;
        double UZTWC;
        double UZFWC;
        double LZTWC;
        double LZFSC;
        double LZFPC;

        // 蒸发
        // 可变不透水面积的蒸散发量
        double AE1;
        double AE3;
        // 各层蒸发量(上土层张力水、上土层自由水、下土层张力水、水面蒸发量)
        double E1;
        double E2;
        double E3;
        double E4;

        // 汇流

        double PAV;//有效降雨计算
        double ADSUR;// 地面径流计算中间值
        double ARS;//不透水层地面径流计算
        double ROIMP;//不透水产流计算
        double RS;//地面径流
        double RI;//壤中流
        double RGS;// 下层快速水
        double RGP;// 下层慢速水

        // 汇流
        double QI;// 壤中出流计算
        double QS;// 快速出流计算
        double QP;// 慢速出流计算
        double QT;// 河网总入流
        double Q;// 出口流量计算
        private double LF; // 中间值(这里是用于传递方法所需参数，意义不用太深究)

        public double getP() {
            return P;
        }

        public double getEP() {
            return EP;
        }

        public double getAUZTW() {
            return AUZTW;
        }

        public void setAUZTW(double AUZTW) {
            this.AUZTW = AUZTW;
        }

        public double getALZTW() {
            return ALZTW;
        }

        public void setALZTW(double ALZTW) {
            this.ALZTW = ALZTW;
        }

        public double getUZTWC() {
            return UZTWC;
        }

        public void setUZTWC(double UZTWC) {
            this.UZTWC = UZTWC;
        }

        public double getUZFWC() {
            return UZFWC;
        }

        public void setUZFWC(double UZFWC) {
            this.UZFWC = UZFWC;
        }

        public double getLZTWC() {
            return LZTWC;
        }

        public void setLZTWC(double LZTWC) {
            this.LZTWC = LZTWC;
        }

        public double getAE1() {
            return AE1;
        }

        public double getAE3() {
            return AE3;
        }

        public double getPAV() {
            return PAV;
        }

        public double getADSUR() {
            return ADSUR;
        }

        public void setCIState(double[] ciState) {
            this.AUZTW = ciState[0];
            this.AE1 = ciState[1];
            this.PAV = ciState[2];
            this.ALZTW = ciState[3];
            this.AE3 = ciState[4];
            this.ADSUR = ciState[5];
            this.ARS = ciState[6];
            this.ROIMP = ciState[7];
        }

        public void setEs(double[] es) {
            this.E1 = es[0];
            this.E2 = es[1];
            this.E3 = es[2];
            this.E4 = es[3];
        }

        public double[] getEs() {
            return new double[]{E1, E2, E3, E4};
        }

        public double[] getSWstorage() {
            return new double[]{UZTWC, UZFWC, LZTWC, LZFSC, LZFPC};
        }

        public double getROIMP() {
            return ROIMP;
        }

        public void setRS(double rs) {
            this.RS = rs;
        }

        public double getRS() {
            return RS;
        }

        public void setRI(double ri) {
            this.RI = ri;
        }

        public double getRI() {
            return RI;
        }

        public double getLZFSC() {
            return LZFSC;
        }

        public void setLZFSC(double lzfsc) {
            this.LZFSC = lzfsc;
        }

        public double getLZFPC() {
            return LZFPC;
        }

        public void setLZFPC(double lzfpc) {
            this.LZFPC = lzfpc;
        }

        public double getQI() {
            return QI;
        }

        public void setQI(double qi) {
            this.QI = qi;
        }

        public double getQS() {
            return QS;
        }

        public void setQS(double qs) {
            this.QS = qs;
        }

        public double getRGS() {
            return RGS;
        }

        public void setRGS(double rgs) {
            this.RGS = rgs;
        }

        public double getQP() {
            return QP;
        }

        public void setQP(double qp) {
            this.QP = qp;
        }

        public double getRGP() {
            return RGP;
        }

        public void setRGP(double rgp) {
            this.RGP = rgp;
        }

        public double getARS() {
            return ARS;
        }

        public void setARS(double ars) {
            this.ARS = ars;
        }

        public double getQ() {
            return Q;
        }

        public void setQ(double q) {
            this.Q = q;
        }

        public double getQT() {
            return QT;
        }

        public void setQT(double qt) {
            this.QT = qt;
        }

        public void setQs(double[] qs) {
            this.QI = qs[0];
            this.QS = qs[1];
            this.QP = qs[2];
            this.QT = qs[3];
            this.Q = qs[4];
        }

        public void setP(double p) {
            this.P = p;
        }

        public void setEP(double E0) {
            double KC = 0.9;
            this.EP = E0 * KC;
        }

        public void setLF(double lf) {
            this.LF = lf;
        }

        public double getLF() {
            return LF;
        }
    }

    // 各阶段状态参数
    ArrayList<SACstate> states;

    // 土壤容量参数
    double AUZTWM;//用于不透水层计算中上含水层含水容量(此处书中暂未找到对应公式)
    double ALZTWM;//用于不透水层计算中下含水层含水容量(此处书中暂未找到对应公式)
    double UZTWM;// 上层张力水容量
    double UZFWM;// 上层自由水容量
    double LZTWM;// 下层张力水容量
    double LZFSM;// 下层浅层自由水容量
    double LZFPM;// 下层深层自由水容量

    // 模型初始值
    double AUZTWC0;// 用于不透水层计算中上含水层初始含水量
    double ALZTWC0;// 用于不透水层计算中下含水层初始含水量
    double UZTWC0; // 上含水层张力水初始量
    double LZTWC0; // 下含水层张力水初始量
    double UZFWC0; // 上含水层自由水初始量
    double LZFSC0;// 下层浅层自由水初始量
    double LZFPC0;// 下层深层自由水初始量

    double RIVA;  // 河网蒸发涉及参数

    // 产流模型
    double UZK; // 壤中流日出流系数
    private double PCTIM;// 不透水面占全流域百分比
    private double RSERV;// 下层自由水未蒸发部分所占比例
    private double LZSK;// 下层浅层自由水出流系数
    private double LZPK;// 下层深层自由水出流系数
    private double ZPERC;// 下渗系数,下层最干旱时的最大下渗率
    private double PFREE;// 下渗水量直接补给下层自由水比例
    private double f;// 流域面积
    private double ADIMP;// 可变不透水面积
    // 汇流参数(具体我也不清楚)
    private double CI;
    private double CGS;
    private double CGP;
    private double CR;
    private double PAREA;

    public double getUZK() {
        return UZK;
    }

    public void setUZK(double uzk) {
        this.UZK = uzk;
    }

    public double getPCTIM() {
        return PCTIM;
    }

    public void setPCTIM(double pctim) {
        this.PCTIM = pctim;
    }

    public double getRSERV() {
        return RSERV;
    }

    public void setRSERV(double rserv) {
        this.RSERV = rserv;
    }

    public double getLZSK() {
        return LZSK;
    }

    public void setLZSK(double lzsk) {
        this.LZSK = lzsk;
    }

    public double getLZPK() {
        return LZPK;
    }

    public void setLZPK(double lzpk) {
        this.LZPK = lzpk;
    }

    public double getZPERC() {
        return ZPERC;
    }

    public void setZPERC(double zperc) {
        this.ZPERC = zperc;
    }

    public double getPFREE() {
        return PFREE;
    }

    public void setPFREE(double pfree) {
        this.PFREE = pfree;
    }

    public double getF() {
        return f;
    }

    public void setF(double f) {
        this.f = f;
    }

    public double getCI() {
        return CI;
    }

    public void setCI(double ci) {
        this.CI = ci;
    }

    public double getCGS() {
        return CGS;
    }

    public void setCGS(double cgs) {
        this.CGS = cgs;
    }

    public double getCGP() {
        return CGP;
    }

    public void setCGP(int cgp) {
        this.CGP = cgp;
    }

    public double getADIMP() {
        return ADIMP;
    }

    public void setADIMP(double adimp) {
        this.ADIMP = adimp;
    }

    public double getPAREA() {
        return PAREA;
    }

    public void setPAREA(double parea) {
        this.PAREA = parea;
    }

    public double getCR() {
        return CR;
    }

    public void setCR(double cr) {
        this.CR = cr;
    }

    public double getUZTWM() {
        return UZTWM;
    }

    public double getLZTWM() {
        return LZTWM;
    }

    public double[] getWM() {
        return new double[]{UZTWM, UZFWM, LZTWM, LZFSM, LZFPM};
    }

    public double getUZTWC0() {
        return UZTWC0;
    }

    public void setUZTWC0(double UZTWC0) {
        this.UZTWC0 = UZTWC0;
    }

    public double getLZTWC0() {
        return LZTWC0;
    }

    public void setLZTWC0(double LZTWC0) {
        this.LZTWC0 = LZTWC0;
    }

    public double getUZFWC0() {
        return UZFWC0;
    }

    public void setUZFWC0(double UZFWC0) {
        this.UZFWC0 = UZFWC0;
    }

    public double getLZFSC0() {
        return LZFSC0;
    }

    public void setLZFSC0(double LZFSC0) {
        this.LZFSC0 = LZFSC0;
    }

    public double getLZFPC0() {
        return LZFPC0;
    }

    public void setLZFPC0(double LZFPC0) {
        this.LZFPC0 = LZFPC0;
    }

    public double getRIVA() {
        return RIVA;
    }

    public ArrayList<SACstate> getStates() {
        return states;
    }

    public void setStates(ArrayList<SACstate> states) {
        this.states = states;
    }

    public void setUZTWM(double UZTWM) {
        this.UZTWM = UZTWM;
    }

    public void setUZFWM(double UZFWM) {
        this.UZFWM = UZFWM;
    }

    public void setLZTWM(double LZTWM) {
        this.LZTWM = LZTWM;
    }

    public void setLZFSM(double LZFSM) {
        this.LZFSM = LZFSM;
    }

    public void setLZFPM(double LZFPM) {
        this.LZFPM = LZFPM;
    }

    public void setRIVA(double RIVA) {
        this.RIVA = RIVA;
    }

    public void setCGP(double CGP) {
        this.CGP = CGP;
    }

    public double getAUZTWM() {
        return AUZTWM;
    }

    public void setAUZTWM(double AUZTWM) {
        this.AUZTWM = AUZTWM;
    }

    public double getALZTWM() {
        return ALZTWM;
    }

    public void setALZTWM(double ALZTWM) {
        this.ALZTWM = ALZTWM;
    }

    public double getUZFWM() {
        return UZFWM;
    }

    public double getLZFSM() {
        return LZFSM;
    }

    public double getLZFPM() {
        return LZFPM;
    }

    public double getAUZTWC0() {
        return AUZTWC0;
    }

    public void setAUZTWC0(double AUZTWC0) {
        this.AUZTWC0 = AUZTWC0;
    }

    public double getALZTWC0() {
        return ALZTWC0;
    }

    public void setALZTWC0(double ALZTWC0) {
        this.ALZTWC0 = ALZTWC0;
    }

    /**
     * 获取优化参数计算得到的实际结果
     *
     * @param para 待优化参数
     * @return 实际结果(流量)
     */
    @Override
    public double[] paraOptimize(double[] para) {
        String[] paraName = new String[]{"AUZTWM", "ALZTWM", "UZTWM", "UZFWM", "LZTWM", "LZFSM", "LZFPM", "RIVA",
                "UZK", "PCTIM", "RESERV", "LZSK", "LZPK", "ZPERC", "PFREE", "ADIMP", "CI", "CGS", "CGP", "CR", "PAREA"};
        HashMap<String, String> map = new HashMap<String, String>();
        for (int i = 0; i < para.length; i++) {
            map.put(paraName[i], String.valueOf(para[i]));
        }
        BeanRefUtil.setFieldValue(this, map);
        this.SACmain();
        double[] Q = new double[this.getStates().size()];
        for (int i = 0; i < this.getStates().size(); i++) {
            Q[i] = this.getStates().get(i).getQ();
        }
        return Q;
    }

    /**
     * 获得最优参数
     *
     * @param model          SAC模型
     * @param populationSize 种群大小
     * @param evolveTime     迭代次数的
     * @param target         训练目标
     * @param initPara       初始参数
     * @param threshold      误差允许范围
     * @return 最优参数
     */
    public static double[] getBestPara(parentModel model, int populationSize, int evolveTime, double[] target, double[] initPara, double threshold) {
        OptimizerByGA ga = new OptimizerByGA(model, populationSize, evolveTime, target, initPara, threshold);
        return ga.GAmain();
    }

    /**
     * 变动不透水层计算
     */
    public void varInsermeable() {
        ArrayList<SACstate> states = this.getStates();
        double AUZTW, AE1, PAV, ALZTW, AE3, ADSUR, ARS, ROIMP;
        double AUZTWM = this.getAUZTWM();
        double ALZTWM = this.getALZTWM();
        for (int i = 0; i < states.size(); i++) {
            SACstate state = states.get(i);
            SACstate initState = stateInit();
            SACstate state0 = (i == 0) ? initState : states.get(i - 1);
            AUZTW = (i == 0) ? state0.getAUZTW() : Math.min(AUZTWM, state0.getAUZTW() - state0.getAE1() + state0.getP());
            ALZTW = (i == 0) ? state0.getALZTW() : Math.min(ALZTWM, state0.getPAV() - state0.getADSUR() + state0.getALZTW() - state0.getAE3());
            AE1 = Math.min(AUZTWM, state.getEP() * AUZTW / AUZTWM);
            AE3 = (state.getEP() - AE1) * ALZTW / (AUZTWM + ALZTWM);
            PAV = Math.max(0, state.getP() - AUZTWM + AUZTW - AE1);
            ADSUR = PAV * (ALZTW - AE3) / 90;
            ARS = Math.max(0, PAV - ADSUR + ALZTW - AE3 - 90);
            ROIMP = this.getPCTIM() * state.getP();
            state.setCIState(new double[]{AUZTW, AE1, PAV, ALZTW, AE3, ADSUR, ARS, ROIMP});
        }
        this.setStates(states);
    }


    /**
     * 上含水层计算(包括蓄水、蒸发、径流)
     */
    public SACstate upPermeable(SACstate state, SACstate state0) {
        double[] wm = this.getWM(); // 各层土壤最大蓄水量
        double UZTWM = wm[0];
        double UZFWM = wm[1];
        // this.getUZTWM()原代码中是UT的最大值
        double UZTWC0 = state0.getUZTWC();
        double UZFWC0 = state0.getUZFWC();
        double UT = Math.min(20, state0.getUZTWC() - state0.getEs()[0] + state0.getP());
        // this.getUZTWM()原代码中是UF的最大值
        double UF = Math.min(30, state0.getP() + state0.getUZTWC() + state0.getUZFWC()
                - state0.getEs()[0] - state0.getEs()[1] - UT);
        double LF = (1 - this.getUZK()) * UF;
        state.setUZTWC((UF / 20 < LF / 30) ? 20 * (UF + LF) / (UZFWM + UZTWM) : UT);
        state.setUZFWC((UF / 20 < LF / 30) ? 30 * (UT + LF) : LF);
        double[] Es = state.getEs();  // 存储四个蒸发量
        Es[0] = Math.min(state.getUZTWC(), state.getEP() * state.getUZTWC() / this.getUZTWM());
        Es[1] = Math.min(state.getUZFWC(), state.getEP() - Es[0]);
        state.setRS(Math.max(0, state.getP() + UZTWC0 + UZFWC0 - Es[0] - Es[1] - wm[0] - wm[1]));
        state.setRI(UF * this.getUZK());
        state.setEs(Es);
        state.setLF(LF);
        return state;
    }

    /**
     * 下透水层计算
     */
    public SACstate lowPermeable(SACstate state, SACstate state0) {
        double[] wm = this.getWM(); // 各层土壤最大蓄水量
        double LZFPM = wm[4], LZFSM = wm[3], LZTWM = wm[2], UZFWM = wm[1], UZTWM = wm[0];
        double SAVED = this.getRSERV() * (LZFSM + LZFPM);
        double PBASE = LZFSM * this.getLZSK() + LZFPM * this.getLZPK();  // 稳定下渗率
        double LZTWC = state0.getLZTWC(), LZFSC = state0.getLZFSC(), LZFPC = state0.getLZFPC();
        double[] Es = state.getEs();
        Es[2] = (state.getEP() - Es[0] - Es[1]) * (state.getLZTWC() / (this.getUZTWM() + this.getLZTWM()));
        Es[3] = this.getRIVA() * state.getEP();
        state.setEs(Es);
        double LF = state.getLF();
        double LT = LZTWC - Es[2];  // 下层张力蓄量
        double DEFR = 1 - (LZFSC + LZFPC + LT) / (LZFPM + LZFSM + LZTWM); // 缺水率
        double PERC = PBASE * (1 + this.getZPERC() * DEFR * DEFR) * LF / UZFWM; // 下渗率
        // 下渗计算
        double RATE = Math.min(PERC, LZFSM + LZFPM + LZTWM - LZFSC - LZFPC - LT);
        //--- 下层自由水计算 ----//
        double FX = Math.min(LZFSM + LZFPM - LZFSC - LZFPC, Math.max(RATE - LZTWM + LT, RATE * this.getPFREE()));
        // 分配慢速自由水
        double PERCP = Math.min(LZFPM - LZFPC, Math.max(FX - LZFSM + LZFSC,
                (LZFPM / (LZFSM + LZFPM)) * (2 * (1 - (LZFPC / LZFPM)) / ((1 - LZFPC / LZFPM) + 1 - LZFSC / LZFSM)) * FX));
        // 下层快速水计算
        double LS = FX - PERCP + LZFSC;
        // 下层慢速水计算
        double LP = PERCP + LZFPC;
        // 快速地下水计算
        double RGS = LS * this.getLZSK();
        // 慢速地下水计算
        double RGP = LP * this.getLZPK();
        // 下层快速水更新
        LS -= RGS;
        // 下层慢速水更新
        LP -= RGP;
        //--- 下层张力水计算 ---//
        double PERCT = RATE - FX;
        // 更新下层张力水蓄量
        LT += PERCT;
        double RATIO = (LS + LP - SAVED + LT) / (LZFSM + LZFPM - SAVED + LZTWM);
        state.setRGP(RGP);
        state.setRGS(RGS);
        state.setLZTWC(Math.max(LZTWM * RATIO, LT)); //下层快自由水计算
        state.setLZFSC(LT / LZTWM < RATIO ? LS - (LZFSC - LT) : LS); //下层快自由水计算
        state.setLZFPC(LT / LZTWM < RATIO ? LP - Math.max(0, LZFSC - LT) : LP);//下层慢自由水计算
        return state;
    }

    /**
     * 汇流计算
     */
    public SACstate confluence(SACstate state, SACstate state0) {
        double U0 = (1 - 0.01 - 0.01) * this.getF() / (3.6 * 6); // U0为单位转换系数
        double QI = this.getCI() * state0.getQI() + (1.0 - this.getCI()) * state0.getRI() * U0;
        double QS = this.getCGS() * state0.getQS() + 0.2 * state0.getRGS() * U0;
        double QP = this.getCGP() * state0.getQP() + (1 - this.getCGP()) * state0.getRGP() * U0;
        double QT = QI + QS + QP + ((state0.getROIMP() + state.getADSUR()
                + state.getARS()) * this.getADIMP() + state.getRS() * this.getPAREA()) * this.getF() / (3.6 * 6);
        double Q = this.getCR() * state0.getQ() + (1 - this.getCR()) * 0.5 * (state0.getQT() + QT);
        state.setQs(new double[]{QI, QS, QP, QT, Q});
//        System.out.println(Q);
        return state;
    }

    /**
     * 设置模型初始不需优化的参数
     *
     * @param values
     */
    public void setInitValue(double[] values) {
        this.AUZTWC0 = values[0];
        this.ALZTWC0 = values[1];
        this.UZTWC0 = values[2];
        this.UZFWC0 = values[3];
        this.LZTWC0 = values[4];
        this.LZFSC0 = values[5];
        this.LZFPC0 = values[6];
        this.f = values[7];
    }

    /**
     * 模型初始化计算
     */
    public SACstate stateInit() {
        SACstate state0 = this.getStates().get(0);
        state0.setAUZTW(this.getAUZTWC0());
        state0.setALZTW(this.getALZTWC0());
        state0.setUZTWC(this.getUZTWC0());
        state0.setUZFWC(this.getUZFWC0());
        state0.setLZTWC(this.getLZTWC0());
        state0.setLZFPC(this.getLZFPC0());
        state0.setLZFSC(this.getLZFSC0());
        double E3 = (state0.getEP() - 0.45 - 0) * (state0.getLZTWC() / (this.getUZTWM() + this.getLZTWM()));
        double E4 = this.getRIVA() * state0.getEP();
        state0.setEs(new double[]{0.45, 0, E3, E4});
        state0.setQs(new double[]{0, 0, 0, 0, 0});
        return state0;
    }

    public SACmodel SACmain() {
        // 变量创建
        ArrayList<SACstate> states = this.getStates();
        SACstate initState = this.stateInit();
        int n = states.size();
        this.varInsermeable();
        for (int i = 0; i < n; i++) {
            SACstate state = states.get(i);
            SACstate state0 = (i == 0) ? initState : states.get(i - 1);
            state = this.upPermeable(state, state0);
            state = this.lowPermeable(state, state0);
            state = this.confluence(state, state0);
            states.set(i, state);
        }
        this.setStates(states);
        return this;
    }

    public static SACmodel getExcelForSAC(String ExcelName) {
        SACmodel model = new SACmodel();
        HSSFWorkbook wb = null;
        ArrayList<SACstate> states = new ArrayList<SACstate>();
        try {
            FileInputStream fis = new FileInputStream(new File(ExcelName));
            wb = new HSSFWorkbook(fis);
            fis.close();
            ArrayList<Double> tempData = new ArrayList<Double>();
            HSSFSheet sheet = wb.getSheet("Sheet2");
            double lastRowNum = sheet.getLastRowNum();
            for (int i = 1; i < lastRowNum + 1; i++) {
                SACstate state = new SACstate();
                HSSFRow row = sheet.getRow(i);
                state.setP(row.getCell(0).getNumericCellValue());
                state.setEP(row.getCell(1).getNumericCellValue());
                states.add(state);
            }
            model.setStates(states);
        } catch (IOException e) {
            e.printStackTrace();
        }
        return model;
    }
}
