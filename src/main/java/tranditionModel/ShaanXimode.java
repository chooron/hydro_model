package tranditionModel;


import optimization.OptimizerByGA;
import tranditionModel.parent.parentModel;

import java.util.ArrayList;


/**
 * 陕北模型
 * @author jing xin
 */
public class ShaanXimode implements parentModel {
    //---------常规变量-----------
    private int type; // 下渗曲线选取
    private ArrayList<Double> ps; // 降雨过程
    private ArrayList<Double> es; //蒸发过程
    private double dt;  // 时间步长
    private double bx; // 下渗能力分布曲线指数
    private double area; //流域面积
    private double winit; //初始土壤含水
    double[] ws; //土壤含水量
    double[] rs; //径流量
    private double[] UH;

    //------霍顿下渗曲线参数-------
    private double f0; // 流域平均最大下渗能力
    private double fc; // 流域平均稳定下渗能力
    private double k; // 下渗能力衰弱系数
    private double err_limit; // 土壤含水计算的允许误差

    //------菲利普下渗曲线参数------
    private double A; // 菲利普下渗曲线参数
    private double B;


    //参数优化
    double[] target;

    public double[] getTarget() {
        return target;
    }

    public void setTarget(double[] target) {
        this.target = target;
    }

    public void setBx(double bx) {
        this.bx = bx;
    }

    public double[] getRs() {
        return rs;
    }

    // 土壤含水量变化量计算
    public ShaanXimode soilWCount() {
        ArrayList<Double> ps = this.ps;
        ArrayList<Double> es = this.es;
        double[] rs = new double[ps.size()];
        double[] ws = new double[ps.size()];
        double ft = 0.0;
        ws[0] = (this.winit);
        for (int i = 0; i < ps.size(); i++) {
            double rs2;
            double ws2;
            if (ps.get(i) - es.get(i) < 0) {
                rs2 = (0.0);
                ws2 = ws[i] - (ps.get(i) - es.get(i));
            } else {
                if (this.type == 0) {
                    ft = HortonInfil(this, ws[i]);
                } else if (this.type == 1) {
                    ft = PhillipsInfli(this, ws[i]);
                }
                // 计算流域该时段的最大点下渗能力
                double fmm = ft * (1 + this.bx);
                if (ps.get(i) - es.get(i) > fmm) {
                    // 全流域产流
                    rs2 = ps.get(i) - es.get(i) - fmm;
                    ws2 = ws[i] + ft;
                } else {
                    rs2 = ps.get(i) - es.get(i) - ft + ft * Math.pow(1 - (ps.get(i) - es.get(i)) / fmm, this.bx + 1);
                    ws2 = ws[i] + ps.get(i) - es.get(i) - rs2;
                }
            }
            rs[i] = rs2;
            if(i+1<ps.size()) ws[i + 1] = ws2;
        }
        this.setResult(ws, rs);
        return this;
    }

    public double[] UnitHydrograph(double[] RS, double[] RIMP) {
        int n = RS.length;
        double[] result = new double[n];
        for (int i = 0; i < n; ++i) {
            double rs = RS[i];
            double rimp = RIMP[i];
            for (int j = 0; j < UH.length && i + j < n; ++j)
                // TODO: 2020/12/22 这里与原计算存在倍比关系，目前找不到问题（5）
                result[i + j] += (rs + rimp) * UH[j] * area / (3.6 * dt)*5;
        }
        return result;
    }

    private void setResult(double[] ws, double[] rs) {
        this.ws = ws;
        this.rs = rs;
    }

    public ShaanXimode modelInit(int type, ArrayList<Double> ps, ArrayList<Double> es, double dt, double area, double winit, double[] UH) {
        this.type = type;
        this.ps = ps;
        this.es = es;
        this.dt = dt;
        this.area = area;
        this.winit = winit;
        this.UH = UH;
        return this;
    }

    public ShaanXimode H_init(double f0, double fc, double k, double err_limit) {
        this.f0 = f0;
        this.fc = fc;
        this.k = k;
        this.err_limit = err_limit;
        return this;
    }

    public ShaanXimode F_init(double A, double B) {
        this.A = A;
        this.B = B;
        return this;
    }

    /**
     * 菲利普斯下渗曲线计算
     *
     * @param shaanXimode 陕北模型
     * @param w1          时段初土壤含水
     * @return 下渗
     */
    private static double PhillipsInfli(ShaanXimode shaanXimode, double w1) {
        return shaanXimode.B * shaanXimode.B * (1 + Math.sqrt(1 + shaanXimode.A * w1
                / (shaanXimode.B * shaanXimode.B))) / w1 + shaanXimode.A;
    }

    /**
     * 霍顿下渗曲线计算
     *
     * @param shaanXimode 陕北模型
     * @param w1          时段初土壤含水
     * @return 下渗
     */
    private static double HortonInfil(ShaanXimode shaanXimode, double w1) {
        // 由经验初设模型参数
        double T = w1 / shaanXimode.f0;

        double w1Bycount = shaanXimode.fc * T + 1 / shaanXimode.k * (shaanXimode.f0 -
                shaanXimode.fc) * (1 - Math.exp(-shaanXimode.k * T));
        double fBycount = 0;
        while (Math.abs(w1 - w1Bycount) > shaanXimode.err_limit) {
            fBycount = shaanXimode.fc + (shaanXimode.f0 - shaanXimode.fc) * Math.exp(-shaanXimode.k * T);
            T = T + (w1 - w1Bycount) / fBycount;
            w1Bycount = shaanXimode.fc * T + 1 / shaanXimode.k * (shaanXimode.f0 -
                    shaanXimode.fc) * (1 - Math.exp(-shaanXimode.k * T));
        }
        return fBycount;
    }

    /**
     * 返回优化后的参数
     *
     * @param para 初始参数
     * @return 优化后的参数
     */
    @Override
    public double[] paraOptimize(double[] para) {
//        this.bx = para[0];
        this.H_init(para[0], para[1], para[2], 0.05);
        this.soilWCount();
        return this.UnitHydrograph(this.getRs(), new double[this.getRs().length]);
    }

    public double[] getBestPara(parentModel model, int populationSize, int evolveTime, double[] target, double[] initPara, double threshold) {
        OptimizerByGA ga = new OptimizerByGA(model, populationSize, evolveTime, target, initPara, threshold);
        return ga.GAmain();
    }

}
