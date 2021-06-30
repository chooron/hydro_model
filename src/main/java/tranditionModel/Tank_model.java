package tranditionModel;

import optimization.OptimizerByGA;
import tranditionModel.parent.parentModel;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * @author jing xin
 * @cite https://blog.csdn.net/qq_37948866/article/details/103145972
 * 水箱模型
 */
public class Tank_model implements parentModel {
    ArrayList<state> states = new ArrayList<state>(); // 时变状态变量
    int n1; // 单孔水箱个数
    int n2; // 双孔水箱个数
    double[] tankHs; // 各水箱初始水深
    double[] d_factors; // 各水箱下渗参数
    double[] r1_factors;// 双孔水箱出流1参数
    double[] r2_factors; // 双孔水箱出流2参数
    double[] r_factors; // 单孔水箱出流参数
    double[] thres1; // 双孔水箱出流1孔高
    double[] thres2; // 双孔水箱出流2孔高
    double[] thres; // 单孔水箱出流孔高

    ArrayList<single_tank> single_tanks; // 未知个数的单孔水箱
    ArrayList<double_tank> double_tanks; // 未知个数的双孔水箱

    /**
     * 输入双孔和单孔水箱个数来构造对象
     *
     * @param i 单孔个数
     * @param j 双孔个数
     */
    public Tank_model(int i, int j) {
        this.n1 = i;
        this.n2 = j;
        tankHs = new double[i + j]; // 各水箱初始水深
        d_factors = new double[i + j]; // 各水箱下渗参数
        r1_factors = new double[j]; // 双孔水箱出流1参数
        r2_factors = new double[j]; // 双孔水箱出流2参数
        r_factors = new double[i]; // 单孔水箱出流参数
        thres1 = new double[j]; // 双孔水箱出流1孔高
        thres2 = new double[j]; // 双孔水箱出流2孔高
        thres = new double[i]; // 单孔水箱出流孔高

        single_tanks = new ArrayList<single_tank>(i); // 未知个数的单孔水箱
        double_tanks = new ArrayList<double_tank>(j); // 未知个数的双孔水箱
    }

    public ArrayList<state> getStates() {
        return states;
    }

    public double[] getTankHs() {
        return tankHs;
    }

    public void setTankHs(double[] tankHs) {
        this.tankHs = tankHs;
    }

    public double[] getD_factors() {
        return d_factors;
    }

    public void setD_factors(double[] d_factors) {
        this.d_factors = d_factors;
    }

    public double[] getR1_factors() {
        return r1_factors;
    }

    public void setR1_factors(double[] r1_factors) {
        this.r1_factors = r1_factors;
    }

    public double[] getR2_factors() {
        return r2_factors;
    }

    public void setR2_factors(double[] r2_factors) {
        this.r2_factors = r2_factors;
    }

    public double[] getR_factors() {
        return r_factors;
    }

    public void setR_factors(double[] r_factors) {
        this.r_factors = r_factors;
    }

    public double[] getThres1() {
        return thres1;
    }

    public void setThres1(double[] thres1) {
        this.thres1 = thres1;
    }

    public double[] getThres2() {
        return thres2;
    }

    public void setThres2(double[] thres2) {
        this.thres2 = thres2;
    }

    public double[] getThres() {
        return thres;
    }

    public void setThres(double[] thres) {
        this.thres = thres;
    }

    public ArrayList<single_tank> getSingle_tanks() {
        single_tank[] single_tanks = new single_tank[tankHs.length - thres1.length];
//        ArrayList<single_tank> single_tanks = new ArrayList<single_tank>(tankHs.length);
        for (int i = thres1.length; i < tankHs.length; i++) {
            single_tank single_tank = new single_tank(tankHs[i], d_factors[i], r_factors[i - 1], thres[i - 1]);
            single_tanks[i - thres1.length] = single_tank;
        }
        return new ArrayList<single_tank>(Arrays.asList(single_tanks));
    }

    public ArrayList<double_tank> getDouble_tanks() {
        double_tank[] double_tanks = new double_tank[thres1.length];
        for (int i = 0; i < thres1.length; i++) {
            double_tank double_tank = new double_tank(tankHs[i], d_factors[i], r1_factors[i], r2_factors[i], thres1[i], thres2[i]);
            double_tanks[i] = double_tank;
        }
        return new ArrayList<double_tank>(Arrays.asList(double_tanks));
    }

    /**
     * 获得最优参数
     *
     * @param model          Tank模型
     * @param populationSize 种群大小
     * @param evolveTime     迭代次数的
     * @param target         训练目标
     * @param initPara       初始参数
     * @param threshold      误差允许范围
     * @return 最优参数
     */
    public double[] getBestPara(parentModel model, int populationSize, int evolveTime, double[] target, double[] initPara, double threshold) {
        OptimizerByGA ga = new OptimizerByGA(model, populationSize, evolveTime, target, initPara, threshold);
        return ga.GAmain();
    }

    @Override
    public double[] paraOptimize(double[] para) {
        return this.tank_model_Main(para);
    }

    /**
     * 各时段的状态变量
     */
    public static class state {
        public String time;
        public double P;
        public double R;
    }

    /**
     * 单孔水箱
     */
    static class single_tank {
        double z_end;//水箱储存量
        double d_factor;//单侧孔水箱的下渗系数
        double r_factor;//单侧孔水箱的流出系数
        double thres;//单侧孔水箱的孔高

        /**
         * 构造函数
         */
        public single_tank(double z_end, double d_factor, double r_factor, double thres) {
            this.z_end = z_end;
            this.d_factor = d_factor;
            this.r_factor = r_factor;
            this.thres = thres;
        }


        /**
         * 单孔水箱计算
         *
         * @return 下时段水深, 下渗, 出流
         */
        public double[] single_tank_count(double P) {
            double z_start = this.z_end + P;
            double out_down = z_start * this.d_factor; // 下渗量输出
            double out_right = (z_start <= thres) ? 0 : (z_start - thres) * r_factor;// 出流量输出
            double z_end = z_start - out_down - out_right;
            return new double[]{z_end, out_down, out_right};
        }
    }

    /**
     * 双侧孔水箱
     */
    static class double_tank {
        double z_end;//水箱储存量
        double d_factor;//双侧孔水箱的下渗系数
        double r1_factor;//双侧孔水箱的流出系数
        double r2_factor;//双侧孔水箱的流出系数
        double thres_1;//双侧孔水箱的孔高
        double thres_2;//双侧孔水箱的孔高

        public double_tank(double z_end, double d_factor, double r1_factor, double r2_factor, double thres_1, double thres_2) {
            this.z_end = z_end;
            this.d_factor = d_factor;
            this.r1_factor = r1_factor;
            this.r2_factor = r2_factor;
            this.thres_1 = thres_1;
            this.thres_2 = thres_2;
        }

        /**
         * 双孔水箱计算
         *
         * @return 下时段水深, 下渗, 出流1, 出流2
         */
        public double[] double_tank_count(double P) {
            double z_start = this.z_end + P;
            double out_down = z_start * d_factor;
            double out_right_1 = (z_start <= this.thres_1) ? 0 : (z_start - this.thres_1) * this.r1_factor;
            double out_right_2 = (z_start <= this.thres_2) ? 0 : (z_start - this.thres_2) * this.r2_factor;
            double z_end = z_start - out_down - out_right_1 - out_right_2;
            return new double[]{z_end, out_down, out_right_1 + out_right_2};
        }
    }

    /**
     * 计算主程序
     *
     * @param para = 水箱初始高度(i+j)+水箱下渗能力(i+j)+双孔高度1(j)+双孔高度2(j)+双孔出流1(j)+双孔出流2(j)+单孔高度(i)+单孔出流(i)
     * @return 径流过程
     */
    public double[] tank_model_Main(double[] para) {
        para[2 * (n1 + n2)-1] = 0;
        this.setTankHs(Arrays.copyOfRange(para, 0, n1 + n2));
        this.setD_factors(Arrays.copyOfRange(para, n1 + n2, 2 * (n1 + n2)));
        this.setThres1(Arrays.copyOfRange(para, 2 * (n1 + n2), 2 * (n1 + n2) + n2));
        this.setThres2(Arrays.copyOfRange(para, 2 * (n1 + n2) + n2, 2 * (n1 + n2) + 2 * n2));
        this.setR1_factors(Arrays.copyOfRange(para, 2 * (n1 + n2) + n2 * 2, 2 * (n1 + n2) + 3 * n2));
        this.setR2_factors(Arrays.copyOfRange(para, 2 * (n1 + n2) + n2 * 3, 2 * (n1 + n2) + 4 * n2));
        this.setThres(Arrays.copyOfRange(para, 2 * (n1 + n2) + n2 * 4, 2 * (n1 + n2) + 4 * n2 + n1));
        this.setR_factors(Arrays.copyOfRange(para, 2 * (n1 + n2) + n2 * 4 + n1, 2 * (n1 + n2) + 4 * n2 + n1 * 2));
        ArrayList<state> states = this.getStates();
        ArrayList<single_tank> STs = this.getSingle_tanks();
        ArrayList<double_tank> DTs = this.getDouble_tanks();
        double[][][] result = new double[states.size()][3][STs.size() + DTs.size()];
        double[] Rs = new double[states.size()];
        for (int i = 0; i < states.size(); i++) {
            for (int j = 0; j < DTs.size(); j++) {
                double p = (j == 0) ? states.get(i).P : result[i][1][j - 1];
                result[i][0][j] = DTs.get(j).double_tank_count(p)[0];
                result[i][1][j] = DTs.get(j).double_tank_count(p)[1];
                result[i][2][j] = DTs.get(j).double_tank_count(p)[2];
                DTs.get(j).z_end = result[i][0][j];
            }
            for (int j = DTs.size(); j < STs.size() + DTs.size(); j++) {
                double p = result[i][1][j - 1];
                result[i][0][j] = STs.get(j - DTs.size()).single_tank_count(p)[0];
                result[i][1][j] = STs.get(j - DTs.size()).single_tank_count(p)[1];
                result[i][2][j] = STs.get(j - DTs.size()).single_tank_count(p)[2];
                STs.get(j - DTs.size()).z_end = result[i][0][j];
            }
            for (int j = 0; j < STs.size() + DTs.size(); j++) {
                Rs[i] += result[i][2][j];
            }
        }
        return Rs;
    }
}
