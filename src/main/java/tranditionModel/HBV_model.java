package tranditionModel;

import optimization.OptimizerByGA;
import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.dom4j.Document;
import org.dom4j.DocumentException;
import org.dom4j.Element;
import org.dom4j.io.SAXReader;
import tranditionModel.parent.parentModel;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;


/**
 * HBV水文模型
 *
 * @author jing xin
 * @Citation: AghaKouchak A., Habib E., 2010, Application of a Conceptual Hydrologic
 * Model in Teaching Hydrologic Processes, International Journal of Engineering Education, 26(4), 963-973.
 * <p>
 * AghaKouchak A., Nakhjiri N., and Habib E., 2012, An educational model for ensemble streamflow
 * simulation and uncertainty analysis, Hydrology and Earth System Sciences Discussions, 9, 7297-7315, doi:10.5194/hessd-9-7297-2012.
 */
public class HBV_model implements parentModel {
    /*待优化参数*/
    double DD;
    double FC;     // 田间持水量
    double Beta;
    double C;
    double L;
    double K0;
    double K1;
    double K2;
    double Kp;
    double PWP;

    /*恒定参数*/
    public double area; // 流域面积
    public double T_t; // 融雪温度
    public double snow0;  // 初始雪量
    public double soil0; // 初始土壤含水
    public double S10; // S_1初始值
    public double S20; // S_2初始值

    public double[] qo; // 实测流量

    public double[] getQo() {
        return qo;
    }

    public ArrayList<state> states; // 时间序列变量
    public double[] monthly;
    public double[] dpem;

    public static class state {
        public int month;
        public double prec;  // 降雨量
        public double snow;  // 雪量
        public double temp;  // 温度
        public double lwater;
        public double pe;
        public double ea; // 实际蒸发量
        public double soil;
        public double dq;
        public double S1;
        public double S2;
        public double q;
        public double R;
    }


    @Override
    public double[] paraOptimize(double[] para) {
        initPara(para);
        HBV_main();
        return getR();
    }

    public double[] getBestPara(int populationSize, int evolveTime, double[] target, double[] initPara, double threshold) {
        OptimizerByGA ga = new OptimizerByGA(this, populationSize, evolveTime, target, initPara, threshold);
        return ga.GAmain();
    }

    public double[] getR() {
        double[] result = new double[states.size()];
        for (int i = 0; i < states.size(); i++) {
            result[i] = states.get(i).R;
        }
        return result;
    }


    public void initPara(double[] para) {
        this.DD = para[0];
        this.FC = para[1];
        this.Beta = para[2];
        this.C = para[3];
        this.L = para[4];
        this.K0 = para[5];
        this.K1 = para[6];
        this.K2 = para[7];
        this.Kp = para[8];
        this.PWP = para[9];
    }

    /**
     * 模型计算主程序
     */
    public void HBV_main() {
        state state_first = states.get(0);
        state_first.snow = snow0;
        state_first.soil = soil0;
        state_first.S1 = S10;
        state_first.S2 = S20;
        for (int i = 1; i < states.size(); i++) {
            state state0 = states.get(i - 1);
            state state1 = states.get(i);
            if (state1.temp < T_t) {
                state1.snow = state0.snow + state1.prec;
                state1.lwater = 0;
            } else {
                state1.snow = Math.max(state0.snow - DD * (state1.temp - T_t), 0);
                state1.lwater = state1.prec + Math.min(state0.snow, DD * (state1.temp - T_t));
            }
            state1.pe = (1 + C * (state1.temp - monthly[state1.month])) * dpem[state1.month];
            state1.ea = (state0.soil > PWP) ? state1.pe : state1.pe * (state0.soil / PWP);
            state1.dq = state1.lwater * Math.pow(state0.soil / FC, Beta);
            state1.soil = state0.soil + state1.lwater - state1.dq - state1.ea;
            state1.S1 = state0.S1 + state1.dq - Math.max(0, state0.S1 - L) * K0 - state0.S1 * K1 - state0.S1 * Kp;
            state1.S2 = state0.S2 + state0.S1 * Kp - state0.S2 * K2;
            state1.q = Math.max(0, state0.S1 - L) * K0 + state1.S1 * K1 + state1.S2 * K2;
            state1.R = state1.q * area * 1000 / (24 * 3600); // 时段为一天
        }

    }

    public void initModel(String xmlPath) throws DocumentException, IOException {
        SAXReader reader = new SAXReader();
        Document document = reader.read(new File(xmlPath));
        Element rootElement = document.getRootElement();
        Element path = rootElement.element("path");
        String excelPath = path.element("excel").getText();
        String sheet1 = path.element("sheet1").getText();
        String sheet2 = path.element("sheet2").getText();
        FileInputStream fis = new FileInputStream(excelPath);
        HSSFWorkbook wb = new HSSFWorkbook(fis);
        fis.close();
        HSSFSheet dataSheet = wb.getSheet(sheet1);
        int lastRowNum = dataSheet.getLastRowNum();
        ArrayList<state> states = new ArrayList<>();
        double[] qo = new double[lastRowNum];
        for (int i = 0; i < lastRowNum; i++) {
            HSSFRow row = dataSheet.getRow(i+1);
            state state = new state();
            state.month = (int) row.getCell(1).getNumericCellValue();
            state.temp = row.getCell(2).getNumericCellValue();
            state.prec = row.getCell(3).getNumericCellValue();
            states.add(state);
            qo[i] = row.getCell(4).getNumericCellValue();
        }
        this.states = states;
        this.qo = qo;
        HSSFSheet avgSheet = wb.getSheet(sheet2);
        int lastRowNum2 = avgSheet.getLastRowNum();
        double[] monthly = new double[lastRowNum2];
        double[] dpem = new double[lastRowNum2];
        for (int i = 0; i < lastRowNum2; i++) {
            HSSFRow row = avgSheet.getRow(i + 1);
            monthly[i] = row.getCell(0).getNumericCellValue();
            dpem[i] = row.getCell(2).getNumericCellValue();
        }
        this.monthly = monthly;
        this.dpem = dpem;
        this.snow0 = Double.parseDouble(path.element("snow0").getText());
        this.soil0 = Double.parseDouble(path.element("soil0").getText());
        this.S10 = Double.parseDouble(path.element("S10").getText());
        this.S20 = Double.parseDouble(path.element("S20").getText());
        this.area = Double.parseDouble(path.element("area").getText());
        wb.close();
    }

}
