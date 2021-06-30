package tranditionModel;


import org.dom4j.DocumentException;
import org.junit.Test;
import tranditionModel.privateUtil.BeanRefUtil;
import util.Container;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * 遗传算法优化传统概念性水文模型的预测测试
 *
 * @author jing xin
 */
public class forecastTest {
    @Test
    public void XajmodelTest() {
        XajModel model = new XajModel();
        double[] target = new double[]{145.04165680890637, 460.7913823034037, 842.5582508682467, 1090.3560174219692, 1076.7379455811042,
                647.5281847355375, 321.697017419185, 157.11330197158256, 176.84446052687775, 186.5792050831855, 255.29466401569954,
                282.6181087245075, 939.727763378969, 2136.1018863118343, 1627.4300758813029, 539.215289601611, 229.19268550806882,
                107.46327766328407, 84.87582590867348, 126.07169689399987, 155.3758564515583, 106.53659209868633, 111.12027207117487, 88.15867583419072};
        double[] initPara = new double[]{20, 70, 78, 0.5, 0.1, 0.3, 0, 20, 2, 0.6, 0.2, 0.8, 0.8};
        double[] P = new double[]{10, 24.1, 20.4, 18.3, 10.1, 5.5, 0.6, 3.1, 1.9, 4.6, 5, 4.8, 36.2, 29, 6, 3.6, 0.4, 0, 0.5, 3.8, 0, 1.8, 0.2, 0.3};
        double[] EI = new double[24];
        ArrayList<XajModel.state> states = new ArrayList<>();
        for (int i = 0; i < P.length; i++) {
            XajModel.state state = new XajModel.state();
            state.setP(P[i]);
            state.setE(EI[i]);
            states.add(state);
        }
        double[] others = new double[]{0, 70, 80, 20, 2, 40, 20};
        model.setStates(states);
        model.setOthers(others);
        double[] bestPara = XajModel.getBestPara(model, 50, 1000, target, initPara, 0.0001);
        model = model.SetSoilWaterStorageParam(bestPara[0], bestPara[1], bestPara[2]);
        model = model.SetEvapotranspirationParam(bestPara[3], bestPara[4]);
        model = model.SetRunoffGenerationParam(bestPara[5], bestPara[6]);
        model = model.SetSourcePartitionParam(bestPara[7], bestPara[8], bestPara[9], bestPara[10]);
        model = model.SetRunoffConcentrationParam(bestPara[11], bestPara[12], 537);

        model = model.ComputeRunoffGeneration(P, EI, others[0], others[1], others[2]);
        model = model.ComputeSourcePartition(others[3], others[4]);
        model = model.ComputeRunoffConcentration(others[5], others[6], others[4]);
        double[] result = model.GetResult();
        System.out.println(Arrays.toString(result));
    }

    @Test
    public void SACmodelTest() {
        double[] initPara = new double[]{17, 78, 24, 24, 80, 30, 110, 0.03, 0.6,
                0.02, 0.4, 0.1, 0.03, 12, 0.1, 0.008, 0.8, 0.7, 0.9, 0.32, 0.89};
        SACmodel model = SACmodel.getExcelForSAC("./xmlFortest/SACdata.xls");
        model.setInitValue(new double[]{20, 90, 20, 20, 90, 10, 80, 100});
        double[] target = new double[]{49.7328373015873, 176.004162754117, 312.745328469911, 302.980670685312,
                178.410972683550, 98.8088489395449, 60.9882809492616, 40.8524087453631, 29.0043042578440, 21.4424173684383,
                16.3570247833621, 12.9813931920615, 10.8170429646819, 9.41138383680922, 8.45871058027875, 7.77506248300261,
                7.25583558448034, 6.84297516586878, 6.50414774147636, 5.2};
        double[] bestPara = SACmodel.getBestPara(model, 20, 1000, target, initPara, 0.0001);
        String[] paraName = new String[]{"AUZTWM", "ALZTWM", "UZTWM", "UZFWM", "LZTWM", "LZFSM", "LZFPM", "RIVA",
                "UZK", "PCTIM", "RESERV", "LZSK", "LZPK", "ZPERC", "PFREE", "ADIMP", "CI", "CGS", "CGP", "CR", "PAREA"};
        HashMap<String, String> map = new HashMap<String, String>();
        for (int i = 0; i < bestPara.length; i++) {
            map.put(paraName[i], String.valueOf(bestPara[i]));
        }
        BeanRefUtil.setFieldValue(this, map);
        model.SACmain();
        double[] Q = new double[model.getStates().size()];
        for (int i = 0; i < model.getStates().size(); i++) {
            Q[i] = model.getStates().get(i).getQ();
        }
        System.out.println(Arrays.toString(Q));
    }

    @Test
    public void TankModelTest() {
        double[] Ps = new double[]{5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 2, 1, 0};
        double[] time = new double[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
        double[] target = new double[]{0.477, 0.563, 0.744, 1.18, 1.966, 3.315, 4.269, 5.271, 5.859, 4.943, 4.365, 3.97, 3.758};
        Tank_model model = new Tank_model(3, 1);
        ArrayList<Tank_model.state> states = model.getStates();
        for (int i = 0; i < time.length; i++) {
            Tank_model.state state = new Tank_model.state();
            state.P = Ps[i];
            state.time = String.valueOf(time[i]);
            state.R = target[i];
            states.add(state);
        }
//        double[] initPara = new double[]{0, 0.45, 0.23, 0.16, 0.2, 0.34, 0.05, 0.03, 0.02, 0.008, 30, 15, 10};
        // 水箱初始高度(i+j)+水箱下渗能力(i+j)+双孔高度1(j)+双孔高度2(j)+双孔出流1(j)+双孔出流2(j)+单孔高度(i)+单孔出流(i)
        double[] initPara = new double[]{0, 0, 0, 0, 0.2, 0.1, 0.05, 0, 45, 15, 0.3, 0.06, 10, 15, 15, 0.06, 0.03, 0.005};
        double[] bestPara = model.getBestPara(model, 20, 1000, target, initPara, 0.0001);
        double[] Rs = model.paraOptimize(bestPara);
        System.out.println(Arrays.toString(Rs));
    }

    @Test
    public void ShaanXimodelTest() {
        ShaanXimode model = new ShaanXimode();
        int type = 0;
        double[] tempps = new double[]{0.0, 3.5, 5.7, 5.7, 5.7, 4.3, 2.5};
        ArrayList<Double> ps = Container.toArrayList(tempps);
        for (int i = 0; i < 19; i++) {
            ps.add(0.0);
        }
        ArrayList<Double> es = new ArrayList<Double>();
        for (int i = 0; i < ps.size(); i++) {
            es.add(0.0);
        }
        double[] UH = new double[]{0, 0, 0, 0.02, 0.03, 0.08, 0.12, 0.15, 0.12, 0.1, 0.08, 0.07, 0.06, 0.05, 0.04, 0.04, 0.02, 0.01, 0.01, 0.01};
        double dt = 0.2;
        double area = 187;
        double winit = 15;
        double[] initPara = new double[]{150, 20, 7.1};
        double[] target = new double[]{0, 0, 0, 0, 1, 4.5, 12.5, 29.6, 54.3, 83.2, 106, 116, 109.1, 93.6, 77.6, 65.2, 54.8, 45.2, 35.5, 26.4, 18.7, 12.7, 8.4, 5.6, 2.8, 0.8};
        model.setTarget(target);
        model = model.modelInit(type, ps, es, dt, area, winit, UH);
        model.setBx(1);
        double[] bestPara = model.getBestPara(model, 20, 1000, target, initPara, 0.0001);
        model = model.H_init(bestPara[0], bestPara[1], bestPara[2], 0.05);
        model = model.soilWCount();
        double[] res = model.UnitHydrograph(model.getRs(), new double[model.getRs().length]);
        System.out.println(Arrays.toString(res));
    }

    @Test
    public void HBVmodelTest() throws IOException, DocumentException {
        HBV_model model = new HBV_model();
        double[] initPara = new double[]{5, 150, 4, 0.04, 3.5, 0.12, 0.05, 0.02, 0.02, 150};
        model.initModel("./xmlFortest/HBV/HBVmodel.xml");
        double[] bestPara = model.getBestPara(30, 1000, model.getQo(), initPara, 0.0001);
        model.initPara(bestPara);
        model.HBV_main();
        double[] Rs = model.getR();
        System.out.println(Arrays.toString(Rs));
    }
}
