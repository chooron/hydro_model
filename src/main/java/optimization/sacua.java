package optimization;


import java.util.*;

/**
 * 全局最优算法
 * SCE-UA 算法的基本思路是将基于确定性的复合型搜索技术和自然界中的生物竞争进化原理相结合.
 * 算法的关键部分为竞争的复合型进化算法( CCE) .在CCE 中, 每个复合型的顶点都是潜在的父辈,、
 * 都有可能参与产生下一代群体的计算.每个子复合型的作用如同一对父辈.随机方式在构建子复合型中的应用,
 * 使得在可行域中的搜索更加彻底.用SCE-UA 算法求解最小化问题的具体步骤如下：
 * step1: 初始化，假定待优化问题是n维问题，选取参与进化的复合型个数p，和每个复合型所包含的顶点数目m，得到样本点个数维s=p*m
 * step2：产生样本点，在可行域内随机产生s个样本点，分别计算每个点的函数值
 * step3：样本点根据函数值升序排序，并根据大小重新编号
 * step3：划分复合型整体，将样本集划分成p个复合型，每个复合型有m点
 * step4：复合型进化
 * step5：复合型混合
 * step6：收敛性判断
 */
public class sacua {
    double xo;//待优化参数的初始值
    double b1;// 待优化参数的下限
    double bu;// 待优化参数的上限
    double maxn;//进化过程中函数最大调用次数
    double kstop;//集中之前最大演化数
    double pcento;//允许在kstop循环前收敛的百分比变化
    double peps;
    double iseed;//生成随机数种子
    double iniflg;//初始参数矩阵的标志
    int ngs;//复合型个数
    int npg;//复合型顶点数目
    int nps;//简单型
    int nspl;//每个复合型进化的迭代数
    int mings;//进化过程复合型的最小数目
    boolean compileFlag;

    double[] BESTX;
    double BESTFX;

    /**
     * 默认构造
     */
    public sacua() {
    }

    /**
     * 自定义构造
     */
    public void compile(double maxn, double kstop, double pcento, double peps, double iseed, double iniflg) {
        this.compileFlag = true;
        this.maxn = maxn;
        this.kstop = kstop;
        this.pcento = pcento;
        this.peps = peps;
        this.iseed = iseed;
        this.iniflg = iniflg;
    }

    /**
     * 主优化函数
     */
//    public double[] optimizeMain(double[] x0, double[][] bounds) {
//        if (!compileFlag) {
//            compile(1000.0, 10.0, 0.1, 0.001, -1, 0);
//        }
//        int nopt = x0.length;  // 点的维度
//        this.npg = nopt * 2 + 1;
//        this.nps = nopt + 1;
//        nspl = npg;
//        mings = ngs;
//        int npt = npg * ngs;
//        // 上下界之差
//        double[] boundDiff = new double[bounds.length];
//        for (int i = 0; i < bounds.length; i++) boundDiff[i] = bounds[i][1] - bounds[i][0];
//        // 初始化种群
//        Random random = new Random((long) iseed);
//        double[][] x = new double[npt][nopt];
//        double[] rand = RandomGenerator.EvenDTVector(nopt, (long) iseed);
//        for (int i = 0; i < npt; i++) {
//            for (int j = 0; j < nopt; j++) {
//                x[i][j] = (bounds[j][0] + rand[j]) * boundDiff[j];
//            }
//        }
//        // 是否包含x0
//        if (iniflg == 1) System.arraycopy(x0, 0, x[0], 0, nopt);
//        // 计算函数值
//        double nloop = 0;
//        double[] xf = new double[npt];
//        Map<Integer, Double> xfi = new TreeMap<Integer, Double>();
//        for (int i = 0; i < npt; i++) xf[i] = xajModel(x[i]);
//        double f0 = xf[0];
//        // 根据函数值排序
//        double[][] xfBySort = sortWithIndex(xf);
//        int[] sortIndex = new int[xfBySort.length];
//        for (int i = 0; i < xfBySort.length; i++) sortIndex[i] = (int) xfBySort[i][0];
//        x = getByListIndex(x, sortIndex);
//        // 记录最好最差的点
//        double[] bestx = x[0];
//        double[] worstx = x[npt - 1];
//        double bestfx = xfBySort[0][1];
//        double worstfx = xfBySort[npt - 1][1];
//        BESTX = bestx;
//        BESTFX = bestfx;
//        // 计算标准误差
//        double xnstd = ListUtil.variance(xfBySort[0]);  // 未使用
//        // 计算收敛指标
//
//        return null;
//    }

    public static void main(String[] args) {
        double[] test = new double[]{1, 2, 4, 2, 6, 8, 10, 2, 4, 2};
        double[][] doubles = sortWithIndex(test);
        System.out.println(Arrays.deepToString(doubles));
    }

    public static double[][] sortWithIndex(double[] xfi) {
        Comparator<Map.Entry<Integer, Double>> vauleComparator = new Comparator<Map.Entry<Integer, Double>>() {
            @Override
            public int compare(Map.Entry<Integer, Double> o1, Map.Entry<Integer, Double> o2) {
                return (int) (o1.getValue() - o2.getValue());
            }
        };
        TreeMap<Integer, Double> map = new TreeMap<Integer, Double>();
        for (int i = 0; i < xfi.length; i++) {
            map.put(i, xfi[i]);
        }
        List<Map.Entry<Integer, Double>> list = new ArrayList<Map.Entry<Integer, Double>>(map.entrySet());
        Collections.sort(list, vauleComparator);
        double[][] result = new double[list.size()][2];
        for (int i = 0; i < list.size(); i++) {
            result[i][0] = list.get(i).getKey();
            result[i][1] = list.get(i).getValue();
        }
        return result;
    }

    public static double[][] getByListIndex(double[][] origin, int[] indexList) {
        double[][] result = new double[origin.length][origin[0].length];
        for (int i = 0; i < origin.length; i++) {
            result[i] = origin[indexList[i]];
        }
        return result;
    }

    public static double[] getGrng(double[][] x, double[] bound){
        return null;
    }
}
