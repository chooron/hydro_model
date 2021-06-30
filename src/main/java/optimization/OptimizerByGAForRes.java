package optimization;


import optimization.ga.*;
import tranditionModel.parent.parentModel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * @author jing xin
 * 面向概念型水文模型的参数优化，适用于XAJ、SAC、ShaanXi和Tank四个模型
 */
public class OptimizerByGAForRes {
    parentModel model;
    int populationSize;
    int evolveTime;
    double[] target;
    double[] initPara;
    double threshold;

    /**
     * 概念型水文模型参数优化对象
     *
     * @param model          概念模型
     * @param populationSize 种群大小
     * @param evolveTime     迭代周期
     * @param target         实际结果
     * @param initPara       初始参数
     * @param threshold      终止条件
     */
    public OptimizerByGAForRes(parentModel model, int populationSize, int evolveTime, double[] target, double[] initPara, double threshold) {
        this.model = model;
        this.populationSize = populationSize;
        this.evolveTime = evolveTime;
        this.target = target;
        this.initPara = initPara;
        this.threshold = threshold;
    }

    public double[] GAmain() {
        Population<MyVector> population = createInitialPopulation(populationSize);
        Fitness<MyVector, Double> fitness = new MyVectorFitness(target,model);
        GeneticAlgorithm<MyVector, Double> ga = new GeneticAlgorithm<MyVector, Double>(population, fitness);
        addListener(ga);
        ga.evolve(evolveTime);
        MyVector best = ga.getBest();
        return best.getVector();
    }

    private Population<MyVector> createInitialPopulation(int populationSize) {
        Population<MyVector> population = new Population<MyVector>();
        MyVector base = new MyVector(initPara);
        for (int i = 0; i < populationSize; i++) {
            MyVector chr = base.mutate();
            population.addChromosome(chr);
        }
        return population;
    }

    /**
     * After each iteration Genetic algorithm notifies listener
     */
    private void addListener(GeneticAlgorithm<MyVector, Double> ga) {
        // just for pretty print
        System.out.printf("%s\t%s\t%s%n", "iter", "fit", "chromosome");
        List<MyVector> history = new ArrayList<MyVector>();
        // Lets add listener, which prints best chromosome after each iteration
        ga.addIterationListener(new IterartionListener<MyVector, Double>() {
            @Override
            public void update(GeneticAlgorithm<MyVector, Double> ga) {

                MyVector best = ga.getBest();
                double bestFit = ga.fitness(best);
                int iteration = ga.getIteration();

                // Listener prints best achieved solution
                System.out.printf("%s\t%s\t%s%n", iteration, bestFit, best);

                // If fitness is satisfying - we can stop Genetic algorithm
                if (bestFit < threshold) {
                    ga.terminate();
                }
            }
        });
    }

    public static class MyVector implements Chromosome<MyVector>, Cloneable {
        private static final Random random = new Random();

        public double[] vector; //

        public MyVector(int l) {
            vector = new double[l];
        }

        public MyVector(double[] vector) {
            this.vector = vector;
        }

        @Override
        public List<MyVector> crossover(MyVector other) {
            MyVector thisClone = this.clone();
            MyVector otherClone = other.clone();
            // 单点交叉
            int index = random.nextInt(this.vector.length - 1);
            for (int i = index; i < this.vector.length; i++) {
                double tmp = thisClone.vector[i];
                thisClone.vector[i] = otherClone.vector[i];
                otherClone.vector[i] = tmp;
            }
            return Arrays.asList(thisClone, otherClone);
        }

        @Override
        public MyVector mutate() {
            MyVector result = this.clone();

            // just select random element of vector
            // and increase or decrease it on small value
            int index = random.nextInt(this.vector.length);
//            double mutationValue = (random.nextDouble() - random.nextDouble());
            result.vector[index] = (double) (random.nextInt(440)/(440-229)+229)/100;
            result.vector[index] = Math.max(0, result.vector[index]);
            return result;
        }

        @Override
        protected MyVector clone() {
            MyVector clone = new MyVector(this.vector.length);
            System.arraycopy(this.vector, 0, clone.vector, 0, this.vector.length);
            return clone;
        }

        public double[] getVector() {
            return this.vector;
        }

        @Override
        public String toString() {
            return Arrays.toString(this.vector);
        }
    }

    /**
     * Fitness function, which calculates difference between chromosomes vector
     * and target vector
     */
    public static class MyVectorFitness implements Fitness<MyVector, Double> {

        private double[] target;
        public parentModel model;

        public MyVectorFitness(double[] target, parentModel model) {
            this.target = target;
            this.model = model;
        }

        public MyVectorFitness(double[] target) {
            this.target = target;
        }

        @Override
        public Double calculate(MyVector chromosome) {
            double delta = 0;
            double[] v = chromosome.getVector();
            double[] N = model.paraOptimize(v);
//            System.out.println(Arrays.toString(q));
            for (int i = 0; i < N.length; i++) {
                delta += N[i];
            }
            return 1/(delta);
        }


//        private double sqr(double x) {
//            return x * x;
//        }
    }
}


