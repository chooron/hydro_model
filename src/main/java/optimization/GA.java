package optimization;

import optimization.ga.*;
import tranditionModel.XajModel;
import optimization.ga.Population;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class GA {
    public static void main(String[] args) {
        int populationSize = 50;
        Population<MyVector> population = createInitialPopulation(populationSize);
        double[] target = new double[]{145.04165680890637, 460.7913823034037, 842.5582508682467, 1090.3560174219692, 1076.7379455811042,
                647.5281847355375, 321.697017419185, 157.11330197158256, 176.84446052687775, 186.5792050831855, 255.29466401569954,
                282.6181087245075, 939.727763378969, 2136.1018863118343, 1627.4300758813029, 539.215289601611, 229.19268550806882,
                107.46327766328407, 84.87582590867348, 126.07169689399987, 155.3758564515583, 106.53659209868633, 111.12027207117487, 88.15867583419072};
        Fitness<MyVector, Double> fitness = new MyVectorFitness(target);

        GeneticAlgorithm<MyVector, Double> ga = new GeneticAlgorithm<MyVector, Double>(population, fitness);

        addListener(ga);

        ga.evolve(1000);
    }

    private static Population<MyVector> createInitialPopulation(int populationSize) {

        Population<MyVector> population = new Population<MyVector>();
        MyVector base = new MyVector();
        for (int i = 0; i < populationSize; i++) {
            // each member of initial population
            // is mutated clone of base chromosome
            MyVector chr = base.mutate();
            population.addChromosome(chr);
        }
        return population;
    }
    /**
     * After each iteration Genetic algorithm notifies listener
     */
    private static void addListener(GeneticAlgorithm<MyVector, Double> ga) {
        // just for pretty print
        System.out.println(String.format("%s\t%s\t%s", "iter", "fit", "chromosome"));

        // Lets add listener, which prints best chromosome after each iteration
        ga.addIterationListener(new IterartionListener<MyVector, Double>() {

            private final double threshold = 1e-5;

            @Override
            public void update(GeneticAlgorithm<MyVector, Double> ga) {

                MyVector best = ga.getBest();
                double bestFit = ga.fitness(best);
                int iteration = ga.getIteration();

                // Listener prints best achieved solution
                System.out.println(String.format("%s\t%s\t%s", iteration, bestFit, best));

                // If fitness is satisfying - we can stop Genetic algorithm
                if (bestFit < this.threshold) {
                    ga.terminate();
                }
            }
        });
    }

    public static class MyVector implements Chromosome<MyVector>,Cloneable{
        private static final Random random = new Random();

        private final double[] vector = new double[]{20,70,78,0.5,0.1,0.3,0,20,2,0.6,0.2,0.8,0.8};
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
            double mutationValue = (random.nextDouble()- random.nextDouble())*0.5;
            result.vector[index] += mutationValue;
            result.vector[index] =Math.max(0,result.vector[index]);
            return result;
        }
        @Override
        protected MyVector clone() {
            MyVector clone = new MyVector();
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

        public MyVectorFitness(double[] target) {
            this.target = target;
        }

        @Override
        public Double calculate(MyVector chromosome) {
            double delta = 0;
            double[] v = chromosome.getVector();
            double[] q = new XajModel().paraOptimize(v);
//            System.out.println(Arrays.toString(q));
            for (int i = 0; i < q.length; i++) {
                delta += this.sqr(q[i] - this.target[i]);
            }
            return delta;
        }

        private double sqr(double x) {
            return x * x;
        }
    }
}
