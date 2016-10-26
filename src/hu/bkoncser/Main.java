package hu.bkoncser;

import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * Created by bkoncser on 2016-10-26.
 */
public class Main {

    public static void main(String[] args){


        NormalDistribution nb = new NormalDistribution();




        System.out.println("working...  " + nb.probability(5.0f));
    }

}
