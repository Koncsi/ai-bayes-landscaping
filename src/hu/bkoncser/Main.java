package hu.bkoncser;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.BlockRealMatrix;

/**
 * Created by bkoncser on 2016-10-26.
 */
public class Main {

    public static void main(String[] args){

        MyMatrix myMatrix = new MyMatrix();

        myMatrix.parseInputToMatrix();

        myMatrix.printMatrix();

    }

}
