package hu.bkoncser;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.BlockRealMatrix;

import static javafx.scene.input.KeyCode.H;

/**
 * Created by bkoncser on 2016-10-26.
 */
public class Main {

    public static void main(String[] args){

        MyMatrix H = new MyMatrix();

        int L = H.parseInputToMatrix();
        int I = H.getRowDimension();
        int J = H.getColumnDimension();
        H.printMatrix();

        MyMatrix U = new MyMatrix(L, I);
        U.setParamAlpha(1.0);
        U.generate();

        MyMatrix V = new MyMatrix(L,J);
        V.setParamAlpha(1.0);
        V.generate();


        //for(int i = 0; i<100; i++) {
           // System.out.println(i);
            U.refresh(V, H, true);
            V.refresh(U, H, false);
        //}


        BlockRealMatrix Hopt = U.getMatrix().transpose();
       // new MyMatrix(Hopt).printMatrix();
        System.out.println("---------------------");
        Hopt = Hopt.multiply(V.getMatrix());

        new MyMatrix(Hopt).printMatrix();

    }

}
