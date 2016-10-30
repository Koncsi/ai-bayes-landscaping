package hu.bkoncser;

import org.apache.commons.math3.linear.BlockRealMatrix;
import java.util.ArrayList;

/**
 * Created by bkoncser on 2016-10-26.
 */
public class Main {

    public static void main(String[] args){


        MyMatrix H = new MyMatrix();

        int L = H.parseInputToMatrix();
        int I = H.getRowDimension();
        int J = H.getColumnDimension();

        MyMatrix U = new MyMatrix(L, I,false);
        U.setParamAlpha(1.0);
        U.generate();

        MyMatrix V = new MyMatrix(L,J,false);
        V.setParamAlpha(1.0);
        V.generate();

        ArrayList<MyMatrix> resultsU = new ArrayList<>();
        ArrayList<MyMatrix> resultsV = new ArrayList<>();
        for(int i = 0; i<120; i++) {
            U.refresh(V, H, true);
            V.refresh(U, H, false);
            if(i > 30){
                resultsU.add(U);
                resultsV.add(V);
            }
        }

        double[][] uResult = new double[U.getRowDimension()][U.getColumnDimension()];
        for(int i = 0; i < U.getRowDimension(); i++){
            for(int j = 0; j < U.getColumnDimension(); j++){
                double val = 0.0;
                for(int k = 0; k < resultsU.size(); k++)
                    val += resultsU.get(k).getMatrix().getData()[i][j];
                uResult[i][j] = val/resultsU.size();
            }
        }

        double[][] vResult = new double[V.getRowDimension()][V.getColumnDimension()];
        for(int i = 0; i < V.getRowDimension(); i++){
            for(int j = 0; j < V.getColumnDimension(); j++){
                double val = 0.0;
                for(int k = 0; k < resultsV.size(); k++)
                    val += resultsV.get(k).getMatrix().getData()[i][j];
                vResult[i][j] = val/resultsV.size();
            }
        }

        new MyMatrix(new BlockRealMatrix(uResult).transpose()).printMatrix();
        System.out.println();
        new MyMatrix(new BlockRealMatrix(vResult).transpose()).printMatrix();
    }

}
