package hu.bkoncser;

import org.apache.commons.math3.linear.BlockRealMatrix;

import java.io.BufferedReader;
import java.io.InputStreamReader;

/**
 * Created by Koncsi on 26/10/2016.
 */
public class MyMatrix {

    private BlockRealMatrix matrix;

    private Double paramBeta;

    public MyMatrix(){
        matrix = null;
        paramBeta = null;
    }

    void parseInputToMatrix(){

        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        String[] inputParameters;
        String[] inputValues;

        try {

            String input = reader.readLine();
            inputParameters = input.split(",");

            double [][] values = new double[Integer.parseInt(inputParameters[0])][Integer.parseInt(inputParameters[1])];

            for(int i = 0; i < Integer.parseInt(inputParameters[0]); i++){
                input = reader.readLine();
                inputValues = input.split(",");
                for(int j = 0; j < Integer.parseInt(inputParameters[1]); j++){
                    values[i][j] = Double.parseDouble(inputValues[j]);
                }
            }

            matrix = new BlockRealMatrix(values);
            paramBeta = Double.parseDouble(inputParameters[3]);

        }catch (Exception e){
            e.printStackTrace();
        }
    }

    void printMatrix(){

        double [][] values = matrix.getData();
        for(int i = 0; i< matrix.getRowDimension(); i++){
            for(int j =0; j<matrix.getColumnDimension(); j++){
                System.out.print(values[i][j] + ",");
            }
            System.out.println();
        }

    }

    int getRowDimension(){
        return matrix.getRowDimension();
    }

    int getColumnDimension(){
        return matrix.getColumnDimension();
    }

    Double getParamBeta(){
        return paramBeta;
    }

}
