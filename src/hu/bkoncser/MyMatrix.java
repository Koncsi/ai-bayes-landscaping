package hu.bkoncser;

import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;

/**
 * Created by Koncsi on 26/10/2016.
 */
public class MyMatrix {

    private BlockRealMatrix matrix;

    private Double paramBeta;


    private Double paramAlpha;

    public MyMatrix(){
        matrix = null;
        paramBeta = null;
    }

    public MyMatrix(BlockRealMatrix mtx) {
        matrix = mtx;
        paramBeta = 0.0;
    }

    public MyMatrix(int i, int l) {
        double[][] identity = new double[i][l];
        for(int a = 0; a < i; a++){
            for (int b = 0; b < l; b++){
                if(b==a)
                    identity[a][b] = 1.0;
                else
                    identity[a][b] = 0.0;
            }
        }
        matrix = new BlockRealMatrix(identity);
    }

    int parseInputToMatrix(){

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

            return Integer.parseInt(inputParameters[2]);

        }catch (Exception e){
            e.printStackTrace();
        }
       return  0;
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

    MyMatrix multiply(MyMatrix other){
        return new MyMatrix(this.matrix.multiply(other.getMatrix()));
    }

    MyMatrix scalarMulri(double d){
        return new MyMatrix((BlockRealMatrix) this.matrix.scalarMultiply(d));
    }

    BlockRealMatrix getMatrix(){
        return  matrix;
    }

    double[][] getDoubles(){
        return matrix.getData();
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

    public Double getParamAlpha() {
        return paramAlpha;
    }

    public void setParamAlpha(Double paramAlpha) {
        this.paramAlpha = paramAlpha;
    }

    public void generate() {
        MyMatrix I = new MyMatrix(this.getRowDimension(), this.getRowDimension());

        double[] nulls = new double[I.getRowDimension()];
        for(int i = 0; i < I.getRowDimension();i++)
            nulls[i] = 0.0;


        MultivariateNormalDistribution normalDist = new MultivariateNormalDistribution(nulls, I.scalarMulri(1.0/this.getParamAlpha()).getDoubles());

        for(int i = 0; i < matrix.getColumnDimension(); i++)
            matrix.setColumn(i, normalDist.sample());

    }

    public void refresh(MyMatrix uv, MyMatrix h, boolean forU) {
        MultivariateNormalDistribution normalDist;
        for(int index = 0; index < this.getColumnDimension(); index++){
            MyMatrix labda = calculateLabda(h.getParamBeta(),uv);
            double[] pszi = calculatePszi(labda,h,uv,index,forU);
            normalDist = new MultivariateNormalDistribution(pszi,invert(labda).getData());
            matrix.setColumn(index, normalDist.sample());
        }
    }

    BlockRealMatrix invert(MyMatrix mtx){
        RealMatrix test = new LUDecomposition(mtx.getMatrix()).getSolver().getInverse();
        return new BlockRealMatrix(test.getData());
    }

    private double[] calculatePszi(MyMatrix labda, MyMatrix h, MyMatrix uv, int index, boolean forU) {

        BlockRealMatrix mtx = invert(labda);

        mtx.scalarMultiply(h.getParamBeta());
        double hij;
        BlockRealMatrix resultMtx;
        if(forU) {
            hij = h.getMatrix().getRow(index)[0];

            resultMtx = uv.getMatrix().getColumnMatrix(0);
            resultMtx = (BlockRealMatrix) resultMtx.scalarMultiply(hij);
            for (int i = 1; i < uv.getColumnDimension(); i++) {
                hij = h.getMatrix().getRow(index)[i];
                BlockRealMatrix temp = uv.getMatrix().getColumnMatrix(0);
                resultMtx = resultMtx.add(temp.scalarMultiply(hij));
            }
        }
        else{
            hij = h.getMatrix().getRow(0)[index];

            resultMtx = uv.getMatrix().getColumnMatrix(0);
            resultMtx = (BlockRealMatrix) resultMtx.scalarMultiply(hij);
            for (int i = 1; i < uv.getColumnDimension(); i++) {
                hij = h.getMatrix().getRow(i)[index];
                BlockRealMatrix temp = uv.getMatrix().getColumnMatrix(0);
                resultMtx = resultMtx.add(temp.scalarMultiply(hij));
            }
        }
        double[] ret = mtx.preMultiply(resultMtx.getColumn(0));

        return ret;
    }

    private MyMatrix calculateLabda(Double paramBeta, MyMatrix uv) {
        BlockRealMatrix I = new MyMatrix(uv.getRowDimension(),uv.getRowDimension()).getMatrix();

        BlockRealMatrix column = uv.getMatrix().getColumnMatrix(0);
        BlockRealMatrix mtx = column.multiply(column.transpose());

        for(int i = 1; i<uv.getColumnDimension(); i++){
            column = uv.getMatrix().getColumnMatrix(i);
            BlockRealMatrix tempMtx  = column.multiply(column.transpose());
            mtx = mtx.add(tempMtx);
        }
        mtx = (BlockRealMatrix) mtx.scalarMultiply(paramBeta);
        mtx = mtx.add(I.scalarMultiply(this.getParamAlpha()));
        return new MyMatrix(mtx);
    }
}
