import cern.colt.matrix.tdouble.DoubleMatrix2D;

import java.util.ArrayList;

public class MöbiusVoting {
    double[][] c1;
    double[][] c2;
    double[][] CorrespondenceMatrix;
    int K;
    double epsilon;
    public MöbiusVoting(double[][] c1,double[][] c2,int threshold,double epsilon){
        this.c1=c1;
        this.c2=c2;
        if(threshold!=-1) this.K=threshold;
        else this.K=(c1.length*40)/100;
        this.epsilon=epsilon;
    }
    private int[][] sampleRandomPoint(){
        /* samples 3 point inside each set*/
        throw new Error("to be implemented");
    }
    private DoubleMatrix2D createMobiusTransform(int[][] inputPoints){
        /* creates the möbius transform using the sampled points*/
        throw new Error("to be implemented");
    }
    private double[][] applyMobius(double[][] c,DoubleMatrix2D Mobius){
        /* this methods apply Mobius to the c plane with Mobius Matrix*/
        throw new Error("to be implemented");
    }

    private int[][] findMutuallyNearestNeighbors(double[][] mc1,double[][] mc2){
        /*this methods find mutually nearest Neighbors in the two planes constructed with mobius transform*/
        throw new Error("to be implemented");
    }
    private void updateCorrespondenceMatrix(){
        /*updates correspondence Matrix, if the number of mutually closes pairs >=K:
        * TO DO:
        * 1) compute the deformation energy
        * 2) For each mutually nearest neighbors: updte the correspondence matrix
        * */
        throw new Error("to be implemented");
    }

    public int[][] computeFuzzyCorrespondenceMatrix(){
        /* This function first compute the correspondence matrix using previous algorithm
        * Then it computes the FUZZY correspondence matrix.
        * Gives back the index of the final list of correspondences.
         */
        throw new Error("to be implemented");
    }
}
