package algo;

import Jcg.polyhedron.Vertex;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import com.sun.applet2.preloader.event.ErrorEvent;
import meshmanager.SurfaceMesh;

import java.util.ArrayList;

public class MobiusVoting {
    /* Once the parametrization on the complexe plane is done we use this Voting class to determines which sampled points are corresponding the most.
     * c1 and c2 represents the coordinates from the sampled points in the complex plane.
     * Those coordinates are obtained associating each point with its nearest mid-edge vertex then projecting on the complex plane.
     *
     */
    double[][] c1;
    double[][] c2;
    double[][] transformedc1;
    double[][] transformedc2;
    double[][] CorrespondenceMatrix;

    int K;
    double epsilon;

    public MobiusVoting(double[][] c1, double[][] c2, int threshold, double epsilon) {
        this.c1 = c1;
        this.c2 = c2;
        if (threshold != -1) this.K = threshold;
        else this.K = (c1.length * 10) / 100;
        this.epsilon = epsilon;
    }

    //thereafter we define useful function on complex number represented as double[2]
    private double[] complexProd(double[] a,double[] b){
        double[] c= {a[0]*b[0]-a[1]*b[1],a[1]*b[0]+a[0]*b[1]};
        return c;
    }
    private double[] complexDiff(double[] a,double[] b){
        double[] c={a[0]-b[0],a[1]-b[1]};
        return c;
    }
    private double complexDist(double[] a, double[] b){
        double dist=Math.pow(a[0]-b[0],2)+Math.pow(a[1]-b[1],2);
        return dist;
    }
    private double[] complexDivision(double[] a,double[] b){
        double[] c=this.complexProd(a,new double[]{b[0],-b[1]});
        double e=b[0]*b[0]+b[1]*b[1];
        if(e==0){
            //this will be the case for the last sampled point.... something need to be done to cope with this issue
            return new double[]{1000000,1000000};
        }
        double[] d={c[0]/e,c[1]/e};
        return d;
    }
    private double[] complexAdd(double[] a,double[] b){
        double[] add={a[0]+b[0],a[1]+b[1]};
        return add;
    }

    //our methods for the voting
    private int[][] sampleRandomPoint() {
        /* samples 3 point inside each set*/
        int length=this.c1.length;
        int[][] sampled=new int[2][3];
        for(int i=0;i<2;i++){
            for(int j=0;j<3;j++){
                boolean isValid=false;
                int idx=0;
                while(!isValid) {
                    idx= (int) (Math.random() * length);
                    isValid = true;
                    for (int h = 0; h < j; h++) {
                        if (sampled[i][h] == idx) {
                            isValid = false;
                        }
                    }
                }
                sampled[i][j] =idx;
            }
        }
        return sampled;
    }

    private double[][][] matrixProd(double[][][] mat1,double[][][] mat2){
        double[][][] result=new double[mat1.length][mat2[0].length][2];
        for(int i=0;i<mat1.length;i++){
            for(int j=0;j<mat2[0].length;j++){
                for(int k=0;k<mat1[0].length;k++){
                    result[i][j]=this.complexAdd(result[i][j],this.complexProd(mat1[i][k],mat2[k][j]));
                }
            }
        }
        return result;
    }

    private double[][][][] createMobiusTransform(int[][] inputPoints) {
        /* creates möbius transforms for each mesh using the sampled points
        *  It returns a double list of 4 elements such as a=double[0][1], b=double[0][1], d=double[0][2], d=double[0][3] for m1
        *  if the computed Möbius matrix is of the form (a*z+b)/(c*z+d)
        * We aligned our initial points on 3 constants points, we could choose e(i*pi*2/3*j), j=1,2,3 as used in the paper
        * But we can also choose 0,1,infinity*/
        /*double[] y1={1/2,Math.sqrt(3)/2};
        double[] y2={-1/2,Math.sqrt(3)/2};
        double[] y3={-1,0};*/
        double[][][] localisationMatrix={{{-5.55111512*Math.pow(10,-17),-0.57735027},{-0.5,-0.28867513}},
                                        {{-1.11022302*Math.pow(10,-16),-0.57735027},{-0.5,-0.28867513}}};
        double[] c1z1=this.c1[inputPoints[0][0]];
        double[] c1z2=this.c1[inputPoints[0][1]];
        double[] c1z3=this.c1[inputPoints[0][2]];
        double[] c2z1=this.c2[inputPoints[1][0]];
        double[] c2z2=this.c2[inputPoints[1][1]];
        double[] c2z3=this.c2[inputPoints[1][2]];

        double[][][] mobTransform=new double[2][2][2];
        mobTransform[0][0]=this.complexDiff(c1z2,c1z3);
        mobTransform[0][1]=this.complexDiff(this.complexProd(c1z1,c1z3),this.complexProd(c1z2,c1z1));
        mobTransform[1][0]=this.complexDiff(c1z2,c1z1);
        mobTransform[1][1]=this.complexDiff(this.complexProd(c1z3,c1z1),this.complexProd(c1z3,c1z2));
        double[][][] mobTransform1=this.matrixProd(localisationMatrix,mobTransform);

        mobTransform=new double[2][2][2];
        mobTransform[0][0]=this.complexDiff(c2z2,c2z3);
        mobTransform[0][1]=this.complexDiff(this.complexProd(c2z1,c2z3),this.complexProd(c2z2,c2z1));
        mobTransform[0][0]=this.complexDiff(c2z2,c2z1);
        mobTransform[0][1]=this.complexDiff(this.complexProd(c2z3,c2z1),this.complexProd(c2z3,c2z2));
        double[][][] mobTransform2=this.matrixProd(localisationMatrix,mobTransform);

        double[][][][] mob=new double[2][2][2][2];
        mob[0]=mobTransform1;
        mob[1]=mobTransform2;
        return mob;
    }
    private void applyMobius(double[][][][] Mobius) {
        /* This methods apply Mobius to the c plane with Mobius Matrix
        *
        * */
        this.transformedc1=new double[this.c1.length][2];
        for(int i=0;i<this.c1.length;i++){
            double[] up=this.complexAdd(this.complexProd(this.c1[i],Mobius[0][0][0]),Mobius[0][0][1]);
            double[] down=this.complexAdd(this.complexProd(this.c1[i],Mobius[0][1][0]),Mobius[0][1][1]);
            this.transformedc1[i]=this.complexDivision(up,down);
        }
        this.transformedc2=new double[this.c2.length][2];
        for(int i=0;i<this.c2.length;i++){
            double[] up=this.complexAdd(this.complexProd(this.c2[i],Mobius[1][0][0]),Mobius[1][0][1]);
            double[] down=this.complexAdd(this.complexProd(this.c2[i],Mobius[1][1][0]),Mobius[1][1][1]);
            this.transformedc2[i]=this.complexDivision(up,down);
        }
        /*
        System.out.println("Typical value of c1");
        for(int i=0;i<this.c1.length;i=i+10){
            for(int j=0;j<2;j++){
                System.out.print(" "+transformedc1[i][j]);
            }

        }
        System.out.println("Typical value of c2");
        for(int i=0;i<this.c1.length;i=i+10){
            for(int j=0;j<2;j++){
                System.out.print(" "+transformedc2[i][j]);
            }
        }
        System.out.println("");
        System.out.println("");*/

    }
    private int[] findNearestNeighbor(int idx,boolean isC1){
        /* a greedy algorithm to find nearest neighbors*/
        double maxDist=-1;
        ArrayList<Integer> indexesArray=new ArrayList<>();
        if(isC1){
            //we need to look in transformedc2 to find the nearest neighbor of transformedc1[idx]!
            for(int i=0;i<this.transformedc2.length;i++){
               double dist=this.complexDist(transformedc1[idx],this.transformedc2[i]);
               if(dist<maxDist || maxDist==-1) {
                   indexesArray.clear();
                   maxDist = dist;
                   indexesArray.add(i);
               }
               else if(dist==maxDist){
                   indexesArray.add(i);
               }
            }
        }
        else{
            //we need to look in transformedc1 to find the nearest neighbor of transformedc1[idx]!
            for(int i=0;i<this.transformedc1.length;i++){
                double dist=this.complexDist(transformedc2[idx],this.transformedc1[i]);
                if(dist<maxDist || maxDist==-1) {
                    indexesArray.clear();
                    maxDist = dist;
                    indexesArray.add(i);
                }
                else if(dist==maxDist){
                    indexesArray.add(i);
                }
            }
        }
        int[] indexes=new int[indexesArray.size()];
        for(int i=0;i<indexesArray.size();i++){
            indexes[i]=indexesArray.get(i);
        }
        return indexes;
    }
    private ArrayList<int[]> findMutuallyNearestNeighbors() {
        /*this methods find mutually nearest Neighbors in the two planes constructed with mobius transform
        * Gives back the list of mutually nearest neighbor as an ArrayList of int[2] such as:
        * int[0] is the index of a point in c1
        * int[1] is the index of a point in c2
        * */
        ArrayList<int[]> c1NearestNeighbor=new ArrayList<>();
        ArrayList<int[]> c2NearestNeighbor=new ArrayList<>();
        //build nearest neighbor list
        for(int i=0;i<this.transformedc1.length;i++){
            c1NearestNeighbor.add(this.findNearestNeighbor(i,true));
        }
        for(int i=0;i<this.transformedc2.length;i++){
            c2NearestNeighbor.add(this.findNearestNeighbor(i,false));
        }
        //build mutually nearest closet number list.
        ArrayList<int[]> mutuallyClosest=new ArrayList<>();
        for(int i=0;i<c1NearestNeighbor.size();i++){
            for(int j:c1NearestNeighbor.get(i)){
                for(int h:c2NearestNeighbor.get(j)){
                    if(h==i){
                        mutuallyClosest.add(new int[]{i,j});
                    }
                }
            }
        }
        return mutuallyClosest;
    }
    private void updateCorrespondenceMatrix(ArrayList<int[]> mutuallyClosestPairs) {
        /*updates correspondence Matrix, if the number of mutually closes pairs >=K:
         * TO DO:
         * 1) compute the deformation energy
         * 2) For each mutually nearest neighbors: updte the correspondence matrix
         * */
        if(mutuallyClosestPairs.size()>this.K){
            double E=0;
            for(int i=0;i<mutuallyClosestPairs.size();i++){
                int[] pair=mutuallyClosestPairs.get(i);
                E+=Math.sqrt(this.complexDist(this.transformedc1[pair[0]],this.transformedc2[pair[1]]));
            }
            E=E/mutuallyClosestPairs.size();
            for(int i=0;i<mutuallyClosestPairs.size();i++){
                int[] pair=mutuallyClosestPairs.get(i);
                this.CorrespondenceMatrix[pair[0]][pair[1]]+=1./(this.epsilon+E);
            }
        }
    }
    public double[][] createConfidenceMatrix(int I){
        /* this algorithm creates the confidence Matrix by executing I mobius voting*/
        this.CorrespondenceMatrix=new double[this.c1.length][this.c2.length];
        for(int i=0;i<I;i++){
            int[][] sampled=this.sampleRandomPoint();
            double[][][][] Mobius=this.createMobiusTransform(sampled);
            this.applyMobius(Mobius);
            ArrayList<int[]> mutuallyNearest=this.findMutuallyNearestNeighbors();
            for(int j=0;j<3;j++){
                boolean isFound=false;
                for(int[] arr:mutuallyNearest){
                    if(arr[0]==sampled[0][j] && arr[1]==sampled[1][j]){
                        isFound=true;
                    }
                }
                if(!isFound){
                    throw new Error("did not find sampled");
                }
            }
            this.updateCorrespondenceMatrix(mutuallyNearest);
        }
        return CorrespondenceMatrix;
    }
}
