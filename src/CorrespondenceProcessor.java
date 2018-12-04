import Jcg.polyhedron.Vertex;
import com.sun.org.apache.xpath.internal.operations.Bool;
import javafx.beans.binding.DoubleExpression;
import meshmanager.SurfaceMesh;

import java.beans.VetoableChangeListener;
import java.util.ArrayList;
import java.util.Hashtable;

public class CorrespondenceProcessor {
    /* this is the class where we compute our fuzzy correspondence Matrix after getting the Mobius Voting*/
    SurfaceMesh m1;
    SurfaceMesh m2;
    double threshold;
    public CorrespondenceProcessor(SurfaceMesh m1,SurfaceMesh m2,double threshold){
        this.m1=m1;
        this.m2=m2;
        this.threshold=threshold;
    }
    public static Hashtable<Vertex, Double> geodesicDistances(SurfaceMesh m1, Vertex p){
        /* we need to find the geodesic distance from p to other vertices.
        * To do so we use Djikstra algorithm, in a similar way as we did in the sampler
        * It computes the distances for the whole mesh.
        * */
        ArrayList<Vertex> sources=new ArrayList<>();
        sources.add(p);
        Hashtable<Vertex,Double> distTable=new Hashtable<>();
        distTable.put(p,0.);
        Sampler.executefps(sources,distTable,m1);
        return distTable;
    }
    private double[] computeFeatureVector(ArrayList<Vertex> selectedSet,Vertex p,Hashtable<Vertex,Double> dist){
        /* Computes feature Vector of p with the selected Set.
         * It is defined as double of length that of selectedSet
         * it ith element is the geodesic distance to vertex p
         */
        double[] featureVector=new double[selectedSet.size()];

        double median=0;
        for(int i=0;i<selectedSet.size();i++){
            featureVector[i]=dist.get(selectedSet.get(i));
            median+=dist.get(selectedSet.get(i));
        }
        for(int i=0;i<selectedSet.size();i++){
            featureVector[i]=featureVector[i]/median;
        }
        return featureVector;
   }
    private double computeL2Dist(double[] a){
        /*L2 distance*/
       double dist=0;
       for(double d:a){
           dist+=d*d;
       }
       return Math.sqrt(dist);
   }
    private double computeFeatureDist(double[] a,double[] b){
        /* L2 feature distance from a-b*/
        double up=0;
        for(int i=0;i<a.length;i++){
            up+=a[i]*a[i]-b[i]*b[i];
       }
       return up/(this.computeL2Dist(a)*this.computeL2Dist(b));
   }
    private int[] findMax(Hashtable<Double,int[]> old,double[][] correspondenceMatrix){
        double max=0;
        int[] coord=new int[2];
        for(int i=0;i<correspondenceMatrix.length;i++){
            for(int j=0;j<correspondenceMatrix.length;j++){
                int[] coordInter={i,j};
                boolean contains=true;
                for(int[] coords:old.values()){
                    if(coordInter[0]==coords[0]){
                        contains=false;
                        break;
                    }
                    else if(coordInter[1]==coords[1]){
                        contains=false;
                        break;
                    }
                }
                if(contains){
                    if(max<correspondenceMatrix[i][j]){
                        max=correspondenceMatrix[i][j];
                        coord=coordInter;
                    }
                }
            }
        }
        return coord;
   }
    public ArrayList<Vertex[]> computeFuzzyCorrespondenceMatrix(double[][] correspondenceMatrix) {
        /* This function first compute the correspondence matrix using previous algorithm
         * Then it computes the FUZZY correspondence matrix, with a minimum threshold of confidence
         * We sort the confidence scores.
         * It then look at low scores, in confidence order and:
         * For each pair (zj,wj):
          * 1) computes features vector of zj => aj
         *  2) Look for closest feature vector of aj => bl, if geodesic dist of this vector to its closest feature vector is smaller than delta (function of surface),
         *  then we do the same for wj and if it is also true for wj, then we add (zj,wj) to our list of final correspondence.
         * Gives back vertexes of final list of correspondences.
         */

        //CREATION OF FIRST CORRESPONDENCE MATRIX:
        Hashtable<Double,int[]> fuzzyCorrespondence=new Hashtable<Double, int[]>();
        for(int i=0;i<correspondenceMatrix.length;i++){
            int[] coord=this.findMax(fuzzyCorrespondence,correspondenceMatrix);
            fuzzyCorrespondence.put(correspondenceMatrix[coord[0]][coord[1]],coord);
        }

        //Find vertex array corresponding to points
        Vertex[][] interestingPoints=new Vertex[fuzzyCorrespondence.keySet().size()][2];
        int boucle=0;
        ArrayList<Vertex[]> realCorrespondence=new ArrayList<Vertex[]>();
        ArrayList<Vertex[]> potentialCorrespondence=new ArrayList<Vertex[]>();
        for(double max:fuzzyCorrespondence.keySet()){
            int[] coords=fuzzyCorrespondence.get(max);
            Vertex[] v=new Vertex[2];
            v[0]=this.m1.polyhedron3D.vertices.get(coords[0]);
            v[1]=this.m2.polyhedron3D.vertices.get(coords[1]);
            interestingPoints[boucle][0]=v[0];
            interestingPoints[boucle][1]=v[1];
            if(max>this.threshold){
                realCorrespondence.add(v);
            }
            else{
                potentialCorrespondence.add(v);
            }
            boucle++;
        }
        // INCREASE OF THIS CORRESPONDENCE MATRIX PRECISION
        ArrayList<Vertex> set1=new ArrayList<Vertex>();
        ArrayList<Vertex> set2=new ArrayList<Vertex>();
        for( Vertex[] v2: potentialCorrespondence){
            set1.add(v2[0]);
            set2.add(v2[1]);
        }
        Hashtable<Vertex,double[]> a1=new Hashtable<>();
        Hashtable<Vertex,double[]> a2=new Hashtable<>();

        //we compute feature vectors
        Hashtable<Vertex,Double> dist=new Hashtable<>();
        for(Vertex[] v:potentialCorrespondence){
            Vertex p=v[0];
            dist=this.geodesicDistances(this.m1,p);
            a1.put(p,this.computeFeatureVector(set1,p,dist));

            p=v[1];
            dist=this.geodesicDistances(this.m2,p);
            a2.put(p,this.computeFeatureVector(set2,p,dist));
        }
        //computes gamma:
        double maxRadius=0;
        for(double d:dist.values()) {
            if(maxRadius<d) {
                maxRadius = d;
            }
        }
        double gamma=0.05*maxRadius*2*Math.PI;
        //finds good correspondence pairs in potential correspondence pairs
        for(Vertex[] v:potentialCorrespondence){
            double min=-1;
            for(Vertex v1:set1){
                double distance=this.computeFeatureDist(a1.get(v[0]),a1.get(v1));
                if((min == -1) || (distance < min)){
                    min=distance;
                }
            }
            if(min<gamma) {
                min=-1;
                for (Vertex v2 : set2) {
                    double distance = this.computeFeatureDist(a1.get(v[0]), a1.get(v2));
                    if ((min == -1) || (distance < min)) {
                        min = distance;
                    }
                }
                if(min<gamma){
                    realCorrespondence.add(v);
                }
            }
        }
        return realCorrespondence;
    }
}
