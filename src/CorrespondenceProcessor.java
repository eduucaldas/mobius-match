import Jcg.polyhedron.Vertex;
import meshmanager.SurfaceMesh;

public class CorrespondenceProcessor {
    /* this is the class where we compute our fuzzy correspondence Matrix after getting the Mobius Voting*/

    public CorrespondenceProcessor(SurfaceMesh m1,SurfaceMesh m2){
        throw new Error("to be implemented");
    }

    private double[] computeFeatureVector(Vertex[] selectedSet,Vertex p){
        /* Computes feature Vector of p with the selected Set.
         * It is defined as double of length that of selectedSet
         * it ith element is the geodesic distance to vertex p
         */
       throw new Error("to be implemented");
   }
   private double computeL2Dist(double[] a){
        /*L2 distance*/
       throw new Error("to be implemented");
   }
   private double computeFeatureDist(double[] a,double[] b){
        /* L2 feature distance from a-b*/
       throw new Error("to be implemented");
   }

    public Vertex[][] computeFuzzyCorrespondenceMatrix(double[][] correspondenceMatrix) {
        /* This function first compute the correspondence matrix using previous algorithm
         * Then it computes the FUZZY correspondence matrix, with a minimum threshold of confidence
         * We sort the confidence scores.
         * It then look at low scores, in confidence order and:
         * For each pair (zj,wj):
          * 1) computes features vector of zj => aj [see def of this vector in **** func]
         *  2) Look for closest feature vector of aj => bl, if geodesic dist of this vector to its closest feature vector is smaller than delta (function of surface),
         *  then we do the same for wj and if it is also true for wj, then we add (zj,wj) to our list of final correspondence.
         * Gives back vertexes of final list of correspondences.
         */



        throw new Error("to be implemented");
    }


}
