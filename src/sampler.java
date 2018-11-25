import java.lang.reflect.Array;
import java.util.ArrayList;
import MeshManager.SurfaceMesh;

public class sampler {
    /* The sampler class enables us to do the first step, that is sampling from two meshes
        Starting by sampling local maxima of Gauss Curvature by using a smoothest version of the angle-deficit formula [Desbrun et al. 2002]

        It then uses the Farthest Point Algorithm to take a spread of points, computing geodesic distances with an approximate algorithm based on
        Dijkstra's Shortest path algorithm
     */
    SurfaceMesh M1;
    SurfaceMesh M2;
    int N;
    public sampler(SurfaceMesh M1, SurfaceMesh M2,int N){
        this.M1=M1;
        this.M2=M2;
        if(N!=-1) {
            this.N = N;
        }
        else this.N=100; //the paper average use
    }

    private ArrayList<Point_3> findLocalGaussCurvature(){
        /* need to implement [desbrun et al. 2002] algorithm to find local maxima of Gauss Curvature*/
        throw new Error("need to be completed");
    }

    private ArrayList<Point_3> extendsSampling(ArrayList<Point_3> gaussCurvLocalMaxima){
        /* Starting From local maxima of Gauss Curvature, it spreads to capture N points */
        throw new Error("need to be completed");
    }
    public ArrayList<Point_3> sample(){
        /* sample points */
        throw new Error("need to be completed");
    }
}
