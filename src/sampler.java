import java.util.ArrayList;

import Jcg.geometry.Point_3;
import Jcg.geometry.Vector_3;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Vertex;
import MeshManager.SurfaceMesh;

public class sampler {
    /* The sampler class enables us to do the first step, that is sampling from two meshes.
        Starting by sampling local maxima of Gauss Curvature by using a smoothest version of the angle-deficit formula [Desbrun et al. 2002]
        It then uses the Farthest Point Algorithm to take a spread of points, computing geodesic distances with an approximate algorithm based on
        Dijkstra's Shortest path algorithm
     */
    SurfaceMesh M1;
    SurfaceMesh M2;
    int N;
    double epsilon;

    public sampler(SurfaceMesh M1, SurfaceMesh M2, int N, double epsilon) {
        this.M1 = M1;
        this.M2 = M2;
        if (N != -1) {
            this.N = N;
        } else this.N = 100; //the paper average use
        if (epsilon != -1) {
            this.epsilon = epsilon;
        } else this.epsilon = 0.001;
    }

    private double gaussCurvature(Vertex<Point_3> v) {
        /* computes gaussian curvature of a vertex -- could stay unused!!!*/
        double curvature = 2 * Math.PI;
        double area = 0;
        Halfedge<Point_3> e = v.getHalfedge();
        Halfedge<Point_3> f1 = e.opposite;
        Halfedge<Point_3> f2 = e.next.opposite;
        Point_3 p0 = v.getPoint();
        while (!f2.equals(e)) {
            f2 = f2.opposite;
            Point_3 p1 = f1.getVertex().getPoint();
            Point_3 p2 = f2.getVertex().getPoint();
            Vector_3 v2 = new Vector_3(p0, p2);
            Vector_3 v1 = new Vector_3(p0, p1);
            double pdt = (double) v1.innerProduct(v2);
            double angle = Math.acos(pdt / Math.sqrt((double) v1.squaredLength() * (double) v2.squaredLength()));
            curvature += -angle;
            area += pdt;
            f1 = f2;
            f2 = f2.opposite.next.opposite;
        }
        return curvature / area;
    }

    private Vector_3 gaussCurvatureDerivative(Vertex<Point_3> v) {
        Vector_3 angleGrad = new Vector_3(0, 0, 0);
        Halfedge<Point_3> e = v.getHalfedge();
        Halfedge<Point_3> f1 = e.opposite;
        Halfedge<Point_3> f2 = e.next.opposite;
        Point_3 p0 = v.getPoint();
        do {
            f2 = f2.opposite;
            Point_3 p1 = f1.getVertex().getPoint();
            Point_3 p2 = f2.getVertex().getPoint();
            //our first angle : delta ij
            Vector_3 v10 = new Vector_3(p1, p0);
            Vector_3 v12 = new Vector_3(p1, p2);
            double pdt = (double) v12.innerProduct(v10);
            double crosspdt = (double) v12.innerProduct(v10);
            double cotDelta = pdt / crosspdt;
            //our second angle: gamma ij
            Point_3 p3 = f1.getNext().getVertex().getPoint();
            Vector_3 v13 = new Vector_3(p1, p3);
            pdt = (double) v13.innerProduct(v10);
            crosspdt = (double) v13.innerProduct(v10);
            double cotGamma = pdt / crosspdt;
            //we then add the participation from v10 (that is the edge f1 to the gradient).
            angleGrad.sum(v10.multiplyByScalar((cotGamma + cotDelta) * (double) v10.squaredLength()));
            f1 = f2;
            f2 = f2.opposite.next.opposite;
        } while (!f2.equals(e.next.opposite));
        return angleGrad;
    }

    private ArrayList<int[]> findLocalGaussCurvature(double epsilon) {
        /* need to implement [desbrun et al. 2002] algorithm to find local maxima of Gauss Curvature
         * Our methods gaussCurvatureDerivative can give us for each Vertex the derivative of Gauss Curvature.
         * We will send back an index array with all selected points
         * */
        ArrayList<Integer> indexOfLocalMaximaM1 = new ArrayList<Integer>();
        ArrayList<Integer> indexOfLocalMaximaM2 = new ArrayList<Integer>();
        for (Vertex<Point_3> v : this.M1.polyhedron3D.vertices) {
            Vector_3 g = this.gaussCurvatureDerivative(v);
            if ((double) g.squaredLength() < epsilon) {
                indexOfLocalMaximaM1.add(v.index);
            }
        }
        for (Vertex<Point_3> v : this.M2.polyhedron3D.vertices) {
            Vector_3 g = this.gaussCurvatureDerivative(v);
            if ((double) g.squaredLength() < epsilon) {
                indexOfLocalMaximaM2.add(v.index);
            }
        }
        ArrayList<int[]> indx = new ArrayList<int[]>();
        int[] b = new int[indexOfLocalMaximaM1.size()];
        for (int i = 0; i < indexOfLocalMaximaM1.size(); i++) {
            b[i] = indexOfLocalMaximaM1.get(i);
        }
        int[] c = new int[indexOfLocalMaximaM2.size()];
        for (int i = 0; i < indexOfLocalMaximaM2.size(); i++) {
            c[i] = indexOfLocalMaximaM2.get(i);
        }
        indx.add(b);
        indx.add(c);
        return indx;
    }

    private ArrayList<int[]> extendsSampling(ArrayList<int[]> idxOFGaussCurvLocalMax) {
        /* Starting From local maxima of Gauss Curvature, it spreads to capture N points
         * TO DO:
         * 1) Implement FPS algorithm in case we don't have enough maxima of Gauss Curvature
         * 1.1) This algorithm idea is to take the vertex with max distances from  its min distances to our presampled points
         *      Then it add this vertex to the sampled points and reiterates.
         *      The issue is that we are using geodesic distances, so we need to build a graph from scratch every time we want to compute
         *      distances from one point to other points.
         * 2) use this algorithm starting from the initial sample and extend while we haven't N sampled points...
         * */
        if (idxOFGaussCurvLocalMax.get(0).length >= N && idxOFGaussCurvLocalMax.get(1).length >= N) {
            return idxOFGaussCurvLocalMax;
        } else throw new Error("need to be completed");
    }

    public ArrayList<int[]> sample() {
        /* sample points */
        ArrayList<int[]> sampled = this.extendsSampling(this.findLocalGaussCurvature(this.epsilon));
        return sampled;
    }
}
