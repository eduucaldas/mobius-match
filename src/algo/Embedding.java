package algo;

import Jcg.geometry.Point_3;
import Jcg.mesh.MeshLoader;
import Jcg.polyhedron.Face;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Polyhedron_3;
import Jcg.polyhedron.Vertex;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.SparseDoubleAlgebra;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;
import utils.GeometryUtils;
import utils.MeshUtils;
import utils.PolyGraph;

public class Embedding {
//    private HashMap<Vertex<Point_3>, double[]> fi;
    private final Polyhedron_3<Point_3> p;
    private DoubleMatrix1D u;
//    private HashMap<Vertex<Point_3>, Double> uStar;
    private Face<Point_3> cutFace;
    private static final double[] uCut = new double[]{-1, 1};


    private double computeCoef(Halfedge<Point_3> e){
        return GeometryUtils.cotSupportAngle(e) + GeometryUtils.cotSupportAngle(e.opposite);
    }

    private DoubleMatrix1D computeU() {
        DoubleMatrix2D A = new SparseRCDoubleMatrix2D(p.vertices.size(), p.vertices.size());
        DoubleMatrix1D b = new SparseDoubleMatrix1D(p.vertices.size());
        int idVCut0 = cutFace.getVertexIndices(p)[0];
        int idVCut1 = cutFace.getVertexIndices(p)[1];
        // setting up linear system
        for (Vertex<Point_3> v : p.vertices) {
            if(v.index == idVCut0){
                b.setQuick(v.index, uCut[0]);
                A.setQuick(v.index, v.index, 1);
            }
            else if(v.index == idVCut1){
                b.setQuick(v.index, uCut[1]);
                A.setQuick(v.index, v.index, 1);
            }
            else{
                double accCoef = 0;
                for (Halfedge<Point_3> e : MeshUtils.findOutEdges(v)) {
                    double coef = computeCoef(e);
                    accCoef += coef;
                    A.setQuick(v.index, e.getVertex().index, -coef);
                }
                A.setQuick(v.index, v.index, accCoef);
            }
        }


        // Solving Linear System
        SparseDoubleAlgebra alg = new SparseDoubleAlgebra(0.00001);
        return alg.solve(A, b);
    }

//    private HashMap<Vertex<Point_3>, Double> computeUStar() {
//
//    }
//
//    private HashMap<Vertex<Point_3>, double[]> computeFi(
//            Polyhedron_3<Point_3> p, Face<Point_3> cutFace,
//            HashMap<Vertex<Point_3>, Double> u, HashMap<Vertex<Point_3>, Double> uStar) {
//
//    }

    public Embedding(Polyhedron_3<Point_3> p){
        this.p = p;
        this.cutFace = PolyGraph.findCutFace(p);

        this.u = computeU();
//        this.uStar = computeUStar();
    }

//
//    public double[] map(Vertex<Point_3> v) {
//
//    }

    // Unit Testing
    public static void main(String[] args) {
        Polyhedron_3<Point_3> tetra = MeshLoader.getSurfaceMesh("DATA/shapes-OFF/tetrahedron.off");
        System.out.println(tetra.edgesToString());
        System.out.println(tetra.verticesToString());

        // index vertex
        System.out.println("Unit Test geoEdge");
        for (Vertex<Point_3> v:tetra.vertices) {
            System.out.println(v.index + "," + v);
            for (Vertex<Point_3> w: MeshUtils.findNeighVertices(v)) {

                System.out.println("\t" + w.index + "," + w);
            }
        }
        Polyhedron_3<Point_3> sphere = MeshLoader.getSurfaceMesh("DATA/shapes-OFF/sphere.off");
        Polyhedron_3<Point_3> cow = MeshLoader.getSurfaceMesh("DATA/shapes-OFF/cow.off");
        // ComputeU
        System.out.println("Unit Test ComputeU");
        Embedding fi = new Embedding(sphere);
        System.out.println(fi.u);



    }


}
