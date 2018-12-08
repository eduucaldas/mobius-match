package parametrization;
import Jcg.geometry.Point_3;
import Jcg.geometry.Vector_3;
import Jcg.polyhedron.Face;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Vertex;
import cern.colt.matrix.Norm;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.solver.DefaultDoubleIterationMonitor;
import cern.colt.matrix.tdouble.algo.solver.DoubleCG;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.algo.solver.preconditioner.DoubleDiagonal;
import cern.colt.matrix.tdouble.algo.solver.preconditioner.DoublePreconditioner;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;
import meshmanager.SurfaceMesh;

import java.util.Hashtable;


public class MobiusParametrization {
    /* Within this class we implemented the parametrization algorithm.
    *
    *
    * */
    int[] vertexOrder;
    boolean[] isInner;
    SurfaceMesh m1;
    double[] Sum;
    int nIn;
    double precision;
    double[] U;
    boolean[] isInitialPoints;

    public MobiusParametrization(SurfaceMesh m1,Face cutFace,double precision){
        this.m1=m1;
        this.precision=precision;
        int n=m1.polyhedron3D.vertices.size();
        int nIn = n-3; // matrix size
        Halfedge faceEdge=cutFace.getEdge();
        this.isInner=new boolean[n];
        this.isInitialPoints=new boolean[n];
        for(int i=0;i<n;i++){
            this.isInitialPoints[i]=false;
            this.isInner[i]=false;
        }
        for(int i=0;i<3;i++){
            if(i==0) {
                this.isInitialPoints[faceEdge.vertex.index]=true;
            }
            isInner[faceEdge.vertex.index]=false;
            faceEdge=faceEdge.next;
        }
        int[] vertexOrder = new int[n];
        int counter = 0; // counter for inner vertices
        for (int i = 0; i < n; i++) {
            if (isInner[i] == true) {
                vertexOrder[i] = counter; // store the permutation between vertex indices: from G to G1 (the sub-graph consisting of inner vertices)
                counter++;
            }
        }
        this.vertexOrder=vertexOrder;
        this.isInner=isInner;
        this.nIn=nIn;
    }

    public static double computesFacingAngle(Halfedge f){
        //computes facing angle of a vertex.
        Point_3 uj=(Point_3) f.vertex.getPoint();
        Point_3 ui=(Point_3) f.opposite.vertex.getPoint();
        Point_3 u1=(Point_3) f.next.vertex.getPoint();
        Point_3 u2=(Point_3) f.opposite.next.vertex.getPoint();
        //computes alpha:
        Vector_3 v1=new Vector_3(ui,u1);
        Vector_3 v2=new Vector_3(ui,u2);
        double pdt = (double) v1.innerProduct(v2);
        double crosspdt =Math.sqrt((double)v1.crossProduct(v2).squaredLength());
        double cotalpha = pdt / crosspdt;
        //computes beta:
        v1=new Vector_3(uj,u1);
        v2=new Vector_3(uj,u2);
        pdt = (double) v1.innerProduct(v2);
        crosspdt =Math.sqrt((double)v1.crossProduct(v2).squaredLength());
        double cotbeta = pdt / crosspdt;
        return cotalpha+cotbeta;
    }

    private DoubleMatrix2D createUMatrix() {
        DoubleMatrix2D U; // (sparse) matrix implementation based on Parallel Colt library
        System.out.print("Creating Laplacian matrix from a graph of size " + nIn+3 + " (using Parallel Colt library)...");
        long startTime = System.nanoTime(), endTime; // for evaluating time performances

        U = new SparseDoubleMatrix2D(nIn, nIn); // create a sparse matrix of size nInxnIn
        for (Vertex v : m1.polyhedron3D.vertices) {
            if(isInner[v.index]){
                U.set(vertexOrder[v.index],vertexOrder[v.index],1);
                //First browse of indexes => computes edge bases angles for each neighbor vertex, and compute their sum
                Halfedge e=v.getHalfedge().opposite;
                Halfedge f=e;
                double sum=0;
                do{
                    double facingAngle=MobiusParametrization.computesFacingAngle(f);
                    U.set(vertexOrder[f.opposite.vertex.index],vertexOrder[f.vertex.index],facingAngle);
                    sum=sum+facingAngle;
                    f=f.opposite.next;
                }while(f!=e);
                //second browse: divide by sum!
                Sum[vertexOrder[f.opposite.vertex.index]]=sum;
                do{
                    double val=U.get(vertexOrder[f.opposite.vertex.index],vertexOrder[f.vertex.index]);
                    U.set(vertexOrder[f.opposite.vertex.index],vertexOrder[f.vertex.index],(-1)*val/sum);
                    f=f.opposite.next;
                }while(f!=e);
            }
        }
        endTime = System.nanoTime();
        double duration = (double) (endTime - startTime) / 1000000000.;
        System.out.println("done (" + duration + " seconds)");
        return U;
    }
    private double[] createsRightHandTerm(){
        //creates righ hand term
        double[] B=new double[nIn];
        for (Vertex v : m1.polyhedron3D.vertices) {
            if (isInner[v.index] && !isInitialPoints[v.index]) {
                Halfedge e=v.getHalfedge().opposite;
                Halfedge f=e;
                do{
                    double facingAngle=MobiusParametrization.computesFacingAngle(f);
                    B[vertexOrder[f.vertex.index]]+=facingAngle/Sum[vertexOrder[f.opposite.vertex.index]];
                }while(f!=e);
            }
        }
        return B;
    }
    private double[] solveForU() {
        /* SOLVE for U using the PColt implementation to have sparse matrix
        * */
        DoubleMatrix2D U=this.createUMatrix();
        double[] b=this.createsRightHandTerm();
        double[] result=this.solve(U,b,precision);
        return result;
    }
    private double[] solve(DoubleMatrix2D A, double[] b, double precision) {
        int size = A.columns(); // matrix size
        // set the linear solver
        System.out.print("Setting solver and preconditioner (diagonal)...");
        double[] start = new double[size]; // initial guess for the CG iterator
        DoubleCG itSolver; // iterative solver implemented in Parallel Colt (Conjugate Gradient)
        DefaultDoubleIterationMonitor m;
        DoublePreconditioner preconditioner; // preconditioner
        preconditioner = new DoubleDiagonal(size); // use the diagonal matrix as preconditioner
        preconditioner.setMatrix(A);
        itSolver = new DoubleCG(new DenseDoubleMatrix1D(size));

        itSolver.setPreconditioner(preconditioner);

        m = (DefaultDoubleIterationMonitor) itSolver.getIterationMonitor();
        m.setMaxIterations(size);
        m.setNormType(Norm.Two); // choose euclidean norm
        m.setRelativeTolerance(precision);
        System.out.println("done");

        // running the iterative linear solver
        double[] x = new double[size]; // solution (and initial guess)
        DoubleMatrix1D X = new DenseDoubleMatrix1D(size); // solution
        DoubleMatrix1D B = new DenseDoubleMatrix1D(size);

        for (int i = 0; i < size; i++) {
            B.set(i, b[i]);
            X.set(i, start[i]); // set initial guess
        }
        try {
            itSolver.solve(A, B, X); // compute solution of the linear system
        } catch (IterativeSolverDoubleNotConvergedException e) {
            e.printStackTrace();
        }
        for (int i = 0; i < size; i++) {
            x[i] = X.get(i);
        }

        int verbosity = 1;
        if (verbosity > 0) {
            int iter = m.iterations();
            double residual = m.residual();
            System.err.println("  Conjugate Gradient solved in " + iter + " turns, dist = " + residual);
        }

        return x;
    }

    /*methods for computing the conjugate*/
    private double conjugatePrevDiff(Halfedge p){
        /* staying inside the face of the halfedge
         it computes the value of up*-ub*, the prev halfedge
        * */
        Vertex uj=p.vertex;
        Vertex ui=p.prev.vertex;
        Vertex uk=p.next.vertex;
        Point_3 vj=(Point_3) uj.getPoint();
        Point_3 vk=(Point_3) ui.getPoint();
        Point_3 vi=(Point_3) uk.getPoint();
        //computes thetai
        Vector_3 v1=new Vector_3(vi,vj);
        Vector_3 v2=new Vector_3(vi,vk);
        double pdt = (double) v1.innerProduct(v2);
        double crosspdt =Math.sqrt((double)v1.crossProduct(v2).squaredLength());
        double cotthetai = pdt / crosspdt;
        //computes thetaj
        v1=new Vector_3(vj,vi);
        v2=new Vector_3(vj,vk);
        pdt = (double) v1.innerProduct(v2);
        crosspdt =Math.sqrt((double)v1.crossProduct(v2).squaredLength());
        double cotthetaj = pdt / crosspdt;
        double val=(U[uk.index]-U[ui.index])*cotthetaj;
        val+=(U[uk.index]-U[uj.index])*cotthetai;
        val=val/2;
        return val;
    }
    private void updateSurround(Hashtable<Halfedge,Double> midEdgeMesh,Halfedge p){
        if(!midEdgeMesh.containsKey(p.opposite)){
            midEdgeMesh.put(p.opposite,midEdgeMesh.get(p));
            this.updateSurround(midEdgeMesh,p.opposite);
        }
        if(!midEdgeMesh.containsKey(p.prev)){
            midEdgeMesh.put(p.prev,midEdgeMesh.get(p)-this.conjugatePrevDiff(p));
            this.updateSurround(midEdgeMesh,p.prev);
        }
        if(!midEdgeMesh.containsKey(p.next)) {
            midEdgeMesh.put(p.next, this.conjugatePrevDiff(p.next) + midEdgeMesh.get(p));
            this.updateSurround(midEdgeMesh, p.next);
        }
    }
    private Hashtable<Halfedge, Double> computeConjugate() {
        /*Instead of seeing a mid-edge mesh as a mesh, we can see it as an hastable between halfedge and their value*/
        Hashtable<Halfedge,Double> midEdgeMesh= new Hashtable<>();
        Halfedge e0=m1.polyhedron3D.vertices.get(0).getHalfedge();
        midEdgeMesh.put(e0,0.);
        this.updateSurround(midEdgeMesh,e0);
        return midEdgeMesh;
    }

    private Hashtable<Halfedge,double[]> planarEmbeding(){
        double[] solvedU=this.solveForU();
        /* This creates the final U
        *
        * */
        this.U=new double[this.nIn+3];
        for(Vertex v:m1.polyhedron3D.vertices){
            if(this.isInner[v.index]) {
                U[v.index] = solvedU[this.vertexOrder[v.index]];
            }
            else{
                if(this.isInitialPoints[v.index])
                    U[v.index] =0;
                else
                    U[v.index] =1;
            }

        }
        Hashtable<Halfedge,Double> midEdge=this.computeConjugate();
        Hashtable<Halfedge,double[]> result=new Hashtable<>();
        for(Halfedge e:m1.polyhedron3D.halfedges){
            double[] val={(U[e.vertex.index]+U[e.opposite.vertex.index])/2,midEdge.get(e)};
            result.put(e,val);
        }
        return result;
    }

    public double[][] getProjectionFromSampled(Vertex[] sample){
        Hashtable<Halfedge,double[]> planarEmbed=this.planarEmbeding();
        double[][] projection=new double[sample.length][2];
        for(int i=0;i<sample.length;i++){
            projection[i]=planarEmbed.get(sample[i]);
        }
        return projection;
    }
}
