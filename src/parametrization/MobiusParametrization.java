package parametrization;
import Jcg.geometry.Point_3;
import Jcg.geometry.Vector_3;
import Jcg.mesh.MeshLoader;
import Jcg.polyhedron.Face;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Polyhedron_3;
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
import viewer.MeshViewer;
import Jama.Matrix;

import java.util.ArrayList;
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
    Face cutFace;

    public MobiusParametrization(SurfaceMesh m1,Face cutFace,double precision){
        this.m1=m1;
        this.precision=precision;
        int n=m1.polyhedron3D.vertices.size();
        int nIn = n-2; // matrix size
        Halfedge faceEdge=cutFace.getEdge();
        this.isInner=new boolean[n];
        this.isInitialPoints=new boolean[n];
        for(int i=0;i<n;i++){
            this.isInitialPoints[i]=false;
            this.isInner[i]=true;
        }
        for(int i=0;i<2;i++){
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
            else{
                vertexOrder[i] = -1;
            }
        }
        this.vertexOrder=vertexOrder;
        this.nIn=nIn;
        this.cutFace=cutFace;
    }

    public static double computesBaseAngle(Halfedge f){
        //This function computes the angle surrounding the edge (see Polthier-Pinkall 1993 for graphical images
        //computes facing angle of a vertex.
        Point_3 ui=(Point_3) f.vertex.getPoint();
        Point_3 uj=(Point_3) f.opposite.vertex.getPoint();
        Point_3 u1=(Point_3) f.next.vertex.getPoint();
        Point_3 u2=(Point_3) f.opposite.prev.opposite.vertex.getPoint();
        //computes alpha:
        Vector_3 v1=new Vector_3(u1,ui);
        Vector_3 v2=new Vector_3(u1,uj);
        double pdt = (double) v1.innerProduct(v2);
        double crosspdt =Math.sqrt((double)v1.crossProduct(v2).squaredLength());
        double cotalpha = pdt / crosspdt;
        //computes beta:
        Vector_3 v3=new Vector_3(u2,ui);
        Vector_3 v4=new Vector_3(u2,uj);
        pdt = (double) v3.innerProduct(v4);
        crosspdt =Math.sqrt((double)v3.crossProduct(v4).squaredLength());
        double cotbeta = pdt / crosspdt;
        return cotalpha+cotbeta;
    }

    private DoubleMatrix2D createUMatrix() {
        DoubleMatrix2D U; // (sparse) matrix implementation based on Parallel Colt library
        System.out.print("Creating Laplacian matrix from a graph of size " + nIn+2 + " (using Parallel Colt library)...");
        long startTime = System.nanoTime(), endTime; // for evaluating time performances

        U = new SparseDoubleMatrix2D(nIn, nIn); // create a sparse matrix of size nInxnIn
        this.Sum=new double[nIn];
        for (Vertex v : m1.polyhedron3D.vertices) {
            if(isInner[v.index]){
                //U.set(vertexOrder[v.index],vertexOrder[v.index],1.);
                //First browse of indexes => computes edge bases angles for each neighbor vertex, and compute their sum
                Halfedge e=v.getHalfedge().opposite;
                Halfedge f=e;
                double sum=0;
                do{
                    double facingAngle=MobiusParametrization.computesBaseAngle(f);
                    if(isInner[f.vertex.index])
                        U.set(vertexOrder[f.opposite.vertex.index],vertexOrder[f.vertex.index],(-1)*facingAngle);
                    sum=sum+facingAngle;
                    f=f.opposite.next;
                }while(f!=e);
                //second browse: divide by sum!
                Sum[vertexOrder[f.opposite.vertex.index]]=sum;
                /*do{
                    if(isInner[f.vertex.index]) {
                        double val = U.get(vertexOrder[f.opposite.vertex.index], vertexOrder[f.vertex.index]);
                        U.set(vertexOrder[f.opposite.vertex.index], vertexOrder[f.vertex.index], val / sum);
                    }
                    f=f.opposite.next;
                }while(f!=e);*/
                U.set(vertexOrder[v.index],vertexOrder[v.index],sum);
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
            if (!isInner[v.index]) {
                Halfedge e=v.getHalfedge().opposite;
                Halfedge f=e;
                do{
                    if(isInner[f.vertex.index]) {
                        double facingAngle = MobiusParametrization.computesBaseAngle(f.opposite);
                        if (isInitialPoints[v.index])
                            B[vertexOrder[f.vertex.index]] += facingAngle;// / Sum[vertexOrder[f.vertex.index]];
                        else
                            B[vertexOrder[f.vertex.index]] -= facingAngle;// / Sum[vertexOrder[f.vertex.index]];
                    }
                    f=f.opposite.next;
                }while(f!=e);
            }
        }
        return B;
    }
    private double[] solveForU() {
        /* SOLVE for U using the PColt implementation to have sparse matrix
        * */
        DoubleMatrix2D myMatrix=this.createUMatrix();
        //System.out.println(myMatrix.toString());
        double[] b=this.createsRightHandTerm();
        double[] result=this.solve(myMatrix,b,this.precision);
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
        System.out.print("Found solution: ");
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
    public double[] solve(DoubleMatrix2D A, double[] b) {
        //JamaSolver...
        double[][] a=new double[nIn][nIn];
        for(int i=0;i<nIn;i++){
            for(int j=0;j<nIn;j++)
                a[i][j]=A.get(i,j);
        }
        Matrix matrix=new Matrix(a);
        double[][] B=new double[1][]; // row vector
        B[0]=b; // the transpose gives the column vector
        Jama.Matrix X=matrix.solve(new Jama.Matrix(B).transpose()); // compute column vector Ax
        return X.transpose().getArray()[0]; // return the transpose of Ax
    }
    /*methods for computing the conjugate*/
    private void updateSurround(Hashtable<Halfedge,Double> midEdgeMesh,Halfedge p){
        if(!midEdgeMesh.containsKey(p.prev)){
            double val;
            if(midEdgeMesh.containsKey(p.prev.opposite)){
                val=midEdgeMesh.get(p.prev.opposite);
            }
            else
                val=midEdgeMesh.get(p)-this.conjugatePrevDiff(p.prev);
            midEdgeMesh.put(p.prev,val);

        }
        if(!midEdgeMesh.containsKey(p.next)) {
            double val;
            if(midEdgeMesh.containsKey(p.next.opposite)){
                val=midEdgeMesh.get(p.next.opposite);
            }
            else
                val=this.conjugatePrevDiff(p) + midEdgeMesh.get(p);
            midEdgeMesh.put(p.next,val );
        }
    }
    private double conjugatePrevDiff(Halfedge p){
        /* staying inside the face of the halfedge
         it computes the value of ur*-us* as in the paper figure.
        * */
        Vertex uj=p.vertex;
        Vertex uk=p.prev.vertex;
        Vertex ui=p.next.vertex;
        Point_3 vj=(Point_3) uj.getPoint();
        Point_3 vk=(Point_3) uk.getPoint();
        Point_3 vi=(Point_3) ui.getPoint();
        //computes thetai
        Vector_3 v1=new Vector_3(vi,vk);
        Vector_3 v2=new Vector_3(vi,vj);
        double pdt = (double) v1.innerProduct(v2);
        double crosspdt =Math.sqrt((double)v1.crossProduct(v2).squaredLength());
        double cotthetai = pdt / crosspdt;
        //computes thetak
        v1=new Vector_3(vk,vi);
        v2=new Vector_3(vk,vj);
        pdt = (double) v1.innerProduct(v2);
        crosspdt =Math.sqrt((double)v1.crossProduct(v2).squaredLength());
        double cotthetak = pdt / crosspdt;
        double val=(U[ui.index]-U[uj.index])*cotthetak;
        val+=(U[uk.index]-U[uj.index])*cotthetai;
        val=val/2;
        return val;
    }

    private Hashtable<Halfedge, Double> computeConjugate() {
        /*Instead of seeing a mid-edge mesh as a mesh, we can see it as an hastable between halfedge and their value
        * We need to browse the whole mesh starting from a halfedge.
        **/
        ArrayList<Halfedge> sources=new ArrayList<>();
        Hashtable<Halfedge,Double> midEdgeMesh= new Hashtable<>();
        Halfedge e0=this.m1.polyhedron3D.vertices.get(10).getHalfedge();
        midEdgeMesh.put(e0,0.);
        sources.add(e0);
        while(!sources.isEmpty()){
            Halfedge e=sources.remove(0);
            updateSurround(midEdgeMesh,e);
            if(!midEdgeMesh.containsKey(e.next.opposite)) {
                midEdgeMesh.put(e.next.opposite, midEdgeMesh.get(e.next));
                if(e.next.opposite.face!=cutFace)
                    sources.add(e.next.opposite);
            }
            if(!midEdgeMesh.containsKey(e.prev.opposite)) {
                midEdgeMesh.put(e.prev.opposite, midEdgeMesh.get(e.prev));
                if(e.prev.opposite.face!=cutFace)
                    sources.add(e.prev.opposite);
            }
            if(!midEdgeMesh.containsKey(e.opposite)){
                midEdgeMesh.put(e.opposite,midEdgeMesh.get(e));
                if(e.opposite.face!=cutFace)
                    sources.add(e.opposite);
            }
        }
        return midEdgeMesh;
    }
    public Hashtable<Halfedge,double[]> planarEmbeding(){
        double[] solvedU=this.solveForU();
        /* This creates the final U
        * */
        this.U=new double[this.nIn+2];
        for(Vertex v:m1.polyhedron3D.vertices){
            if(this.isInner[v.index]) {
                U[v.index] = solvedU[this.vertexOrder[v.index]];
            }
            else{
                if(this.isInitialPoints[v.index])
                    U[v.index] =1;
                else
                    U[v.index] =-1;

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
    private Halfedge findClosestMidEdgeVertex(Vertex v){
        Halfedge e=v.getHalfedge();
        Halfedge f=e.opposite;
        Halfedge selectedEdge=f;
        double distMin=-1;
        do{
            double dist=(double)f.getVertex().getPoint().squareDistance(v.getPoint());
            if(distMin==-1 || dist<distMin){
                distMin=dist;
                selectedEdge=f;
            }
            f=f.opposite.next;
        }while(f!=e.opposite);
        return selectedEdge;
    }
    public double[][] getProjectionFromSampled(Vertex[] sample){
        Hashtable<Halfedge,double[]> planarEmbed=this.planarEmbeding();
        double[][] projection=new double[sample.length][2];
        for(int i=0;i<sample.length;i++){
            //we need to find closest mid-edge vertex.
            projection[i]=planarEmbed.get(this.findClosestMidEdgeVertex(sample[i]));
        }
        return projection;
    }
}
