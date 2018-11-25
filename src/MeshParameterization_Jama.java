import java.util.List;

import Jama.Matrix;

import Jcg.geometry.Point_2;
import Jcg.geometry.Point_3;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Polyhedron_3;
import Jcg.polyhedron.Vertex;

/**
 * Implementation of the planar (2D) Tutte barycentric method, based on the resolution of 2 linear systems
 * 
 * @author Luca Castelli Aleardi, Ecole Polytechnique
 * @version 2017
 */
public class MeshParameterization_Jama extends TutteLayout2D {

	/** vertex numbering in the sub-graph consisting of inner vertices (array of size 'n') */
	public int[] vertexOrder;
	
	/**
	 * Initialize the parameterization
	 */	
	public MeshParameterization_Jama(Polyhedron_3<Point_3> mesh, int face) {
		super(mesh, face); // class the constructor of the ancestor class
		
		vertexOrder=new int[n];
		int counter=0; // counter for inner vertices
		for(int i=0;i<n;i++) {
			if(this.isInside[i]==true) {
				vertexOrder[i]=counter; // store the permutation between vertex indices: from G to G1 (the sub-graph consisting of inner vertices)
				counter++;
			}
		}
	}

	/**
	 * Compute the Tutte drawing of a planar graph, by iteratively
	 * computing barycenters of neighboring vertices <p>
	 * <p>
	 * The algorithm stops when a given convergence condition is reached
	 */	
	public void computeLayout(double tolerance) {
		if(this.nInnerVertices==0) // nothing to compute: there are no inner vertices (all vertices belong to the boundary cycle)
			return;
		
		double[][] rightTerm=this.getRightTerm();
		double[] rightTermX=rightTerm[0];
		double[] rightTermY=rightTerm[1];
		
		double[][] solution=null;
		
		solution=this.solveLinearSystem(rightTermX, rightTermY); // an array containing the two vectors of size 'n' (for the x and y coordinates)
				
		for(Vertex u: this.mesh.vertices) {
			if(this.isInside[u.index]==true) { // set coordinated only for inner vertices
				int index=this.vertexOrder[u.index]; // index of vertex 'u' in the sub-graph: an integer in [0..n-k)
				
				this.points[u.index]=new Point_2(solution[0][index], solution[1][index]);
			}
		}
		
	}

	/**
	 * Solve the linear system (A+D)x=b using Jama library <p>
	 * <p>
	 * Remark: this method is not suitable for a large linear system, since Jama matrices are dense and the solver is not fast
	 * 
	 * @param A  the adjacency matrix (for the sub-graph consisting of inner vertices)
	 * @param D  the diagonal matrix (for the original graph)
	 * @return the solution of the linear system (A+D)x=b
	 */	
	public double[][] solveLinearSystem(double[] x, double y[]) {
		System.out.print("Creating Laplacian (dense) Tutte laplacian matrix of size"+this.nInnerVertices+" (using Jama library)...");
    	Matrix TutteMatrix=this.computeJama();
    	System.out.println("done");

    	System.out.print("Solving linear system with Jama...");
    	long startTime=System.nanoTime(), endTime; // for evaluating time performances
    	/*
    	System.out.println(" ");
    	for(int i=0;i<TutteMatrix.getRowDimension();i++){
    		for(int j=0;j<TutteMatrix.getColumnDimension();j++){
    			System.out.print(TutteMatrix.get(i, j)+" ");
    		}
    		System.out.println(" ");
    	}
    	for(int i=0;i<x.length;i++){
    		System.out.println(x[i]+" "+y[i]);
    	}*/
    	double[] solutionX=this.solve(TutteMatrix, x); // solve the linear system for x-coordinates
    	double[] solutionY=this.solve(TutteMatrix, y); // solve the linear system for y-coordinates
		
    	// evaluate time performances
    	endTime=System.nanoTime();
    	double duration=(double)(endTime-startTime)/1000000000.;
    	System.out.println("done ("+duration+" seconds)");
    	
		return new double[][] {solutionX, solutionY};
	}
	
	/**
	 * Solve linear system Ax=b
	 * 
	 * @param b
	 *            right hand side vector
	 * 
	 * @return the vector solution x[]
	 */
	public double[] solve(Matrix A, double[] b) {
		double[][] B=new double[1][]; // row vector
		B[0]=b; // the transpose gives the column vector
		Jama.Matrix X=A.solve(new Jama.Matrix(B).transpose()); // compute column vector Ax
		return X.transpose().getArray()[0]; // return the transpose of Ax
	}

	
	/**
	 * Computes and returns the Tutte laplacian matrix of the graph (Jama matrix)
	 */	
	public Matrix computeJama() {
		double[][] diag=this.getDenseDiagonalMatrix();
		double[][] adj=this.getDenseAdjacencyMatrix();
		Matrix diagMat=new Matrix(diag);
		Matrix adjMat=new Matrix(adj);
		return diagMat.minus(adjMat);
		//throw new Error("To be completed");	}
	}
	/**
	 * Compute and return the diagonal matrix corresponding to the (n-k) inner vertices
	 * 
	 * @return the (n-k)x(n-k) diagonal matrix storing the degrees of the inner vertices in the original (entire) graph G
	 */	
	public double[][] getDenseDiagonalMatrix() {
		double[][] diag=new double[this.nInnerVertices][this.nInnerVertices];
		for(Vertex v:this.mesh.vertices){
			if(this.isInside[v.index]==true) {
				diag[this.vertexOrder[v.index]][this.vertexOrder[v.index]]=this.mesh.vertexDegree(v);
			}
		}
		return diag;
		//throw new Error("To be completed");
	}

	/**
	 * Compute and return the adjacency matrix corresponding to the (n-k) inner vertices
	 * 
	 * @return the (n-k)x(n-k) adjacency matrix corresponding to the inner vertices
	 */	
	public double[][] getDenseAdjacencyMatrix() {
		double[][] adj=new double[this.nInnerVertices][this.nInnerVertices];
		for(Vertex<Point_3> v:this.mesh.vertices){
			if(this.isInside[v.index]){
				Halfedge<Point_3> e= v.getHalfedge();
				if(this.isInside[e.opposite.vertex.index]){
					adj[this.vertexOrder[v.index]][this.vertexOrder[e.opposite.vertex.index]]=1;
				}
				Halfedge<Point_3> f=e.getNext().getOpposite();
				while(!f.equals(e)){
					f=f.getOpposite();
					if(this.isInside[f.vertex.index]){
						adj[this.vertexOrder[v.index]][this.vertexOrder[f.vertex.index]]=1;
					}
					f=f.getOpposite().getNext().getOpposite();
				}
			}
		}
		return adj;
		//throw new Error("To be completed");
	}

	/**
	 * Compute and return the vector right term given by the locations of the boundary vertices (on the outer cycle)
	 * 
	 * @return the (n-k)x(n-k) adjacency matrix corresponding to the inner vertices
	 */	
	public double[][] getRightTerm() {
		double[][] rightTerm=new double[this.nInnerVertices][this.nInnerVertices];
		for(Vertex<Point_3> v:this.mesh.vertices){
			//for each outlier we add its coordinate in the right term at each of its neighbors index.
			if(this.isInside[v.index]==false) {	
				Halfedge e= v.getHalfedge();
				if(this.isInside[e.opposite.vertex.index]){
					rightTerm[0][this.vertexOrder[e.opposite.vertex.index]]+=this.points[v.index].x;
					rightTerm[1][this.vertexOrder[e.opposite.vertex.index]]+=this.points[v.index].y;
				}
				Halfedge f=e.getNext().getOpposite();
				while(!f.equals(e)){
					if(this.isInside[f.opposite.vertex.index]){
						rightTerm[0][this.vertexOrder[f.opposite.vertex.index]]+=this.points[v.index].x;
						rightTerm[1][this.vertexOrder[f.opposite.vertex.index]]+=this.points[v.index].y;
					}
					f=f.getNext().getOpposite();
				}
			}
		}
		return rightTerm;
	}

}
