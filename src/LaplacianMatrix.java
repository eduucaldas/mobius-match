import java.util.Collection;

import Jcg.graph.arraybased.ArrayBasedGraph;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;
//import linalg.JamaMatrices.JamaMatrix;
import Jama.Matrix;
//import mtj.*;
//import linalg.PColt.PColtMatrix;
//import parallelcolt.*;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

public class LaplacianMatrix {
	/**
	 * Computes and returns the laplacian matrix of a graph (Jama matrix)
	 */	
	public static Matrix computeJama(ArrayBasedGraph g) {
		int n=g.sizeVertices();
		
		double[][] m=new double[n][n];
		for(int i=0;i<n;i++) {
			for(int j=0;j<n;j++) {
				if(i==j) m[i][j]=g.degree(i);
				else if(g.adjacent(i,j)==true)
					m[i][j]=-1.;
				else m[i][j]=0.;
			}
	   	}
	   	return new Matrix(m);
	}

	/**
	 * Create the Laplacian matrix of a graph (with n vertices)
	 * (not efficient implementation: it takes nxn time)
	 * 
	 * @param g
	 * 			the input graph
	 */	
	public static MTJSparseMatrix computeMTJ(ArrayBasedGraph g) {
		CompRowMatrix A; // the result
		int n=g.sizeVertices(); // matrix size
		System.out.print("Creating Laplacian (sparse) matrix from a graph of size "+n+" (using MTJ library)...");
		
		int[][] m=new int[n][]; // indices of non-zero entries for each row
		for(int i=0;i<n;i++) { // iterate over all vertices
			Collection<Integer> neighbors=g.neighbors(i); // neighbors of vertex i
			m[i]=new int[neighbors.size()];
			int count=0;
			for(Integer j: neighbors) {
				//System.out.print(" "+j);
				m[i][count]=j;
				count++;
			}
			//System.out.println();
		}
		
		A=new CompRowMatrix(n, n, m); // create a sparse matrix of size nxn

		for(int i=0;i<n;i++) { // iterate over all vertices
			Collection<Integer> neighbors=g.neighbors(i); // neighbors of vertex i
			for(Integer j: neighbors) {
				//System.out.print(" "+j);
				if(i==j) A.set(i, j, neighbors.size());
				else A.set(i, j, -1.);
			}
			//System.out.println();
		}

		System.out.println("done");
		return new MTJSparseMatrix(A);
	}

	/**
	 * Create the Laplacian matrix of a graph (with n vertices)
	 * (not efficient implementation: it takes nxn time)
	 * 
	 * @param g
	 * 			the input graph
	 */	
	public static PColtMatrix computePColt(ArrayBasedGraph g) {
		DoubleMatrix2D A; // (sparse) matrix implementation based on Parallel Colt library
		int n=g.sizeVertices(); // matrix size
		System.out.print("Creating Laplacian matrix from a graph of size "+n+" (using Parallel Colt library)...");
		
		A=new SparseDoubleMatrix2D(n, n); // create a sparse matrix of size nxn

		for(int i=0;i<g.sizeVertices();i++) {
			A.setQuick(i, i, g.degree(i)); // set diagonal entries
			Collection<Integer> neighbors=g.neighbors(i); // neighbors of vertex i
			for(Integer j: neighbors) {
				A.setQuick(i, j, -1.);
				A.setQuick(j, i, -1.);
			}
	   	}
		System.out.println("done");
		
		return new PColtMatrix(A);
	}

}
