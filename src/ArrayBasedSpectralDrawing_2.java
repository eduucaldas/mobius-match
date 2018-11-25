import java.util.*;

import Jcg.geometry.*;
import Jcg.graph.arraybased.*;
import Jcg.graph.arraybased.drawing.*;
import linalg.EigenSolver;
import linalg.JamaMatrices.JamaEigenSolver;
import linalg.JamaMatrices.JamaMatrix;
import linalg.MTJ.MTJSparseEigenSolver;
import linalg.MTJ.MTJSparseMatrix;
import linalg.PColt.PColtEigenSolver;
import linalg.PColt.PColtMatrix;

/**
 * Provides methods for drawing graphs in 2D using spectral based methods
 *
 * @author Luca Castelli Aleardi
 */
public class ArrayBasedSpectralDrawing_2<X extends Point_> extends ArrayBasedGraphDrawing<X>{

	EigenSolver eigensolver; // eigensolver for laplacian matrices
	public int k; // number of eigenvalues to compute

	public ArrayBasedSpectralDrawing_2(ArrayBasedGraph g, int k) {
		this.g=g;
    	this.points=new ArrayList<X>(g.sizeVertices());
    	this.k=k;
    	
    	// use Jama library (not efficient)
    	//JamaMatrix laplacian=LaplacianMatrix.computeJama(g);
    	//this.eigensolver=new JamaEigenSolver(laplacian);

    	// use Parallel Colt library (sparse matrices, iterative eigensolver)
    	PColtMatrix laplacian=LaplacianMatrix.computePColt(g);
    	this.eigensolver=new PColtEigenSolver(laplacian);

    	// use MTJ library (sparse matrices, iterative eigensolver)
    	//MTJSparseMatrix laplacian=LaplacianMatrix.computeMTJ(g);
    	//this.eigensolver=new MTJSparseEigenSolver(laplacian);
}
	
	/**
	 * Return the 2d coordinates of the i-th vertex
	 */
	public Point_2 computeCoordinates_2(int i, double[] eigenvalues, double[][] eigenvectors) {
		double x,y;
		x=eigenvectors[1][i]/Math.sqrt(eigenvalues[1]);
		y=eigenvectors[2][i]/Math.sqrt(eigenvalues[2]);
		return new Point_2(x,y);
	}
	
	public void computeDrawing() {
		this.eigensolver.computeEigenvalueDecomposition(k);
		System.out.println("Laplacian Matrix computed");

		double[] eValues=this.eigensolver.getEigenvalues();
		double[][] eVectors=this.eigensolver.getEigenvectors();
		
		X p;
		for(int i=0;i<this.g.sizeVertices();i++) {
			p=(X)computeCoordinates_2(i, eValues, eVectors);
			this.points.add(p);
		}
	}

}
