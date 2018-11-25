import Jcg.geometry.Point_2;
import Jcg.graph.arraybased.*;
import Jcg.graph.arraybased.drawing.*;

/**
 * This class provides input/output methods for testing (planar) drawing algorithms
 */
public class TestSpectralDrawing {

	/**
	 * Test graph drawing algorithms
	 */
    public static void main (String[] args) throws InterruptedException {
		System.out.println("Testing spectral drawing");
		if(args.length!=2) {
			System.out.println("Error: wrong input parameters (two parameters required)");
			System.out.println("\tfirst parameter: input file");
			System.out.println("\tsecond parameter: k (integer value: number to eigenvalues to compute)");
			System.exit(0);
		}
		
    	String filename=args[0];
    	int k=Integer.parseInt(args[1]);
 
    	ArrayBasedGraph g=null; // the input graph
 		g=ArrayBasedGraphLoader.readGraphFromFile(filename); // load graph from file
		ArrayBasedGraphDrawing<Point_2> d=new ArrayBasedSpectralDrawing_2<Point_2>(g, k);

    	d.computeDrawing();
    	System.out.println("planar representation computed");
    	Thread.sleep(1000);
    	d.draw2D();
    		
    	//ArrayBasedGraphLoader.writeMeshSkeletonToFile("Data/tri_round_cube.off", "round_cube.txt");
    }

}
