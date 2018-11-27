import meshmanager.SurfaceMesh;

public class MidEdgeUniformization {
    /* This class intends to build the mid-edge mesh, then to parametrize it on the complex plane.
     */
    SurfaceMesh originalMesh;
    SurfaceMesh midEdgeMesh;

    public MidEdgeUniformization(SurfaceMesh m) {
        this.originalMesh = m;
        this.computeMidEdgeMesh();
    }

    private void computeMidEdgeMesh() {
        /* create mid edge mesh*/
        throw new Error("to be implemented");
        //this.midEdgeMesh=;
    }

    public double[][] mapToComplexePlane() {
        /* this function maps our mid edge mesh to the complex plane
         *  TO DO:
         *  1) Determines cut face
         *  2) Solve linear system for u given the border constraint resulting from 1)
         *  3) Go through the graph to get u*
         *  4) embeds for each midge edge vertices its coordinates
         * */
        throw new Error("to be implemented");
    }

    private int findCutFace() {
        /* this function find the face having a vertex of minimal geodesic distance measured to all other vertice
         *  this face is cut, then the first vertex X(understand u)value is set to 0
         *  and other vertex values from the face are set to 1
         * */
        throw new Error("to be implemented");
    }

    private double[] solveForU() {
        /* SOLVE FOR U using the PColt implementation to have sparse matrix
         * Will probably need to change some elements in MeshParametrization_PColt to get the right implementation
         *
         * */
        throw new Error("to be implemented");
    }

    private double[] computeConjugate() {
        throw new Error("to be implemented");
    }

}
