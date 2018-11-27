package Parametrization;

import Jcg.geometry.Point_2;
import Jcg.geometry.Point_3;
import Jcg.polyhedron.Face;
import Jcg.polyhedron.Polyhedron_3;
import Jcg.polyhedron.Vertex;

/**
 * Planar parameterization of a surface mesh, based on the Tutte barycentric method
 *
 * @author Luca Castelli Aleardi, Ecole Polytechnique
 * @version 2018
 */
public abstract class TutteLayout2D {
    /**
     * 'number of vertices
     */
    public int n;

    /**
     * input mesh to draw (with 'n' vertices)
     */
    public Polyhedron_3<Point_3> mesh;

    /**
     * says whether a vertex in an inner vertex and thus must be drawn (boolean array of size 'n')
     */
    public boolean[] isInside;

    /**
     * says whether a vertex belongs to the exterior boundary
     */
    public boolean[] isOnBoundary;

    /**
     * 2D positions of the vertices of the mesh (to be computed)
     */
    public Point_2[] points;

    /**
     * number of inner vertices of the mesh
     */
    public int nInnerVertices = 0;

    /**
     * number of the boundary vertices of the graph (on the outer cycle)
     */
    public int nBoundaryVertices = 0;

    /**
     * Initialize the drawing
     * <p>
     * Warning: only works is the boundary face is a triangle
     **/
    public TutteLayout2D(Polyhedron_3<Point_3> mesh, int faceIndex) {
        this.mesh = mesh;
        this.n = this.mesh.sizeOfVertices();

        this.isOnBoundary = new boolean[n];
        this.isInside = new boolean[n];
        this.points = new Point_2[n];

        int index = 0;
        for (Vertex v : mesh.vertices) {
            this.isOnBoundary[index] = false;
            index++;
        }

        // set the vertices on the outer face
        Face f = this.mesh.facets.get(faceIndex); // triangle boundary face (the vertices are fixed)
        Vertex v0 = f.getEdge().getVertex();
        Vertex v1 = f.getEdge().getNext().getVertex();
        Vertex v2 = f.getEdge().getNext().getNext().getVertex();
        this.isOnBoundary[v0.index] = true;
        this.isOnBoundary[v1.index] = true;
        this.isOnBoundary[v2.index] = true;
        System.out.println("Outer face: " + v0.index + ", " + v1.index + ", " + v2.index);

        for (Vertex v : mesh.vertices) {
            this.isInside[v.index] = !this.isOnBoundary[v.index];
            this.points[v.index] = new Point_2();
        }

        Point_2[] boundaryLocations = regularPolygonVertices(3, 1.);
        this.points[v0.index] = boundaryLocations[0];
        this.points[v1.index] = boundaryLocations[1];
        this.points[v2.index] = boundaryLocations[2];

        for (int i = 0; i < this.n; i++) {
            if (isInside[i] == true)
                this.nInnerVertices++;
        }

        for (int i = 0; i < this.n; i++) {
            if (isOnBoundary[i] == true)
                this.nBoundaryVertices++;
        }
        System.out.println("Boundary vertices: " + this.nBoundaryVertices);
        System.out.println("Inner vertices: " + this.nInnerVertices);

        for (int i = 0; i < points.length; i++) {
            if (this.points[i] == null && isOnBoundary[i] == true)
                throw new Error("Boundary vertex location not defined for vertex v" + i);
            else if (this.points[i] == null && isInside[i] == true)
                throw new Error("Inner vertex location not defined for vertex v" + i);
        }


    }

    /**
     * A planar parameterization of a surface mesh (using Tutte parameterization)
     * The algorithm computes the parameterization up to a given numeric tolerance
     *
     * @param tolerance the numeric tolerance (defining the halt condition)
     */
    public abstract void computeLayout(double tolerance);

    /**
     * Return the vertices of a regular polygon (equally spaced on a circle of radius r)
     *
     * @param r the radius of the circle
     * @param k the number of points on the outer cycle
     * @return Point_2[] an array of 2D points, storing the vertices of a regular polygon
     */
    public static Point_2[] regularPolygonVertices(int k, double r) {
        Point_2[] vertices = new Point_2[k];
        double x, y;

        for (int i = 0; i < k; i++) {
            x = r * Math.cos((2. * Math.PI / k) * i);
            y = r * Math.sin((2. * Math.PI / k) * i);
            vertices[i] = new Point_2(x, y);
        }
        return vertices;
    }

    /**
     * Project the vertex of the mesh into the plane
     * Each vertex v is assigned the coordinates (x, y, 0.), where 'x' and 'y' are the planar coordinates
     */
    public void projectVertices() {

        for (Vertex v : this.mesh.vertices) {
            //if(this.isInside[v.index]==true) {
            double x = this.points[v.index].getX().doubleValue();
            double y = this.points[v.index].getY().doubleValue();
            v.setPoint(new Point_3(x, y, 0.));

        }
    }

}
