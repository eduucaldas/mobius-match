package meshmanager;

import Jcg.geometry.Point_3;
import Jcg.mesh.MeshLoader;
import Jcg.polyhedron.Face;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Polyhedron_3;
import Jcg.polyhedron.Vertex;
import viewer.MeshViewer;

import java.util.Collection;
import java.util.Hashtable;

/**
 * Class for rendering a surface triangle mesh (using Processing)
 *
 * @author Luca Castelli Aleardi (INF555, 2012)
 */
public class SurfaceMesh {

    double scaleFactor = 60; // scaling factor: useful for 3d rendering
    MeshViewer view; // Processing 3d frame (where meshes are rendered)
    public Polyhedron_3<Point_3> polyhedron3D; // triangle mesh
    public boolean displaySampled;
    public boolean displayFound;
    public Vertex[] correspondence;
    public Vertex[] sampled;
    public Face cutFace;
    public boolean displayCutFace;

    public double getScaleFactor() {
        return scaleFactor;
    }

    public double zoom = 1;
    public Hashtable<Halfedge, double[]> planarEmbedding;
    public boolean drawSurface;

    /**
     * Create a surface mesh from an OFF file
     */
    public SurfaceMesh(MeshViewer view, String filename) {
        this.view = view;
        this.displayFound = false;
        this.displaySampled = false;
        this.displayCutFace = false;
        this.polyhedron3D = MeshLoader.getSurfaceMesh(filename);
        this.drawSurface = true;

        //System.out.println(polyhedron3D.verticesToString());
        //System.out.println(polyhedron3D.facesToString());
        polyhedron3D.isValid(false);

        this.scaleFactor = this.computeScaleFactor();
    }

    /**
     * Draw a segment between two points
     */
    public void drawSegment(Point_3 p, Point_3 q) {
        float s = (float) this.scaleFactor;
        float x1 = (float) p.getX().doubleValue() * s;
        float y1 = (float) p.getY().doubleValue() * s;
        float z1 = (float) p.getZ().doubleValue() * s;
        float x2 = (float) q.getX().doubleValue() * s;
        float y2 = (float) q.getY().doubleValue() * s;
        float z2 = (float) q.getZ().doubleValue() * s;
        this.view.line(x1, y1, z1, x2, y2, z2);
    }

    /**
     * Draw a vertex (as a small sphere)
     */
    public void drawVertex(Point_3 p) {
        float s = (float) this.scaleFactor;
        float x1 = (float) p.getX().doubleValue() * s;
        float y1 = (float) p.getY().doubleValue() * s;
        float z1 = (float) p.getZ().doubleValue() * s;

        view.translate(x1, y1, z1);
        view.sphere(s / 25f);
        view.translate(-x1, -y1, -z1);
    }

    public void drawVertex(Point_3 p, float multiplyScale) {
        float s = (float) this.scaleFactor;
        float x1 = (float) p.getX().doubleValue() * s;
        float y1 = (float) p.getY().doubleValue() * s;
        float z1 = (float) p.getZ().doubleValue() * s;

        view.translate(x1, y1, z1);
        view.sphere(multiplyScale);
        //view.sphere(s * multiplyScale / 25f);
        view.translate(-x1, -y1, -z1);
    }


    /**
     * Draw a triangle
     */
    public void drawTriangle(Point_3 p, Point_3 q, Point_3 r) {
        float s = (float) this.scaleFactor;
        view.vertex((float) (p.getX().doubleValue() * s), (float) (p.getY().doubleValue() * s), (float) (p.getZ().doubleValue() * s));
        view.vertex((float) (q.getX().doubleValue() * s), (float) (q.getY().doubleValue() * s), (float) (q.getZ().doubleValue() * s));
        view.vertex((float) (r.getX().doubleValue() * s), (float) (r.getY().doubleValue() * s), (float) (r.getZ().doubleValue() * s));
    }

    /**
     * Draw a (triangle or polygonal) face
     */
    public void drawFace(Face<Point_3> f) {
        Halfedge<Point_3> h = f.getEdge();
        Halfedge<Point_3> pEdge = h.getNext();

        Point_3 u = h.getOpposite().getVertex().getPoint();
        //view.noStroke();
        view.fill(200, 200, 200, 255); // color of the triangle

        while (pEdge.getVertex() != h.getOpposite().getVertex()) {
            Point_3 v = pEdge.getOpposite().getVertex().getPoint();
            Point_3 w = pEdge.getVertex().getPoint();
            this.drawTriangle(u, v, w); // draw a triangle face

            pEdge = pEdge.getNext();
        }
    }

    /**
     * Draw the entire mesh
     */
    private void drawPlanarEmbedding() {
        for (Halfedge<Point_3> e : planarEmbedding.keySet()) {
            Point_3 p = new Point_3(planarEmbedding.get(e)[0], planarEmbedding.get(e)[1], 0);
            Point_3 q = new Point_3(planarEmbedding.get(e.next)[0], planarEmbedding.get(e.next)[1], 0);
            this.drawSegment(p, q); // draw edge (p,q)
        }
    }

    public void draw(int type) {
        if (this.drawSurface) {
            if (type == 0) {
                view.beginShape(view.TRIANGLES);
                for (Face<Point_3> f : this.polyhedron3D.facets) {
                    this.drawFace(f);
                }
                view.endShape();
            } else {
                // draw all edges
                view.strokeWeight(1); // line width (for edges)
                view.stroke(4);
                for (Halfedge<Point_3> e : this.polyhedron3D.halfedges) {
                    Point_3 p = e.vertex.getPoint();
                    Point_3 q = e.opposite.vertex.getPoint();

                    this.drawSegment(p, q); // draw edge (p,q)
                }
            }
            view.strokeWeight(1);
            if (this.displaySampled) {
                this.displaySampled();
            }
            if (this.displayFound) {
                this.displayFound();
            }
            if (this.displayCutFace) {
                this.displayCutFace();
            }
        } else {
            this.drawPlanarEmbedding();
        }
    }

    private void displayCutFace() {
        Halfedge e = this.cutFace.getEdge();
        Halfedge f = e.next;
        view.stroke(0, 0, 255);
        while (f != e) {
            this.drawVertex((Point_3) f.vertex.getPoint(), 6);
            f = f.next;
        }

    }

    private void displaySampled() {
        view.stroke(255, 0, 0);
        for (Vertex v : sampled) {
            this.drawVertex((Point_3) v.getPoint(), 4);
        }
    }

    private void displayFound() {
        //view.stroke(0,255,0);

        for (int i = 0; i < correspondence.length; i++) {
            Vertex v = correspondence[i];
            if (i % 2 == 0)
                view.stroke(((i + 1) * 31 + 7) % 256 , ((i + 1) * 31 * 31 + 9) % 256, ((i + 1) * 17 + 10) % 256);
            this.drawVertex((Point_3) v.getPoint(), 8);
        }
    }

    /**
     * Draw the X, Y and Z axis
     */
    public void drawAxis() {
        double s = 1;
        Point_3 p000 = new Point_3(0., 0., 0.);
        Point_3 p100 = new Point_3(s, 0., 0.);
        Point_3 p010 = new Point_3(0., s, 0.);
        Point_3 p011 = new Point_3(0., 0., s);

        drawSegment(p000, p100);
        drawSegment(p000, p010);
        drawSegment(p000, p011);
    }


    /**
     * Return the value after truncation
     */
    public static double round(double x, int precision) {
        return ((int) (x * precision) / (double) precision);
    }

    /**
     * Compute the scale factor (depending on the max distance of the point set)
     */
    public double computeScaleFactor() {
        if (this.polyhedron3D == null || this.polyhedron3D.vertices.size() < 1)
            return 1;
        double maxDistance = 0.;
        Point_3 origin = new Point_3(0., 0., 0.);
        for (Vertex<Point_3> v : this.polyhedron3D.vertices) {
            double distance = Math.sqrt(v.getPoint().squareDistance(origin).doubleValue());
            maxDistance = Math.max(maxDistance, distance);
        }
        return Math.sqrt(3) / maxDistance * 150;
    }

    /**
     * Update the scale factor
     */
    public void updateScaleFactor() {
        if (this.drawSurface)
            this.scaleFactor = this.computeScaleFactor() * this.zoom;
        else {
            Collection<double[]> S = planarEmbedding.values();
            double maxDistance = 0.;
            Point_3 origin = new Point_3(0., 0., 0.);
            for (double[] d : S) {
                Vertex v = new Vertex(new Point_3(d[0], d[1], 0));
                double distance = Math.sqrt(v.getPoint().squareDistance(origin).doubleValue());
                maxDistance = Math.max(maxDistance, distance);
            }
            this.scaleFactor = Math.sqrt(3) / maxDistance * 150 * this.zoom;
        }
    }


}
