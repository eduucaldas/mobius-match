package viewer;

import Jcg.geometry.Point_3;
import Jcg.polyhedron.Face;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Vertex;
import algo.runAlgo;
import meshmanager.SurfaceMesh;
import processing.core.PApplet;
import utils.MidEdgeFromSurface;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Hashtable;


/**
 * A viewer to debug the parametrization!
 */
public class ParametrizationViewer extends MeshViewer {
    String filename1 = "DATA/shapes-OFF/cow.off";
    String filename2 = "DATA/shapes-OFF/cow.off";
    Hashtable<Halfedge, double[]> h1;
    Hashtable<Halfedge, double[]> h2;
    int drawnMesh;
    boolean drawSurface;
    double scaleFactor;
    double zoom = 1;

    public void setup() {
        size(800, 600, P3D);
        ArcBall arcball = new ArcBall(this);
        m1 = new SurfaceMesh(this, filename1);
        m2 = new SurfaceMesh(this, filename2);
        runAlgo runner = new runAlgo(this);
        runner.initializeDebug();
        h1 = runner.executeDebug(m1);
        h2 = runner.executeDebug(m2);
        drawSurface = true;
        drawnMesh = 0;
        this.updateScaleFactor();
    }

    private void drawSurface(SurfaceMesh m0) {
        ArrayList<Halfedge> midEdge = MidEdgeFromSurface.MidEdgeFromSurface(m0);
        for (Halfedge<Point_3> e : midEdge) {
            Point_3 p = e.vertex.getPoint();
            Point_3 q = e.opposite.vertex.getPoint();
            this.drawSegment(p, q); // draw edge (p,q)
        }
        if (this.drawnMesh == 0) {
            Point_3 p = m1.polyhedron3D.vertices.get(10).getPoint();
            Point_3 p2 = m1.polyhedron3D.vertices.get(0).getPoint();
            this.strokeWeight(10);
            this.stroke(255, 0, 0);
            this.drawSegment(p2, p2);
            this.drawSegment(p, p);
        } else {
            Point_3 p = m2.polyhedron3D.vertices.get(10).getPoint();
            Point_3 p2 = m2.polyhedron3D.vertices.get(0).getPoint();
            this.strokeWeight(10);
            this.stroke(255, 0, 0);
            this.drawSegment(p, p);
            this.drawSegment(p2, p2);
        }
    }

    private void drawPlanarEmbedding(Hashtable<Halfedge, double[]> h0) {
        for (Halfedge<Point_3> e : h0.keySet()) {
            Point_3 p = new Point_3(h0.get(e)[0], h0.get(e)[1], 0);
            Point_3 q = new Point_3(h0.get(e.next)[0], h0.get(e.next)[1], 0);
            this.drawSegment(p, q); // draw edge (p,q)
        }
        Face f;
        Point_3 pIniConj;
        if (this.drawnMesh == 0) {
            f = m1.cutFace;
            double[] r = h0.get(m1.polyhedron3D.vertices.get(10).getHalfedge());
            pIniConj = new Point_3(r[0], r[1], 0);
        } else {
            f = m2.cutFace;
            double[] r = h0.get(m2.polyhedron3D.vertices.get(10).getHalfedge());
            pIniConj = new Point_3(r[0], r[1], 0);
        }
        Halfedge h = f.getEdge();
        this.stroke(255, 0, 0);
        this.strokeWeight(10);
        Point_3 p = new Point_3(h0.get(h)[0], h0.get(h)[1], 0);
        Point_3 q = new Point_3(h0.get(h.next)[0], h0.get(h.next)[1], 0);
        this.drawSegment(p, q);
        this.drawSegment(pIniConj, pIniConj);

    }

    public void draw() {
        background(255);
        translate(width / 2.f, height / 2.f, -1 * height / 2.f);
        this.strokeWeight(1);
        stroke(150, 150, 150);
        this.stroke(4);
        this.updateScaleFactor();
        if (this.drawnMesh == 0) {
            if (this.drawSurface)
                this.drawSurface(m1);
            else
                this.drawPlanarEmbedding(h1);
        } else {
            if (this.drawSurface)
                this.drawSurface(m2);
            else
                this.drawPlanarEmbedding(h2);
        }

    }

    public void drawSegment(Point_3 p, Point_3 q) {
        float s = (float) this.scaleFactor;
        float x1 = (float) p.getX().doubleValue() * s;
        float y1 = (float) p.getY().doubleValue() * s;
        float z1 = (float) p.getZ().doubleValue() * s;
        float x2 = (float) q.getX().doubleValue() * s;
        float y2 = (float) q.getY().doubleValue() * s;
        float z2 = (float) q.getZ().doubleValue() * s;
        this.line(x1, y1, z1, x2, y2, z2);
    }

    public void updateScaleFactor() {
        if (this.drawSurface) {
            if (this.drawnMesh == 1)
                this.scaleFactor = m1.computeScaleFactor();
            else
                this.scaleFactor = m2.computeScaleFactor();
        } else {
            Collection<double[]> S;
            if (this.drawnMesh == 1) {
                S = h1.values();
            } else
                S = h2.values();
            double maxDistance = 0.;
            Point_3 origin = new Point_3(0., 0., 0.);
            for (double[] d : S) {
                Vertex v = new Vertex(new Point_3(d[0], d[1], 0));
                double distance = Math.sqrt(v.getPoint().squareDistance(origin).doubleValue());
                maxDistance = Math.max(maxDistance, distance);
            }
            this.scaleFactor = Math.sqrt(3) / maxDistance * 150;
        }
        this.scaleFactor = this.scaleFactor * zoom;
    }

    public void keyPressed() {
        switch (key) {
            case ('m'):
            case ('M'):
                this.drawnMesh = 1 - this.drawnMesh;
                break;
            case ('e'):
            case ('E'):
                this.drawSurface = !this.drawSurface;
                break;
            case ('i'):
            case ('I'):
                this.zoom = zoom * 2;
                break;
            case ('o'):
            case ('O'):
                this.zoom = zoom / 2;
                break;
            default:
                break;
        }
        if (this.drawnMesh == 0)
            this.m1.updateScaleFactor();
        else
            this.m2.updateScaleFactor();
    }

    /**
     * For running the PApplet as Java application
     */
    public static void main(String[] args) {
        //PApplet pa=new MeshViewer();
        //pa.setSize(400, 400);
        PApplet.main(new String[]{"viewer.ParametrizationViewer"});
    }
}
