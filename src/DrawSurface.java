import processing.core.*;

/**
 * @author Luca Castelli Aleardi (INF555, 2014)
 */
public class DrawSurface extends PApplet {

    Surface s = null;
    ArcBall arcball;
    int patchType = 0;
    int n = 7;

    public void setup() {
        size(600, 400, P3D);
        this.arcball = new ArcBall(this);

        initApplet();
    }

    public void initApplet() {
        if (this.patchType == 0) this.s = new Paraboloid(this);
        else if (this.patchType == 1) this.s = new BezierPatchTensorProduct(this);
        else this.s = new BezierPatchDeCasteljau(this);


        this.s.initialize(n, n);
    }

    public void draw() {
        background(0);
        this.lights();

        translate(width / 2.f, height / 2.f, -1 * height / 2.f);
        //this.strokeWeight(1);
        //stroke(255,0,0);

        if (s != null) {
            this.s.drawControlPolygon();
            this.s.plotSurface(0.02, 0.02);
        }
    }

    public void drawSegment(Point_3 p, Point_3 q) {
        line((float) p.getX(), (float) p.getY(), (float) p.getZ(), (float) q.getX(), (float) q.getY(), (float) q.getZ());
    }

    /**
     * Draw a vertex (as a small sphere)
     */
    public void drawVertex(Point_3 p) {
        float x1 = (float) p.getX();
        float y1 = (float) p.getY();
        float z1 = (float) p.getZ();

        this.translate(x1, y1, z1);
        this.sphere(0.02f);
        this.translate(-x1, -y1, -z1);
    }

    public void keyPressed() {
        switch (key) {
            case ('i'):
            case ('I'): {
                zoomIn();
            }
            break;
            case ('o'):
            case ('O'): {
                zoomOut();
            }
            break;
            case ('s'):
            case ('S'): {
                this.patchType = (this.patchType + 1) % 3;
                initApplet();
            }
            break;
        }
    }

    public void zoomIn() {
        System.out.println("zoom in");
        this.arcball.scaleFactor = this.arcball.scaleFactor * 1.1f;
    }

    public void zoomOut() {
        System.out.println("zoom out");
        this.arcball.scaleFactor = this.arcball.scaleFactor * 0.9f;
    }

}
