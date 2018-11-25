/**
 * @author Luca Castelli Aleardi (INF555, 2014)
 */
public abstract class Surface {
    public Point_3[][] points; // control points
    DrawSurface frame; // drawing frame

    public Surface(DrawSurface frame) {
        this.frame = frame;
    }

    public abstract void initialize(int N, int M);

    /**
     * Draw a segment between two points (in the given frame)
     */
    public void drawSegment(Point_3 p, Point_3 q) {
        this.frame.drawSegment(p, q);
    }

    /**
     * Draw the control polygon
     */
    public void drawControlPolygon() {
        //this.frame.stroke(255, 0, 0);
        for (int i = 0; i < points.length; i++) {
            for (int j = 0; j < points[i].length; j++)
                this.frame.drawVertex(points[i][j]);
        }
    }

    /**
     * Evaluate the surface at (u,v)
     * Return point (x(u,v), y(u,v), z(u,v))
     */
    public abstract Point_3 evaluate(double u, double v);

    /**
     * Plot the surface (in the frame), for u, v=0..1, with step dt
     */
    public void plotSurface(double du, double dv) {
        for (double u = 0.; u <= 1. - du; u += du) {
            frame.beginShape(frame.QUAD_STRIP);
            for (double v = 0.; v <= 1.; v += dv) {
                Point_3 p0, p1; // four points defining a quad
                p0 = evaluate(u, v);
                p1 = evaluate(u + du, v);
                frame.noStroke();
                frame.fill(255);
                frame.vertex((float) p0.getX(), (float) p0.getY(), (float) p0.getZ());
                frame.vertex((float) p1.getX(), (float) p1.getY(), (float) p1.getZ());
            }
            frame.endShape();
        }
    }

}
