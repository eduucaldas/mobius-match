import java.util.LinkedList;

/**
 * Abstract class defining methods for implementing Bezier curves
 *
 * @author Luca Castelli Aleardi (INF555, 2014)
 */
public abstract class Curve {

    DrawCurve frame; // drawing frame
    Point_2[] points; // the vertices of the control polygon
    int[] weights;
    LinkedList<Integer> controlWeights = new LinkedList<Integer>();

    public Curve(DrawCurve frame, LinkedList<Point_2> p) {
        this.frame = frame;
        updateInputPoints(p);
    }

    public Curve(DrawCurve frame, Point_2[] points) {
        this.frame = frame;
        this.points = points;
    }

    public void updateInputPoints(LinkedList<Point_2> controlPoints) {
        this.points = new Point_2[controlPoints.size()];
        this.weights = new int[controlPoints.size()];
        while (controlWeights.size() < controlPoints.size())
            controlWeights.add((int) (Math.round(Math.random() * 9) + 1));
        int i = 0;
        for (Point_2 p : controlPoints) {
            this.points[i] = p;
            this.weights[i] = controlWeights.get(i);
            i++;
        }
    }

    /**
     * Draw a segment between two points (in the given frame)
     */
    public void drawSegment(Point_2 p, Point_2 q) {
        this.frame.line((float) p.getX(), (float) p.getY(), (float) q.getX(), (float) q.getY());
    }

    /**
     * Return the string "a[0]+a[1]x+a[2]x^2+...+a[n]x^n"
     * array a[] gives the coefficients of the polynomial expression
     */
    public static String polynomialToString(double[] a) {
        if (a == null || a.length == 0) return "polynome not defined";
        String result = "" + round(a[0], 1000);
        for (int i = 1; i < a.length; i++) {
            String signe = "+";
            if (a[i] < 0) signe = "";
            result = result + " " + signe + round(a[i], 1000) + "x^" + i;
        }
        return result;
    }

    /**
     * Return the value after truncation
     */
    public static double round(double x, int precision) {
        return ((int) (x * precision) / (double) precision);
    }

    /**
     * Draw a circle at given location in the frame
     */
    public void drawPoint(Point_2 p) {
        this.frame.ellipse((float) p.getX(), (float) p.getY(), 5, 5);
    }

    /**
     * Draw the control polygon
     */
    public void drawControlPolygon() {
        this.frame.stroke(0, 0, 255);
        for (int i = 1; i < this.points.length; i++) {
            drawSegment(this.points[i], this.points[i - 1]);
        }
        this.frame.stroke(0, 0, 0);
    }

    /**
     * Draw the control polygon (with a given color)
     */
    public void drawControlPolygon(int r, int g, int b) {
        this.frame.stroke(r, g, b);
        for (int i = 1; i < this.points.length; i++) {
            drawSegment(this.points[i], this.points[i - 1]);
        }
        this.frame.stroke(0, 0, 0);
    }

    /**
     * Draw the control polygon
     */
    public void drawControlPoints() {
        for (int i = 0; i < this.points.length; i++)
            drawPoint(this.points[i]);
    }

    /**
     * Evaluate the curve for parameter t
     * Return point (x(t), y(t))
     */
    public abstract Point_2 evaluate(double t);


    /**
     * Plot the curve (in the frame), for t=0..1, with step dt
     */
    public abstract void plotCurve(double dt);

    /**
     * Perform the rendering of the curve using subdivision approach
     * Perform the subdivision n times
     */
    public abstract void subdivisionRendering(int n);

}
