import java.util.LinkedList;

public class Bezier extends Curve {

    Transformation_2 transformation;

    public Bezier(DrawCurve frame, LinkedList<Point_2> p) {
        super(frame, p);
        transformation = null;
    }

    public Bezier(DrawCurve frame, LinkedList<Point_2> p, Transformation_2 transformation) {
        super(frame, p);
        this.transformation = transformation;
    }

    public Bezier(DrawCurve frame, Point_2[] points, Transformation_2 transformation) {
        super(frame, points);
        this.transformation = transformation;
    }

    public Bezier(DrawCurve frame, Point_2[] points) {
        super(frame, points);
    }
    /**
     * Draw the control polygon
     */
    public void drawControlPolygon() {
        this.drawPolygon(this.points);
    }

    public void  drawPolygon(Point_2[] points) {
        this.frame.stroke(0, 0, 255);
        for (int i = 1; i < points.length; i++) {
            drawSegment(points[i], points[i - 1]);
            //drawSegment(transformation.transform(points[i]), transformation.transform(points[i-1]));
        }
        this.frame.stroke(0, 0, 0);
    }

    public void drawTSegment(Point_2 point1, Point_2 point2, Transformation_2 transformation){
        drawSegment(transformation.transform(point1), transformation.transform(point2));
    }

    public Bezier transformedBezier(){
        Point_2[] tpoints = new Point_2[this.points.length];
        for (int i = 0; i < this.points.length; i++) {
            tpoints[i] = this.transformation.transform(this.points[i]);
        }
        return new Bezier(this.frame, tpoints);
    }

    /**
     * Evaluate the curve for parameter t
     * Return point (x(t), y(t))
     */
    public Point_2 evaluate(double t) {
//        return recursiveDeCasteljau(this.points.length - 1, 0, t);
//        return iterativeDeCasteljau(t);
        return bernsteinBezier(t);
    }

    /**
     * Perform the subdivision (once) of the Bezier curve (with parameter t)
     * Return two Bezier curves (with n control points each)
     */
    public Bezier[] subdivide(double t) {
        int n = this.points.length - 1; // degree and number of edges of the control polygon
        Point_2[] b0 = new Point_2[n + 1]; // first control polygon
        Bezier[] result = new Bezier[2]; // the pair of Bezier curves to return as result

        Point_2[] defPoints = this.points.clone();
        double[] coef = {1 - t, t};
        for (int r = 1; r < this.points.length; r++) {
            b0[r - 1] = defPoints[0];
            for (int i = 0; i < this.points.length - r; i++) {
                defPoints[i] = Point_2.linearCombination(
                        new Point_2[]{defPoints[i], defPoints[i+1]},
                        coef);
            }
        }
        b0[n] = defPoints[0];
        result[0] = new Bezier(this.frame, b0);
        result[1] = new Bezier(this.frame, defPoints);

        return result;
    }

    /**
     * Plot the curve (in the frame), for t=0..1, with step dt
     */
    public void plotCurve(double dt) {
        this.drawControlPolygon();
        this.drawControlPoints();
        this.drawBezier(dt);
        // to be completed TD INF574
        Transformation_2 rotation = new Transformation_2(Math.PI/3);
        Transformation_2 scaling = new Transformation_2(0.9, 0.6);
        Transformation_2 translation = new Transformation_2(new Vector_2(0.2, 0.3));
        Bezier trans = new Bezier(this.frame, this.points, rotation.compose(scaling.compose(translation)));
        // draw by appplying the transformation to b(t)
        trans.transformedBezier();
        // transform control points and then draw
        trans = trans.transformedBezier();
        trans.drawControlPolygon();
        trans.drawControlPoints();
        trans.drawBezier(dt);

    }


    public void drawBezier(double dt) {
        int nPoints = (int)Math.floor(1./dt);
        Point_2 previous = this.points[0];
        Point_2 next;
        for (int i = 1; i < nPoints; i++) {
            next = evaluate(i*dt);
            drawSegment(next, previous);
            previous = next;
        }
        next = evaluate(1.);
        drawSegment(next, previous);
    }

    public void drawTBezier(double dt){
        int nPoints = (int)Math.floor(1./dt);
        Point_2 previous = this.points[0];
        Point_2 next;
        for (int i = 1; i < nPoints; i++) {
            next = evaluate(i*dt);
            drawTSegment(next, previous, this.transformation);
            previous = next;
        }
        next = evaluate(1.);
        drawTSegment(next, previous, this.transformation);
    }


    /**
     * Perform the rendering of the curve using subdivision approach
     * Perform the subdivision n times
     */
    public void subdivisionRendering(int n) {
        if (this.points.length < 3) return;
        if (n == 1) {
            this.drawControlPolygon();
            return;
        }
        for (Bezier b: this.subdivide(0.5)
             ) {
            b.subdivisionRendering(n/2);
        }
    }

    public Point_2 recursiveDeCasteljau(int r, int i, double t) {
        if(r == 0) return this.points[i];
        return Point_2.linearCombination(
                new Point_2[]{recursiveDeCasteljau(r-1, i, t), recursiveDeCasteljau(r-1, i+1, t)},
                new double[]{1.-t, t});
    }

    /**
     * Perform the (iterative) De Casteljau algorithm to evaluate b(t)
     */
    public Point_2 iterativeDeCasteljau(double t) {
        Point_2[] defPoints = this.points.clone();
        double[] coef = {1-t, t};
        for (int r = 1; r < this.points.length; r++) {
            for (int i = 0; i < this.points.length - r; i++) {
                defPoints[i] = Point_2.linearCombination(
                        new Point_2[]{defPoints[i], defPoints[i+1]},
                        coef
                );
            }
        }
        return defPoints[0];
    }

    public Point_2 bernsteinBezier(double t) {
        int n = this.points.length - 1;
        return Point_2.linearCombination(this.points, bernsteinCoef(n, t));
    }

    private static double[] bernsteinCoef(int n, double t){
        double s = 1 - t;
        double powt = 1, pows = 1;
        double[] coef = binCoef(n);
        for (int i = 0; i <= n; i++, powt *= t, pows *= s) {
            coef[i] *= powt;
            coef[n-i] *= pows;
        }
        return coef;

    }
    private static double[] binCoef(int n){
        double[] coef = new double[n+1];
        coef[0] = 1.0;
        for (int i = 0; i < n; i++) {
            coef[i+1] = ((n-i)*coef[i])/(i+1);
        }
        return coef;
    }

}
