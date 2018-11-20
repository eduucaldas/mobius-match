import java.util.LinkedList;

public class RationalBezier extends Curve {

    Transformation_2 transformation;

    public RationalBezier(DrawCurve frame, LinkedList<Point_2> p) {
        super(frame, p);
        transformation = null;
    }

    public RationalBezier(DrawCurve frame, LinkedList<Point_2> p,
                          Transformation_2 transformation) {
        super(frame, p);
        this.transformation = transformation;
    }

    public Point_2 evaluate(double t) {
        return rationalDeCasteljau(this.points.length - 1, 0, t, this.points);
    }

    public void plotCurve(double dt) {
        throw new Error("TD INF574: to be completed");
    }

    // apply the transformation to all points (control points and curve)
    public void plotCurveAffine(double dt) {
        throw new Error("TD INF574: to be completed");
    }

    //apply the transformation only to the control points, and then compute the curve
    public void plotControl(double dt) {
        throw new Error("TD INF574: to be completed");
    }

    public void subdivisionRendering(int n) {
    }

    Point_2 rationalDeCasteljau(int r, int i, double t, Point_2[] pnts) {
        throw new Error("TD INF574: to be completed");
    }
}
