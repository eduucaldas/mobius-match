/**
 * @author Luca Castelli Aleardi (INF555, 2014)
 * <p>
 * Compute a Bezier Patch surface in 3D (with tensor product approach)
 */
public class BezierPatchTensorProduct extends Surface {

    public BezierPatchTensorProduct(DrawSurface frame) {
        super(frame);
    }

    /**
     * Initialize control polygon
     */
    public void initialize(int N, int M) {
        points = new Point_3[N][N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double x = ((double) i) / N;
                double y = ((double) j) / N;
                double z = (Math.pow(1.8 * (0.2f - x), 3)) * Math.pow(1.8 * (0.5f - y), 3);
                points[i][j] = new Point_3(x, y, z);
                //System.out.print(" "+points[i][j]);
            }
            System.out.println();
        }
    }

    public Point_3 evaluate(double u, double v) {
        Point_3[] bezierInter = new Point_3[this.points.length];

        int n = points[0].length - 1;
        double[] coef = bernsteinCoef(n, u);

        for (int i = 0; i < this.points.length; i++) {
            bezierInter[i] = Point_3.linearCombination(points[i], coef);
        }

        return bernsteinBezier(v, bezierInter);
    }

    public Point_3 bernsteinBezier(double t, Point_3[] points) {
        int n = points.length - 1;
        return Point_3.linearCombination(points, bernsteinCoef(n, t));
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

    /**
     * Perform the (iterative) De Casteljau algorithm (to compute a 3D curve)
     */
    public Point_3 iterativeDeCasteljau(Point_3[] controlPolygon, double t) {
        throw new Error("To be completed: TD3");
    }

    /**
     * Compute the tensor product of two Bezier curves B(u), B(v)
     */
    public Point_3 tensorProduct(int m, int n, double u, double v) {
        throw new Error("To be completed: TD3");
    }

}
