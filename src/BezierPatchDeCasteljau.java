/**
 * @author Luca Castelli Aleardi (INF555, 2014)
 * <p>
 * Compute a Bezier Patch surface in 3D (iterative De Casteljau)
 */
public class BezierPatchDeCasteljau extends Surface {

    public BezierPatchDeCasteljau(DrawSurface frame) {
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
        throw new Error("To be completed: TD3");
    }


}
