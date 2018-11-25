import processing.core.PVector;

/**
 * This example illustrates the use of class Surface
 * We render the hyperbolic paraboloid y=uz, for u,v:[0,1]x[0,1]
 * We apply a scaling in order to obtain a correct view
 *
 * @author Luca Castelli Aleardi (INF555, 2012)
 */
public class Paraboloid extends Surface {

    public Paraboloid(DrawSurface frame) {
        super(frame);
    }

    /**
     * This method is not necessary
     */
    public void initialize(int N, int M) {
        int n = N;
        points = new Point_3[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double x = ((double) i) / n;
                double z = ((double) j) / n;
                double y = x * z;
                points[i][j] = new Point_3(x, y, z);
                System.out.print(" " + points[i][j]);
            }
            System.out.println();
        }
    }

    public Point_3 evaluate(double u, double v) {
        return new Point_3(u, (u * v), v);
    }

}
