package utils;

import Jcg.geometry.Point_3;
import Jcg.geometry.Vector_3;
import Jcg.polyhedron.Halfedge;

public class GeometryUtils {
    public static Point_3 mid(Point_3 p, Point_3 q){
        return Point_3.linearCombination(
                new Point_3[]{p, q},
                new Number[]{0.5, 0.5});
    }

    public static Vector_3 geoEdge(Halfedge<Point_3> e){
        return new Vector_3(e.opposite.getVertex().getPoint(), e.getVertex().getPoint());
    }

    public static double normEdge(Halfedge<Point_3> e){
        return Math.sqrt(geoEdge(e).squaredLength().doubleValue());
    }
}
