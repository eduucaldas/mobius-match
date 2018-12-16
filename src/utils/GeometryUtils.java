package utils;

import Jcg.geometry.Point_3;
import Jcg.geometry.Vector_3;
import Jcg.mesh.MeshLoader;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Polyhedron_3;

public class GeometryUtils {
    public static Point_3 mid(Point_3 p, Point_3 q){
        return Point_3.linearCombination(
                new Point_3[]{p, q},
                new Number[]{0.5, 0.5});
    }

    public static Vector_3 geoEdge(Halfedge<Point_3> e){
        return new Vector_3(e.opposite.getVertex().getPoint(), e.getVertex().getPoint());
    }

    public static double norm(Vector_3 v) {
        return Math.sqrt(v.squaredLength().doubleValue());
    }

    public static double normEdge(Halfedge<Point_3> e){
        return norm(geoEdge(e));
    }

    public static Vector_3 project(Vector_3 v, Vector_3 intoProject){
        return intoProject.multiplyByScalar(
                        v.innerProduct(intoProject).doubleValue()/intoProject.squaredLength().doubleValue());
    }

    public static double cotSupportAngle(Halfedge<Point_3> e){
        Vector_3 vHyp = GeometryUtils.geoEdge(e.next.opposite);
        Vector_3 vToProject = GeometryUtils.geoEdge(e.next.next);
        Vector_3 vCathAdj = project(vHyp, vToProject);

        return norm(vCathAdj)/norm(vCathAdj.difference(vHyp));
    }


    // Unit Testing
    public static void main(String[] args) {
        // mid
        System.out.println("Unit Test mid");
        Point_3 a = new Point_3(0, 0, 0);
        Point_3 b = new Point_3(0, 0, 3);
        Point_3 c = new Point_3(1, 1, 1);
        System.out.println(mid(a, c));

        // geoEdge
        System.out.println("Unit Test geoEdge");
        Polyhedron_3<Point_3> p = MeshLoader.getSurfaceMesh("DATA/shapes-OFF/tetrahedron.off");
        System.out.println(p.edgesToString());
        System.out.println(p.verticesToString());
        Vector_3 v1 = geoEdge(p.halfedges.get(0));
        System.out.println(v1);

        // project
        System.out.println("Unit Test geoEdge");
        Vector_3 v2 = new Vector_3(a, c);
        System.out.println(project(v1, v2));

        // cotSupportAngle
        System.out.println(cotSupportAngle(p.halfedges.get(0)));



    }
}
