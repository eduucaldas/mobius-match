package utils;

import Jcg.geometry.Point_3;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Vertex;

import java.util.List;

public class MidEdgeMesh {
    public static Point_3 midVertexFromEdge(Halfedge<Point_3> e) {
        return GeometryUtils.mid(
                e.getVertex().getPoint(),
                e.opposite.getVertex().getPoint());
    }

    public static Point_3 findClosestMidNeighbor(Vertex v) {
        List<Halfedge<Point_3>> outE = MeshUtils.findOutEdges(v);
        Halfedge<Point_3> minEdge = outE.get(0);
        double minNorm = Double.POSITIVE_INFINITY;
        double normE;
        for (Halfedge e:outE) {
            normE = GeometryUtils.normEdge(e);
            if(minNorm > normE){
                minNorm = normE;
                minEdge = e;
            }
        }

        return midVertexFromEdge(minEdge);
    }
}
