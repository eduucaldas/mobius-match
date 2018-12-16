package utils;

import Jcg.geometry.Point_3;
import Jcg.mesh.MeshLoader;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Polyhedron_3;
import Jcg.polyhedron.Vertex;

import java.util.ArrayList;
import java.util.List;

public class MeshUtils {
    public static List<Halfedge<Point_3>> findOutEdges(Vertex<Point_3> v) {
        ArrayList<Halfedge<Point_3>> outE = new ArrayList<>();
        Halfedge<Point_3> begin = v.getHalfedge().opposite;
        outE.add(begin);
        for (Halfedge<Point_3> e = begin.prev.opposite; e != begin; e = e.prev.opposite) {
            outE.add(e);
        }

        return outE;
    }

    public static List<Vertex<Point_3>> findNeighVertices(Vertex<Point_3> v) {
        List<Halfedge<Point_3>> outE = findOutEdges(v);
        ArrayList<Vertex<Point_3>> neigh = new ArrayList<>(outE.size());
        for (Halfedge<Point_3> e:outE) {
            neigh.add(e.getVertex());
        }
        return neigh;
    }


    // Unit Testing
    public static void main(String[] args) {
        Polyhedron_3<Point_3> p = MeshLoader.getSurfaceMesh("DATA/shapes-OFF/cube.off");
        System.out.println(p.edgesToString());
        System.out.println(p.verticesToString());
        List<Halfedge<Point_3>> neigh = findOutEdges(p.vertices.get(0));
        for (Halfedge edge:neigh) {
            System.out.println(edge.toString() + " " + edge.getVertex().toString());
        }

    }
}
