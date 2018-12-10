package utils;

import Jcg.geometry.Point_;
import Jcg.geometry.Point_3;
import Jcg.polyhedron.Face;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Polyhedron_3;
import Jcg.polyhedron.Vertex;
import meshmanager.SurfaceMesh;

import java.util.ArrayList;

public class MidEdgeFromSurface {
    public static Vertex getMidVertex(Halfedge e){
        Point_3 p1=(Point_3)e.vertex.getPoint();
        Point_3 p2=(Point_3)e.opposite.vertex.getPoint();
        Point_ p=new Point_3((p1.x+p2.x)/2,(p1.y+p2.y)/2,(p1.z+p2.z)/2);
        return new Vertex(p);
    }
    public static ArrayList<Halfedge> MidEdgeFromSurface(SurfaceMesh m1){
        Polyhedron_3 midEdge=new Polyhedron_3(m1.polyhedron3D.vertices.size(),m1.polyhedron3D.facets.size(),m1.polyhedron3D.halfedges.size());
        ArrayList<Halfedge> midEdgeHalfEdges=new ArrayList<>();
        for(Face f:m1.polyhedron3D.facets){
            Halfedge e=f.getEdge();
            Vertex v0=MidEdgeFromSurface.getMidVertex(e);
            Vertex v1=MidEdgeFromSurface.getMidVertex(e.next);
            Vertex v2=MidEdgeFromSurface.getMidVertex(e.prev);
            Halfedge e0=new Halfedge();
            Halfedge e0Opposite=new Halfedge();
            Halfedge e1=new Halfedge();
            Halfedge e1Opposite=new Halfedge();
            Halfedge e2=new Halfedge();
            Halfedge e2Opposite=new Halfedge();
            e0.vertex=v0;
            e0Opposite.vertex=v1;
            e0.opposite=e0Opposite;
            e1.vertex=v1;
            e1Opposite.vertex=v2;
            e1.opposite=e1Opposite;
            e2.vertex=v2;
            e2Opposite.vertex=v0;
            e2.opposite=e2Opposite;
            e0.prev=e1;
            e0.next=e2;
            e2.prev=e0;
            e2.next=e1;
            e1.next=e0;
            e1.prev=e2;
            midEdgeHalfEdges.add(e0);
            midEdgeHalfEdges.add(e1);
            midEdgeHalfEdges.add(e2);
        }
        return midEdgeHalfEdges;
    }
}
