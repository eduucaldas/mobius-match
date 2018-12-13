package test;

import Jcg.geometry.Point_3;
import Jcg.geometry.Vector_3;
import Jcg.mesh.MeshLoader;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Polyhedron_3;
import Jcg.polyhedron.Vertex;
import meshmanager.SurfaceMesh;

import java.util.ArrayList;
import java.util.Hashtable;

public class MeshBrowsingTest {
    static String filename="OFF/cow.off";
    private static void updateSurround(Hashtable<Halfedge,Double> midEdgeMesh,Halfedge p){
        /*to minimize computation time, if we see that opposite halfedge is in the hash table we use its value*/
        if(!midEdgeMesh.containsKey(p))
            midEdgeMesh.put(p,midEdgeMesh.get(p.opposite));
        if(!midEdgeMesh.containsKey(p.prev)){
            double val;
            if(midEdgeMesh.containsKey(p.prev.opposite)){
                val=midEdgeMesh.get(p.prev.opposite);
            }
            else
                val=0;
            midEdgeMesh.put(p.prev,val);

        }
        if(!midEdgeMesh.containsKey(p.next)) {
            double val;
            if(midEdgeMesh.containsKey(p.next.opposite)){
                val=midEdgeMesh.get(p.next.opposite);
            }
            else
                val=0;
            midEdgeMesh.put(p.next,val);
        }
    }
    public static Hashtable<Halfedge,Double> getMidEdgeMesh(Polyhedron_3 f){
        ArrayList<Halfedge> sources=new ArrayList<>();
        Hashtable<Halfedge,Double> midEdgeMesh= new Hashtable<>();
        Halfedge e0=((Vertex)f.vertices.get(0)).getHalfedge();
        midEdgeMesh.put(e0,0.);
        sources.add(e0);
        while(!sources.isEmpty()){
            Halfedge e=sources.remove(0);
            updateSurround(midEdgeMesh,e);
            if(!midEdgeMesh.containsKey(e.next.opposite)) {
                midEdgeMesh.put(e.next.opposite, midEdgeMesh.get(e.next));
                sources.add(e.next.opposite);
            }
            if(!midEdgeMesh.containsKey(e.prev.opposite)) {
                midEdgeMesh.put(e.prev.opposite, midEdgeMesh.get(e.prev));
                sources.add(e.prev.opposite);
            }
        }
        return midEdgeMesh;
    }

    //We initially try a browsing algorithm:
    public static Hashtable<Halfedge, Double> getMidEdgeMesh2(Polyhedron_3 f){
        Hashtable<Halfedge,Double> midEdgeMesh= new Hashtable<>();
        Halfedge e0=((Vertex)f.vertices.get(0)).getHalfedge();
        midEdgeMesh.put(e0,0.);
        updateSurround(midEdgeMesh,e0);
        boolean browse=true;
        Halfedge e=e0.prev.opposite;
        while(browse){
            //We compute values for the face
            updateSurround(midEdgeMesh,e);
            //we change the face we are looking at.
            e=e.prev.opposite;
            if(midEdgeMesh.containsKey(e)){
                e=e.opposite.next.next.opposite;
                if(midEdgeMesh.containsKey(e)){
                    browse=false;
                }
            }
        }
        return midEdgeMesh;
    }

    public static void main(String[] args){
        /*Polyhedron_3 f= MeshLoader.getSurfaceMesh(filename);
        Hashtable<Halfedge,Double> midEdge=getMidEdgeMesh(f);
        System.out.println("Mid edge number of halfedge: "+midEdge.size());
        System.out.println("Initial number of halfedge: "+f.halfedges.size());*/

    }


}
