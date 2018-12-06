package test;

import Jcg.polyhedron.Polyhedron_3;
import meshmanager.SurfaceMesh;
import viewer.MeshViewer;
import Jcg.polyhedron.Vertex;
import algo.Sampler;

public class SamplerTest {
    String name="OFF/cow.off";
    private void test(){
        MeshViewer viewer=new MeshViewer();
        SurfaceMesh m1=new SurfaceMesh(viewer,this.name);
        Sampler s=new Sampler(-1,-1);
        Vertex[] vTbl=s.sample(m1);
        Polyhedron_3 p=new Polyhedron_3();
        for(Vertex v:vTbl) {
            System.out.println(v.getPoint().toString());
        }
    }
    public static void main(String[] args) {
       SamplerTest st=new SamplerTest();
       st.test();
    }
}
