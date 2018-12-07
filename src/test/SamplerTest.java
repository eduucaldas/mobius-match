package test;

import Jcg.geometry.Point_3;
import Jcg.polyhedron.Polyhedron_3;
import meshmanager.SurfaceMesh;
import viewer.MeshViewer;
import Jcg.polyhedron.Vertex;
import algo.Sampler;

public class SamplerTest {
    SurfaceMesh m1;
    MeshViewer view;
    Vertex[] vTbl;
    public SamplerTest(SurfaceMesh m1,MeshViewer view){
        this.m1=m1;
        this.view=view;
        Sampler s=new Sampler(-1,-1);
        this.vTbl=s.sample(m1);
        System.out.println("Scale factor is "+m1.getScaleFactor());
        System.out.println("Number of sampled points: "+this.vTbl.length);
    }
    public void test(){
        view.stroke(255,0,0);
        for(Vertex v:vTbl) {
            m1.drawVertex((Point_3) v.getPoint(),4);
        }
        view.strokeWeight(1);
    }
}
