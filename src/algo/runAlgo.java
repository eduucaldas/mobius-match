package algo;

import Jcg.polyhedron.Face;
import Jcg.polyhedron.Vertex;
import algo.CorrespondenceProcessor;
import algo.MobiusVoting;
import algo.Sampler;
import javafx.scene.shape.Mesh;
import meshmanager.SurfaceMesh;
import parametrization.MobiusParametrization;
import processing.core.PApplet;
import utils.PolyGraph;
import viewer.MeshViewer;

import java.util.ArrayList;

public class runAlgo {
    MeshViewer mV;
    double[][] correspondences;

    public runAlgo(MeshViewer mv){
        this.mV=mv;
    }
    public void executeAlgorithm(){
        System.out.println("Parametrization of first mesh");
        SurfaceMesh m1= mV.m1;
        Sampler sampler1=new Sampler(-1,-1);
        System.out.println("Sampling mesh 1");
        Vertex[] sampled1=sampler1.sample(m1);
        m1.sampled=sampled1;
        m1.displaySampled=true;
        System.out.println("Finished sampling mesh 1");

        SurfaceMesh m2= mV.m2;
        Sampler sampler2=new Sampler(-1,-1);
        System.out.println("Sampling mesh 1");
        Vertex[] sampled2=sampler2.sample(m2);
        m2.sampled=sampled2;
        m2.displaySampled=true;
        System.out.println("Finished sampling mesh 2");

        System.out.println("Parametrizing 1:");
        System.out.println("finding cut face");
        Face cut1=PolyGraph.findCutFace(m1.polyhedron3D);
        System.out.println("found cut face");
        MobiusParametrization mp1=new MobiusParametrization(m1,cut1,0.001);
        System.out.println("finding projection sampled");
        double[][] c1=mp1.getProjectionFromSampled(sampled1);
        System.out.println("Ended parametrization of 1");

        System.out.println("Parametrizing 2:");
        System.out.println("finding cut face");
        Face cut2=PolyGraph.findCutFace(m2.polyhedron3D);
        MobiusParametrization mp2=new MobiusParametrization(m2,cut2,0.001);
        System.out.println("finding projection sampled");
        double[][] c2=mp2.getProjectionFromSampled(sampled1);
        System.out.println("Ended parametrization of 2");

        MobiusVoting mv=new MobiusVoting(c1,c2,-1,0.001);
        System.out.println("Computing correspondences matrix");
        correspondences=mv.createConfidenceMatrix(100);
        System.out.println("Computing fuzzy matrix and giving corresponding vertex");
        CorrespondenceProcessor cp=new CorrespondenceProcessor(m1,m2,0.001);
        ArrayList<Vertex[]> aV=cp.computeFuzzyCorrespondenceMatrix(correspondences);
        Vertex[] v1=new Vertex[aV.size()];
        Vertex[] v2=new Vertex[aV.size()];
        for(int i=0;i<aV.size();i++){
            v1[i]= aV.get(i)[0];
            v2[i]=aV.get(i)[1];
        }
        /*indicate to SurfaceMesh the vertex corresponding to correspondences*/
        m1.correspondence=v1;
        m2.correspondence=v2;
        m1.displayFound=true;
        m2.displayFound=true;
    }
}
