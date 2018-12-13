package algo;

import Jcg.geometry.Point_3;
import Jcg.mesh.MeshBuilder;
import Jcg.mesh.MeshLoader;
import Jcg.mesh.SharedVertexRepresentation;
import Jcg.polyhedron.Face;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Vertex;
import meshmanager.SurfaceMesh;
import parametrization.MobiusParametrization;
import utils.PolyGraph;
import utils.SVRUtils;
import viewer.MeshViewer;
import viewer.ParametrizationViewer;

import java.util.ArrayList;
import java.util.Hashtable;

public class runAlgo {
    MeshViewer mV;
    double[][] correspondences;

    SurfaceMesh m1;
    SurfaceMesh m2;
    Vertex[] sampled1;
    Vertex[] sampled2;
    double[][] c1;
    double[][] c2;

    public runAlgo(MeshViewer mv){
        this.mV=mv;
    }

    private void sample(){
        System.out.println("Parametrization of first mesh");
        m1= mV.m1;
        Sampler sampler1=new Sampler(-1,-1);
        System.out.println("Sampling mesh 1");
        sampled1=sampler1.sample(m1);
        m1.sampled=sampled1;
        m1.displaySampled=true;
        System.out.println("Finished sampling mesh 1");

        m2= mV.m2;
        Sampler sampler2=new Sampler(-1,-1);
        System.out.println("Sampling mesh 2");
        sampled2=sampler2.sample(m2);
        m2.sampled=sampled2;
        m2.displaySampled=true;
        System.out.println("Finished sampling mesh 2");
    }
    private void findCutFace(){
        System.out.println("finding cut face for 1");
        Face cut1=PolyGraph.findCutFace(m1.polyhedron3D);
        m1.cutFace=cut1;
        m1.displayCutFace=true;
        System.out.println("found cut face 1" +cut1.getEdge().vertex.index);

        System.out.println("finding cut face for 2");
        Face cut2=PolyGraph.findCutFace(m2.polyhedron3D);
        m2.cutFace=cut2;
        m2.displayCutFace=true;
        System.out.println("found cut face 2" +cut2.getEdge().vertex.index);

    }
    private Hashtable<Halfedge,double[]> findPlanarEmbedingForDebug(SurfaceMesh m1){
        MobiusParametrization mp1=new MobiusParametrization(m1,m1.cutFace,0.0001);
        System.out.println("Computing parametrization");
        return mp1.planarEmbeding();
    }
    private void parametrize(){
        MobiusParametrization mp1=new MobiusParametrization(m1,m1.cutFace,0.001);
        System.out.println("finding projection sampled");
        c1=mp1.getProjectionFromSampled(sampled1);
        System.out.println("Ended parametrization of 1");
        //Face cut2=m2.polyhedron3D.facets.get(0);
        MobiusParametrization mp2=new MobiusParametrization(m2,m2.cutFace,0.001);
        System.out.println("finding projection sampled");
        c2=mp2.getProjectionFromSampled(sampled2);
        System.out.println("Ended parametrization of 2");
    }
    private void mobiusVote(){
        MobiusVoting mv=new MobiusVoting(c1,c2,-1,0.001);
        System.out.println("Computing correspondences matrix");
        correspondences=mv.createConfidenceMatrix(100);
        /*for(int i=0;i<correspondences.length;i++){
            for(int j=0;j<correspondences[i].length;j++){
                System.out.print(correspondences[i][j]+" ");
            }
            System.out.println("");
        }*/
        System.out.println("Computing fuzzy matrix and giving corresponding vertex");
    }
    private void establishCorrespondenceProcessor(){
        CorrespondenceProcessor cp=new CorrespondenceProcessor(m1,m2,0.7);
        ArrayList<Vertex[]> aV=cp.computeFuzzyCorrespondenceMatrix(correspondences);
        Vertex[] v1=new Vertex[aV.size()];
        Vertex[] v2=new Vertex[aV.size()];
        for(int i=0;i<aV.size();i++){
            v1[i]= aV.get(i)[0];
            v2[i]=aV.get(i)[1];
        }
        System.out.println("Found : "+v1.length+" coresponding vertex in v1 and "+v2.length+" corresponding vertex in v2");
        /*indicate to SurfaceMesh the vertex corresponding to correspondences*/
        m1.correspondence=v1;
        m2.correspondence=v2;
        m1.displayFound=true;
        m2.displayFound=true;
    }
    public void executeAlgorithm(){
        this.sample();
        this.findCutFace();
        this.parametrize();
        this.mobiusVote();
        this.establishCorrespondenceProcessor();
    }
    public void initializeDebug(){
        this.sample();
        this.findCutFace();
    }
    public Hashtable<Halfedge,double[]> executeDebug(SurfaceMesh m0){
        Hashtable<Halfedge,double[]> planarEmbed1=this.findPlanarEmbedingForDebug(m0);
        return planarEmbed1;
    }
}
