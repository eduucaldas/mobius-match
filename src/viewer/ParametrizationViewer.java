package viewer;


import algo.runAlgo;
import processing.core.*;
import Jcg.geometry.*;
import Jcg.polyhedron.*;
import meshmanager.*;
import test.SamplerTest;

/**
 * A simple 3d viewer for visualizing surface meshes (based on Processing)
 *
 * @author Luca Castelli Aleardi (INF574, 2018)
 *
 * We modified this MeshViewer class to fit to our projects need.
 */
public class ParametrizationViewer extends MeshViewer {

    public SurfaceMesh m1;
    public SurfaceMesh m2;
    /*
    public SurfaceMesh m1; // 3d surface mesh renderer
    public SurfaceMesh m2;
    int linAlgLibrary = 0;
    int numLinearAlgLibraries = 2;
    int renderType = 1; // choice of type of rendering
    int renderModes = 2; // number of rendering modes
    Polyhedron_3<Point_3> mesh;
    double tolerance = 0.00001;
    public SamplerTest sampler;
    String filename;
    public runAlgo vm;
    int drawnMesh=0;*/

    String filename1="OFF/parametrize1.off";
    String filename2="OFF/parametrize2.off";

    public void setup() {
        size(800, 600, P3D);
        ArcBall arcball = new ArcBall(this);
        // open mesh from input file
        this.m1 = new SurfaceMesh(this, filename1);
        //this.mesh = this.m1.polyhedron3D;
        this.m2 = new SurfaceMesh(this,filename2);
        //this.sampler=new SamplerTest(this.renderer,this);
        int index = 0;
        for (Vertex v : m1.polyhedron3D.vertices) {
            v.index = index;
            index++;
        }
        index=0;
        for(Vertex v:m2.polyhedron3D.vertices){
            v.index=index;
            index++;
        }
        //vm=new runAlgo(this);
    }
    public void draw() {
        background(255);
//        this.lights();
        /*
        directionalLight(101, 204, 255, -1, 0, 0);
        directionalLight(51, 102, 126, 0, -1, 0);
        directionalLight(51, 102, 126, 0, 0, -1);
        directionalLight(102, 50, 126, 1, 0, 0);
        directionalLight(51, 50, 102, 0, 1, 0);
        directionalLight(51, 50, 102, 0, 0, 1);
        */
        translate(width / 2.f, height / 2.f, -1 * height / 2.f);
        this.strokeWeight(1);
        stroke(150, 150, 150);
        if(this.drawnMesh==0)
            this.m1.draw(renderType);
        else
            this.m2.draw(renderType);
    }
    public void keyPressed() {
        switch (key) {
            case('m'):
            case('M'):
                this.drawnMesh=1-this.drawnMesh;
                break;
            case('s'):
            case('S'):
                this.m1.displaySampled=!this.m1.displaySampled;
                this.m2.displaySampled=!this.m2.displaySampled;
            default:
                break;
        }
        if(this.drawnMesh==0)
            this.m1.updateScaleFactor();
        else
            this.m2.updateScaleFactor();
    }
    /**
     * For running the PApplet as Java application
     */
    public static void main() {
        //PApplet pa=new MeshViewer();
        //pa.setSize(400, 400);
        PApplet.main(new String[]{"viewer.ParametrizationViewer"});
    }
}
