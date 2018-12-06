package viewer;

import processing.core.*;
import test.*;
import Jcg.geometry.*;
import Jcg.polyhedron.*;
import meshmanager.*;
import parametrization.*;
import test.SamplerTest;

/**
 * A simple 3d viewer for visualizing surface meshes (based on Processing)
 *
 * @author Luca Castelli Aleardi (INF574, 2018)
 */
public class MeshViewer extends PApplet {

    SurfaceMesh renderer; // 3d surface mesh renderer
    int linAlgLibrary = 0;
    int numLinearAlgLibraries = 2;
    int renderType = 1; // choice of type of rendering
    int renderModes = 2; // number of rendering modes
    Polyhedron_3<Point_3> mesh;
    double tolerance = 0.00001;
    public SamplerTest sampler;

    //String filename="OFF/sphere.off";
    //String filename="OFF/cube.off";
    //String filename="OFF/star.off";
    //String filename="OFF/horse1.off";
    String filename = "OFF/tri_triceratops.off";

    public void setup() {
        size(800, 600, P3D);
        ArcBall arcball = new ArcBall(this);

        // open mesh from input file
        this.renderer = new SurfaceMesh(this, filename);
        this.mesh = this.renderer.polyhedron3D;
        this.sampler=new SamplerTest(this.renderer,this);

        int index = 0;
        for (Vertex v : mesh.vertices) {
            v.index = index;
            index++;
        }
    }

    public void draw() {
        background(255);
//        this.lights();
        directionalLight(101, 204, 255, -1, 0, 0);
        directionalLight(51, 102, 126, 0, -1, 0);
        directionalLight(51, 102, 126, 0, 0, -1);
        directionalLight(102, 50, 126, 1, 0, 0);
        directionalLight(51, 50, 102, 0, 1, 0);
        directionalLight(51, 50, 102, 0, 0, 1);

        translate(width / 2.f, height / 2.f, -1 * height / 2.f);
        this.strokeWeight(1);
        stroke(150, 150, 150);

        this.renderer.draw(renderType);
    }

    public void keyPressed() {
        switch (key) {
            case ('p'):
            case ('P'):
                this.computeParameterization(0);
                break;
            case ('r'):
            case ('R'):
                this.renderType = (this.renderType + 1) % this.renderModes;
                break;
            case ('l'):
            case ('L'):
                this.changeLALibrary();
                break;
            default:
                break;
        }

        this.renderer.updateScaleFactor();
    }

    public void changeLALibrary() {
        System.out.print("Changing linear algebra library: ");
        this.linAlgLibrary = (this.linAlgLibrary + 1) % this.numLinearAlgLibraries;

        if (this.linAlgLibrary == 0)
            System.out.println("Jama");
        else if (this.linAlgLibrary == 1)
            System.out.println("PColt");
        else
            System.out.println("Jama");
    }

    public void computeParameterization(int face) {
        TutteLayout2D layout;

        if (this.linAlgLibrary == 0)
            layout = new MeshParameterization_Jama(this.mesh, face);
        else if (this.linAlgLibrary == 1)
            layout = new MeshParameterization_PColt(this.mesh, face);
        else
            layout = new MeshParameterization_Jama(this.mesh, face);

        layout.computeLayout(tolerance);
        layout.projectVertices();
    }

    /**
     * For running the PApplet as Java application
     */
    public static void main(String[] args) {
        //PApplet pa=new MeshViewer();
        //pa.setSize(400, 400);
        PApplet.main(new String[]{"viewer.MeshViewer"});
    }

}
