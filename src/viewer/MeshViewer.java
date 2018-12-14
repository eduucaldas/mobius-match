package viewer;

import Jcg.polyhedron.Vertex;
import algo.runAlgo;
import meshmanager.SurfaceMesh;
import processing.core.PApplet;

/**
 * A simple 3d viewer for visualizing surface meshes (based on Processing)
 *
 * @author Luca Castelli Aleardi (INF574, 2018)
 * <p>
 * We modified this MeshViewer class to fit to our projects need.
 */
public class MeshViewer extends PApplet {

    public SurfaceMesh m1; // 3d surface mesh renderer
    public SurfaceMesh m2;
    int renderType = 1; // choice of type of rendering
    public runAlgo vm;
    int drawnMesh = 0;
    boolean algoIsDone;

    String filename1 = "DATA/non-rigid-world-OFF/cat1.off";
    String filename2 = "DATA/non-rigid-world-OFF/cat2.off";

    //String filename2="DATA/shapes-OFF/tri_triceratops.off";
    public void setup() {
        size(800, 600, P3D);
        ArcBall arcball = new ArcBall(this);

        // open mesh from input file
        this.m1 = new SurfaceMesh(this, filename1);
        //this.mesh = this.m1.polyhedron3D;
        this.m2 = new SurfaceMesh(this, filename2);
        //this.sampler=new SamplerTest(this.renderer,this);
        int index = 0;
        for (Vertex v : m1.polyhedron3D.vertices) {
            v.index = index;
            index++;
        }
        index = 0;
        for (Vertex v : m2.polyhedron3D.vertices) {
            v.index = index;
            index++;
        }
        vm = new runAlgo(this);
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
        if (this.drawnMesh == 0)
            this.m1.draw(renderType);
        else
            this.m2.draw(renderType);
    }

    public void keyPressed() {
        switch (key) {
            case ('e'):
            case ('E'):
                //We execute our algorithm if the e or E button is pressed.
                this.vm.executeAlgorithm();
                this.algoIsDone = true;
                break;
            case ('m'):
            case ('M'):
                this.drawnMesh = 1 - this.drawnMesh;
                break;
            case ('s'):
            case ('S'):
                this.m1.displaySampled = !this.m1.displaySampled;
                this.m2.displaySampled = !this.m2.displaySampled;
                break;
            case ('i'):
            case ('I'):
                this.m1.zoom = this.m1.zoom * 2;
                this.m2.zoom = this.m2.zoom * 2;
                break;
            case ('o'):
            case ('O'):
                this.m1.zoom = this.m1.zoom / 2;
                this.m2.zoom = this.m2.zoom / 2;
                break;
            case ('p'):
            case ('P'):
                if (this.algoIsDone) {
                    this.m1.drawSurface = !this.m1.drawSurface;
                    this.m2.drawSurface = !this.m2.drawSurface;
                }
            default:
                break;
        }
        if (this.drawnMesh == 0)
            this.m1.updateScaleFactor();
        else
            this.m2.updateScaleFactor();
    }

    public static void main(String[] args) {
        PApplet.main(new String[]{"viewer.MeshViewer"});
    }
}
