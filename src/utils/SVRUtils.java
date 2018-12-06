package utils;

import Jcg.geometry.Point_3;
import Jcg.geometry.Triangle_3;
import Jcg.mesh.MeshBuilder;
import Jcg.mesh.MeshLoader;
import Jcg.mesh.SharedVertexRepresentation;
import Jcg.polyhedron.Polyhedron_3;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.DoubleFunction;

public class SVRUtils {
    /*Only works with triangle meshes*/
    public static void writeMidOFF(String filename) {
        SharedVertexRepresentation svr = new SharedVertexRepresentation(filename);
        svr2off(midEdge(svr), appendMidToFilename(filename));
    }

    public static SharedVertexRepresentation midEdge(SharedVertexRepresentation svr){
        int[] face;
        Triangle_3 t;
        Point_3[] nFace = new Point_3[3];
        ArrayList<Triangle_3> lt = new ArrayList<>();

        for (int k = 0; k < svr.faces.length; k++) {
            if (svr.faceDegrees[k] != 3)
                throw new IllegalArgumentException("computeMidEdgeMesh only works for triangle meshes");
            face = svr.faces[k];
            for (int i = 0, j = 1; i < face.length; i++, j = (j + 1) % svr.faceDegrees[k]) {
                nFace[i] = Point_3.linearCombination(
                        new Point_3[]{svr.points[face[i]], svr.points[face[j]]},
                        new Number[]{0.5, 0.5});
            }
            t = new Triangle_3(nFace[0], nFace[1], nFace[2]);
            lt.add(t);
        }
        return new SharedVertexRepresentation(lt);
    }
    public static void svr2off(SharedVertexRepresentation svr, String filename){
        try (
                FileWriter fw = new FileWriter(filename);
                BufferedWriter bw = new BufferedWriter(fw)
        ){
            bw.write("OFF");
            bw.newLine();
            bw.write(svr.points.length + " " + svr.faces.length + " " + svr.sizeHalfedges);
            bw.newLine();

            for (Point_3 p: svr.points) {
                bw.write(point2str(p));
                bw.newLine();
            }

            bw.newLine();

            for (int[] face: svr.faces){
                bw.write(face2str(face));
                bw.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /* TODO: smarter generation of name */
    private static String appendMidToFilename(String filename) {
        return "OFF/mid.off";

    }

    private static String face2str(int[] face){
        /* Use reduce here */
        String[] sf = Arrays.stream(face)
                            .mapToObj(Integer::toString)
                            .toArray(String[]::new);
        return face.length + "  " + String.join(" ", sf);
    }

    private static String point2str(Point_3 p){
        DoubleFunction<String> format = x -> String.format("%.6f", x);
        return format.apply(p.x) + " " + format.apply(p.y) + " " + format.apply(p.z);
    }


    //Unit Testing
    public static void testOff2Svr2Off(String filename){
        SharedVertexRepresentation svr = new SharedVertexRepresentation(filename);
        svr2off(svr, "OFF/test.off");
    }// At the end should make a diff of test.off and filename

    public static void main(String[] args) {
        writeMidOFF("OFF/tri_triceratops.off");
    }
//
//    public double[][] mapToComplexePlane() {
//        /* this function maps our mid edge mesh to the complex plane
//         *  TO DO:
//         *  1) Determines cut face
//         *  2) Solve linear system for u given the border constraint resulting from 1)
//         *  3) Go through the graph to get u*
//         *  4) embeds for each midge edge vertices its coordinates
//         * */
//        throw new Error("to be implemented");
//    }
//
//    private int findCutFace() {
//        /* this function find the face having a vertex of minimal geodesic distance measured to all other vertice
//         *  this face is cut, then the first vertex X(understand u)value is set to 0
//         *  and other vertex values from the face are set to 1
//         * */
//        throw new Error("to be implemented");
//    }
//
//    private double[] solveForU() {
//        /* SOLVE FOR U using the PColt implementation to have sparse matrix
//         * Will probably need to change some elements in MeshParametrization_PColt to get the right implementation
//         *
//         * */
//        throw new Error("to be implemented");
//    }
//
//    private double[] computeConjugate() {
//        throw new Error("to be implemented");
//    }

}
