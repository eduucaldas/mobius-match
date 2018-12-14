package utils;

import Jcg.geometry.Point_3;
import Jcg.geometry.Triangle_3;
import Jcg.mesh.SharedVertexRepresentation;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;


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
        //DoubleFunction<String> format = x -> String.format("%.6f", x);
        //return format.apply(p.x) + " " + format.apply(p.y) + " " + format.apply(p.z);
        return p.x+" "+p.y+" "+p.z;
    }


    //Unit Testing
    public static void testOff2Svr2Off(String filename){
        SharedVertexRepresentation svr = new SharedVertexRepresentation(filename);
        svr2off(svr, "OFF/test.off");
    }// At the end should make a diff of test.off and filename

    public static void main(String[] args) {
        writeMidOFF("OFF/tri_triceratops.off");
    }

}
