package utils;

import Jcg.geometry.Point_3;
import Jcg.mesh.MeshLoader;
import Jcg.polyhedron.Face;
import Jcg.polyhedron.Polyhedron_3;

import java.util.Arrays;

public class PolyGraph {
    /* Pseudocode for finding the findCutFace
     *   For every face:
     *       getVertexIndices i,j,k
     *       w[a,b] <- dist(V[a], V[b]), for i,j,k
     *
     *   floyd-warshall using w -> compute the wanted vertex index x
     *
     *   V[x] -(getHalfedge)> -(getFace)> wanted result
     * */

    /* For the moment we`zre using a functional approach */

    private static double[][] computeEdges(Polyhedron_3<Point_3> poly) {
        int nVert = poly.vertices.size();
        double[][] edge = new double[nVert][nVert];
        /*  Remark: As the graph from a mesh is undirected
                    we work with the upper triangle of the edgeWeight matrix
                    and then we reflect it
        */

        // Initialization
        for (int i = 0; i < nVert; i++) {
            for (int j = i + 1; j < nVert; j++) {
                edge[i][j] = Double.POSITIVE_INFINITY;
            }
            edge[i][i] = 0.;
        }

        Point_3[] vertexCoords;
        int[] vertexIds;

        // Compute edges for every face (upper triangle)
        for (Face face : poly.facets) {
            vertexIds = face.getVertexIndices(poly);
            vertexCoords = Arrays.stream(vertexIds)
                    .mapToObj(i -> poly.vertices.get(i).getPoint())
                    .toArray(Point_3[]::new);

            for (int i = 0, j = 1; i < face.degree(); i++, j = (j + 1) % face.degree()) {
                edge[vertexIds[i]][vertexIds[j]] = vertexCoords[j].distanceFrom(vertexCoords[i]).doubleValue();
            }
        }

        // reflection of the edgeWeight matrix
        for (int i = 0; i < nVert; i++) {
            for (int j = i + 1; j < nVert; j++) {
                edge[j][i] = edge[i][j];
            }
        }

        return edge;
    }

    private static void unitTestComputeEdges(String filename) {
        // They should be equal
        double[][] edgeWeightsViaPolyhedron = computeEdges(MeshLoader.getSurfaceMesh(filename));
        // First let`s just compare their sum, should be enough
        double ew;
        for (int i = 0; i < edgeWeightsViaPolyhedron.length; i++) {
            for (int j = 0; j < edgeWeightsViaPolyhedron[0].length; j++) {
                ew = edgeWeightsViaPolyhedron[i][j];
                System.out.print(ew + " ");
            }
            System.out.println();
        }
    }

    private static double[][] floydWarshall(double[][] edges) {
        int nVert = edges.length;
        double[][] dist = new double[nVert][nVert];
        for (int k = 0; k < nVert; k++) {
            for (int i = 0; i < nVert; i++) {
                for (int j = 0; j < nVert; j++) {
                    if (dist[i][j] > dist[i][k] + dist[k][j])
                        dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
        return dist;
    }

    private static int indexOfMinGeodesicAvgVertex(double[][] dist) {
        int indexOfMin = 0;
        int nVert = dist.length;
        double geodesicSum;
        double minGeodesicSum = Double.POSITIVE_INFINITY;

        for (int i = 0; i < nVert; i++) {
            geodesicSum = 0;
            for (int j = 0; j < nVert; j++) {
                geodesicSum += dist[i][j];
            }
            if (geodesicSum < minGeodesicSum) {
                minGeodesicSum = geodesicSum;
                indexOfMin = i;
            }
        }
        return indexOfMin;
    }

    public static Face<Point_3> findCutFace(Polyhedron_3<Point_3> poly) {
        double[][] edges = computeEdges(poly);

        double[][] geodesicDistances = floydWarshall(edges);

        int cutVertexIndex = indexOfMinGeodesicAvgVertex(geodesicDistances);

        return poly.vertices.get(cutVertexIndex).getHalfedge().getFace();

    }

    public static void main(String[] args) {
        String uTComputeEdgesFilename = "DATA/shapes-OFF/tetrahedron.off";
        unitTestComputeEdges(uTComputeEdgesFilename);

        String findCutFaceFilename = "DATA/shapes-OFF/tri_triceratops.off";
        Polyhedron_3<Point_3> triceratopsPoly = MeshLoader.getSurfaceMesh(findCutFaceFilename);
        Face<Point_3> cutFace = findCutFace(triceratopsPoly);
        System.out.println(cutFace.toString());

    }


}
