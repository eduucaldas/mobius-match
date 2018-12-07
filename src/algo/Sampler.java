package algo;

import java.util.*;

import Jcg.geometry.Point_3;
import Jcg.geometry.Vector_3;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Vertex;
import meshmanager.SurfaceMesh;

public class Sampler {
    /* The algo.Sampler class enables us to do the first step, that is sampling from two meshes.
        Starting by sampling local maxima of Gauss Curvature by using a smoothest version of the angle-deficit formula [Desbrun et al. 2002]
        It then uses the Farthest Point Algorithm to take a spread of points, computing geodesic distances with an approximate algorithm based on
        Dijkstra's Shortest path algorithm.
     */
    int N;
    double epsilon;
    public Sampler(int N,double epsilon){
        if(N!=-1) {
            this.N = N;
        }
        else this.N=100; //the paper average use
        if(epsilon!=-1){
            this.epsilon=epsilon;
        }
        else this.epsilon=0.001;
    }
    private double gaussCurvature(Vertex<Point_3> v){
        /* computes gaussian curvature of a vertex -- could stay unused!!!*/
        double curvature=2*Math.PI;
        double area=0;
        Halfedge<Point_3> e=v.getHalfedge();
        Halfedge<Point_3> f1=e.opposite;
        Halfedge<Point_3> f2=e.next.opposite;
        Point_3 p0=v.getPoint();
        while(!f2.equals(e)) {
            f2 = f2.opposite;
            Point_3 p1 = f1.getVertex().getPoint();
            Point_3 p2 = f2.getVertex().getPoint();
            Vector_3 v2 = new Vector_3(p0, p2);
            Vector_3 v1 = new Vector_3(p0, p1);
            double pdt = (double) v1.innerProduct(v2);
            double angle = Math.acos(pdt / Math.sqrt((double) v1.squaredLength() * (double) v2.squaredLength()));
            curvature += -angle;
            area += pdt;
            f1 = f2;
            f2 = f2.opposite.next.opposite;
        }
        return curvature/area;
    }
    private static Vector_3 gaussCurvatureDerivative(Vertex<Point_3> v){
        Vector_3 angleGrad=new Vector_3(0,0,0);
        Halfedge<Point_3> e=v.getHalfedge();
        Halfedge<Point_3> f1=e.opposite;
        Halfedge<Point_3> f2=e.next.opposite;
        Point_3 p0=v.getPoint();
        do {
            f2 = f2.opposite;
            Point_3 p1 = f1.getVertex().getPoint();
            Point_3 p2 = f2.getVertex().getPoint();
            //our first angle : delta ij
            Vector_3 v10 = new Vector_3(p1,p0);
            Vector_3 v12 = new Vector_3(p1, p2);
            double pdt = (double) v12.innerProduct(v10);
            double crosspdt =Math.sqrt((double)v12.crossProduct(v10).squaredLength());
            double cotDelta=pdt/crosspdt;
            //our second angle: gamma ij
            Point_3 p3=f1.getNext().getVertex().getPoint();
            Vector_3 v13 = new Vector_3(p1,p3);
            pdt = (double) v13.innerProduct(v10);
            crosspdt =Math.sqrt((double)v13.crossProduct(v10).squaredLength());
            double cotGamma=pdt/crosspdt;
            //we then add the participation from v10 (that is the edge f1 to the gradient).
            angleGrad=angleGrad.sum(v10.multiplyByScalar((cotGamma+cotDelta)/(double)v10.squaredLength()));
            f1 = f2;
            f2 = f2.opposite.next.opposite;
        } while(!f2.equals(e.next.opposite));
        return angleGrad;
    }
    private static Vertex[] findLocalGaussCurvature(double epsilon,SurfaceMesh M1){
        /*
        * Our methods gaussCurvatureDerivative can give us for each Vertex the derivative of Gauss Curvature.
        * We will send back an array with all selected points
        * */
        ArrayList<Vertex> indexOfLocalMaximaM1=new ArrayList<Vertex>();
        for(Vertex<Point_3> v:M1.polyhedron3D.vertices){
            Vector_3 g=Sampler.gaussCurvatureDerivative(v);
            if((double)g.squaredLength()<epsilon) {
                System.out.print(g.x+" "+g.y+" "+g.z);
                indexOfLocalMaximaM1.add(v);
            }
        }
        Vertex[] v=new Vertex[indexOfLocalMaximaM1.size()];
        for(int i=0;i<indexOfLocalMaximaM1.size();i++){
            v[i]=indexOfLocalMaximaM1.get(i);
        }
        return v;
    }
    public static Vertex[] findNeighbour(Vertex v,SurfaceMesh m1){
        /* gives back all neighbours vertex*/
        Vertex[] neighbors=new Vertex[m1.polyhedron3D.vertexDegree(v)];
        Halfedge e=v.getHalfedge();
        Halfedge f=e.next.opposite;
        int i=0;
        while(!f.equals(e)){
            neighbors[i]=f.opposite.vertex;
            i++;
            f=f.next.opposite;
        }
        neighbors[i]=f.opposite.vertex;
        return neighbors;
    }
    public static Hashtable<Vertex,Double> executefps(ArrayList<Vertex> sources, Hashtable<Vertex,Double> distTable,SurfaceMesh m1){
        boolean continuer=true;
        while(continuer) {
            ArrayList<Vertex> newSource = new ArrayList<>();
            for (Vertex source : sources) {
                Vertex[] neighbors = Sampler.findNeighbour(source, m1);
                for (Vertex neighbor : neighbors) {
                    double distance = (double) neighbor.getPoint().squareDistance(source.getPoint());
                    if (distTable.get(neighbor) == -1.) {
                        distTable.replace(neighbor, distTable.get(source) + distance);
                        newSource.add(neighbor);
                    } else if (distTable.get(neighbor) > distTable.get(source) + distance) {
                        distTable.replace(neighbor, distTable.get(source) + distance);
                        newSource.add(neighbor);
                    }
                }
            }
            if (newSource.isEmpty()) {
                continuer = false;
            } else {
                sources = newSource;
            }
        }
        return distTable;
    }
    private static Vertex fps(ArrayList<Vertex> initialSample,SurfaceMesh m1){
        /* Realise farthest point search: from sourced points it looks for the geodesically most distant points
        * It repeats the following algorithm:
        * Initialisation:
        *     each source vertex is put in the hash, with dist = 0.
        *     all other vertex are put in the Hash, with dist= -1.(infinity)
        *     try: 1) For each source points look for Reacheable points.
        *                   if its dist is -1: then set new real dist and add this point to the new sources table.
        *                   else if vertexCurrentDist<sources dist+dist(sources-vertex) do nothing
        *                        else set new real min dist and add this point to the new sources table.
        *    if sources is not empty:
        *          4) Update the sources list.
        *     if not: we have explore the whole graph, thus finish by looking at the max dist in visited sources.
        *
        *     We could also implement the algorithm the following way:
        *     each source vertex is put in the hash, with dist = 0, and in a SortedHash relative to distance
        *     all other vertex are put in the Hash, with dist= -1.(infinity)
        *     While there are elements in the sortedHash:
        *       We take the first (minimal distance) element, and expand it to neighbours, the same way as before.
        *
        * */
        Hashtable<Vertex,Double> distTable=new Hashtable<>();
        for(Vertex s:m1.polyhedron3D.vertices){
            distTable.put(s,-1.);
        }
        if(initialSample.size()>0) {
            for (Vertex s : initialSample) {
                distTable.replace(s, 0.);
            }
        }
        else{
            distTable.replace(m1.polyhedron3D.vertices.get(0),0.);
        }
        ArrayList<Vertex> sources=initialSample;
        distTable=Sampler.executefps(sources,distTable,m1);
        //each point's dist should be at the minimum distances from each initial sources.
        double max=0;
        Vertex vmax=null;
        for(Vertex v:distTable.keySet()){
            if(distTable.get(v)>max){
                max=distTable.get(v);
                vmax=v;
            }
        }
        return vmax;
    }
    private Vertex[] extendsSampling(Vertex[] GaussCurvLocalMax,SurfaceMesh m1){
        /* Starting From local maxima of Gauss Curvature, it spreads to capture N points
        * TO DO:
        * 1) Implement FPS algorithm in case we don't have enough maxima of Gauss Curvature
        * 1.1) This algorithm idea is to take the vertex with max distances from  its min distances to our presampled points
        *      Then it add this vertex to the sampled points and reiterates.
        *      The issue is that we are using geodesic distances, so we need to build a graph from scratch every time we want to compute
        *      distances from one point to other points.
        * 2) use this algorithm starting from the initial sample and extend while we haven't N sampled points...
        * */
        ArrayList<Vertex> v=new ArrayList<>();
        if(GaussCurvLocalMax.length>=this.N){
            Vertex[] v2=new Vertex[this.N];
            for(int i=0;i<N;i++){
                v2[i]=GaussCurvLocalMax[i];
            }
            return v2;
        }
        else{
            System.out.println(" ");
            System.out.println(" number of Gauss Max: "+GaussCurvLocalMax.length);
            for(int i=0;i<GaussCurvLocalMax.length;i++){
                v.add(GaussCurvLocalMax[i]);
            }
            for(int i=0;i<Math.min(this.N,m1.polyhedron3D.vertices.size())-GaussCurvLocalMax.length;i++){
                v.add(Sampler.fps(v,m1));
            }
        }
        Vertex[] v2=new Vertex[v.size()];
        for(int i=0;i<v.size();i++){
            v2[i]=v.get(i);
        }
        return v2;
    }
    public Vertex[] sample(SurfaceMesh m1){
        /* sample points */
        Vertex[] GaussCurvLocalMax;
        this.epsilon/=10;
        do {
            this.epsilon*=10;
            GaussCurvLocalMax = Sampler.findLocalGaussCurvature(this.epsilon, m1);
        }while(GaussCurvLocalMax.length==0);
        Vertex[] sampled=this.extendsSampling(GaussCurvLocalMax,m1);
        return sampled;
    }
}
