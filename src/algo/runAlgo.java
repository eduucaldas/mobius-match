package algo;

import Jcg.polyhedron.Face;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Vertex;
import meshmanager.SurfaceMesh;
import parametrization.MobiusParametrization;
import utils.PolyGraph;
import viewer.MeshViewer;

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

    public runAlgo(MeshViewer mv) {
        this.mV = mv;
    }

    private void sample(int N, double epsilon) {
        /* Sample the mesh
                - N defines the number of sample to take
                - epsilon is used to discover local maximum of Gauss Curvature.
                    Vertex are considered local maxima of Gauss Curvature when their gradient value is below epsilon
                    While we don't have find any local maxima of Gauss Curvature we multiply epsilon by 10 at each epoch....
                Giving N=-1 and epsilon=-1 calls for default values (100 and 0.001 respectively)
         */
        System.out.println("Parametrization of first mesh");
        m1 = mV.m1;
        Sampler sampler1 = new Sampler(N, epsilon);
        System.out.println("Sampling mesh 1");
        sampled1 = sampler1.sample(m1);
        m1.sampled = sampled1;
        m1.displaySampled = true;
        System.out.println("Finished sampling mesh 1");

        m2 = mV.m2;
        Sampler sampler2 = new Sampler(N, epsilon);
        System.out.println("Sampling mesh 2");
        sampled2 = sampler2.sample(m2);
        m2.sampled = sampled2;
        m2.displaySampled = true;
        System.out.println("Finished sampling mesh 2");
    }

    private void findCutFace() {
        /*  Finds cut faces
         */
        System.out.println("finding cut face for 1");
        Face cut1 = PolyGraph.findCutFace(m1.polyhedron3D);
        m1.cutFace = cut1;
        m1.displayCutFace = true;
        System.out.println("found cut face 1" + cut1.getEdge().vertex.index);

        System.out.println("finding cut face for 2");
        Face cut2 = PolyGraph.findCutFace(m2.polyhedron3D);
        m2.cutFace = cut2;
        m2.displayCutFace = true;
        System.out.println("found cut face 2" + cut2.getEdge().vertex.index);

    }

    private Hashtable<Halfedge, double[]> findPlanarEmbedingForDebug(SurfaceMesh m1) {
        MobiusParametrization mp1 = new MobiusParametrization(m1, m1.cutFace, 0.0001);
        System.out.println("Computing parametrization");
        return mp1.planarEmbeding();
    }

    private void parametrize(double precision) {
        /* Parametrizing the mesh,
            - the precision parameter is used in the linear solver as the precision to reach.
         */
        MobiusParametrization mp1 = new MobiusParametrization(m1, m1.cutFace, precision);
        System.out.println("finding projection sampled");
        c1 = mp1.getProjectionFromSampled(sampled1);
        this.m1.planarEmbedding = mp1.planarEmbed;
        System.out.println("Ended parametrization of 1");
        //Face cut2=m2.polyhedron3D.facets.get(0);
        MobiusParametrization mp2 = new MobiusParametrization(m2, m2.cutFace, precision);
        System.out.println("finding projection sampled");
        c2 = mp2.getProjectionFromSampled(sampled2);
        this.m2.planarEmbedding = mp2.planarEmbed;
        System.out.println("Ended parametrization of 2");
    }

    private void mobiusVote(double epsilon, int numberOfEpoch) {
        /*Launch Mobius Vote,
            - The epsilon parameter is used in the computation of correspondence Matrix value
            - The second parameter is the number of epoch to run. The paper recommends to take
         */
        MobiusVoting mv = new MobiusVoting(c1, c2, -1, epsilon);
        System.out.println("Computing correspondences matrix");

        correspondences = mv.createConfidenceMatrix(numberOfEpoch);
    }

    private void establishCorrespondenceProcessor() {
        /* Computes fuzzy Matrix
            - the threshold parameter is the parameter to do discriminate first values.
         */
        System.out.println("Computing fuzzy matrix and giving corresponding vertex");
        CorrespondenceProcessor cp = new CorrespondenceProcessor(m1, m2, 0.97);
        ArrayList<Vertex[]> aV = cp.computeFuzzyCorrespondenceMatrix(correspondences);
        Vertex[] v1 = new Vertex[aV.size()];
        Vertex[] v2 = new Vertex[aV.size()];
        for (int i = 0; i < aV.size(); i++) {
            v1[i] = aV.get(i)[0];
            v2[i] = aV.get(i)[1];
        }
        /*indicate to SurfaceMesh the vertex corresponding to correspondences*/
        m1.correspondence = v1;
        m2.correspondence = v2;
        m1.displayFound = true;
        m2.displayFound = true;
    }

    public void executeAlgorithm() {
        this.sample(-1, -1);
        this.findCutFace();
        this.parametrize(0.001);
        //For more accurate solution use NumberOfEpoch=10*(int)Math.pow(this.c1.length,3)
        this.mobiusVote(0.001, 10 * (int) Math.pow(this.c1.length, 2));
        this.establishCorrespondenceProcessor();
    }

    public void initializeDebug() {
        this.sample(-1, -1);
        this.findCutFace();
    }

    public Hashtable<Halfedge, double[]> executeDebug(SurfaceMesh m0) {
        Hashtable<Halfedge, double[]> planarEmbed1 = this.findPlanarEmbedingForDebug(m0);
        return planarEmbed1;
    }
}
