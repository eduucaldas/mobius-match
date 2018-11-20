import processing.core.*;

import java.util.*;
import java.applet.*;
import java.awt.*;
import java.awt.event.*;

/**
 * A class for drawing 2D curves and handling mouse events
 *
 * @author Luca Castelli Aleardi (INF555, 2012)
 */
public class DrawCurve extends PApplet implements ActionListener {

    int maxCurves = 10;
    int currentCurve = -1; // no curve scheme defined (at the beginning)
    int curveType = 0;

    Curve[] schemes = new Curve[maxCurves]; // curve
    LinkedList<Point_2>[] points = new LinkedList[maxCurves]; // the input points to interpolate
    final static double epsilon = 15.;


    public Point_2 selectedPoint = null; // point selected with mouse events

    public void setup() {
        initButton();
        size(500, 500); // size of the window
        for (int i = 0; i < maxCurves; i++) {
            this.points[i] = new LinkedList<Point_2>();
        }
    }

    public void initButton() {
        button1 = new Button("New Bezier");
        add(button1);
        button1.addActionListener(this);

			/*button2 = new Button("Reset");
			add(button2);
			button2.addActionListener(this);*/

        button3 = new Button("Rational Bezier");
        add(button3);
        button3.addActionListener(this);
    }

    public void actionPerformed(ActionEvent e) {
        if (e.getSource() == button1) {
            curveType = 0;
            currentCurve++;
            this.schemes[currentCurve] = new Bezier(this, this.points[currentCurve]);

            System.out.println("New Bezier");
        }
        if (e.getSource() == button3) {
            curveType = 1;
            currentCurve++;
            this.schemes[currentCurve] = new RationalBezier(this, this.points[currentCurve]);

            System.out.println("New rational Bezier");
        } else
            System.out.println("Button 2 was pressed");
    }

    Button button1, button2, button3;

    public void draw() {
        background(220);

        if (points.length < 1) return; // no points to interpolate
        if (currentCurve < 0) return;
        // choose the rendering method
        for (int i = 0; i < schemes.length; i++) {
            if (points[i].size() > 1){
                schemes[i].plotCurve(1. / 100.);
                schemes[i].subdivisionRendering(4); // render the curve by recursive subdivision
            }
        }
    }

    public void removePoint(int x, int y) {
        int curve = this.findCurve(x, y);
        int index = findPoint(x, y, curve);
        if (index >= 0 && index < this.points[curve].size())
            this.points[curve].remove(index);
    }

    public int findPoint(int x, int y, int curve) {
        if (curve < 0 || curve >= this.schemes.length)
            return -1;
        Point_2 p = new Point_2(x, y);

        int index = 0;
        boolean found = false;
        for (Point_2 q : this.points[curve]) {
            if (q.squareDistance(p) < epsilon) {
                found = true;
                break;
            }
            index++;
        }
        if (found == true)
            return index;
        else
            return -1;
    }

    public int findCurve(int x, int y) {
        for (int i = 0; i < this.schemes.length; i++) {
            int result = findPoint(x, y, i);
            if (result >= 0)
                return i;
        }
        return -1;
    }

    public Point_2 selectPoint(int x, int y) {
        Point_2 p = new Point_2(x, y);

        for (int i = 0; i < currentCurve; i++) {
            for (Point_2 q : this.points[i]) {
                if (q.squareDistance(p) < epsilon) {
                    return q;
                }
            }
        }
        return null;
    }

    public void mouseClicked() {
        if (currentCurve < 0) return;
        Point_2 p = new Point_2(mouseX, mouseY);

        if (mouseButton == LEFT && this.selectedPoint == null) {
            this.points[currentCurve].add(p);
        } else if (mouseButton == RIGHT)
            removePoint(mouseX, mouseY);

        this.schemes[currentCurve].updateInputPoints(this.points[currentCurve]);
    }

    public void mousePressed() {
        this.selectedPoint = selectPoint(mouseX, mouseY);
    }

    public void mouseReleased() {

        if (this.selectedPoint != null) {
            this.selectedPoint.setX(mouseX);
            this.selectedPoint.setY(mouseY);
        }
    }

    public void mouseDragged() {
    }

}
