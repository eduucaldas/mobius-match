/*

  Adapted into Processing library 5th Feb 2006 Tom Carden
  from "simple arcball use template" 9.16.03 Simon Greenwold
   
  Copyright (c) 2003 Simon Greenwold
  Copyright (c) 2006 Tom Carden

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General
  Public License along with this library; if not, write to the
  Free Software Foundation, Inc., 59 Temple Place, Suite 330,
  Boston, MA  02111-1307  USA

*/


import processing.core.PApplet;

import java.awt.event.MouseEvent;

public class ArcBall {

    PApplet parent;
    public float scaleFactor = 150.f;

    float center_x, center_y, center_z, radius;
    Vec3 v_down, v_drag;
    Quaternion q_now, q_down, q_drag;
    Vec3[] axisSet;
    int axis;

    /**
     * defaults to radius of min(width/2,height/2) and center_z of -radius
     */
    public ArcBall(PApplet parent) {
        this(parent.g.width / 2.0f, parent.g.height / 2.0f, -PApplet.min(parent.g.width / 2.0f, parent.g.height / 2.0f), PApplet.min(parent.g.width / 2.0f, parent.g.height / 2.0f), parent);
    }

    public ArcBall(float center_x, float center_y, float center_z, float radius, PApplet parent) {

        this.parent = parent;

        parent.registerMouseEvent(this);
        parent.registerPre(this);

        this.center_x = center_x;
        this.center_y = center_y;
        this.center_z = center_z;
        this.radius = radius;

        v_down = new Vec3();
        v_drag = new Vec3();

        q_now = new Quaternion();
        q_down = new Quaternion();
        q_drag = new Quaternion();

        axisSet = new Vec3[]{
                new Vec3(1.0f, 0.0f, 0.0f), new Vec3(0.0f, 1.0f, 0.0f), new Vec3(0.0f, 0.0f, 1.0f)};
        axis = -1;  // no constraints...
    }

    public void mouseEvent(MouseEvent event) {
        int id = event.getID();
        if (id == MouseEvent.MOUSE_DRAGGED) {
            mouseDragged();
        } else if (id == MouseEvent.MOUSE_PRESSED) {
            mousePressed();
        }
    }

    public void mousePressed() {
        v_down = mouse_to_sphere(parent.mouseX, parent.mouseY);
        q_down.set(q_now);
        q_drag.reset();
    }

    public void mouseDragged() {
        v_drag = mouse_to_sphere(parent.mouseX, parent.mouseY);
        Vec3 cross = Vec3.cross(v_down, v_drag);
        Vector_3 crossVector = new Vector_3(cross.x, cross.y, cross.z);
        //q_drag.set(Vec3.dot(v_down, v_drag), Vec3.cross(v_down, v_drag));
        q_drag.set((float) Vec3.dot(v_down, v_drag), crossVector);
    }

    public void pre() {
        parent.translate(center_x, center_y, center_z);
        q_now = Quaternion.multiply(q_drag, q_down);
        applyQuat2Matrix(q_now);
        parent.scale(this.scaleFactor);
        parent.translate(-center_x, -center_y, -center_z);
    }

    Vec3 mouse_to_sphere(float x, float y) {
        Vec3 v = new Vec3();
        v.x = (x - center_x) / radius;
        v.y = (y - center_y) / radius;

        float mag = v.x * v.x + v.y * v.y;
        if (mag > 1.0f) {
            v.normalize();
        } else {
            v.z = PApplet.sqrt(1.0f - mag);
        }

        return (axis == -1) ? v : constrain_vector(v, axisSet[axis]);
    }

    Vec3 constrain_vector(Vec3 vector, Vec3 axis) {
        Vec3 res = new Vec3();
        res.sub(vector, Vec3.mul(axis, Vec3.dot(axis, vector)));
        res.normalize();
        return res;
    }

    void applyQuat2Matrix(Quaternion q) {
        // instead of transforming q into a matrix and applying it...

        float[] aa = q.getValue();
        parent.rotate(aa[0], aa[1], aa[2], aa[3]);
    }

}


