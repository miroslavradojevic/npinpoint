package com.braincadet.npinpoint;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.plugin.PlugIn;

import java.awt.*;
import java.io.File;

public class Viewer implements PlugIn {

    public void run(String s) {

//        IJ.open();
//        ImagePlus img = IJ.getImage();

        String img_path = "";
        String det_path = ""; // img.getOriginalFileInfo().directory;
        float r = 3;
        boolean show_dirs = false;

        GenericDialog gd = new GenericDialog("Select Detection");
        gd.addStringField("image path", img_path, 60);
        gd.addStringField("detection path [.det]", det_path, 60);
        gd.addNumericField("radius", r, 0);
        gd.addCheckbox("directions", show_dirs);
        gd.showDialog();

        if (gd.wasCanceled()) return;

        img_path = gd.getNextString();
        det_path = gd.getNextString();
        r = (float) gd.getNextNumber();
        show_dirs = gd.getNextBoolean();

        // read image
        ImagePlus img = new ImagePlus(img_path);
        if (img==null) {
            IJ.log(img_path+"  was not read");
            return;
        }

		// read detection
        File f_det = new File(det_path);

        if (!f_det.exists()) {
            IJ.log("file does not exist");
            return;
        }

        if (!Tools.getFileExtension(det_path).equals("det")) {
            IJ.log("file needs to be .det");
            return;
        }

        // read csv file (.det) at specified location
        ReadDET det_reader = new ReadDET(det_path);

        det_reader.print();

        // loop through read values and extract an Overlay for the viewer

        Overlay ov_read = new Overlay();

        for (int i = 0; i < det_reader.x.size(); i++) { // loop read detections

            float dx,dy;

            float x = det_reader.x.get(i);
            float y = det_reader.y.get(i);
            if (r==-1) r = det_reader.r.get(i); // will change the size of the detection overlay circles

            boolean body_added = false;

            Color clr=null;

            if (det_reader.t.get(i).equals("END")) {
//				clr = Color.YELLOW;
                clr = new Color(1, 1, 0, 0.6f);
                body_added = true;
            }
            else if (det_reader.t.get(i).equals("BIF")) {
//				clr = Color.RED;
                clr = new Color(1, 0, 0, 0.6f);
                body_added = true;
            }
            else if (det_reader.t.get(i).equals("CROSS")){
//				clr = Color.GREEN;
                //clr = new Color(0, 1, 0, 0.6f);
                clr = new Color(1, 0, 0, 0.6f);
                body_added = true;
            }

            OvalRoi regroi = new OvalRoi(x-r, y-r, 2*r, 2*r);
            regroi.setStrokeWidth(1);
            regroi.setStrokeColor(clr);
            regroi.setFillColor(clr);
//			if (det_reader.t.get(i).equals("END"))
            ov_read.add(regroi);

            if (show_dirs && body_added) { // add directions if they are !NaN

                float scale_direction = 2f; // just for visualization
                int nr_directions = det_reader.v.get(i).length;
                for (int j = 0; j < nr_directions; j++) {

                    dx = det_reader.v.get(i)[j][0];
                    dy = det_reader.v.get(i)[j][1];

                    if (!Float.isNaN(dx) && !Float.isNaN(dy)) {
                        Line l = new Line(x+r*dx, y+r*dy, x+scale_direction*r*dx, y+scale_direction*r*dy);
                        l.setStrokeWidth(1f);
                        l.setStrokeColor(clr);
                        l.setFillColor(clr);
                        ov_read.add(l);
                    }

                }


            }

        }

        img.setOverlay(ov_read);
        img.show();

    }
}
