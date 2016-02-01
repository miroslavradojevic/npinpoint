package com.braincadet.npinpoint; /** modified code of the ImageJ plugin "Find Connected Regions" under the terms of the GNU General Public License
    Copyright 2006, 2007 Mark Longair
    modified by Miroslav Radojevic 15-01-16
    find connected regions with the same value in 8 bit images */
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;

import java.util.Collections;
import java.util.Iterator;
import java.util.ArrayList;
import java.util.Random;

import ij.measure.Calibration;
import ij.process.FloatProcessor;
import java.awt.image.ColorModel;
import ij.measure.ResultsTable;

import java.awt.Polygon;

public class Conn {

    ImagePlus 	img;

    int 		width;
    int 		height;
    int 		depth;

    ArrayList<Region> 	results;
    boolean 	saveLocations;

    public Conn(ImagePlus img, boolean saveLoc){

        this.img = img;

        width 	= 	img.getWidth();
        height 	= 	img.getHeight();
        depth 	= 	img.getStackSize();

        this.saveLocations = saveLoc;

        results 	= new ArrayList<Region>();

    }

    public void set(ImagePlus img, boolean saveLoc){

        this.img = img;

        width 	= 	img.getWidth();
        height 	= 	img.getHeight();
        depth 	= 	img.getStackSize();

        this.saveLocations = saveLoc;

        results 	= new ArrayList<Region>();

    }

    boolean pleaseStop = false;

    public void cancel() {
        pleaseStop = true;
    }

    /* An inner class to make the results list sortable. */
    private class Region implements Comparable {

        Region(int value, String materialName, int points, boolean sameValue) {
            byteImage = true;
            this.value = value;
            this.materialName = materialName;
            this.points = points;
            this.sameValue = sameValue;
            this.locations = new ArrayList<int[]>();
        }

        Region(int points, boolean sameValue) {
            byteImage = false;
            this.points = points;
            this.sameValue = sameValue;
            this.locations = new ArrayList<int[]>();
        }

        boolean byteImage;
        int points;
        String materialName;
        int value;
        boolean sameValue;

        ArrayList<int[]> 	locations;

        public void addLocations(ArrayList<int[]> loc_list){
            locations.clear();
            for (int i = 0; i < loc_list.size(); i++) {
                locations.add((int[])loc_list.get(i));
            }
        }

        public int[][] returnLocations(){
            int[][] out = new int[locations.size()][3];
            for (int i = 0; i < locations.size(); i++) {
                out[i] = (int[])locations.get(i);
            }
            return out;
        }

        public int compareTo(Object otherRegion) {
            Region o = (Region) otherRegion;
            return (points < o.points) ? -1 : ((points > o.points) ? 1 : 0);
        }

        @Override
        public String toString() {
            if (byteImage) {
                String materialBit = "";
                if (materialName != null) {
                    materialBit = " (" + materialName + ")";
                }
                return "Region of value " + value + materialBit + " containing " + points + " points";
            } else {
                return "Region containing " + points + " points";
            }
        }

        public void addRow( ResultsTable rt ) {
            rt.incrementCounter();
            if(byteImage) {
                if(sameValue)
                    rt.addValue("Value in Region",value);
                rt.addValue("Points In Region",points);
                if(materialName!=null)
                    rt.addLabel("Material Name",materialName);
            } else {
                rt.addValue("Points in Region",points);
            }
        }

    }

    private static final byte IN_QUEUE = 1;

    private static final byte ADDED = 2;

    public void run(String ignored) {

        boolean diagonal = true;//"Allow_diagonal connections?"
        boolean display = false;//"Display_image_for_each region?"
        boolean showResults = false;//"Display_results table?"
        boolean mustHaveSameValue = true;//"Regions_must have the same value?"
        boolean startFromPointROI = false;//"Start_from_point selection?"
        boolean autoSubtract = false;//"Autosubtract discovered regions from original image?"
        double valuesOverDouble = 0.0;//"Regions_for_values_over: "
        double minimumPointsInRegionDouble = 1;//"Minimum_number_of_points in a region"
        int stopAfterNumberOfRegions = -1;//"Stop_after this number of regions are found: "


        //ImageCalculator iCalc = new ImageCalculator();

        ImagePlus imagePlus = img;//IJ.getImage();
        if (imagePlus == null) {
            System.out.println("No image to operate on.");
            return;
        }

        int type = imagePlus.getType();

        if (!(ImagePlus.GRAY8 == type || ImagePlus.COLOR_256 == type || ImagePlus.GRAY32 == type)) {
            IJ.error("The image must be either 8 bit or 32 bit for this plugin.");
            return;
        }

        boolean byteImage = false;
        if (ImagePlus.GRAY8 == type || ImagePlus.COLOR_256 == type) {
            byteImage = true;
        }

        if (!byteImage && mustHaveSameValue) {
            IJ.error("You can only specify that each region must have the same value for 8 bit images.");
            return;
        }

        boolean startAtMaxValue = !mustHaveSameValue;

        int point_roi_x = -1;
        int point_roi_y = -1;
        int point_roi_z = -1;

        if( startFromPointROI ) {

            Roi roi = imagePlus.getRoi();
            if (roi == null) {
                IJ.error("There's no point selected in the image.");
                return;
            }
            if (roi.getType() != Roi.POINT) {
                IJ.error("There's a selection in the image, but it's not a point selection.");
                return;
            }
            Polygon p = roi.getPolygon();
            if(p.npoints > 1) {
                IJ.error("You can only have one point selected.");
                return;
            }

            point_roi_x = p.xpoints[0];
            point_roi_y = p.ypoints[0];
            point_roi_z = imagePlus.getCurrentSlice()-1;

            System.out.println("Fetched ROI with co-ordinates: "+p.xpoints[0]+", "+p.ypoints[0]);
        }



        if (width * height * depth > Integer.MAX_VALUE) {
            IJ.error("This stack is too large for this plugin (must have less than " + Integer.MAX_VALUE + " points.");
            return;
        }

//        String[] materialList = null;

//        AmiraParameters parameters = null;
//        if (AmiraParameters.isAmiraLabelfield(imagePlus)) {
//            parameters = new AmiraParameters(imagePlus);
//            materialList = parameters.getMaterialList();
//        }

        ImageStack stack = imagePlus.getStack();

        byte[][] sliceDataBytes = null;
        float[][] sliceDataFloats = null;

        if (byteImage) {
            sliceDataBytes = new byte[depth][];
            for (int z = 0; z < depth; ++z) {
                ByteProcessor bp = (ByteProcessor) stack.getProcessor(z+1);
                sliceDataBytes[z] = (byte[]) bp.getPixelsCopy();
            }
        } else {
            sliceDataFloats = new float[depth][];
            for (int z = 0; z < depth; ++z) {
                FloatProcessor bp = (FloatProcessor) stack.getProcessor(z+1);
                sliceDataFloats[z] = (float[]) bp.getPixelsCopy();
            }
        }

        // Preserve the calibration and colour lookup tables
        // for generating new images of each individual
        // region.
        Calibration calibration = imagePlus.getCalibration();

        ColorModel cm = null;
        if (ImagePlus.COLOR_256 == type) {
            cm = stack.getColorModel();
        }

        ResultsTable rt=ResultsTable.getResultsTable();
        rt.reset();

//		CancelDialog cancelDialog=new CancelDialog(this);
//		cancelDialog.show();

        boolean firstTime = true;

        while (true) {

            if( pleaseStop )
                break;

			/* Find one pixel that's above the minimum, or
			   find the maximum in the case where we're
			   not insisting that all regions are made up
			   of the same color.  These are set in all
			   cases... */

            int initial_x = -1;
            int initial_y = -1;
            int initial_z = -1;

            int foundValueInt = -1;
            float foundValueFloat = Float.MIN_VALUE;
            int maxValueInt = -1;
            float maxValueFloat = Float.MIN_VALUE;

            if (firstTime && startFromPointROI ) {

                initial_x = point_roi_x;
                initial_y = point_roi_y;
                initial_z = point_roi_z;

                if(byteImage)
                    foundValueInt = sliceDataBytes[initial_z][initial_y * width + initial_x] & 0xFF;
                else
                    foundValueFloat = sliceDataFloats[initial_z][initial_y * width + initial_x];

            } else if (byteImage && startAtMaxValue) {

                for (int z = 0; z < depth; ++z) {
                    for (int y = 0; y < height; ++y) {
                        for (int x = 0; x < width; ++x) {
                            int value = sliceDataBytes[z][y * width + x] & 0xFF;
                            if (value > maxValueInt) {
                                initial_x = x;
                                initial_y = y;
                                initial_z = z;
                                maxValueInt = value;
                            }
                        }
                    }
                }

                foundValueInt = maxValueInt;

				/* If the maximum value is below the
				   level we care about, we're done. */

                if (foundValueInt < valuesOverDouble) {
                    break;
                }

            } else if (byteImage && !startAtMaxValue) {

                // Just finding some point in the a region...
                for (int z = 0; z < depth && foundValueInt == -1; ++z) {
                    for (int y = 0; y < height && foundValueInt == -1; ++y) {
                        for (int x = 0; x < width; ++x) {
                            int value = sliceDataBytes[z][y * width + x] & 0xFF;
                            if (value > valuesOverDouble) {

                                initial_x = x;
                                initial_y = y;
                                initial_z = z;
                                foundValueInt = value;
                                break;
                            }
                        }
                    }
                }

                if (foundValueInt == -1) {
                    break;
                }

            } else {

                // This must be a 32 bit image and we're starting at the maximum
                assert (!byteImage && startAtMaxValue);

                for (int z = 0; z < depth; ++z) {
                    for (int y = 0; y < height; ++y) {
                        for (int x = 0; x < width; ++x) {
                            float value = sliceDataFloats[z][y * width + x];
                            if (value > valuesOverDouble) {
                                initial_x = x;
                                initial_y = y;
                                initial_z = z;
                                maxValueFloat = value;
                            }
                        }
                    }
                }

                foundValueFloat = maxValueFloat;

                if (foundValueFloat == Float.MIN_VALUE) {
                    break;

                    // If the maximum value is below the level we
                    // care about, we're done.
                }
                if (foundValueFloat < valuesOverDouble) {
                    break;
                }

            }

            firstTime = false;

            int vint = foundValueInt;
            //float vfloat = foundValueFloat;

            String materialName = null;
//            if (materialList != null) {
//                materialName = materialList[vint];
//            }
            int pointsInQueue = 0;
            int queueArrayLength = 1024;
            int[] queue = new int[queueArrayLength];

            byte[] pointState = new byte[depth * width * height];
            int i = width * (initial_z * height + initial_y) + initial_x;
            pointState[i] = IN_QUEUE;
            queue[pointsInQueue++] = i;

            int pointsInThisRegion = 0;

            while (pointsInQueue > 0) {

                if(pleaseStop)
                    break;

                int nextIndex = queue[--pointsInQueue];

                int currentPointStateIndex = nextIndex;
                int pz = nextIndex / (width * height);
                int currentSliceIndex = nextIndex % (width * height);
                int py = currentSliceIndex / width;
                int px = currentSliceIndex % width;

                pointState[currentPointStateIndex] = ADDED;

                if (byteImage) {
                    sliceDataBytes[pz][currentSliceIndex] = 0;
                } else {
                    sliceDataFloats[pz][currentSliceIndex] = Float.MIN_VALUE;
                }
                ++pointsInThisRegion;

                int x_unchecked_min = px - 1;
                int y_unchecked_min = py - 1;
                int z_unchecked_min = pz - 1;

                int x_unchecked_max = px + 1;
                int y_unchecked_max = py + 1;
                int z_unchecked_max = pz + 1;

                int x_min = (x_unchecked_min < 0) ? 0 : x_unchecked_min;
                int y_min = (y_unchecked_min < 0) ? 0 : y_unchecked_min;
                int z_min = (z_unchecked_min < 0) ? 0 : z_unchecked_min;

                int x_max = (x_unchecked_max >= width) ? width - 1 : x_unchecked_max;
                int y_max = (y_unchecked_max >= height) ? height - 1 : y_unchecked_max;
                int z_max = (z_unchecked_max >= depth) ? depth - 1 : z_unchecked_max;

                for (int z = z_min; z <= z_max; ++z) {
                    for (int y = y_min; y <= y_max; ++y) {
                        for (int x = x_min; x <= x_max; ++x) {

                            // If we're not including diagonals,
                            // skip those points.
                            if ((!diagonal) && (x == x_unchecked_min || x == x_unchecked_max) && (y == y_unchecked_min || y == y_unchecked_max) && (z == z_unchecked_min || z == z_unchecked_max)) {
                                continue;
                            }
                            int newSliceIndex = y * width + x;
                            int newPointStateIndex = width * (z * height + y) + x;

                            if (byteImage) {

                                int neighbourValue = sliceDataBytes[z][newSliceIndex] & 0xFF;

                                if (mustHaveSameValue) {
                                    if (neighbourValue != vint) {
                                        continue;
                                    }
                                } else {
                                    if (neighbourValue <= valuesOverDouble) {
                                        continue;
                                    }
                                }
                            } else {

                                float neighbourValue = sliceDataFloats[z][newSliceIndex];

                                if (neighbourValue <= valuesOverDouble) {
                                    continue;
                                }
                            }

                            if (0 == pointState[newPointStateIndex]) {
                                pointState[newPointStateIndex] = IN_QUEUE;
                                if (pointsInQueue == queueArrayLength) {
                                    int newArrayLength = (int) (queueArrayLength * 1.2);
                                    int[] newArray = new int[newArrayLength];
                                    System.arraycopy(queue, 0, newArray, 0, pointsInQueue);
                                    queue = newArray;
                                    queueArrayLength = newArrayLength;
                                }
                                queue[pointsInQueue++] = newPointStateIndex;
                            }
                        }
                    }
                }
            }

            if(pleaseStop)
                break;

            // So now pointState should have no IN_QUEUE
            // status points...
            Region region;
            if (byteImage) {
                region = new Region(vint, materialName, pointsInThisRegion, mustHaveSameValue );
            } else {
                region = new Region(pointsInThisRegion, mustHaveSameValue);
            }
            if (pointsInThisRegion < minimumPointsInRegionDouble) {
                // System.out.println("Too few points - only " + pointsInThisRegion);
                continue;
            }

            byte replacementValue;
            if (byteImage) {
                replacementValue = (byte) ( (cm == null) ? 255 : vint );
            } else {
                replacementValue = (byte) 255;
            }

            if (display || autoSubtract) {

                ImageStack newStack = new ImageStack(width, height);
                for (int z = 0; z < depth; ++z) {
                    byte[] sliceBytes = new byte[width * height];
                    for (int y = 0; y < height; ++y) {
                        for (int x = 0; x < width; ++x) {

                            byte status = pointState[width * (z * height + y) + x];

                            if (status == IN_QUEUE) {
                                IJ.log("BUG: point " + x + "," + y + "," + z + " is still marked as IN_QUEUE");
                            }

                            if (status == ADDED) {
                                sliceBytes[y * width + x] = replacementValue;
                            }
                        }
                    }
                    ByteProcessor bp = new ByteProcessor(width, height);
                    bp.setPixels(sliceBytes);
                    newStack.addSlice("", bp);
                }

                if (ImagePlus.COLOR_256 == type) {
                    if (cm != null) {
                        newStack.setColorModel(cm);
                    }
                }

                ImagePlus newImagePlus = new ImagePlus(region.toString(), newStack);

                if (calibration != null) {
                    newImagePlus.setCalibration(calibration);
                }

//                if (parameters != null) {
//                    parameters.setParameters(newImagePlus, true);
//                }

                if (autoSubtract) {
                    //iCalc.calculate("Subtract stack", imagePlus, newImagePlus);
                }

                if (display) {
                    newImagePlus.show();
                } else {
                    newImagePlus.changes = false;
                    newImagePlus.close();
                }
            }

            if(saveLocations){

                // make a location list for this region
                ArrayList<int[]> loc_list = new ArrayList<int[]>();

                for (int z = 0; z < depth; ++z) {
                    //byte[] sliceBytes = new byte[width * height];
                    for (int y = 0; y < height; ++y) {
                        for (int x = 0; x < width; ++x) {
                            byte status = pointState[width * (z * height + y) + x];

                            if (status == IN_QUEUE) {
                                System.out.println("BUG: point " + x + "," + y + "," + z + " is still marked as IN_QUEUE");
                            }

                            if (status == ADDED) {
                                loc_list.add(new int[]{y, x, z});
                                //labels[z*(height*width)+y*width+x] = results.size();
                                //((ArrayList)locations.get(index)).add(new int[]{y, x, z});
                            }
                        }
                    }
                }

                region.addLocations(loc_list);

            }

            // add it wit or without the locations
            results.add(region);

            if ( (stopAfterNumberOfRegions > 0) && (results.size() >= stopAfterNumberOfRegions) ) {
                break;
            }

        } // while

        Collections.sort(results, Collections.reverseOrder());

//		cancelDialog.dispose();

        for (Iterator<Region> it = results.iterator(); it.hasNext();) {
            Region r = it.next();
//			System.out.println(r.toString());
//			int[][] l = r.returnLocations();
//			for (int i = 0; i < l.length; i++) {
//				System.out.print("| "+l[i][0]+" , "+l[i][1]+" , "+l[i][2]+" | ");
//			}

            if( showResults ) {
                r.addRow(rt);
            }
        }

        if( showResults )
            rt.show("Results");

    }

    private byte[] randomColorGen(){
        // returns r,g,b in bytes
        Random generator = new Random();

        int randomIndex = generator.nextInt(256);

        byte[][] jet256 = new byte[][]{
                {(byte)   0, (byte)   0, (byte) 131},  {(byte)   0, (byte)   0, (byte) 135},
                {(byte)   0, (byte)   0, (byte) 139},  {(byte)   0, (byte)   0, (byte) 143},
                {(byte)   0, (byte)   0, (byte) 147},  {(byte)   0, (byte)   0, (byte) 151},
                {(byte)   0, (byte)   0, (byte) 155},  {(byte)   0, (byte)   0, (byte) 159},
                {(byte)   0, (byte)   0, (byte) 163},  {(byte)   0, (byte)   0, (byte) 167},
                {(byte)   0, (byte)   0, (byte) 171},  {(byte)   0, (byte)   0, (byte) 175},
                {(byte)   0, (byte)   0, (byte) 179},  {(byte)   0, (byte)   0, (byte) 183},
                {(byte)   0, (byte)   0, (byte) 187},  {(byte)   0, (byte)   0, (byte) 191},
                {(byte)   0, (byte)   0, (byte) 195},  {(byte)   0, (byte)   0, (byte) 199},
                {(byte)   0, (byte)   0, (byte) 203},  {(byte)   0, (byte)   0, (byte) 207},
                {(byte)   0, (byte)   0, (byte) 211},  {(byte)   0, (byte)   0, (byte) 215},
                {(byte)   0, (byte)   0, (byte) 219},  {(byte)   0, (byte)   0, (byte) 223},
                {(byte)   0, (byte)   0, (byte) 227},  {(byte)   0, (byte)   0, (byte) 231},
                {(byte)   0, (byte)   0, (byte) 235},  {(byte)   0, (byte)   0, (byte) 239},
                {(byte)   0, (byte)   0, (byte) 243},  {(byte)   0, (byte)   0, (byte) 247},
                {(byte)   0, (byte)   0, (byte) 251},  {(byte)   0, (byte)   0, (byte) 255},
                {(byte)   0, (byte)   4, (byte) 255},  {(byte)   0, (byte)   8, (byte) 255},
                {(byte)   0, (byte)  12, (byte) 255},  {(byte)   0, (byte)  16, (byte) 255},
                {(byte)   0, (byte)  20, (byte) 255},  {(byte)   0, (byte)  24, (byte) 255},
                {(byte)   0, (byte)  28, (byte) 255},  {(byte)   0, (byte)  32, (byte) 255},
                {(byte)   0, (byte)  36, (byte) 255},  {(byte)   0, (byte)  40, (byte) 255},
                {(byte)   0, (byte)  44, (byte) 255},  {(byte)   0, (byte)  48, (byte) 255},
                {(byte)   0, (byte)  52, (byte) 255},  {(byte)   0, (byte)  56, (byte) 255},
                {(byte)   0, (byte)  60, (byte) 255},  {(byte)   0, (byte)  64, (byte) 255},
                {(byte)   0, (byte)  68, (byte) 255},  {(byte)   0, (byte)  72, (byte) 255},
                {(byte)   0, (byte)  76, (byte) 255},  {(byte)   0, (byte)  80, (byte) 255},
                {(byte)   0, (byte)  84, (byte) 255},  {(byte)   0, (byte)  88, (byte) 255},
                {(byte)   0, (byte)  92, (byte) 255},  {(byte)   0, (byte)  96, (byte) 255},
                {(byte)   0, (byte) 100, (byte) 255},  {(byte)   0, (byte) 104, (byte) 255},
                {(byte)   0, (byte) 108, (byte) 255},  {(byte)   0, (byte) 112, (byte) 255},
                {(byte)   0, (byte) 116, (byte) 255},  {(byte)   0, (byte) 120, (byte) 255},
                {(byte)   0, (byte) 124, (byte) 255},  {(byte)   0, (byte) 128, (byte) 255},
                {(byte)   0, (byte) 131, (byte) 255},  {(byte)   0, (byte) 135, (byte) 255},
                {(byte)   0, (byte) 139, (byte) 255},  {(byte)   0, (byte) 143, (byte) 255},
                {(byte)   0, (byte) 147, (byte) 255},  {(byte)   0, (byte) 151, (byte) 255},
                {(byte)   0, (byte) 155, (byte) 255},  {(byte)   0, (byte) 159, (byte) 255},
                {(byte)   0, (byte) 163, (byte) 255},  {(byte)   0, (byte) 167, (byte) 255},
                {(byte)   0, (byte) 171, (byte) 255},  {(byte)   0, (byte) 175, (byte) 255},
                {(byte)   0, (byte) 179, (byte) 255},  {(byte)   0, (byte) 183, (byte) 255},
                {(byte)   0, (byte) 187, (byte) 255},  {(byte)   0, (byte) 191, (byte) 255},
                {(byte)   0, (byte) 195, (byte) 255},  {(byte)   0, (byte) 199, (byte) 255},
                {(byte)   0, (byte) 203, (byte) 255},  {(byte)   0, (byte) 207, (byte) 255},
                {(byte)   0, (byte) 211, (byte) 255},  {(byte)   0, (byte) 215, (byte) 255},
                {(byte)   0, (byte) 219, (byte) 255},  {(byte)   0, (byte) 223, (byte) 255},
                {(byte)   0, (byte) 227, (byte) 255},  {(byte)   0, (byte) 231, (byte) 255},
                {(byte)   0, (byte) 235, (byte) 255},  {(byte)   0, (byte) 239, (byte) 255},
                {(byte)   0, (byte) 243, (byte) 255},  {(byte)   0, (byte) 247, (byte) 255},
                {(byte)   0, (byte) 251, (byte) 255},  {(byte)   0, (byte) 255, (byte) 255},
                {(byte)   4, (byte) 255, (byte) 251},  {(byte)   8, (byte) 255, (byte) 247},
                {(byte)  12, (byte) 255, (byte) 243},  {(byte)  16, (byte) 255, (byte) 239},
                {(byte)  20, (byte) 255, (byte) 235},  {(byte)  24, (byte) 255, (byte) 231},
                {(byte)  28, (byte) 255, (byte) 227},  {(byte)  32, (byte) 255, (byte) 223},
                {(byte)  36, (byte) 255, (byte) 219},  {(byte)  40, (byte) 255, (byte) 215},
                {(byte)  44, (byte) 255, (byte) 211},  {(byte)  48, (byte) 255, (byte) 207},
                {(byte)  52, (byte) 255, (byte) 203},  {(byte)  56, (byte) 255, (byte) 199},
                {(byte)  60, (byte) 255, (byte) 195},  {(byte)  64, (byte) 255, (byte) 191},
                {(byte)  68, (byte) 255, (byte) 187},  {(byte)  72, (byte) 255, (byte) 183},
                {(byte)  76, (byte) 255, (byte) 179},  {(byte)  80, (byte) 255, (byte) 175},
                {(byte)  84, (byte) 255, (byte) 171},  {(byte)  88, (byte) 255, (byte) 167},
                {(byte)  92, (byte) 255, (byte) 163},  {(byte)  96, (byte) 255, (byte) 159},
                {(byte) 100, (byte) 255, (byte) 155},  {(byte) 104, (byte) 255, (byte) 151},
                {(byte) 108, (byte) 255, (byte) 147},  {(byte) 112, (byte) 255, (byte) 143},
                {(byte) 116, (byte) 255, (byte) 139},  {(byte) 120, (byte) 255, (byte) 135},
                {(byte) 124, (byte) 255, (byte) 131},  {(byte) 128, (byte) 255, (byte) 128},
                {(byte) 131, (byte) 255, (byte) 124},  {(byte) 135, (byte) 255, (byte) 120},
                {(byte) 139, (byte) 255, (byte) 116},  {(byte) 143, (byte) 255, (byte) 112},
                {(byte) 147, (byte) 255, (byte) 108},  {(byte) 151, (byte) 255, (byte) 104},
                {(byte) 155, (byte) 255, (byte) 100},  {(byte) 159, (byte) 255, (byte)  96},
                {(byte) 163, (byte) 255, (byte)  92},  {(byte) 167, (byte) 255, (byte)  88},
                {(byte) 171, (byte) 255, (byte)  84},  {(byte) 175, (byte) 255, (byte)  80},
                {(byte) 179, (byte) 255, (byte)  76},  {(byte) 183, (byte) 255, (byte)  72},
                {(byte) 187, (byte) 255, (byte)  68},  {(byte) 191, (byte) 255, (byte)  64},
                {(byte) 195, (byte) 255, (byte)  60},  {(byte) 199, (byte) 255, (byte)  56},
                {(byte) 203, (byte) 255, (byte)  52},  {(byte) 207, (byte) 255, (byte)  48},
                {(byte) 211, (byte) 255, (byte)  44},  {(byte) 215, (byte) 255, (byte)  40},
                {(byte) 219, (byte) 255, (byte)  36},  {(byte) 223, (byte) 255, (byte)  32},
                {(byte) 227, (byte) 255, (byte)  28},  {(byte) 231, (byte) 255, (byte)  24},
                {(byte) 235, (byte) 255, (byte)  20},  {(byte) 239, (byte) 255, (byte)  16},
                {(byte) 243, (byte) 255, (byte)  12},  {(byte) 247, (byte) 255, (byte)   8},
                {(byte) 251, (byte) 255, (byte)   4},  {(byte) 255, (byte) 255, (byte)   0},
                {(byte) 255, (byte) 251, (byte)   0},  {(byte) 255, (byte) 247, (byte)   0},
                {(byte) 255, (byte) 243, (byte)   0},  {(byte) 255, (byte) 239, (byte)   0},
                {(byte) 255, (byte) 235, (byte)   0},  {(byte) 255, (byte) 231, (byte)   0},
                {(byte) 255, (byte) 227, (byte)   0},  {(byte) 255, (byte) 223, (byte)   0},
                {(byte) 255, (byte) 219, (byte)   0},  {(byte) 255, (byte) 215, (byte)   0},
                {(byte) 255, (byte) 211, (byte)   0},  {(byte) 255, (byte) 207, (byte)   0},
                {(byte) 255, (byte) 203, (byte)   0},  {(byte) 255, (byte) 199, (byte)   0},
                {(byte) 255, (byte) 195, (byte)   0},  {(byte) 255, (byte) 191, (byte)   0},
                {(byte) 255, (byte) 187, (byte)   0},  {(byte) 255, (byte) 183, (byte)   0},
                {(byte) 255, (byte) 179, (byte)   0},  {(byte) 255, (byte) 175, (byte)   0},
                {(byte) 255, (byte) 171, (byte)   0},  {(byte) 255, (byte) 167, (byte)   0},
                {(byte) 255, (byte) 163, (byte)   0},  {(byte) 255, (byte) 159, (byte)   0},
                {(byte) 255, (byte) 155, (byte)   0},  {(byte) 255, (byte) 151, (byte)   0},
                {(byte) 255, (byte) 147, (byte)   0},  {(byte) 255, (byte) 143, (byte)   0},
                {(byte) 255, (byte) 139, (byte)   0},  {(byte) 255, (byte) 135, (byte)   0},
                {(byte) 255, (byte) 131, (byte)   0},  {(byte) 255, (byte) 128, (byte)   0},
                {(byte) 255, (byte) 124, (byte)   0},  {(byte) 255, (byte) 120, (byte)   0},
                {(byte) 255, (byte) 116, (byte)   0},  {(byte) 255, (byte) 112, (byte)   0},
                {(byte) 255, (byte) 108, (byte)   0},  {(byte) 255, (byte) 104, (byte)   0},
                {(byte) 255, (byte) 100, (byte)   0},  {(byte) 255, (byte)  96, (byte)   0},
                {(byte) 255, (byte)  92, (byte)   0},  {(byte) 255, (byte)  88, (byte)   0},
                {(byte) 255, (byte)  84, (byte)   0},  {(byte) 255, (byte)  80, (byte)   0},
                {(byte) 255, (byte)  76, (byte)   0},  {(byte) 255, (byte)  72, (byte)   0},
                {(byte) 255, (byte)  68, (byte)   0},  {(byte) 255, (byte)  64, (byte)   0},
                {(byte) 255, (byte)  60, (byte)   0},  {(byte) 255, (byte)  56, (byte)   0},
                {(byte) 255, (byte)  52, (byte)   0},  {(byte) 255, (byte)  48, (byte)   0},
                {(byte) 255, (byte)  44, (byte)   0},  {(byte) 255, (byte)  40, (byte)   0},
                {(byte) 255, (byte)  36, (byte)   0},  {(byte) 255, (byte)  32, (byte)   0},
                {(byte) 255, (byte)  28, (byte)   0},  {(byte) 255, (byte)  24, (byte)   0},
                {(byte) 255, (byte)  20, (byte)   0},  {(byte) 255, (byte)  16, (byte)   0},
                {(byte) 255, (byte)  12, (byte)   0},  {(byte) 255, (byte)   8, (byte)   0},
                {(byte) 255, (byte)   4, (byte)   0},  {(byte) 255, (byte)   0, (byte)   0},
                {(byte) 251, (byte)   0, (byte)   0},  {(byte) 247, (byte)   0, (byte)   0},
                {(byte) 243, (byte)   0, (byte)   0},  {(byte) 239, (byte)   0, (byte)   0},
                {(byte) 235, (byte)   0, (byte)   0},  {(byte) 231, (byte)   0, (byte)   0},
                {(byte) 227, (byte)   0, (byte)   0},  {(byte) 223, (byte)   0, (byte)   0},
                {(byte) 219, (byte)   0, (byte)   0},  {(byte) 215, (byte)   0, (byte)   0},
                {(byte) 211, (byte)   0, (byte)   0},  {(byte) 207, (byte)   0, (byte)   0},
                {(byte) 203, (byte)   0, (byte)   0},  {(byte) 199, (byte)   0, (byte)   0},
                {(byte) 195, (byte)   0, (byte)   0},  {(byte) 191, (byte)   0, (byte)   0},
                {(byte) 187, (byte)   0, (byte)   0},  {(byte) 183, (byte)   0, (byte)   0},
                {(byte) 179, (byte)   0, (byte)   0},  {(byte) 175, (byte)   0, (byte)   0},
                {(byte) 171, (byte)   0, (byte)   0},  {(byte) 167, (byte)   0, (byte)   0},
                {(byte) 163, (byte)   0, (byte)   0},  {(byte) 159, (byte)   0, (byte)   0},
                {(byte) 155, (byte)   0, (byte)   0},  {(byte) 151, (byte)   0, (byte)   0},
                {(byte) 147, (byte)   0, (byte)   0},  {(byte) 143, (byte)   0, (byte)   0},
                {(byte) 139, (byte)   0, (byte)   0},  {(byte) 135, (byte)   0, (byte)   0},
                {(byte) 131, (byte)   0, (byte)   0},  {(byte) 128, (byte)   0, (byte)   0}
        };

        return jet256[randomIndex];
    }

    public boolean connected(int[] pt1, int[] pt2){
        // TODO: check that they're within image
        int r1 = getRank(pt1);
        int r2 = getRank(pt2);
        return (r1==r2) && (r1!=-1);

    }

    public boolean belongsToBiggestRegion(int[] pt){
        // TODO: check that they're within image
        int r = getRank(pt);
        return (r==0);

    }

    public int getRank(int[] pt){
        int rank = -1;
        for (Iterator<Region> it = results.iterator(); it.hasNext();) {
            Region r = it.next();
            rank++;
            int[][] locs = r.returnLocations();
            for (int i = 0; i < locs.length; i++) {
                if(equal(locs[i], pt)){
                    return rank;
                }
            }
        }
        return -1;

    }

    public int getNrConnectedRegions(){
        // return clusters here
        return results.size();
    }

    public ArrayList<ArrayList<int[]>> getConnectedRegions()
    {
        ArrayList<ArrayList<int[]>> regs = new ArrayList<ArrayList<int[]>>(results.size());

        for (int i=0; i<results.size(); i++) {
            ArrayList<int[]> locs = new ArrayList<int[]>(results.get(i).locations.size());
            for (int j=0; j<results.get(i).locations.size(); j++) {
                locs.add(new int[]{results.get(i).locations.get(j)[0], results.get(i).locations.get(j)[1]});
            }
            regs.add(locs);
        }

        return regs;
    }

    public ArrayList<ArrayList<int[]>> getConnectedRegions3D_XYZ()
    {

        ArrayList<ArrayList<int[]>> regs = new ArrayList<ArrayList<int[]>>(results.size());

        for (int i=0; i<results.size(); i++) {
            ArrayList<int[]> locs = new ArrayList<int[]>(results.get(i).locations.size());
            for (int j=0; j<results.get(i).locations.size(); j++) {
                locs.add(new int[]{
                        results.get(i).locations.get(j)[1],
                        results.get(i).locations.get(j)[0],
                        results.get(i).locations.get(j)[2]
                });
            }
            regs.add(locs);
        }

        return regs;


    }

    public ImagePlus showLabels(){

        byte[][] bytes_red = new byte[height][width * height];
        byte[][] bytes_grn = new byte[height][width * height];
        byte[][] bytes_blu = new byte[height][width * height];

        for (Iterator<Region> it = results.iterator(); it.hasNext();) {
            Region r = it.next();
            byte[] color = randomColorGen();
            int[][] locs = r.returnLocations();
            for (int i = 0; i < locs.length; i++) {
                int row = locs[i][0];
                int col = locs[i][1];
                int lay = locs[i][2];
                bytes_red[lay][row * width + col] = color[0];
                bytes_grn[lay][row * width + col] = color[1];
                bytes_blu[lay][row * width + col] = color[2];
            }
        }

        ImageStack newStack = new ImageStack(width, height);
        for (int z = 0; z < depth; ++z) {
            ColorProcessor cp = new ColorProcessor(width, height);
            cp.setRGB(bytes_red[z], bytes_grn[z], bytes_blu[z]);
            newStack.addSlice("", cp);
        }
        ImagePlus outImage = new ImagePlus(img.getShortTitle(), newStack);
        return outImage;

    }

    private boolean equal(int[] pt1, int[] pt2){
        return (pt1[0]==pt2[0] && pt1[1]==pt2[1] && pt1[2]==pt2[2]);
    }

}